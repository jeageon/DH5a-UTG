from __future__ import annotations

from typing import Any, Optional
from io import StringIO

from Bio import SeqIO

from ..config import (
    NCBI_EFETCH,
    ENSEMBL_OVERLAP,
    ENSEMBL_OVERLAP_CHUNK_BP,
    ENSEMBL_OVERLAP_MAX_BP,
    DEFAULT_FEATURES,
    FeatureScanOptions,
)
from ..models.data_schemas import GenomicCoordinates, NegativeFeature
from ..utils import seq_utils
from ..utils.coord_utils import build_chunks, ensembl_to_relative
from ..utils.exceptions import ToolError
from ..utils.feature_utils import dedupe_features, merge_by_type
from ..utils.seq_utils import scan_ambiguous, scan_extreme_gc_windows, scan_homopolymers
from ..utils.api_client import ApiClient


def _first_non_null(*values: Any) -> Any:
    for value in values:
        if value is not None:
            return value
    return None


def _to_float(value: Any) -> Optional[float]:
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


class FeatureScanner:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def scan(
        self,
        coordinates: GenomicCoordinates,
        full_sequence: str,
        requested_features: Optional[list[str]] = None,
        options: Optional[FeatureScanOptions] = None,
    ) -> tuple[list[NegativeFeature], list[str]]:
        requested_features = requested_features or list(DEFAULT_FEATURES)
        options = options or FeatureScanOptions()
        warnings: list[str] = []
        seq_len = len(full_sequence)
        requested = set(requested_features)

        collected: list[NegativeFeature] = []
        if coordinates.coordinate_source == "ensembl":
            collected.extend(self._scan_overlap(coordinates, requested, seq_len, warnings, options))
        else:
            warnings.append("NCBI sequence source does not support Ensembl overlap lookup; skipped repeat/variant-based features")
            collected.extend(self._scan_ncbi_annotations(coordinates, requested, seq_len, warnings))
        collected.extend(self._scan_internal(
            full_sequence,
            requested,
            options,
            warnings,
        ))

        deduped = dedupe_features(collected)
        merge_gaps = {
            "extreme_gc": options.gc_step,
            "homopolymer": 0,
            "ambiguous": 0,
            "annotation": 0,
            "repeat": 0,
            "simple": 0,
            "variation": 0,
            "structural_variation": 0,
        }
        normalized = merge_by_type(deduped, merge_gaps=merge_gaps)
        return normalized, warnings

    def _scan_overlap(
        self,
        coordinates: GenomicCoordinates,
        requested: set[str],
        seq_len: int,
        warnings: list[str],
        options: FeatureScanOptions,
    ) -> list[NegativeFeature]:
        del warnings
        if not requested.intersection({"repeat", "simple", "variation", "structural_variation"}):
            return []
        features: list[NegativeFeature] = []
        region_start = coordinates.ext_start_1based
        region_end = coordinates.ext_end_1based
        region_len = region_end - region_start + 1
        if region_len <= ENSEMBL_OVERLAP_MAX_BP:
            chunks = [(region_start, region_end)]
        else:
            chunks = [(chunk.start, chunk.end) for chunk in build_chunks(region_start, region_end, ENSEMBL_OVERLAP_CHUNK_BP)]

        feature_types = requested.intersection({"repeat", "simple", "variation", "structural_variation"})
        for start, end in chunks:
            region = f"{coordinates.seq_region_name}:{start}..{end}:{coordinates.strand}"
            for ftype in sorted(feature_types):
                try:
                    params = {"feature": ftype}
                    url = ENSEMBL_OVERLAP.format(species=coordinates.species, region=region)
                    resp = self.api.get(
                        url,
                        headers={"Accept": "application/json"},
                        params=params,
                    )
                except ToolError:
                    continue
                if not isinstance(resp.json_obj, list):
                    continue
                for item in resp.json_obj:
                    feature = self._to_negative_feature(item, ftype, coordinates, seq_len, options.maf_threshold)
                    if feature is None:
                        continue
                    features.append(feature)
        return features

    def _scan_ncbi_annotations(
        self,
        coordinates: GenomicCoordinates,
        requested: set[str],
        seq_len: int,
        warnings: list[str],
    ) -> list[NegativeFeature]:
        if "annotation" not in requested or not coordinates.ncbi_accession:
            return []
        if not coordinates.ext_end_1based or coordinates.ext_end_1based < coordinates.ext_start_1based:
            return []
        if not coordinates.ncbi_accession:
            return []

        params = {
            "db": "nuccore",
            "id": coordinates.ncbi_accession,
            "rettype": "gbwithparts",
            "retmode": "text",
            "seq_start": coordinates.ext_start_1based,
            "seq_stop": coordinates.ext_end_1based,
        }
        if coordinates.strand == -1:
            params["strand"] = 2
        try:
            resp = self.api.get(
                NCBI_EFETCH,
                headers={"Accept": "text/plain"},
                params=params,
                disable_cache=True,
            )
        except ToolError as exc:
            warnings.append(f"NCBI annotation fetch failed: {exc}")
            return []

        try:
            records = list(SeqIO.parse(StringIO(resp.text or ""), "genbank"))
        except Exception as exc:
            warnings.append(f"failed to parse NCBI GenBank annotation: {exc}")
            return []
        if not records:
            return []
        record = records[0]

        results: list[NegativeFeature] = []
        for item in record.features:
            if not item.type or str(item.type).lower() == "source":
                continue
            try:
                raw_start = int(item.location.start)
                raw_end = int(item.location.end)
            except Exception:
                continue
            if raw_start < 0 or raw_end <= raw_start:
                continue

            # NCBI may return region-relative or absolute feature coordinates depending on service behavior.
            # Normalize to region-relative when needed.
            if raw_end > (coordinates.ext_end_1based - coordinates.ext_start_1based + 1) and raw_start >= coordinates.ext_start_1based:
                rel_start = raw_start - coordinates.ext_start_1based + 1
                rel_end = raw_end - coordinates.ext_start_1based + 1
            else:
                rel_start = raw_start
                rel_end = raw_end

            if rel_start < 0:
                rel_start = 0
            if rel_end > seq_len:
                rel_end = seq_len
            if rel_start >= rel_end:
                continue

            qualifiers = {str(k).lower(): v for k, v in item.qualifiers.items()}
            attrs: dict[str, Any] = {}
            for key, value in qualifiers.items():
                if isinstance(value, list) and value:
                    attrs[key] = value[0]
                elif value:
                    attrs[key] = value

            feature_name = _first_non_null(
                attrs.get("gene"),
                attrs.get("locus_tag"),
                attrs.get("product"),
                attrs.get("note"),
                item.type,
            )
            description = f"{item.type}: {feature_name}" if feature_name else str(item.type)

            results.append(
                NegativeFeature(
                    feature_type="annotation",
                    start=rel_start,
                    end=rel_end,
                    description=description,
                    source="ncbi_gb",
                    attributes=attrs,
                )
            )
        return results

    def _to_negative_feature(
        self,
        item: dict[str, Any],
        feature_type: str,
        coordinates: GenomicCoordinates,
        seq_len: int,
        maf_threshold: float,
    ) -> Optional[NegativeFeature]:
        if not isinstance(item, dict):
            return None

        start1 = _to_int(_first_non_null(item.get("start"), item.get("seq_region_start")), 0)
        end1 = _to_int(_first_non_null(item.get("end"), item.get("seq_region_end")), 0)
        if start1 <= 0 or end1 <= 0:
            return None

        maf = None
        if feature_type in {"variation", "structural_variation"}:
            raw_maf = _first_non_null(
                item.get("minor_allele_frequency"),
                item.get("minor_allele_frequency"),
                item.get("MAF"),
                item.get("maf"),
            )
            maf = _to_float(raw_maf)
            if maf is not None and maf < maf_threshold:
                return None

        rel_start, rel_end = ensembl_to_relative(start1, end1, coordinates.ext_start_1based, seq_len)
        if rel_start >= rel_end:
            return None

        fid = _first_non_null(item.get("id"), item.get("variant_accession"), item.get("variation_name"), "unknown")
        desc = {
            "repeat": f"repeat_region: {fid}",
            "simple": "simple_repeat",
            "variation": f"variant {fid}",
            "structural_variation": f"structural variation {fid}",
        }.get(feature_type, feature_type)
        if feature_type in {"variation", "structural_variation"} and maf is not None:
            desc += f", MAF={maf}"

        attrs = {}
        if feature_type == "variation":
            alleles = _first_non_null(item.get("alleles"), item.get("variant_alleles"), item.get("alleleString"))
            if alleles is not None:
                attrs["alleles"] = alleles
            consequence = _first_non_null(item.get("most_severe_consequence"), item.get("consequence_types"))
            if consequence is not None:
                attrs["consequence"] = consequence
            if fid != "unknown":
                attrs["id"] = fid
        return NegativeFeature(
            feature_type=feature_type,
            start=rel_start,
            end=rel_end,
            description=desc,
            source="ensembl_overlap",
            score=maf,
            strand=_to_int(item.get("strand")),
            attributes=attrs,
        )

    def _scan_internal(
        self,
        full_sequence: str,
        requested: set[str],
        options: FeatureScanOptions,
        warnings: list[str],
    ) -> list[NegativeFeature]:
        del warnings
        results: list[NegativeFeature] = []
        seq_len = len(full_sequence)
        if not requested:
            return results

        if "extreme_gc" in requested:
            windows = scan_extreme_gc_windows(
                full_sequence,
                window_size=options.gc_window,
                step=options.gc_step,
                gc_min=options.gc_min,
                gc_max=options.gc_max,
            )
            merged = seq_utils.merge_intervals_with_gap(windows, gap=options.gc_step)
            for start, end, gc in merged:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="extreme_gc",
                        start=start,
                        end=min(end, seq_len),
                        description=f"Extreme GC window(s): GC<{options.gc_min}% or GC>{options.gc_max}%",
                        source="internal_gc",
                        score=gc,
                    )
                )

        if "homopolymer" in requested:
            hits = scan_homopolymers(full_sequence, at_run=options.homopolymer_at, gc_run=options.homopolymer_gc)
            for base, start, end in hits:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="homopolymer",
                        start=start,
                        end=end,
                        description=f"Homopolymer run: {base}x{end-start}",
                        source="internal_regex",
                        score=float(end - start),
                    )
                )

        if "ambiguous" in requested:
            blocks = scan_ambiguous(full_sequence)
            for start, end in blocks:
                if start >= end or end > seq_len:
                    continue
                results.append(
                    NegativeFeature(
                        feature_type="ambiguous",
                        start=start,
                        end=end,
                        description="Ambiguous base(s) present",
                        source="internal_regex",
                    )
                )
        return results
