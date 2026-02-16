from __future__ import annotations

from dataclasses import dataclass
from io import StringIO
import re
from typing import Any, Optional

from Bio import SeqIO

from ..config import (
    DH5A_ACCESSION,
    DH5A_TAXID,
    ENSEMBL_LOOKUP,
    EBI_COORDINATES_URL,
    NCBI_EFETCH,
    NCBI_ESEARCH,
    NCBI_ESUMMARY,
    UNIPROT_ENTRY_URL,
    UNIPROT_IDMAP_RESULTS,
    UNIPROT_IDMAP_RUN,
    UNIPROT_IDMAP_STATUS,
)
from ..models.data_schemas import GenomicCoordinates
from ..utils.api_client import ApiClient
from ..utils.coord_utils import apply_flank
from ..utils.exceptions import NoMappingError, ToolError


@dataclass
class ResolverResult:
    coordinates: GenomicCoordinates
    warnings: list[str]


_ENSEMBL_GENE_ID_RE = re.compile(r"^ENS\w*G\d+$", re.IGNORECASE)
_UNIPROT_ACCESSION_RE = re.compile(r"^[A-Z]{1}[0-9][A-Z0-9]{4,11}$", re.IGNORECASE)


def _normalize_ensembl_gene_id(raw: Any) -> Optional[str]:
    if not isinstance(raw, str):
        return None
    cleaned = raw.strip()
    if not cleaned:
        return None
    cleaned = cleaned.split(".", 1)[0]
    return cleaned if _ENSEMBL_GENE_ID_RE.match(cleaned) else None


def _normalize_accession(raw: Optional[str]) -> str:
    if not raw:
        return ""
    return str(raw).strip().split(".", 1)[0].upper()


def _extract_gene_from_ref_value(ref: Any) -> Optional[str]:
    if not isinstance(ref, dict):
        return None

    direct = ref.get("id")
    normalized = _normalize_ensembl_gene_id(direct)
    if normalized:
        return normalized

    for prop in ref.get("properties", []) or []:
        if not isinstance(prop, dict):
            continue
        key = str(prop.get("key") or "").lower()
        if key == "geneid":
            normalized = _normalize_ensembl_gene_id(prop.get("value"))
            if normalized:
                return normalized
    return None


def _extract_gene_from_note(note: Any) -> Optional[str]:
    if note is None:
        return None
    text = str(note)
    matches = re.findall(
        r"\b(gene|locus_tag|old_locus_tag|gene_name)\s*[:=]\s*([A-Za-z0-9_\-\.]+)",
        text,
        flags=re.IGNORECASE,
    )
    if not matches:
        return None
    for key, value in matches:
        if key.lower() == "gene":
            return value
    for key, value in matches:
        if key.lower() in {"gene_name", "locus_tag", "old_locus_tag"}:
            return value
    return None


class CoordinateResolver:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def resolve(
        self,
        query: str,
        flank_bp: int,
        flank_mode: str = "genomic",
        query_type: str = "auto",
        taxid_filter: Optional[int] = None,
        ncbi_accession: str = DH5A_ACCESSION,
    ) -> ResolverResult:
        query = self._coerce_str(query)
        if not query:
            raise NoMappingError("Query is empty")

        taxid_filter = DH5A_TAXID if taxid_filter is None else taxid_filter
        ncbi_accession = self._coerce_str(ncbi_accession) or DH5A_ACCESSION
        requested_type = self._normalize_query_type(query_type, query)
        warnings: list[str] = []

        if requested_type == "gene_name":
            ncbi_lookup, ncbi_warnings = self._resolve_ncbi_gene_by_name(
                gene_name=query,
                taxid=taxid_filter,
                ncbi_accession=ncbi_accession,
            )
            warnings.extend(ncbi_warnings)
            if ncbi_lookup:
                return self._build_ncbi_result(
                    query=query,
                    ncbi_lookup=ncbi_lookup,
                    flank_bp=flank_bp,
                    flank_mode=flank_mode,
                    taxid_filter=taxid_filter,
                    warnings=warnings,
                    query_type="gene_name",
                )
            detail = "; ".join(warnings) if warnings else "no valid mapping found"
            raise NoMappingError(f"No NCBI mapping found for gene {query}: {detail}")

        entry_for_routing = self._fetch_uniprot_entry(query)
        use_ncbi_first = self._is_bacterial_entry(entry_for_routing)
        if query_type != "uniprot_id" and requested_type == "auto":
            use_ncbi_first = True

        if use_ncbi_first:
            warnings.append("microbial mode: attempting NCBI-first resolution")

        if use_ncbi_first:
            ncbi_lookup, ncbi_warnings = self._resolve_ncbi_gene(
                query,
                uniprot_entry=entry_for_routing,
                taxid_filter=taxid_filter,
                ncbi_accession=ncbi_accession,
            )
            warnings.extend(ncbi_warnings)
            if ncbi_lookup:
                return self._build_ncbi_result(
                    query=query,
                    ncbi_lookup=ncbi_lookup,
                    flank_bp=flank_bp,
                    flank_mode=flank_mode,
                    taxid_filter=taxid_filter,
                    warnings=warnings,
                    query_type="uniprot_id",
                )
            if _normalize_accession(ncbi_accession):
                detail = "; ".join(warnings) if warnings else "no valid NCBI mapping found"
                raise NoMappingError(
                    f"No NCBI mapping found for '{query}' on preferred accession {_normalize_accession(ncbi_accession)}: {detail}"
                )

        ensembl_gene_id, fallback_warnings = self._resolve_ensembl_gene(query)
        warnings.extend(fallback_warnings)
        if ensembl_gene_id:
            lookup = self._lookup_ensembl_gene(ensembl_gene_id)
            if lookup:
                return self._build_ensembl_result(
                    query=query,
                    lookup=lookup,
                    flank_bp=flank_bp,
                    flank_mode=flank_mode,
                    taxid_filter=taxid_filter,
                    warnings=warnings,
                )

        if not use_ncbi_first:
            ncbi_lookup, ncbi_warnings = self._resolve_ncbi_gene(
                query,
                uniprot_entry=entry_for_routing,
                taxid_filter=taxid_filter,
                ncbi_accession=ncbi_accession,
            )
            warnings.extend(ncbi_warnings)
            if ncbi_lookup:
                return self._build_ncbi_result(
                    query=query,
                    ncbi_lookup=ncbi_lookup,
                    flank_bp=flank_bp,
                    flank_mode=flank_mode,
                    taxid_filter=taxid_filter,
                    warnings=warnings,
                    query_type="uniprot_id",
                )

        detail = "; ".join(warnings) if warnings else "no valid mapping found"
        raise NoMappingError(f"No mapping found for {query}: {detail}")

    def _normalize_query_type(self, query_type: str, query: str) -> str:
        normalized = self._coerce_str(query_type).lower()
        if normalized in {"gene_name", "uniprot_id"}:
            return normalized
        if _UNIPROT_ACCESSION_RE.match(query):
            return "uniprot_id"
        return "gene_name"

    def _is_bacterial_entry(self, entry: Optional[dict[str, Any]]) -> bool:
        if not isinstance(entry, dict):
            return False
        organism = entry.get("organism") or {}
        if not isinstance(organism, dict):
            return False
        lineages = self._collect_lineages(organism)
        microbial_markers = {"bacteria", "archaea", "viral", "viruses"}
        if any(marker in lineages for marker in microbial_markers):
            return True
        organism_name = self._coerce_str(organism.get("scientificName") or organism.get("taxon-scientific-name"))
        return "bacteria" in organism_name.lower() or "archaea" in organism_name.lower()

    def _collect_lineages(self, organism: dict[str, Any]) -> set[str]:
        values: list[Any] = []
        for key in ("lineage", "lineages"):
            value = organism.get(key)
            if value is None:
                continue
            if isinstance(value, list):
                values.extend(value)
            else:
                values.append(value)
        normalized: set[str] = set()
        for item in values:
            if isinstance(item, dict):
                for key in ("scientificName", "name", "value", "taxon"):
                    candidate = self._coerce_str(item.get(key))
                    if candidate:
                        normalized.add(candidate.lower())
                continue
            candidate = self._coerce_str(item)
            if candidate:
                normalized.add(candidate.lower())
        return normalized

    def _build_ensembl_result(
        self,
        query: str,
        lookup: dict[str, Any],
        flank_bp: int,
        flank_mode: str,
        taxid_filter: Optional[int],
        warnings: list[str],
    ) -> ResolverResult:
        gene_start = lookup["gene_start_1based"]
        gene_end = lookup["gene_end_1based"]
        if gene_start is None or gene_end is None:
            raise NoMappingError(f"Invalid coordinate span for {lookup.get('ensembl_gene_id', query)}")
        if gene_start > gene_end:
            warnings.append(
                f"coordinate span was reordered for {lookup.get('ensembl_gene_id', query)}: {gene_start}>{gene_end}"
            )
            gene_start, gene_end = gene_end, gene_start
        if taxid_filter is not None and lookup.get("taxid") != taxid_filter:
            warnings.append(f"taxid mismatch: got {lookup.get('taxid')}, requested {taxid_filter}")

        strand = lookup["strand"]
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)
        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=query,
                query_type="uniprot_id",
                query_gene=self._coerce_str(lookup.get("display_name")),
                ensembl_gene_id=lookup["ensembl_gene_id"],
                species=lookup["species"],
                assembly_name=lookup["assembly_name"],
                seq_region_name=lookup["seq_region_name"],
                gene_start_1based=gene_start,
                gene_end_1based=gene_end,
                strand=strand,
                taxid=lookup.get("taxid"),
                ext_start_1based=ext_start,
                ext_end_1based=ext_end,
            ),
            warnings=warnings,
        )

    def _build_ncbi_result(
        self,
        query: str,
        ncbi_lookup: dict[str, Any],
        flank_bp: int,
        flank_mode: str,
        taxid_filter: Optional[int],
        warnings: list[str],
        query_type: str = "gene_name",
    ) -> ResolverResult:
        gene_start = ncbi_lookup["gene_start_1based"]
        gene_end = ncbi_lookup["gene_end_1based"]
        if gene_start is None or gene_end is None:
            raise NoMappingError(f"Invalid coordinate span for {query}")
        if gene_start > gene_end:
            warnings.append(f"coordinate span was reordered for {query}: {gene_start}>{gene_end}")
            gene_start, gene_end = gene_end, gene_start
        if taxid_filter is not None and ncbi_lookup.get("taxid") is not None and ncbi_lookup.get("taxid") != taxid_filter:
            warnings.append(f"taxid mismatch: got {ncbi_lookup.get('taxid')}, requested {taxid_filter}")

        strand = ncbi_lookup["strand"]
        genome_len = ncbi_lookup.get("ncbi_genome_length")
        if isinstance(genome_len, int) and genome_len > 0:
            warnings.append(f"NCBI genome length used for clamp: {genome_len} bp")
        ext_start, ext_end = apply_flank(gene_start, gene_end, flank_bp, flank_mode, strand)
        if isinstance(genome_len, int) and genome_len > 0 and ext_end > genome_len:
            ext_end = genome_len
            if ext_end < ext_start:
                ext_start = max(1, genome_len - (gene_end - gene_start + 1))
            warnings.append("NCBI region end was clamped to genome length")

        return ResolverResult(
            coordinates=GenomicCoordinates(
                uniprot_id=query,
                query_type=query_type,
                query_gene=self._coerce_str(ncbi_lookup.get("query_name")),
                ensembl_gene_id=ncbi_lookup["ensembl_gene_id"],
                coordinate_source="ncbi",
                ncbi_accession=ncbi_lookup["ncbi_accession"],
                ncbi_genome_length=genome_len if isinstance(genome_len, int) else None,
                species=ncbi_lookup["species"],
                assembly_name=ncbi_lookup["assembly_name"],
                seq_region_name=ncbi_lookup["seq_region_name"],
                gene_start_1based=gene_start,
                gene_end_1based=gene_end,
                strand=strand,
                display_name=self._coerce_str(ncbi_lookup.get("display_name")),
                taxid=ncbi_lookup.get("taxid"),
                ext_start_1based=ext_start,
                ext_end_1based=ext_end,
            ),
            warnings=warnings,
        )

    def _resolve_ensembl_gene(self, query: str) -> tuple[Optional[str], list[str]]:
        warnings: list[str] = []
        try:
            response = self.api.get(
                EBI_COORDINATES_URL.format(accession=query),
                headers={"Accept": "application/json"},
            )
            payload = response.json_obj or []
            if isinstance(payload, dict):
                payload = [payload]
            if not isinstance(payload, list):
                raise ToolError("Invalid EBI response structure")
            candidates: list[str] = []
            for item in payload:
                gene_id = self._extract_ensembl_gene_id(item)
                if gene_id:
                    candidates.append(gene_id)
            if candidates:
                if len(candidates) > 1:
                    warnings.append("multiple EBI coordinate candidates found; selected first valid one")
                return candidates[0], warnings
        except ToolError as exc:
            warnings.append(f"EBI coordinates lookup failed for {query}: {exc}")

        try:
            return self._fallback_uniprot_mapping(query)
        except ToolError as exc:
            warnings.append(f"UniProt mapping lookup failed for {query}: {exc}")
            return None, warnings

    def _extract_ensembl_gene_id(self, entry: dict[str, Any]) -> Optional[str]:
        direct = entry.get("ensemblGeneId") or entry.get("ensembl_gene_id")
        normalized = _normalize_ensembl_gene_id(direct)
        if normalized:
            return normalized
        for ref in entry.get("crossReferences", []) or []:
            db_name = str(ref.get("dbDisplayName", "")).lower()
            if "ensembl" not in db_name:
                continue
            normalized = _extract_gene_from_ref_value(ref)
            if normalized:
                return normalized
        return None

    def _parse_coords(self, entry: dict[str, Any]) -> dict[str, Any]:
        loc = entry.get("genomicLocation") or entry.get("genomic_location") or {}
        start = loc.get("start") or entry.get("genomicStart") or entry.get("start")
        end = loc.get("end") or entry.get("genomicEnd") or entry.get("end")
        chrom = (
            loc.get("chromosome")
            or loc.get("seqRegion")
            or loc.get("seqRegionName")
            or entry.get("seqRegion")
            or entry.get("seq_region_name")
        )
        return {
            "start": _to_int(start),
            "end": _to_int(end),
            "chrom": str(chrom) if chrom is not None else None,
            "taxid": _to_int(entry.get("organism", {}).get("taxid") or entry.get("taxId") or entry.get("taxid")),
            "strand": _to_int(entry.get("strand")) or 1,
        }

    def _fallback_uniprot_mapping(self, query: str) -> tuple[Optional[str], list[str]]:
        warnings: list[str] = []
        payload = {"from": "UniProtKB_AC-ID", "to": "Ensembl", "ids": query}
        run = self.api.post(
            UNIPROT_IDMAP_RUN,
            headers={"Content-Type": "application/x-www-form-urlencoded"},
            data=payload,
        )
        if not isinstance(run.json_obj, dict):
            warnings.append("UniProt mapping response structure is unexpected")
            fallback, _ = self._fallback_uniprot_crossrefs(query)
            return fallback, warnings
        job_id = run.json_obj.get("jobId") or run.json_obj.get("job_id")
        if not job_id:
            fallback, fallback_sources = self._fallback_uniprot_crossrefs(query)
            if not fallback and fallback_sources:
                warnings.append(
                    "mapping API did not return jobId; available crossrefs include: "
                    f"{', '.join(sorted(set(fallback_sources)))}"
                )
            if fallback:
                warnings.append("used UniProt cross-reference fallback mapping")
            return fallback, warnings

        status = _poll_uniprot_status(self.api, str(job_id))
        if not status:
            warnings.append("UniProt mapping job did not complete")
            fallback, fallback_sources = self._fallback_uniprot_crossrefs(query)
            if not fallback and fallback_sources:
                warnings.append(
                    "mapping job not completed; crossrefs available for "
                    f"{', '.join(sorted(set(fallback_sources)))}"
                )
            if fallback:
                warnings.append("used UniProt cross-reference fallback mapping")
            return fallback, warnings

        results = self.api.get(
            UNIPROT_IDMAP_RESULTS.format(job_id=job_id),
            headers={"Accept": "application/json"},
        )
        mapped_gene = self._extract_gene_from_mapping(results.json_obj)
        if mapped_gene:
            return mapped_gene, warnings
        fallback, fallback_sources = self._fallback_uniprot_crossrefs(query)
        if not fallback and fallback_sources:
            warnings.append(
                "mapping results lacked Ensembl gene ID; crossrefs include "
                f"{', '.join(sorted(set(fallback_sources)))}"
            )
        if fallback:
            warnings.append("used UniProt cross-reference fallback mapping")
        return fallback, warnings

    def _fallback_uniprot_crossrefs(self, query: str) -> tuple[Optional[str], list[str]]:
        try:
            response = self.api.get(
                UNIPROT_ENTRY_URL.format(accession=query),
                headers={"Accept": "application/json"},
            )
        except ToolError:
            return None, []

        if not isinstance(response.json_obj, dict):
            return None, []

        refs = response.json_obj.get("uniProtKBCrossReferences", [])
        if not isinstance(refs, list):
            return None, []

        ensembl_sources: list[str] = []
        for ref in refs:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            if "ensembl" not in db_name:
                continue
            ensembl_sources.append(db_name)
            normalized = _extract_gene_from_ref_value(ref)
            if normalized:
                return normalized, ensembl_sources
        return None, ensembl_sources

    def _resolve_ncbi_gene_by_name(
        self,
        gene_name: str,
        taxid: Optional[int],
        ncbi_accession: str,
    ) -> tuple[Optional[dict[str, Any]], list[str]]:
        warnings: list[str] = []
        query = self._coerce_str(gene_name)
        if not query:
            return None, ["empty gene name"]
        preferred_terms: list[str] = [f"{query}[Gene Name]"]
        if taxid is not None:
            preferred_terms.insert(0, f"{query}[Gene Name] AND {taxid}[Taxonomy ID]")
            preferred_terms.insert(1, f"{query}[All Fields] AND {taxid}[Taxonomy ID]")
        if ncbi_accession:
            preferred_terms.insert(0, f"{query}[Gene Name] AND {ncbi_accession}[Accession]")

        gene_ids: list[str] = []
        for term in preferred_terms:
            try:
                ids = self._ncbi_esearch_gene_ids(term)
            except ToolError:
                continue
            for gid in ids:
                if gid not in gene_ids:
                    gene_ids.append(gid)
            if gene_ids:
                break

        if not gene_ids:
            return None, ["NCBI gene search returned no ID"]

        try:
            summaries = self._ncbi_gene_summaries(gene_ids)
        except ToolError:
            return None, ["NCBI gene summary lookup failed"]

        chosen = self._choose_best_ncbi_gene_summary(
            summaries=summaries,
            aliases=[query],
            accessions=[ncbi_accession],
            prefer_accession_match=bool(ncbi_accession),
        )
        if not chosen:
            direct = self._resolve_ncbi_gene_from_accession_features(
                aliases=[query],
                accessions=[ncbi_accession],
                species_name="Escherichia coli",
            )
            if direct:
                return direct, [*warnings, "NCBI fallback: resolved via NCBI nuccore feature scan"]
            return None, ["NCBI gene summaries had no usable coordinates"]

        ncbi_coordinates = self._ncbi_summary_to_coordinates(
            chosen,
            [query],
            self._coerce_str(chosen.get("organism", {}).get("commonname")),
            accessions=[ncbi_accession],
        )
        if not ncbi_coordinates:
            direct = self._resolve_ncbi_gene_from_accession_features(
                aliases=[query],
                accessions=[ncbi_accession],
                species_name="Escherichia coli",
            )
            if direct:
                return direct, [*warnings, "NCBI fallback: resolved via NCBI nuccore feature scan"]
            return None, ["NCBI summary lacked usable genomic coordinates"]
        ncbi_coordinates["query_name"] = query
        return ncbi_coordinates, [*warnings, "NCBI fallback: resolved via gene name search"]

    def _resolve_ncbi_gene(
        self,
        query: str,
        uniprot_entry: Optional[dict[str, Any]] = None,
        taxid_filter: Optional[int] = None,
        ncbi_accession: str = DH5A_ACCESSION,
    ) -> tuple[Optional[dict[str, Any]], list[str]]:
        warnings: list[str] = []
        entry = uniprot_entry if isinstance(uniprot_entry, dict) else self._fetch_uniprot_entry(query)
        if not isinstance(entry, dict):
            return None, ["unable to read UniProt entry for NCBI fallback"]

        organism = entry.get("organism") or {}
        taxid = _to_int(organism.get("taxonId"))
        if taxid_filter is not None:
            taxid = taxid_filter
        organism_name = self._coerce_str(organism.get("scientificName"))
        gene_aliases = self._collect_ncbi_gene_aliases(entry)
        if not gene_aliases:
            return None, ["no gene alias found for NCBI fallback"]

        primary_ncbi_accessions = [self._coerce_str(ncbi_accession)] if self._coerce_str(ncbi_accession) else []
        ncbi_accessions = self._collect_ncbi_nucleotide_accessions(entry)
        if ncbi_accessions or primary_ncbi_accessions:
            combined_accessions = list(dict.fromkeys([*ncbi_accessions, *primary_ncbi_accessions]))
            if combined_accessions:
                warnings.append(
                    "NCBI fallback: collected RefSeq/EMBL accession hints: "
                    f"{', '.join(sorted(set(combined_accessions)))}"
                )
            ncbi_accessions = combined_accessions

        matching_accessions = primary_ncbi_accessions or ncbi_accessions

        summaries = self._collect_ncbi_gene_summaries(gene_aliases, taxid, organism_name)
        if not summaries:
            return None, ["NCBI fallback: no NCBI gene record found for candidate identifiers"]

        chosen = self._choose_best_ncbi_gene_summary(
            summaries=summaries,
            aliases=gene_aliases,
            accessions=matching_accessions,
            prefer_accession_match=bool(matching_accessions),
        )
        if not chosen:
            direct = self._resolve_ncbi_gene_from_accession_features(
                aliases=gene_aliases,
                accessions=matching_accessions,
                species_name=organism_name,
            )
            if direct:
                return direct, warnings + ["NCBI fallback: resolved via NCBI nuccore feature scan"]
            return None, ["NCBI fallback: gene records lacked usable genomic location"]

        ncbi_coordinates = self._ncbi_summary_to_coordinates(
            chosen,
            gene_aliases,
            organism_name,
            accessions=matching_accessions,
        )
        if not ncbi_coordinates:
            direct = self._resolve_ncbi_gene_from_accession_features(
                aliases=gene_aliases,
                accessions=matching_accessions,
                species_name=organism_name,
            )
            if direct:
                return direct, warnings + ["NCBI fallback: resolved via NCBI nuccore feature scan"]
            return None, ["NCBI fallback: failed to extract genomic coordinates from chosen NCBI summary"]
        return ncbi_coordinates, warnings + ["NCBI fallback: resolved via NCBI gene summary"]

    def _resolve_ncbi_gene_from_accession_features(
        self,
        aliases: list[str],
        accessions: list[str],
        species_name: Optional[str] = None,
    ) -> Optional[dict[str, Any]]:
        preferred_aliases: list[str] = []
        alias_name_map: dict[str, str] = {}
        for alias in aliases:
            value = self._coerce_str(alias)
            normalized = value.lower()
            if value and normalized not in alias_name_map:
                alias_name_map[normalized] = value
            if value and normalized not in preferred_aliases:
                preferred_aliases.append(normalized)
        if not preferred_aliases:
            return None

        preferred_accessions = [self._coerce_str(acc) for acc in accessions if self._coerce_str(acc)]
        if not preferred_accessions:
            return None

        feature_priority = {"gene": 0, "cds": 1, "mrna": 2, "rrna": 3, "trna": 3}
        for accession in preferred_accessions:
            record = self._fetch_ncbi_record_by_accession(accession)
            if record is None:
                continue

            selected: Optional[tuple[int, int, int, int, str, int]] = None
            genome_len = len(record.seq) if getattr(record, "seq", None) is not None else None
            for item in record.features:
                ftype = str(item.type).lower() if item.type else "misc_feature"
                if ftype == "source":
                    continue
                feature_aliases = self._collect_feature_aliases(item)
                matched_alias: Optional[str] = None
                for wanted in preferred_aliases:
                    if wanted in feature_aliases:
                        matched_alias = wanted
                        break
                if matched_alias is None:
                    continue
                try:
                    feature_start = int(item.location.start)
                    feature_end = int(item.location.end)
                except Exception:
                    continue
                if feature_end <= feature_start:
                    continue

                start_1based = feature_start + 1
                end_1based = feature_end
                strand = item.location.strand
                if strand not in (-1, 1):
                    strand = 1
                length = end_1based - start_1based + 1
                priority = feature_priority.get(ftype, 5)
                score = (priority, -length)
                if selected is None or score < selected[:2]:
                    selected = (score[0], score[1], start_1based, end_1based, matched_alias, strand)

            if selected is not None:
                _, _, start_1based, end_1based, matched_alias, strand = selected
                accession_norm = _normalize_accession(accession)
                display_alias = alias_name_map.get(matched_alias or "", matched_alias or "unknown")
                return {
                    "query_name": display_alias,
                    "ensembl_gene_id": display_alias,
                    "ncbi_accession": accession_norm or self._coerce_str(accession),
                    "species": species_name or "Escherichia coli DH5alpha",
                    "assembly_name": f"NCBI {species_name or 'Escherichia coli DH5alpha'}",
                    "seq_region_name": accession_norm or self._coerce_str(accession),
                    "gene_start_1based": start_1based,
                    "gene_end_1based": end_1based,
                    "strand": int(strand) if strand in (-1, 1) else 1,
                    "display_name": display_alias,
                    "taxid": DH5A_TAXID,
                    "ncbi_genome_length": genome_len,
                }
        return None

    def _collect_feature_aliases(self, feature: Any) -> set[str]:
        aliases: set[str] = set()
        qualifiers = getattr(feature, "qualifiers", {}) or {}
        if not isinstance(qualifiers, dict):
            return aliases

        for key in (
            "gene",
            "gene_synonym",
            "locus_tag",
            "old_locus_tag",
            "protein_id",
            "standard_name",
        ):
            raw = qualifiers.get(key)
            if raw is None:
                continue
            values = raw if isinstance(raw, list) else [raw]
            for value in values:
                value_text = self._coerce_str(value)
                if value_text:
                    aliases.add(value_text.lower())

        note = qualifiers.get("note")
        if note is not None:
            values = note if isinstance(note, list) else [note]
            for value in values:
                extracted = _extract_gene_from_note(value)
                if extracted:
                    aliases.add(extracted.lower())

        return aliases

    def _fetch_ncbi_record_by_accession(self, accession: str) -> Any:
        response = self.api.get(
            NCBI_EFETCH,
            headers={"Accept": "text/plain"},
            params={
                "db": "nuccore",
                "id": self._coerce_str(accession),
                "rettype": "gbwithparts",
                "retmode": "text",
            },
        )
        if not response.text:
            return None
        try:
            records = list(SeqIO.parse(StringIO(response.text), "genbank"))
        except Exception:
            return None
        if not records:
            return None
        return records[0]

    def _fetch_uniprot_entry(self, query: str) -> Optional[dict[str, Any]]:
        response = self.api.get(
            UNIPROT_ENTRY_URL.format(accession=query),
            headers={"Accept": "application/json"},
        )
        if not isinstance(response.json_obj, dict):
            return None
        return response.json_obj

    def _collect_ncbi_gene_aliases(self, entry: dict[str, Any]) -> list[str]:
        aliases: list[str] = []

        for gene in entry.get("genes", []) or []:
            if not isinstance(gene, dict):
                continue

            gene_name = self._coerce_str(gene.get("geneName", {}).get("value"))
            if gene_name and gene_name not in aliases:
                aliases.append(gene_name)

            for key in ("orderedLocusNames", "orfNames"):
                raw = gene.get(key)
                if not isinstance(raw, list):
                    raw = [raw]
                for item in raw:
                    if isinstance(item, dict):
                        value = self._coerce_str(item.get("value"))
                    else:
                        value = self._coerce_str(item)
                    if value and value not in aliases:
                        aliases.append(value)

            for key in ("synonyms", "name"):
                raw = gene.get(key)
                if raw is None:
                    continue
                if isinstance(raw, list):
                    iterable = raw
                else:
                    iterable = [raw]
                for item in iterable:
                    if isinstance(item, dict):
                        value = self._coerce_str(item.get("value"))
                    else:
                        value = self._coerce_str(item)
                    if value and value not in aliases:
                        aliases.append(value)

        for ref in entry.get("uniProtKBCrossReferences", []) or []:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            if "ensemblbacteria" not in db_name:
                continue
            properties = self._coerce_properties(ref)
            for key in ("geneid", "genesymbol", "proteinid"):
                value = self._coerce_str(properties.get(key))
                if value and value not in aliases:
                    aliases.append(value)

        return aliases

    def _collect_ncbi_nucleotide_accessions(self, entry: dict[str, Any]) -> list[str]:
        accessions: list[str] = []
        for ref in entry.get("uniProtKBCrossReferences", []) or []:
            if not isinstance(ref, dict):
                continue
            db_name = str(ref.get("database") or ref.get("dbDisplayName") or "").lower()
            props = self._coerce_properties(ref)
            if "refseq" in db_name and props.get("nucleotidesequenceid"):
                accessions.append(props["nucleotidesequenceid"])
            elif db_name == "embl":
                molecule_type = (props.get("moleculetype") or "").lower()
                accession = self._coerce_str(ref.get("id"))
                if accession and accession not in accessions and "genomic_dna" in molecule_type.replace(" ", "_"):
                    accessions.append(accession)
        return accessions

    def _collect_ncbi_gene_summaries(
        self,
        aliases: list[str],
        taxid: Optional[int],
        organism_name: str,
    ) -> list[dict[str, Any]]:
        gene_ids: list[str] = []
        for alias in aliases:
            for term in self._build_ncbi_gene_terms(alias, taxid, organism_name):
                try:
                    ids = self._ncbi_esearch_gene_ids(term)
                except ToolError:
                    continue
                for gene_id in ids:
                    if gene_id not in gene_ids:
                        gene_ids.append(gene_id)
                if ids:
                    break
            if gene_ids:
                break

        if not gene_ids:
            return []

        return self._ncbi_gene_summaries(gene_ids)

    def _build_ncbi_gene_terms(self, alias: str, taxid: Optional[int], organism_name: str) -> list[str]:
        token = self._coerce_str(alias)
        if not token:
            return []
        terms: list[str] = []
        if taxid is not None:
            terms.append(f"{token}[Gene Name] AND {taxid}[Taxonomy ID]")
            terms.append(f"{token}[All Fields] AND {taxid}[Taxonomy ID]")
        if organism_name:
            genus_species = self._extract_genus_species(organism_name)
            if genus_species:
                terms.append(f'{token}[Gene Name] AND "{genus_species}"[Organism]')
        terms.append(f"{token}[Gene Name]")
        return terms

    def _ncbi_esearch_gene_ids(self, term: str) -> list[str]:
        response = self.api.get(
            NCBI_ESEARCH,
            headers={"Accept": "application/json"},
            params={"db": "gene", "term": term, "retmode": "json", "retmax": 50},
        )
        if not isinstance(response.json_obj, dict):
            return []
        result = response.json_obj.get("esearchresult")
        if not isinstance(result, dict):
            return []
        raw_ids = result.get("idlist", [])
        if not isinstance(raw_ids, list):
            return []
        ids: list[str] = []
        for raw_id in raw_ids:
            id_text = self._coerce_str(raw_id)
            if id_text and id_text not in ids:
                ids.append(id_text)
        return ids

    def _ncbi_gene_summaries(self, gene_ids: list[str]) -> list[dict[str, Any]]:
        response = self.api.get(
            NCBI_ESUMMARY,
            headers={"Accept": "application/json"},
            params={"db": "gene", "id": ",".join(gene_ids), "retmode": "json"},
        )
        if not isinstance(response.json_obj, dict):
            return []
        result = response.json_obj.get("result")
        if not isinstance(result, dict):
            return []
        uid_list = result.get("uids")
        if not isinstance(uid_list, list):
            return []
        summaries: list[dict[str, Any]] = []
        for uid in uid_list:
            record = result.get(uid)
            if isinstance(record, dict):
                record.setdefault("uid", uid)
                summaries.append(record)
        return summaries

    def _choose_best_ncbi_gene_summary(
        self,
        summaries: list[dict[str, Any]],
        aliases: list[str],
        accessions: list[str],
        prefer_accession_match: bool = False,
    ) -> Optional[dict[str, Any]]:
        target_accessions = {_normalize_accession(acc) for acc in accessions if self._coerce_str(acc)}
        alias_set = {self._coerce_str(alias).lower() for alias in aliases if self._coerce_str(alias)}
        normalized_accessions = {acc for acc in target_accessions if acc}

        for summary in summaries:
            genomic_info = self._coerce_ncbi_gene_genomic_info(
                summary,
                preferred_accessions=normalized_accessions,
                require_preferred=bool(prefer_accession_match and normalized_accessions),
            )
            if not genomic_info:
                continue
            if not normalized_accessions or _normalize_accession(genomic_info.get("chraccver")) in normalized_accessions:
                return summary

        if prefer_accession_match and normalized_accessions:
            return None

        for summary in summaries:
            genomic_info = self._coerce_ncbi_gene_genomic_info(summary)
            if not genomic_info:
                continue
            summary_aliases = self._coerce_ncbi_summary_aliases(summary)
            if summary_aliases.intersection(alias_set):
                return summary

        for summary in summaries:
            if self._coerce_ncbi_gene_genomic_info(summary):
                return summary
        return None

    def _coerce_ncbi_summary_aliases(self, summary: dict[str, Any]) -> set[str]:
        aliases: set[str] = set()
        for key in ("name", "nomenclaturesymbol", "nomenclaturename", "otherdesignations"):
            value = self._coerce_str(summary.get(key))
            if value:
                aliases.add(value.lower())
        other_aliases = self._coerce_str(summary.get("otheraliases"))
        if other_aliases:
            for item in other_aliases.split(","):
                alias = self._coerce_str(item)
                if alias:
                    aliases.add(alias.lower())
        return aliases

    def _coerce_ncbi_gene_genomic_info(
        self,
        summary: dict[str, Any],
        preferred_accessions: Optional[set[str]] = None,
        require_preferred: bool = False,
    ) -> Optional[dict[str, Any]]:
        genomicinfo = summary.get("genomicinfo")
        if not isinstance(genomicinfo, list):
            return None
        preferred = {acc for acc in (preferred_accessions or set()) if acc}
        if preferred:
            for item in genomicinfo:
                if not isinstance(item, dict):
                    continue
                chr_start = _to_int(item.get("chrstart"))
                chr_stop = _to_int(item.get("chrstop"))
                chr_accver = self._coerce_str(item.get("chraccver"))
                if chr_start is None or chr_stop is None or not chr_accver:
                    continue
                if _normalize_accession(chr_accver) in preferred:
                    return item
            if require_preferred:
                return None
        for item in genomicinfo:
            if not isinstance(item, dict):
                continue
            chr_start = _to_int(item.get("chrstart"))
            chr_stop = _to_int(item.get("chrstop"))
            chr_accver = self._coerce_str(item.get("chraccver"))
            if chr_start is not None and chr_stop is not None and chr_accver:
                return item
        return None

    def _ncbi_summary_to_coordinates(
        self,
        summary: dict[str, Any],
        aliases: list[str],
        organism_name: str,
        accessions: Optional[list[str]] = None,
    ) -> Optional[dict[str, Any]]:
        preferred_accessions = {_normalize_accession(acc) for acc in (accessions or []) if acc}
        genomic_info = self._coerce_ncbi_gene_genomic_info(
            summary,
            preferred_accessions=preferred_accessions,
            require_preferred=bool(preferred_accessions),
        )
        if not genomic_info:
            return None

        raw_start = _to_int(genomic_info.get("chrstart"))
        raw_end = _to_int(genomic_info.get("chrstop"))
        if raw_start is None or raw_end is None:
            return None

        start = min(raw_start, raw_end)
        end = max(raw_start, raw_end)
        strand = 1 if raw_start <= raw_end else -1
        organism = summary.get("organism") or {}
        ncbi_accession = self._coerce_str(genomic_info.get("chraccver"))
        species_name = self._coerce_str(organism.get("scientificname")) or organism_name or "unknown"
        genome_len = _to_int(genomic_info.get("chrlength")) or _to_int(genomic_info.get("chrlen"))

        return {
            "query_name": self._coerce_str(aliases[0]) if aliases else None,
            "ensembl_gene_id": self._coerce_str(summary.get("name")) or (self._coerce_str(aliases[0]) if aliases else "ncbi_gene"),
            "ncbi_accession": ncbi_accession or "unknown",
            "species": species_name,
            "assembly_name": f"NCBI {species_name}",
            "seq_region_name": ncbi_accession or "unknown",
            "gene_start_1based": start,
            "gene_end_1based": end,
            "strand": strand,
            "display_name": self._coerce_str(summary.get("name")) or self._coerce_str(summary.get("nomenclaturesymbol")),
            "taxid": _to_int(organism.get("taxid")),
            "ncbi_genome_length": genome_len,
        }

    def _coerce_str(self, value: Any) -> str:
        if value is None:
            return ""
        if isinstance(value, str):
            return value.strip()
        return str(value).strip()

    def _coerce_properties(self, entry_ref: dict[str, Any]) -> dict[str, str]:
        properties: dict[str, str] = {}
        for prop in entry_ref.get("properties", []) or []:
            if not isinstance(prop, dict):
                continue
            key = self._coerce_str(prop.get("key")).lower().replace(" ", "")
            value = self._coerce_str(prop.get("value"))
            if key and value:
                properties[key] = value
        return properties

    def _extract_genus_species(self, organism_name: str) -> str:
        base = organism_name.split("(")[0].strip()
        parts = base.split()
        if len(parts) >= 2:
            return f"{parts[0]} {parts[1]}"
        return ""

    def _lookup_ensembl_gene(self, ensembl_gene_id: str, _depth: int = 0) -> Optional[dict[str, Any]]:
        ensembl_gene_id = _normalize_ensembl_gene_id(ensembl_gene_id) or ensembl_gene_id
        resp = self.api.get(
            ENSEMBL_LOOKUP.format(ensembl_id=ensembl_gene_id),
            headers={"Accept": "application/json"},
            params={"expand": 0},
        )
        if not isinstance(resp.json_obj, dict):
            return None
        data = resp.json_obj
        if _depth > 3:
            return None
        object_type = str(data.get("object_type") or "").lower()
        if object_type and object_type != "gene":
            parent = data.get("Parent") or data.get("parent")
            if isinstance(parent, list):
                for item in parent:
                    if isinstance(item, str):
                        parent_id = _normalize_ensembl_gene_id(item)
                        if parent_id and parent_id != ensembl_gene_id:
                            return self._lookup_ensembl_gene(parent_id, _depth + 1)
                    if isinstance(item, dict):
                        parent_id = _normalize_ensembl_gene_id(item.get("id"))
                        if parent_id and parent_id != ensembl_gene_id:
                            return self._lookup_ensembl_gene(parent_id, _depth + 1)
            elif isinstance(parent, str):
                parent_id = _normalize_ensembl_gene_id(parent)
                if parent_id and parent_id != ensembl_gene_id:
                    return self._lookup_ensembl_gene(parent_id, _depth + 1)
            elif isinstance(parent, dict):
                parent_id = _normalize_ensembl_gene_id(parent.get("id"))
                if parent_id and parent_id != ensembl_gene_id:
                    return self._lookup_ensembl_gene(parent_id, _depth + 1)

        species = data.get("species")
        if isinstance(species, dict):
            species = species.get("name") or species.get("display_name") or species.get("scientific_name")
        return {
            "ensembl_gene_id": ensembl_gene_id,
            "species": str(species or "unknown"),
            "assembly_name": str(data.get("assembly_name", "unknown")),
            "seq_region_name": str(data.get("seq_region_name", data.get("seq_region", ""))),
            "gene_start_1based": _to_int(data.get("start")),
            "gene_end_1based": _to_int(data.get("end")),
            "strand": _to_int(data.get("strand"), default=1),
            "display_name": data.get("display_name") or data.get("external_name"),
            "taxid": _to_int(data.get("taxonomy_id") or data.get("taxid")),
        }

    def _extract_gene_from_mapping(self, mapping_payload: Any) -> Optional[str]:
        if not isinstance(mapping_payload, dict):
            return None
        entries = mapping_payload.get("results") or []
        if isinstance(mapping_payload.get("to"), str) and _normalize_ensembl_gene_id(mapping_payload.get("to")):
            return _normalize_ensembl_gene_id(mapping_payload.get("to"))
        if not isinstance(entries, list):
            entries = []
        for item in entries:
            to = self._normalize_mapping_value(item.get("to"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("toPrimaryAccession"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("to_id"))
            if to:
                return to
            to = self._normalize_mapping_value(item.get("toSecondary"))
            if to:
                return to
        return None

    def _normalize_mapping_value(self, raw: Any) -> Optional[str]:
        return _normalize_ensembl_gene_id(raw)


def _to_int(value: Any, default: Optional[int] = None) -> Optional[int]:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _poll_uniprot_status(api: ApiClient, job_id: str, max_attempts: int = 20, delay_seconds: float = 1.0) -> bool:
    import time

    for _ in range(max_attempts):
        status = api.get(UNIPROT_IDMAP_STATUS.format(job_id=job_id), headers={"Accept": "application/json"})
        data = status.json_obj or {}
        if isinstance(data, dict):
            if "results" in data:
                results = data.get("results")
                if results:
                    return True
                failed_ids = data.get("failedIds")
                if failed_ids:
                    return False
            if data.get("jobStatus") in {"FINISHED", "COMPLETED", "FINISHED_WITH_WARNINGS"}:
                return bool(data.get("results"))
            if data.get("jobStatus") == "ERROR":
                return False
        time.sleep(delay_seconds)
    return False
