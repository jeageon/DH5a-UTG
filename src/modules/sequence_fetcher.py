from __future__ import annotations

from ..config import (
    ENSEMBL_SEQUENCE_REGION,
    NCBI_EFETCH,
    ENSEMBL_SEQUENCE_MAX_BP,
    ENSEMBL_SEQUENCE_SAFETY_BP,
)
from ..models.data_schemas import GenomicCoordinates
from ..utils.api_client import ApiClient
from ..utils.coord_utils import clamp_region_length, build_region_string
from ..utils.exceptions import SequenceLengthMismatchError


class SequenceFetcher:
    def __init__(self, api_client: ApiClient) -> None:
        self.api = api_client

    def fetch(
        self,
        coordinates: GenomicCoordinates,
        mask: str = "soft",
        strict_region: bool = False,
    ) -> tuple[str, list[str]]:
        del strict_region
        warnings: list[str] = []
        use_ensembl = coordinates.coordinate_source == "ensembl"
        if not use_ensembl and mask != "none":
            warnings.append("NCBI sequence source does not support masking; using mask=none")

        region_start = coordinates.ext_start_1based
        region_end = coordinates.ext_end_1based
        if not use_ensembl and coordinates.ncbi_genome_length:
            max_end = max(region_start, coordinates.ncbi_genome_length)
            if region_end > max_end:
                warnings.append(
                    f"clamped end coordinate to genome length ({coordinates.ncbi_genome_length}) for NCBI fetch"
                )
                region_end = max_end

        region_len = region_end - region_start + 1
        expected = region_len
        if region_len > ENSEMBL_SEQUENCE_MAX_BP:
            region_start, region_end, clamped = clamp_region_length(
                region_start,
                region_end,
                ENSEMBL_SEQUENCE_SAFETY_BP,
            )
            if clamped:
                warnings.append(
                    "sequence span exceeded API limit; flank automatically reduced for safety"
                )
                expected = region_end - region_start + 1
                region_len = expected

        if use_ensembl:
            seq, used = self._fetch_region(
                species=coordinates.species,
                chr_name=coordinates.seq_region_name,
                start_1based=region_start,
                end_1based=region_end,
                strand=coordinates.strand,
                mask=mask,
            )
        else:
            seq, used = self._fetch_region_ncbi(
                accession=coordinates.seq_region_name,
                start_1based=region_start,
                end_1based=region_end,
                strand=coordinates.strand,
                mask=mask,
            )
            if not used and len(seq) == expected:
                used = ["NCBI sequence fetched from nuccore"]
        provider = "Ensembl" if use_ensembl else "NCBI"
        if len(seq) != expected:
            if use_ensembl and mask != "none":
                seq_fallback, fallback_warnings = self._fetch_region(
                    species=coordinates.species,
                    chr_name=coordinates.seq_region_name,
                    start_1based=region_start,
                    end_1based=region_end,
                    strand=coordinates.strand,
                    mask="none",
                )
                used.extend(fallback_warnings)
                if len(seq_fallback) == expected:
                    seq = seq_fallback
                else:
                    raise SequenceLengthMismatchError(
                        f"Expected {expected} nt but received {len(seq)} nt with mask={mask} "
                        f"and {len(seq_fallback)} nt with mask=none from Ensembl for "
                        f"region={build_region_string(coordinates.seq_region_name, region_start, region_end, coordinates.strand)}"
                    )
            else:
                raise SequenceLengthMismatchError(
                    f"Expected {expected} nt but received {len(seq)} nt from {provider} for "
                    f"region={build_region_string(coordinates.seq_region_name, region_start, region_end, coordinates.strand)}"
                )
        if len(seq) != expected:
            raise SequenceLengthMismatchError(
                f"Expected {expected} nt but received {len(seq)} nt from {provider}"
            )
        if coordinates.ext_start_1based != region_start or coordinates.ext_end_1based != region_end:
            used.append("sequence region was internally clamped by max request size")
        return seq.upper(), used

    def _fetch_region(
        self,
        species: str,
        chr_name: str,
        start_1based: int,
        end_1based: int,
        strand: int,
        mask: str,
    ) -> tuple[str, list[str]]:
        region = build_region_string(chr_name, start_1based, end_1based, strand)
        params = {}
        if mask and mask != "none":
            params["mask"] = mask
        resp = self.api.get(
            ENSEMBL_SEQUENCE_REGION.format(species=species, region=region),
            headers={"Accept": "text/plain"},
            params=params,
            disable_cache=True,
        )
        seq, warnings = self._normalize_fasta_or_plain(resp.text or "")
        return seq, warnings

    def _fetch_region_ncbi(
        self,
        accession: str,
        start_1based: int,
        end_1based: int,
        strand: int,
        mask: str,
    ) -> tuple[str, list[str]]:
        del mask
        params = {
            "db": "nuccore",
            "id": accession,
            "rettype": "fasta",
            "retmode": "text",
            "seq_start": start_1based,
            "seq_stop": end_1based,
        }
        if strand == -1:
            params["strand"] = 2
        resp = self.api.get(
            NCBI_EFETCH,
            headers={"Accept": "text/plain"},
            params=params,
            disable_cache=True,
        )
        seq, warnings = self._normalize_fasta_or_plain(resp.text or "")
        if not seq:
            return "", ["empty FASTA from NCBI efetch"]
        return seq, warnings + [f"fetched sequence from NCBI nuccore accession {accession}"]

    def _normalize_fasta_or_plain(self, text: str) -> tuple[str, list[str]]:
        lines = [line.strip() for line in text.splitlines() if line.strip()]
        if not lines:
            return "", []
        if lines[0].startswith(">"):
            body = "".join(lines[1:])
            return body, ["input treated as FASTA; header removed"]
        return "".join(lines), []
