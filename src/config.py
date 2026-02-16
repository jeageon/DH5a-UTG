from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


DEFAULT_FEATURES = [
    "annotation",
    "repeat",
    "simple",
    "variation",
    "structural_variation",
    "extreme_gc",
    "homopolymer",
    "ambiguous",
]

DEFAULT_TIMEOUT = 20.0
DEFAULT_RETRIES = 5
DEFAULT_CACHE_TTL_HOURS = 24
DEFAULT_FLANK = 10_000
DH5A_ACCESSION = "CP076470"
DH5A_TAXID = 511145
DH5A_NAME = "Escherichia coli DH5alpha"

EBI_COORDINATES_URL = "https://www.ebi.ac.uk/proteins/api/coordinates/{accession}"
UNIPROT_IDMAP_RUN = "https://rest.uniprot.org/idmapping/run"
UNIPROT_IDMAP_STATUS = "https://rest.uniprot.org/idmapping/status/{job_id}"
UNIPROT_IDMAP_RESULTS = "https://rest.uniprot.org/idmapping/results/{job_id}"
UNIPROT_ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
ENSEMBL_LOOKUP = "https://rest.ensembl.org/lookup/id/{ensembl_id}"
ENSEMBL_SEQUENCE_REGION = "https://rest.ensembl.org/sequence/region/{species}/{region}"
ENSEMBL_SEQUENCE_ID = "https://rest.ensembl.org/sequence/id/{ensembl_id}"
ENSEMBL_OVERLAP = "https://rest.ensembl.org/overlap/region/{species}/{region}"

OUTPUT_FILE_SUFFIX = ".negfeatures.gb"
OUTPUT_DIR = Path("data/output")
CACHE_DIR = Path("data/cache")
OUTPUT_FASTA_SUFFIX = ".target_region.fasta"

ENSEMBL_OVERLAP_MAX_BP = 5_000_000
ENSEMBL_OVERLAP_CHUNK_BP = 4_500_000
ENSEMBL_SEQUENCE_MAX_BP = 10_000_000
ENSEMBL_SEQUENCE_SAFETY_BP = 9_500_000

USER_AGENT = "UTG/1.0.0 (+https://github.com/)"


@dataclass(frozen=True)
class FeatureScanOptions:
    maf_threshold: float = 0.01
    gc_window: int = 50
    gc_step: int = 10
    gc_min: float = 30.0
    gc_max: float = 70.0
    homopolymer_at: int = 5
    homopolymer_gc: int = 4


GENBANK_FEATURE_MAP = {
    "annotation": "misc_feature",
    "gene": "gene",
    "cds": "CDS",
    "trna": "tRNA",
    "rrna": "rRNA",
    "nc_rna": "ncRNA",
    "rep_origin": "rep_origin",
    "repeat_region": "repeat_region",
    "repeat": "repeat_region",
    "simple": "repeat_region",
    "variation": "variation",
    "structural_variation": "misc_feature",
    "extreme_gc": "misc_feature",
    "homopolymer": "misc_feature",
    "ambiguous": "misc_feature",
}
