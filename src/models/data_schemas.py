from __future__ import annotations

from typing import Any
from typing import Optional

from pydantic import BaseModel, Field


class GenomicCoordinates(BaseModel):
    uniprot_id: str
    query_type: str = "uniprot"
    query_gene: Optional[str] = None
    ensembl_gene_id: str
    coordinate_source: str = "ensembl"
    ncbi_accession: Optional[str] = None
    ncbi_genome_length: Optional[int] = None
    ncbi_annotations: list[dict[str, Any]] = Field(default_factory=list)
    species: str
    assembly_name: str
    seq_region_name: str
    gene_start_1based: int
    gene_end_1based: int
    strand: int
    display_name: Optional[str] = None
    taxid: Optional[int] = None

    ext_start_1based: int
    ext_end_1based: int


class NegativeFeature(BaseModel):
    feature_type: str
    start: int
    end: int
    description: str
    source: str
    score: Optional[float] = None
    strand: Optional[int] = None
    attributes: dict[str, Any] = Field(default_factory=dict)


class SequenceRecordBundle(BaseModel):
    coordinates: GenomicCoordinates
    full_sequence: str
    features: list[NegativeFeature]
    metadata: dict[str, Any]
