from __future__ import annotations

from datetime import date, datetime
from pathlib import Path
from typing import Optional

import json
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from ..config import GENBANK_FEATURE_MAP, OUTPUT_FILE_SUFFIX
from ..models.data_schemas import SequenceRecordBundle


def _flatten_qualifier_value(value):
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        return [str(v) for v in value]
    return [str(value)]


def _as_location(start: int, end: int, strand: Optional[int]):
    if start < 0:
        start = 0
    if end < start:
        end = start
    if strand == 0:
        strand = None
    return FeatureLocation(start, end, strand=strand)


def _feature_qualifiers(feature):
    qualifiers: dict[str, list[str]] = {}
    qualifiers["note"] = _flatten_qualifier_value(feature.description)
    gene_name = (
        feature.attributes.get("gene_name")
        or feature.attributes.get("gene")
        or feature.attributes.get("locus_tag")
        or feature.attributes.get("old_locus_tag")
        or feature.attributes.get("protein_id")
    )
    if isinstance(gene_name, str):
        labels = {item.lower() for item in ["cds", "gene", "mrna", "exon", "misc_feature", "misc", "trna", "rrna"]}
        if not gene_name or str(gene_name).strip().lower() in labels:
            qualifiers["label"] = [feature.feature_type]
        else:
            qualifiers["label"] = [str(gene_name)]
    else:
        qualifiers["label"] = [feature.feature_type]
    qualifiers["source"] = [feature.source]
    if feature.score is not None:
        qualifiers["score"] = [f"{feature.score}"]
    for key, value in feature.attributes.items():
        if key == "id":
            qualifiers.setdefault("db_xref", []).append(str(value))
        else:
            qualifiers[key] = _flatten_qualifier_value(value)
    return qualifiers


def _build_record(bundle: SequenceRecordBundle) -> SeqRecord:
    seq = bundle.full_sequence.upper()
    coords = bundle.coordinates
    source_label = "ENSEMBL" if coords.coordinate_source == "ensembl" else "NCBI"
    display = coords.query_gene or coords.uniprot_id
    ncbi_accession = coords.ncbi_accession or coords.seq_region_name
    if ncbi_accession:
        source_id_line = f"{source_label} reference {ncbi_accession}"
    else:
        source_id_line = source_label

    source_db_xrefs: list[str] = []
    if coords.uniprot_id:
        source_db_xrefs.append(f"UniProtKB:{coords.uniprot_id}")
    if coords.ncbi_accession:
        source_db_xrefs.append(f"NCBI_nuccore:{coords.ncbi_accession}")
    if coords.ensembl_gene_id:
        source_db_xrefs.append(f"Ensembl:{coords.ensembl_gene_id}")

    record = SeqRecord(
        Seq(seq),
        id=display,
        name=display,
        description=(
            f"DH5a-UTG target region for {display} "
            f"({coords.assembly_name} {coords.seq_region_name}:{coords.ext_start_1based}-{coords.ext_end_1based}, "
            f"strand={coords.strand}, source={source_id_line})"
        ),
    )
    record.annotations["molecule_type"] = "DNA"
    record.annotations["organism"] = coords.species
    record.annotations["taxonomy"] = [source_label]
    record.annotations["data_file_division"] = "UNC"
    record.annotations["date"] = date.today().strftime("%d-%b-%Y").upper()

    source_feature = SeqFeature(
        _as_location(0, len(seq), 1),
        type="source",
        qualifiers={
            "organism": [coords.species],
            "db_xref": source_db_xrefs,
            "note": [
                f"extracted with ±{coords.ext_end_1based - coords.ext_start_1based + 1} bp flank",
                f"original genomic: {coords.assembly_name} {coords.seq_region_name}:{coords.ext_start_1based}-{coords.ext_end_1based}",
            ],
        },
    )
    record.features.append(source_feature)

    gene_start = max(0, coords.gene_start_1based - coords.ext_start_1based)
    gene_end = max(gene_start, coords.gene_end_1based - coords.ext_start_1based + 1)
    gene_feature = SeqFeature(
        _as_location(gene_start, gene_end, coords.strand),
        type="gene",
        qualifiers={
            "gene": [coords.display_name or coords.query_gene or coords.ensembl_gene_id or coords.uniprot_id],
            "db_xref": source_db_xrefs,
            "note": [f"target span from {source_label}"],
        },
    )
    record.features.append(gene_feature)

    for feature in sorted(
        bundle.features,
        key=lambda item: (item.start, -(item.end - item.start), item.feature_type),
    ):
        feature_type = feature.feature_type
        seq_feature_type = GENBANK_FEATURE_MAP.get(feature_type, "misc_feature")
        qualifiers = _feature_qualifiers(feature)
        rec_feature = SeqFeature(
            _as_location(feature.start, feature.end, feature.strand),
            type=seq_feature_type,
            qualifiers=qualifiers,
        )
        record.features.append(rec_feature)

    return record


def _feature_counts(features):
    counts: dict[str, int] = {}
    for feature in features:
        counts[feature.feature_type] = counts.get(feature.feature_type, 0) + 1
    return counts


def output_paths(
    outdir: Path,
    record_id: str,
    assembly: str,
    chr_name: str,
    ext_start: int,
    ext_end: int,
) -> tuple[Path, Path]:
    safe_chr = "".join(ch if ch.isalnum() else "_" for ch in chr_name)
    safe_asm = "".join(ch if ch.isalnum() else "_" for ch in assembly)
    safe_id = "".join(ch if ch.isalnum() or ch in "._-" else "_" for ch in record_id)
    base = f"{safe_id}.{safe_asm}.{safe_chr}_{ext_start}_{ext_end}"
    gb_path = outdir / f"{base}{OUTPUT_FILE_SUFFIX}"
    json_path = outdir / f"{base}.metadata.json"
    return gb_path, json_path


def write_outputs(
    bundle: SequenceRecordBundle,
    outdir: Path,
    write_metadata_json: bool = True,
) -> tuple[Path, Optional[Path]]:
    outdir.mkdir(parents=True, exist_ok=True)
    record_id = bundle.coordinates.query_gene or bundle.coordinates.uniprot_id
    gb_path, json_path = output_paths(
        outdir,
        record_id,
        bundle.coordinates.assembly_name,
        bundle.coordinates.seq_region_name,
        bundle.coordinates.ext_start_1based,
        bundle.coordinates.ext_end_1based,
    )
    record = _build_record(bundle)
    SeqIO.write(record, gb_path, "genbank")

    if write_metadata_json:
        metadata = dict(bundle.metadata)
        metadata.setdefault("run_timestamp", datetime.utcnow().isoformat() + "Z")
        metadata.setdefault("feature_counts", _feature_counts(bundle.features))
        json_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    else:
        json_path = None
    return gb_path, json_path
