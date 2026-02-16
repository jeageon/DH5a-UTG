from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import click

from .config import (
    CACHE_DIR,
    DH5A_ACCESSION,
    DH5A_TAXID,
    DEFAULT_FEATURES,
    DEFAULT_FLANK,
    DEFAULT_CACHE_TTL_HOURS,
    DEFAULT_RETRIES,
    DEFAULT_TIMEOUT,
    OUTPUT_DIR,
    FeatureScanOptions,
)
from .models.data_schemas import SequenceRecordBundle
from .modules.coordinate_resolver import CoordinateResolver
from .modules.feature_scanner import FeatureScanner
from .modules.output_generator import write_outputs
from .modules.sequence_fetcher import SequenceFetcher
from .utils.api_client import ApiClient
from .utils.exceptions import UTGError


def _parse_features(features_csv: str) -> list[str]:
    items = [item.strip() for item in features_csv.split(",") if item.strip()]
    return items if items else list(DEFAULT_FEATURES)


@click.command()
@click.argument("query")
@click.option("--outdir", default=str(OUTPUT_DIR), type=click.Path(file_okay=False, path_type=Path), help="output directory")
@click.option("--flank", default=DEFAULT_FLANK, type=int, help="flanking length in bp")
@click.option(
    "--flank-mode",
    type=click.Choice(["genomic", "strand_relative"], case_sensitive=False),
    default="genomic",
    help="flank expansion mode",
)
@click.option(
    "--mask",
    type=click.Choice(["none", "soft", "hard"]),
    default="soft",
    help="Masking mode (only used by Ensembl paths)",
)
@click.option(
    "--query-type",
    type=click.Choice(["auto", "gene_name", "uniprot_id"], case_sensitive=False),
    default="auto",
    help="auto detects query type; gene_name is recommended for DH5a",
)
@click.option("--taxid", default=DH5A_TAXID, type=int, help=f"NCBI taxid filter (default: {DH5A_TAXID})")
@click.option("--ncbi-accession", default=DH5A_ACCESSION, help=f"Preferred NCBI accession (default: {DH5A_ACCESSION})")
@click.option(
    "--features",
    default=",".join(DEFAULT_FEATURES),
    help="comma separated feature types",
)
@click.option("--maf-threshold", default=0.01, type=float)
@click.option("--gc-window", default=50, type=int)
@click.option("--gc-step", default=10, type=int)
@click.option("--gc-min", default=30.0, type=float)
@click.option("--gc-max", default=70.0, type=float)
@click.option("--homopolymer-at", default=5, type=int)
@click.option("--homopolymer-gc", default=4, type=int)
@click.option("--tandem-repeat-min-motif", default=2, type=int)
@click.option("--tandem-repeat-max-motif", default=6, type=int)
@click.option("--tandem-repeat-min-copies", default=3, type=int)
@click.option("--low-complexity-window", default=30, type=int)
@click.option("--low-complexity-step", default=10, type=int)
@click.option("--low-complexity-max-entropy", default=1.2, type=float)
@click.option("--timeout", default=DEFAULT_TIMEOUT, type=float)
@click.option("--retries", default=DEFAULT_RETRIES, type=int)
@click.option("--cache", type=click.Choice(["on", "off"]), default="on")
@click.option("--cache-ttl-hours", default=DEFAULT_CACHE_TTL_HOURS, type=int)
@click.option("--offline", is_flag=True, default=False)
@click.option("--debug", is_flag=True, default=False)
@click.option("--write-metadata-json", is_flag=True, default=True)
def cli(
    query: str,
    outdir: Path,
    flank: int,
    flank_mode: str,
    mask: str,
    query_type: str,
    taxid: int,
    ncbi_accession: str,
    features: str,
    maf_threshold: float,
    gc_window: int,
    gc_step: int,
    gc_min: float,
    gc_max: float,
    homopolymer_at: int,
    homopolymer_gc: int,
    tandem_repeat_min_motif: int,
    tandem_repeat_max_motif: int,
    tandem_repeat_min_copies: int,
    low_complexity_window: int,
    low_complexity_step: int,
    low_complexity_max_entropy: float,
    timeout: float,
    retries: int,
    cache: str,
    cache_ttl_hours: int,
    offline: bool,
    debug: bool,
    write_metadata_json: bool,
):
    del debug
    try:
        selected_features = _parse_features(features)
        feature_options = FeatureScanOptions(
            maf_threshold=maf_threshold,
            gc_window=gc_window,
            gc_step=gc_step,
            gc_min=gc_min,
            gc_max=gc_max,
            homopolymer_at=homopolymer_at,
            homopolymer_gc=homopolymer_gc,
            tandem_repeat_min_motif=tandem_repeat_min_motif,
            tandem_repeat_max_motif=tandem_repeat_max_motif,
            tandem_repeat_min_copies=tandem_repeat_min_copies,
            low_complexity_window=low_complexity_window,
            low_complexity_step=low_complexity_step,
            low_complexity_max_entropy=low_complexity_max_entropy,
        )

        cache_enabled = cache == "on"
        api = ApiClient(
            timeout=timeout,
            retries=retries,
            cache_enabled=cache_enabled,
            cache_path=str(CACHE_DIR),
            ttl_hours=cache_ttl_hours,
            offline=offline,
        )

        resolver = CoordinateResolver(api)
        resolver_result = resolver.resolve(
            query=query,
            flank_bp=flank,
            flank_mode=flank_mode,
            taxid_filter=taxid,
            query_type=query_type.lower(),
            ncbi_accession=ncbi_accession,
        )
        coordinates = resolver_result.coordinates

        fetcher = SequenceFetcher(api)
        sequence, fetch_warnings = fetcher.fetch(coordinates=coordinates, mask=mask)

        scanner = FeatureScanner(api)
        detected_features, scan_warnings = scanner.scan(
            coordinates=coordinates,
            full_sequence=sequence,
            requested_features=selected_features,
            options=feature_options,
        )

        metadata = {
            "query": query,
            "query_type": coordinates.query_type,
            "query_gene": coordinates.query_gene,
            "uniprot_id": coordinates.uniprot_id,
            "ensembl_gene_id": coordinates.ensembl_gene_id,
            "coordinate_source": coordinates.coordinate_source,
            "ncbi_accession": coordinates.ncbi_accession,
            "ncbi_genome_length": coordinates.ncbi_genome_length,
            "organism": coordinates.species,
            "assembly": coordinates.assembly_name,
            "region": f"{coordinates.seq_region_name}:{coordinates.ext_start_1based}-{coordinates.ext_end_1based}:{coordinates.strand}",
            "flank_bp": flank,
            "flank_mode": flank_mode,
            "mask": mask,
            "ncbi_accession_preference": ncbi_accession,
            "api_cache": cache_enabled,
            "options": asdict(feature_options),
            "warnings": [*resolver_result.warnings, *fetch_warnings, *scan_warnings],
        }

        bundle = SequenceRecordBundle(
            coordinates=coordinates,
            full_sequence=sequence,
            features=detected_features,
            metadata=metadata,
        )
        gb_path, metadata_path = write_outputs(
            bundle=bundle,
            outdir=outdir,
            write_metadata_json=write_metadata_json,
        )
        click.echo(f"GenBank: {gb_path}")
        if metadata_path:
            click.echo(f"Metadata: {metadata_path}")
    except UTGError as exc:
        raise click.ClickException(str(exc))
    except Exception as exc:
        raise click.ClickException(f"Unexpected error: {exc}") from exc


if __name__ == "__main__":
    cli()
