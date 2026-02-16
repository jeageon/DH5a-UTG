from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Optional
import json

import streamlit as st

from src.config import (
    CACHE_DIR,
    DEFAULT_CACHE_TTL_HOURS,
    DEFAULT_FEATURES,
    DEFAULT_FLANK,
    DH5A_ACCESSION,
    DH5A_TAXID,
    DEFAULT_TIMEOUT,
    DEFAULT_RETRIES,
    FeatureScanOptions,
)
from src.models.data_schemas import SequenceRecordBundle
from src.modules.coordinate_resolver import CoordinateResolver
from src.modules.feature_scanner import FeatureScanner
from src.modules.output_generator import write_outputs
from src.modules.sequence_fetcher import SequenceFetcher
from src.utils.api_client import ApiClient
from src.utils.exceptions import NoMappingError, SequenceLengthMismatchError, ToolError, UTGError


st.set_page_config(page_title="DH5a-UTG", page_icon="🧬", layout="wide")


def _run_pipeline(
    query: str,
    flank: int,
    flank_mode: str,
    mask: str,
    query_type: str,
    taxid: int,
    ncbi_accession: str,
    features_csv: str,
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
    cache_ttl_hours: int,
    api_cache: bool,
    palindrome_min_len: int,
    palindrome_max_len: int,
    hairpin_min_arm: int,
    hairpin_max_arm: int,
    hairpin_max_spacer: int,
    offline: bool = False,
) -> tuple[Path, Optional[Path], SequenceRecordBundle, TemporaryDirectory]:
    selected_features = [item for item in features_csv.split(",") if item]
    if not selected_features:
        selected_features = list(DEFAULT_FEATURES)

    feature_options = FeatureScanOptions(
        maf_threshold=maf_threshold,
        gc_window=gc_window,
        gc_step=gc_step,
        gc_min=gc_min,
        gc_max=gc_max,
        homopolymer_at=homopolymer_at,
        homopolymer_gc=homopolymer_gc,
        palindrome_min_len=palindrome_min_len,
        palindrome_max_len=palindrome_max_len,
        hairpin_min_arm=hairpin_min_arm,
        hairpin_max_arm=hairpin_max_arm,
        hairpin_max_spacer=hairpin_max_spacer,
        tandem_repeat_min_motif=tandem_repeat_min_motif,
        tandem_repeat_max_motif=tandem_repeat_max_motif,
        tandem_repeat_min_copies=tandem_repeat_min_copies,
        low_complexity_window=low_complexity_window,
        low_complexity_step=low_complexity_step,
        low_complexity_max_entropy=low_complexity_max_entropy,
    )

    api = ApiClient(
        timeout=timeout,
        retries=retries,
        cache_enabled=api_cache,
        cache_path=str(CACHE_DIR),
        ttl_hours=cache_ttl_hours,
        offline=offline,
    )

    resolver = CoordinateResolver(api)
    resolver_result = resolver.resolve(
        query=query,
        flank_bp=flank,
        flank_mode=flank_mode,
        query_type=query_type,
        taxid_filter=taxid,
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
    warnings = [*resolver_result.warnings, *fetch_warnings, *scan_warnings]

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
        "region": (
            f"{coordinates.seq_region_name}:{coordinates.ext_start_1based}-"
            f"{coordinates.ext_end_1based}:{coordinates.strand}"
        ),
        "flank_bp": flank,
        "flank_mode": flank_mode,
        "mask": mask,
        "ncbi_accession_preference": ncbi_accession,
        "api_cache": api_cache,
        "options": asdict(feature_options),
        "warnings": warnings,
    }

    bundle = SequenceRecordBundle(
        coordinates=coordinates,
        full_sequence=sequence,
        features=detected_features,
        metadata=metadata,
    )

    tmp_dir = TemporaryDirectory()
    gb_path, json_path = write_outputs(
        bundle=bundle,
        outdir=Path(tmp_dir.name),
        write_metadata_json=True,
    )
    # keep files alive outside this function
    return gb_path, json_path, bundle, tmp_dir


def _to_feature_rows(features: SequenceRecordBundle) -> list[dict[str, object]]:
    items: list[dict[str, object]] = []
    generic_labels = {"cds", "gene", "exon", "intron", "mrna", "misc_feature", "misc"}

    def _choose_display_name(feature: NegativeFeature) -> str:
        for value in (
            feature.attributes.get("gene_name"),
            feature.attributes.get("gene"),
            feature.attributes.get("locus_tag"),
            feature.attributes.get("old_locus_tag"),
            feature.attributes.get("protein_id"),
            feature.attributes.get("product"),
            feature.attributes.get("annotation_type"),
        ):
            if value is None:
                continue
            text = str(value).strip()
            if not text or text.lower() in generic_labels:
                continue
            return text
        return feature.feature_type

    for idx, feature in enumerate(features.features, start=1):
        gene_name = _choose_display_name(feature)
        product = feature.attributes.get("product", "")
        annotation_type = feature.attributes.get("annotation_type", "")
        product_display = product or feature.description
        items.append(
            {
                "idx": idx,
                "type": feature.feature_type,
                "start": feature.start,
                "end": feature.end,
                "source": feature.source,
                "score": feature.score,
                "annotation_type": annotation_type,
                "gene_name": gene_name,
                "product": product_display,
                "description": feature.description,
            }
        )
    return items


def _maybe_float(value: float) -> float:
    if value is None:
        return 0.0
    return float(value)


def main() -> None:
    st.title("DH5a-UTG")
    st.caption("NCBI GenBank accession 기준으로 주변 gDNA 추출 + 간섭 feature 표시 도구")

    with st.expander("실행 설정", expanded=True):
        query = st.text_input("유전자명 또는 UniProt ID", value="", help="예: lacZ, xylR, P04637")
        query_type = st.selectbox(
            "입력 유형",
            options=["auto", "gene_name", "uniprot_id"],
            index=1,
        )
        flank = st.number_input("양쪽 flank (bp)", min_value=0, max_value=50000, value=DEFAULT_FLANK, step=500)
        flank_mode = st.selectbox("flank 계산 방식", options=["genomic", "strand_relative"], index=0)
        mask = st.selectbox("mask", options=["none", "soft", "hard"], index=0)
        taxid = st.number_input("NCBI taxid", min_value=1, value=DH5A_TAXID, step=1)
        ncbi_accession = st.text_input(
            "NCBI nuccore accession",
            value=DH5A_ACCESSION,
            help="예: CP076470, NC_000913.3. 값 변경 시 해당 GenBank 염색체/플라스미드 기준으로 분석 대상을 전환합니다.",
        )
        st.markdown("Feature")
        selected_features = st.multiselect(
            "표시할 feature",
            options=DEFAULT_FEATURES,
            default=DEFAULT_FEATURES,
        )
        feature_csv = ",".join(selected_features)

        st.markdown("스캔 옵션")
        maf_threshold = st.slider("variation MAF 임계값", min_value=0.0, max_value=1.0, value=0.01, step=0.01)
        gc_window = st.number_input("GC 슬라이딩 창", min_value=10, max_value=500, value=50)
        gc_step = st.number_input("GC step", min_value=1, max_value=100, value=10)
        gc_min = st.slider("GC 최소값", min_value=0.0, max_value=100.0, value=30.0, step=1.0)
        gc_max = st.slider("GC 최대값", min_value=0.0, max_value=100.0, value=70.0, step=1.0)
        homopolymer_at = st.number_input("Homopolymer 길이", min_value=3, max_value=50, value=5)
        homopolymer_gc = st.number_input("Homopolymer GC", min_value=2, max_value=20, value=4)
        st.markdown("Primer/ARM 방해 구조 스캔")
        palindrome_min_len = st.number_input("완전 palindrome 최소 길이", min_value=4, max_value=30, value=8, step=1)
        palindrome_max_len = st.number_input("완전 palindrome 최대 길이", min_value=4, max_value=30, value=14, step=1)
        hairpin_min_arm = st.number_input("inverted repeat arm 최소 길이", min_value=4, max_value=30, value=8, step=1)
        hairpin_max_arm = st.number_input("inverted repeat arm 최대 길이", min_value=4, max_value=30, value=12, step=1)
        hairpin_max_spacer = st.number_input("헤어핀 spacer 최대 길이", min_value=0, max_value=200, value=20, step=1)
        tandem_repeat_min_motif = st.number_input("Tandem repeat motif 최소 길이", min_value=1, max_value=12, value=2, step=1)
        tandem_repeat_max_motif = st.number_input("Tandem repeat motif 최대 길이", min_value=1, max_value=20, value=6, step=1)
        tandem_repeat_min_copies = st.number_input("Tandem repeat 최소 반복 횟수", min_value=2, max_value=20, value=3, step=1)
        low_complexity_window = st.number_input("저복잡도 창 길이", min_value=10, max_value=200, value=30, step=5)
        low_complexity_step = st.number_input("저복잡도 sliding step", min_value=1, max_value=100, value=10, step=1)
        low_complexity_max_entropy = st.number_input("저복잡도 최대 엔트로피", min_value=0.1, max_value=2.0, value=1.2, step=0.1)

    with st.expander("네트워크/캐시 설정", expanded=False):
        timeout = st.number_input("API timeout", min_value=1.0, max_value=120.0, value=DEFAULT_TIMEOUT, step=1.0)
        retries = st.number_input("재시도 횟수", min_value=0, max_value=20, value=DEFAULT_RETRIES, step=1)
        cache_ttl_hours = st.number_input("캐시 TTL(시간)", min_value=1, value=DEFAULT_CACHE_TTL_HOURS, step=1)
        api_cache = st.toggle("캐시 사용", value=True)
        offline = st.toggle("오프라인(캐시 전용)", value=False)

    run_clicked = st.button("Run", type="primary")

    if not run_clicked:
        return

    query_value = query.strip()
    if not query_value:
        st.warning("유전자명 또는 UniProt ID를 입력하세요.")
        return
    if taxid <= 0:
        st.warning("taxid는 1 이상의 정수여야 합니다.")
        return
    if not ncbi_accession.strip():
        st.warning("NCBI accession을 입력하세요.")
        return

    with st.spinner("처리 중..."):
        try:
            result = _run_pipeline(
                query=query_value,
                flank=int(flank),
                flank_mode=flank_mode,
                mask=mask,
                query_type=query_type,
                taxid=int(taxid),
                ncbi_accession=ncbi_accession.strip(),
                features_csv=feature_csv,
                maf_threshold=_maybe_float(maf_threshold),
                gc_window=int(gc_window),
                gc_step=int(gc_step),
                gc_min=_maybe_float(gc_min),
                gc_max=_maybe_float(gc_max),
                homopolymer_at=int(homopolymer_at),
                homopolymer_gc=int(homopolymer_gc),
                timeout=float(timeout),
                retries=int(retries),
                cache_ttl_hours=int(cache_ttl_hours),
                api_cache=api_cache,
                palindrome_min_len=int(palindrome_min_len),
                palindrome_max_len=int(palindrome_max_len),
                hairpin_min_arm=int(hairpin_min_arm),
                hairpin_max_arm=int(hairpin_max_arm),
                hairpin_max_spacer=int(hairpin_max_spacer),
                tandem_repeat_min_motif=int(tandem_repeat_min_motif),
                tandem_repeat_max_motif=int(tandem_repeat_max_motif),
                tandem_repeat_min_copies=int(tandem_repeat_min_copies),
                low_complexity_window=int(low_complexity_window),
                low_complexity_step=int(low_complexity_step),
                low_complexity_max_entropy=_maybe_float(low_complexity_max_entropy),
                offline=offline,
            )
            gb_path, metadata_path, bundle, tmp_dir = result
        except (NoMappingError, SequenceLengthMismatchError, ToolError, UTGError) as exc:
            st.error(f"{type(exc).__name__}: {exc}")
            return
        except Exception as exc:
            st.error(f"예상치 못한 오류: {exc}")
            return

        try:
            coord = bundle.coordinates
            st.success("분석 완료")
            st.subheader("요약")
            st.write(
                {
                    "Query": query,
                    "Query type": coord.query_type,
                    "Gene": coord.query_gene,
                    "좌표원천": coord.coordinate_source,
                    "종": coord.species,
                    "서열 범위": f"{coord.seq_region_name}:{coord.ext_start_1based}-{coord.ext_end_1based}({coord.strand})",
                    "총 길이": len(bundle.full_sequence),
                }
            )

            warnings = bundle.metadata.get("warnings") or []
            if warnings:
                st.subheader("Warnings")
                for item in warnings:
                    st.warning(item)

            st.subheader("특성 분포")
            feature_rows = _to_feature_rows(bundle)
            if feature_rows:
                st.dataframe(feature_rows, hide_index=True)
            else:
                st.info("negative feature가 없습니다.")

            metadata = bundle.metadata.copy()
            metadata["run_summary"] = {
                "sequence_length": len(bundle.full_sequence),
                "feature_count": len(bundle.features),
            }
            metadata_bytes = json.dumps(metadata, ensure_ascii=False, indent=2).encode("utf-8")
            gb_bytes = gb_path.read_bytes()

            st.subheader("결과 다운로드")
            st.download_button(
                "GenBank 다운로드",
                data=gb_bytes,
                file_name=gb_path.name,
                mime="application/octet-stream",
            )
            st.download_button(
                "Metadata JSON 다운로드",
                data=metadata_bytes,
                file_name=Path(metadata_path).name if metadata_path else "metadata.json",
                mime="application/json",
            )

        finally:
            tmp_dir.cleanup()


if __name__ == "__main__":
    main()
