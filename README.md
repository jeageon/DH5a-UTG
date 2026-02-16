# DH5a-UTG

`DH5a-UTG`는 미생물(기본: *E. coli* DH5α, NCBI nuccore `CP076470`)의
특정 유전자를 중심으로 ±10kb 구간의 gDNA를 추출하고,
동일 구간에 존재하는 간섭 요소를 `negative feature`로 표시해 PCR 삽입 설계를 돕는 도구입니다.
프라이머/동형상동 arm 설계를 위한 보조로, hairpin 유발 구조나 반복서열도 함께 표시합니다.

기본 전략은 NCBI 우선 해석입니다.
- 입력이 UniProtID(Pxxxx)인 경우: 해당 단백질/유전자와 연결된 NCBI 유전자 좌표를 우선 탐색
- 입력이 유전자명인 경우: NCBI Gene 기반으로 좌표 조회

## 설치

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## 실행

```bash
python -m src.main lacZ --query-type gene_name
python -m src.main P04637 --query-type uniprot_id
```

## WebUI 실행

```bash
python3 -m pip install -r requirements.txt
streamlit run src/webui.py
```

- 기본값은 DH5a(E. coli K-12 계열, taxid 511145 / CP076470) 기반입니다.
- 입력창에 유전자명 또는 UniProt ID를 입력하면 간섭 feature 분석 결과를 즉시 확인할 수 있습니다.

기본값:
- `--flank`: 유전자 좌우 확장 bp(기본 10000)
- `--taxid`: `511145` (E. coli DH5α)
- `--ncbi-accession`: `CP076470`
- `--query-type`: `auto|gene_name|uniprot_id`
- `--flank-mode`: `genomic` 또는 `strand_relative`
- `--features`: `annotation,low_complexity,palindrome,inverted_repeat,tandem_repeat,repeat,simple,variation,structural_variation,extreme_gc,homopolymer,ambiguous`
- `--mask`: `none|soft|hard` (Ensembl용)
- `--maf-threshold`: 변이 MAF 임계값
- `--gc-window`, `--gc-step`, `--gc-min`, `--gc-max`
- `--homopolymer-at`, `--homopolymer-gc`
- `--tandem-repeat-min-motif`, `--tandem-repeat-max-motif`, `--tandem-repeat-min-copies`
- `--low-complexity-window`, `--low-complexity-step`, `--low-complexity-max-entropy`
- `--offline`: 캐시만 사용

출력 파일명:
`{Query}.{assembly}.{chr}_{extStart}_{extEnd}.negfeatures.gb`

동일 basename의 `...metadata.json`도 생성됩니다.

실행 중 유전자 단위 간섭(예: 기존 주석/기능 요소)까지 포함하려면 `annotation`을 기본값에 추가해 두었으며,
필요 시 `--features annotation,repeat,...` 형태로 원하는 항목만 지정할 수 있습니다.
