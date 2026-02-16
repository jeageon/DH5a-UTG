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

### Windows에서 바로 실행(클릭형)

Windows 사용자라면 다음 배치 파일을 더블클릭해서 바로 실행할 수 있습니다.

- `run-UTG-CLI.bat` : 기존 `utg.py` 기반 CLI 실행
- `run-UTG-WebUI.bat` : WebUI 실행(UTG 모드)
- `run-DH5aUTG-CLI.bat` : `src/main.py` 기반 CLI 실행
- `run-DH5aUTG-WebUI.bat` : `src/webui.py` 기반 WebUI 실행

바탕화면 바로가기도 PowerShell로 한 번에 만들 수 있습니다.

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass -Force
powershell -ExecutionPolicy Bypass -File .\create-windows-shortcuts.ps1
```

바로가기는 기본적으로 바탕화면에 다음을 생성합니다.

- `UTG-CLI`
- `UTG-WebUI`
- `DH5aUTG-CLI`
- `DH5aUTG-WebUI`

- 기본값은 DH5a(E. coli K-12 계열, taxid 511145 / CP076470) 기반입니다.
- 입력창에 유전자명 또는 UniProt ID를 입력하면 간섭 feature 분석 결과를 즉시 확인할 수 있습니다.

기본값:
- `--flank`: 유전자 좌우 확장 bp(기본 10000)
- `--taxid`: `511145` (E. coli DH5α)
- `--ncbi-accession`: `CP076470` (원하면 `--ncbi-accession NC_123456` 등으로 엑세션을 바꿔 다른 염기서열 기준으로 즉시 전환)
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

### 엑세션 전환 예시

```bash
python -m src.main lacZ --query-type gene_name --ncbi-accession CP076470 --taxid 511145
python -m src.main lacZ --query-type gene_name --ncbi-accession NC_000913.3 --taxid 511145
```

WebUI에서는 `NCBI nuccore accession` 입력창에 원하는 GenBank accession을 입력한 뒤 실행하면,
해당 염색체/Plasmid 기준으로 분석 대상 유전체가 자동 전환됩니다.

실행 중 유전자 단위 간섭(예: 기존 주석/기능 요소)까지 포함하려면 `annotation`을 기본값에 추가해 두었으며,
필요 시 `--features annotation,repeat,...` 형태로 원하는 항목만 지정할 수 있습니다.

### Windows 설치형 EXE 생성

`PyInstaller`로 실행 파일(.exe)도 만들 수 있습니다.  
아래는 `WebUI` 단일 실행 파일을 만들기 위한 기본 명령입니다.

```powershell
powershell -ExecutionPolicy Bypass -File .\packaging\windows\build-windows-exe.ps1
```

기본 동작은 빌드가 끝나면 바탕화면에 바로가기 아이콘을 자동으로 생성합니다.

```powershell
powershell -ExecutionPolicy Bypass -File .\packaging\windows\build-windows-exe.ps1 -CreateDesktopShortcuts:$false
```

위처럼 지정하면 바탕화면 바로가지를 만들지 않습니다.

필요시 개발 의존성은 아래로 설치할 수 있습니다.

```powershell
pip install -r requirements-dev.txt
```

생성 결과:
- `dist\DH5a-UTG-WebUI.exe`

`CLI`도 같이 만들려면:

```powershell
powershell -ExecutionPolicy Bypass -File .\packaging\windows\build-windows-exe.ps1 -BuildCLI
```

생성 결과:
- `dist\DH5a-UTG-WebUI.exe`
- `dist\DH5a-UTG-CLI.exe`

설치형 exe가 필요하면 Inno Setup이 설치된 Windows에서 아래처럼 실행합니다.

```powershell
powershell -ExecutionPolicy Bypass -File .\packaging\windows\build-windows-exe.ps1 -BuildCLI -BuildInstaller
```

생성 결과:
- `packaging/windows/installer/DH5a-UTG-Setup.exe`

주의:
- 설치형 exe는 별도 의존성 패키징이 더 크고 빌드 시간이 오래 걸릴 수 있습니다.

### GitHub 릴리스에 EXE 업로드하기

이 저장소는 Windows 빌드 GitHub Actions를 제공합니다.  
태그(`v*`)를 푸시하면 자동으로 `win` 환경에서 exe를 빌드해 `Artifacts`에 업로드합니다.

```bash
git tag v1.0.0
git push origin v1.0.0
```

또는 GitHub Actions 화면에서 `Build and Release Windows EXE`를 수동 실행하고
`create_release`를 `true`로 지정하면 `Release`에 실행 파일을 등록합니다.

- `DH5a-UTG-WebUI.exe`
- `DH5a-UTG-CLI.exe`
- `DH5a-UTG-windows-executables.zip`
