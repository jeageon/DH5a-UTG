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
release에서 zip 파일을 C:\Program Files 에 압축풀기

```
## WebUI 실행

```bash
C:\Program Files\DH5a-UTG-1.0.9\DH5a-UTG-1.0.9 폴더내에
run-DH5aUTG-WebUI.bat 실행

*필요시 해당 파일로 바로가기 만들기
```

