# Bulk RNA-seq 분석 리포트 Workflow

`<sub>`[EN](README.md) · **한국어** `</sub>`

**nf-core/rnaseq**의 salmon gene-count 출력을 입력으로 받아, **MultiQC**, **DESeq2**, **GSEA(MsigDB)**, **ORA(ClusterProfiler)**, **TFEA(CollecTRI)**, **PROGENy pathway scoring**, **L2S2 (LINCS L1000) cMap** 결과를 HTML 리포트로 생성하는 Snakemake 파이프라인

## Workflow
```
nf-core/rnaseq ─► salmon counts + multiqc_data
                        │
                        ▼
   ┌─ Snakemake ──────────────────────────────────────┐
   │                                                  │
   │  MultiQC summary            (QC)                 │
   │  DESeq2 VST → PCA · corr    (exploratory)        │
   │  DESeq2 + apeglm            (differential)       │
   │     ├─► fgsea (MSigDB)      (GSEA)               │
   │     ├─► clusterProfiler     (ORA)                │
   │     ├─► decoupleR + CollecTRI   (TF activity)    │
   │     └─► L2S2 / LINCS L1000      (drug reposition)│
   │  PROGENy                    (pathway activity)   │
   │                                                  │
   └──────────────────────┬───────────────────────────┘
                          │
                          ▼
                 Quarto report (HTML)
```
---

## Quick start

```bash
# 1. Snakemake 설치 (Python venv)
python3 -m venv .venv
.venv/bin/pip install snakemake pandas

# 2. dry-run 확인
.venv/bin/snakemake \
  --snakefile workflow/Snakefile \
  --configfile tests/test_data/config_test.yaml \
  --use-conda -n

# 3. 실행
.venv/bin/snakemake \
  --snakefile workflow/Snakefile \
  --configfile tests/test_data/config_test.yaml \
  --use-conda -c1
```

HTML Report는 `results/report/report.html`에 생성됨.

### 초기 설정

아래 명령어를 실행 후 샘플에 대한 정보 입력

```bash
.venv/bin/python workflow/scripts/init_project.py
```

 `config/config.yaml`, `config/samples.tsv`, `config/contrasts.tsv` 자동 생성

---

## Docker

Snakemake + conda/mamba가 포함된 독립 이미지. per-rule R/Python env는 첫 실행 시
`.snakemake/conda/`에 빌드되어 호스트에 캐시됩니다.

```bash
# 빌드 (linux/amd64; Apple Silicon에서는 자동으로 linux/arm64 빌드)
docker build -t bulk-rnaseq:latest .

# Fixture 실행
docker run --rm \
    -u "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -v "$PWD":/project \
    bulk-rnaseq:latest \
    --configfile tests/test_data/config_test.yaml -c1

# 실제 프로젝트 실행 (호스트에 config.yaml + counts + multiqc_data 준비)
docker run --rm \
    -u "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -v "$PWD":/project \
    bulk-rnaseq:latest \
    --configfile config/config.yaml -c1
```

`-u`와 `-e HOME`은 필수:

- `-u "$(id -u):$(id -g)"` — `.snakemake/`, `results/`, `logs/` 산출물이 호스트 사용자
  소유로 작성되도록 강제. 누락 시 컨테이너 내부 root 소유 파일이 호스트에 남음.
- `-e HOME=/tmp` — conda/mamba가 캐시를 쓸 writable HOME 필요. 이미지 기본
  `mambauser` HOME은 사용자 오버라이드 시 쓰기 권한 없음.

---

## 인터랙티브 탐색

`Dockerfile.jupyter`는 JupyterLab + R kernel(IRkernel) + 파이프라인과 동일한
R 스택을 얹은 컨테이너를 빌드합니다. 최종 리포트는 여전히 `bulk-rnaseq`
(Snakemake → Quarto)가 생성하고, 이 이미지는 플롯 수정·추가 그림·임시 질의 등
파이프라인 산출물에 대한 탐색 전용입니다.

```bash
docker build -f Dockerfile.jupyter -t bulk-rnaseq-jupyter:latest .

docker run --rm \
    -u "$(id -u):$(id -g)" \
    -e HOME=/tmp \
    -v "$PWD":/project \
    -p 8888:8888 \
    bulk-rnaseq-jupyter:latest
```

시작 시 출력되는 `http://127.0.0.1:8888/lab?token=…` URL을 열고
`notebooks/explore.ipynb`를 실행. 템플릿은 리포트 각 섹션 — QC, PCA, volcano,
enrichment, TF / pathway activity, drug repositioning — 을 그대로 미러링하며
대응되는 `results/**/*.{csv,tsv,rds}` 파일을 읽어 cutoff·라벨·플롯 코드를
Snakemake 재실행 없이 수정 가능.

---

## 입력 파일

| 파일                     | 용도                                                                                                                |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------- |
| `config/config.yaml`   | 전역 설정: 경로, DE·enrichment cutoff, TF/pathway 도구, drug repositioning 서비스.                                 |
| `config/samples.tsv`   | 컬럼:`sample, condition, replicate, batch`.                                                                       |
| `config/contrasts.tsv` | 컬럼:`contrast_id, factor, numerator, denominator, description`.                                                  |
| `<counts>.tsv`         | nf-core/rnaseq의 `salmon.merged.gene_counts_length_scaled.tsv` (1열=`gene_id` Ensembl, 2열=`gene_name` HGNC). |
| `multiqc_data/`        | nf-core/rnaseq MultiQC 출력 디렉토리.                                                                               |

---

## 리포트 섹션

| 단계                                 | 방법                                                       | 주 산출물                                                              |
| ------------------------------------ | ---------------------------------------------------------- | ---------------------------------------------------------------------- |
| **QC**                         | MultiQC로 FastQC / STAR / Salmon / RSeQC metric 집계.      | 샘플별 QC 표, library-size·mapping-rate plot.                         |
| **탐색 분석**                  | VST 변환 후 top-500 variable gene PCA, sample correlation. | PCA, scree, dendrogram, correlation heatmap.                           |
| **차등발현**                   | DESeq2 Wald test + apeglm LFC shrinkage.                   | DEG 요약, volcano, MA, top-30 DEG heatmap, 전체 결과.                 |
| **Gene-set enrichment (GSEA)** | pre-ranked GSEA (ranking matric: Wald stat).               | MSigDB H / C2:CP (Reactome, WikiPathways, PID, BioCarta) / C2:CGP / C6 |
| **Over-representation (ORA)**  | `clusterProfiler::enricher()` + KEGG live REST.          | DB(GOBP, KEGG, Reactome, Hallmark)별 top-10 up/down                    |
| **TF 활성도**                  | decoupler + ULM + CollecTRI                               | Top-30 TF + 전체 score                                                 |
| **Pathway 활성도**             | decoupler + PROGENy                                        | z-scored heatmap by sample + treated−control delta (Wilcoxon).        |
| **Drug repositioning**         | Up/down DEG signature로 L2S2 paired query.                 | Ranked perturbagen                                                     |
| **Audit trail**                | Config snapshot, MD5, session info                        | 재현성 정보                                                            |

---

## 저장소 구조

```
Bulk-RNAseq/
├── config/
│   ├── config.yaml              # 전역 설정
│   ├── samples.tsv              # 샘플 메타데이터
│   └── contrasts.tsv            # DE 비교 정의
├── workflow/
│   ├── Snakefile                # 파이프라인 진입점
│   ├── rules/
│   │   ├── qc.smk
│   │   ├── de.smk
│   │   ├── enrichment.smk
│   │   └── report.smk
│   ├── envs/                    # 단계별 conda env yaml
│   └── scripts/                 # R / Python 구현
├── report/
│   ├── template.qmd             # 파라미터화된 Quarto 리포트
│   ├── sections/                # 단계별 partial
│   └── assets/                  # CSS, JS
├── tests/
│   └── test_data/               # 6샘플 × 10유전자 fixture
├── results/                     # 파이프라인 출력
├── README.md
└── README.ko.md
```

---

## Reference

### 정량화·카운트·QC

- **Salmon** — Patro R. et al. *Salmon provides fast and bias-aware quantification of transcript expression.* Nat Methods 14, 417–419 (2017). https://doi.org/10.1038/nmeth.4197
- **tximport** — Soneson C., Love M.I., Robinson M.D. *Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences.* F1000Research 4:1521 (2015). https://doi.org/10.12688/f1000research.7563.2
- **MultiQC** — Ewels P. et al. *MultiQC: summarize analysis results for multiple tools and samples in a single report.* Bioinformatics 32(19), 3047–3048 (2016). https://doi.org/10.1093/bioinformatics/btw354

### 차등발현

- **DESeq2** — Love M.I., Huber W., Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biol. 15, 550 (2014). https://doi.org/10.1186/s13059-014-0550-8
- **apeglm** — Zhu A., Ibrahim J.G., Love M.I. *Heavy-tailed prior distributions for sequence count data: removing the noise and preserving large differences.* Bioinformatics 35(12), 2084–2092 (2019). https://doi.org/10.1093/bioinformatics/bty895

### Gene-set 분석

- **GSEA(Method)** — Subramanian A. et al. *Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.* PNAS 102(43), 15545–15550 (2005). https://doi.org/10.1073/pnas.0506580102
- **fgsea** — Korotkevich G. et al. *Fast gene set enrichment analysis.* bioRxiv (2021). https://doi.org/10.1101/060012
- **clusterProfiler 4.0** — Wu T. et al. *clusterProfiler 4.0: a universal enrichment tool for interpreting omics data.* The Innovation 2(3), 100141 (2021). https://doi.org/10.1016/j.xinn.2021.100141
- **MSigDB / Hallmark** — Liberzon A. et al. *The Molecular Signatures Database (MSigDB) Hallmark Gene Set Collection.* Cell Systems 1(6), 417–425 (2015). https://doi.org/10.1016/j.cels.2015.12.004

### 조절·pathway 활성도 추정

- **decoupleR / decoupler-py** — Badia-i-Mompel P. et al. *decoupleR: ensemble of computational methods to infer biological activities from omics data.* Bioinformatics Advances 2(1), vbac016 (2022). https://doi.org/10.1093/bioadv/vbac016
- **CollecTRI** — Müller-Dott S. et al. *Expanding the coverage of regulons from high-confidence prior knowledge for accurate estimation of transcription factor activities.* Nucleic Acids Research gkad841 (2023). https://doi.org/10.1093/nar/gkad841
- **PROGENy** — Schubert M. et al. *Perturbation-response genes reveal signaling footprints in cancer gene expression.* Nat Commun 9, 20 (2018). https://doi.org/10.1038/s41467-017-02391-6

### Drug repositioning / connectivity

- **L1000 Connectivity Map** — Subramanian A. et al. *A next generation Connectivity Map: L1000 platform and the first 1,000,000 profiles.* Cell 171(6), 1437–1452 (2017). https://doi.org/10.1016/j.cell.2017.10.049
- **L2S2** — L1000 signature 검색 (Ma'ayan Lab). https://l2s2.maayanlab.cloud
