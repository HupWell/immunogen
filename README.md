# ImmunoGen MVP

本项目用于将 BioDriver 给出的候选突变肽与 HLA 分型，转为可交付的 mRNA 疫苗设计基础结果。

## 1. 你会得到什么

运行完成后会在 `results/<run_id>/` 生成：

- `peptide_mhc_ranking.csv`
- `selected_peptides.csv`
- `mrna_vaccine.fasta`
- `mrna_design.json`
- `qc_metrics.json`
- `REPORT.md`
- `feasibility_report.json`
- `FEASIBILITY.md`

## 2. 环境准备（Windows + Conda）

```bash
conda create -n immunogen python=3.10 -y
conda activate immunogen
pip install pandas numpy biopython mhcflurry viennarna
set PYTHONUTF8=1
mhcflurry-downloads fetch
```

说明：
- 若你已安装好环境，可直接 `conda activate immunogen`。
- 在 Windows 下建议每次运行前执行 `set PYTHONUTF8=1`，避免编码问题。

## 3. 输入目录结构

将输入文件放到：

`deliveries/<run_id>/to_immunogen/`

必须包含：

- `neoantigen_candidates.csv`
- `hla_typing.json`
- `meta.json`

## 4. 一键运行（推荐）

在项目根目录执行：

```bash
conda activate immunogen
set PYTHONUTF8=1
python scripts/run_all.py --run_id R001 --top_n 10 --min_dissimilarity 0.1 --linker AAY --poly_a_len 120
```

## 5. 分步运行（调试用）

```bash
python scripts/validate_input.py --run_id R001
python scripts/predict_mhc_ranking.py --run_id R001
python scripts/select_top_peptides.py --run_id R001 --top_n 10 --min_dissimilarity 0.1
python scripts/build_multivalent_mrna.py --run_id R001 --linker AAY --poly_a_len 120 --codon_mode optimized
python scripts/run_qc_and_report.py --run_id R001
python scripts/prepare_simhub_delivery.py --run_id R001 --top_k 3
python scripts/validate_feasibility.py --run_id R001 --top_n 10
```

## 6. 当前版本边界（MVP）

- `rank_score` 目前以 MHC 亲和力为主，后续可加入 immunogenicity、VAF、WT 差异加权。
- 免疫原性的 DeepImmuno/PRIME/Repitope 当前使用可复现代理评分；如有真实模型可替换。
- mRNA 支持 `basic/optimized/lineardesign` 三种模式；`lineardesign` 需提供本地可执行文件。
- `figures` 当前可生成真实热图；二级结构图在未安装 RNAfold 时使用回退 profile。
