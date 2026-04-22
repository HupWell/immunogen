# ImmunoGen（mRNA 疫苗 / 免疫基因）

**负责人：** 梁心恬  
**模块定位：** 基于病人特异性新抗原（Neoantigen）设计多价 mRNA 疫苗。  
**一句话目标：** 从 BioDriver 的候选新抗原 + HLA 分型，输出候选 mRNA 全序列（含 UTR、密码子优化、多价串联设计），并按要求准备 Simulation Hub 交付包。

**上游：** `deliveries/<run_id>/to_immunogen/`（`neoantigen_candidates.csv`、`hla_typing.json`、`meta.json`）  
**下游：** `deliveries/<run_id>/to_simhub/<case_id>/`（peptide-MHC，`molecule_type: peptide_mhc`）

---

## 0. 一分钟导读

本仓库实现：读入候选突变肽与 HLA → 打分排序 → 选 Top-N → linker 串联成多价 ORF → 拼接 UTR / polyA → 二级结构（MFE）与图表 → **自证文档（POSITIVE_CONTROL / SELF_CHECK / REPORT）** → **按 SimHub 契约生成交付包（无 SDF）**。

**第一阶段 smoke：** 跑通 1 个 `run_id` → Top-10（可调）→ 1 条多价 mRNA → 至少 1 套 SimHub 目录（`complex.pdb` 等）。  
**交付节奏：** 先完成 `results/<run_id>/` 下自证三件套，再提交 `to_simhub/`。

---

## 1. 你会得到什么

运行 `run_all.py` 后，在 `results/<run_id>/`：

- `peptide_mhc_ranking.csv`
- `selected_peptides.csv`
- `mrna_vaccine.fasta`
- `mrna_design.json`
- `qc_metrics.json`
- `REPORT.md`
- **`POSITIVE_CONTROL.md`**（已知强免疫原 neoantigen 自证）
- **`SELF_CHECK.md`**（双工具记录、阈值、局限）
- `figures/binding_affinity_heatmap.png`
- `figures/mrna_secondary_structure.png`（标题含 MFE）
- `feasibility_report.json`、`FEASIBILITY.md`

在 `deliveries/<run_id>/to_simhub/<case_id>/`（**契约分支 B，禁止 SDF**）：

- `complex.pdb` — Chain **M**（MHC α 粗迹）、**B**（β2m 代理）、**P**（肽段）；生产请换 AlphaFold-Multimer / PANDORA
- `hla_allele.txt` — 可选，当前为 Top1 对应等位基因一行
- `meta.json` — `molecule_type: "peptide_mhc"`
- `selected_for_md.csv` — 送 MD 的 Top-K 列表

---

## 2. 环境准备（Windows + Conda）

```bash
conda create -n immunogen python=3.10 -y
conda activate immunogen
pip install pandas numpy biopython mhcflurry matplotlib seaborn ViennaRNA
set PYTHONUTF8=1
mhcflurry-downloads fetch
```

说明：

- **`ViennaRNA`（pip）**：用于在无 `RNAfold` 可执行文件时计算 MFE 与二级结构（与任务书 RNAfold 目标一致）。
- Windows 下建议每次运行前执行 `set PYTHONUTF8=1`，避免编码问题。

---

## 3. 输入目录

`deliveries/<run_id>/to_immunogen/` 下必须包含：

- `neoantigen_candidates.csv`（列：`mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf`）
- `hla_typing.json`（`HLA-A` / `HLA-B` / `HLA-C` 数组）
- `meta.json`（建议含 `case_id`，用于 `to_simhub` 子目录名）

---

## 4. 一键运行（推荐）

```bash
conda activate immunogen
set PYTHONUTF8=1
python scripts/run_all.py --run_id R_public_001 --top_n 20 --top_k_md 5 --feasibility_top_n 20
```

---

## 5. 分步运行（调试用）

```bash
python scripts/validate_input.py --run_id R001
python scripts/predict_mhc_ranking.py --run_id R001
python scripts/select_top_peptides.py --run_id R001 --top_n 10 --min_dissimilarity 0.1
python scripts/build_multivalent_mrna.py --run_id R001 --linker AAY --poly_a_len 120 --codon_mode optimized
python scripts/run_qc_and_report.py --run_id R001
python scripts/prepare_self_certification.py --run_id R001
python scripts/prepare_simhub_delivery.py --run_id R001 --top_k 3
python scripts/validate_feasibility.py --run_id R001 --top_n 10
```

---

## 6. 公开数据基准

- `data/public/`：pVACtools 示例 TSV + 文献案例合并流程见 `scripts/merge_public_sources.py`、`scripts/prepare_public_dataset.py`。
- 推荐验证 `run_id`：`R_public_001`。

---

## 7. 风险与边界（任务书对齐）

- MHC-II 预测准确率低于 MHC-I；当前含代理分，报告见 `SELF_CHECK.md`。
- AI / 代理免疫原性不确定性大，结论不替代湿实验（ELISpot / 四聚体等）。
- **禁止**对多肽使用 SDF + AM1-BCC 等小分子流程；SimHub 走 **AMBER14SB 纯蛋白力场**。
- 本阶段不做 LNP 详细设计。

---

## 8. 相关文档

- `FINAL_CHECKLIST.md` — 交付前自检
- `RELEASE_NOTES.md` — 版本与里程碑
- `docs/TODO.md` — **后续完善清单**（NetMHCIIpan、免疫原性真工具、AF/PANDORA、信号肽/LNP 等）
