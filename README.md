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
- `hla_typing.json` — **必备** `HLA-A` / `HLA-B` / `HLA-C`；**可选** II 类：`HLA-DRB1`、`HLA-DQA1`、`HLA-DQB1`、`HLA-DPA1`、`HLA-DPB1`（供 NetMHCIIpan 等；OptiType 仅 I 类，II 常由 HLA-HD 等补充）。完整说明与上游映射见 **`docs/hla_typing.md`**，含 II 示例 `data/examples/hla_typing.class_ii.example.json`。
- `meta.json`（建议含 `case_id`，用于 `to_simhub` 子目录名）

**NetMHCIIpan（MHC-II，可选）**：不装也能跑（`auto` / `proxy`）。要真 II 类预测：conda 管 Python、工具单独下载，**最简三步**见 **`docs/netmhciipan_setup.md`**（情况 A / 情况 B）。

---

## 4. 一键运行（推荐）

```bash
conda activate immunogen
set PYTHONUTF8=1
python scripts/run_all.py --run_id R_public_001 --top_n 20 --top_k_md 5 --feasibility_top_n 20 --mhc2_backend auto --backend_mhc1_netmhcpan auto --backend_mhc1_bigmhc auto --wi_deepimmuno 1 --wi_prime 1 --wi_repitope 1 --signal_peptide_preset mitd --prepare_structure_inputs --structure_backend coarse
```

### 4.1 分段执行（新增 `--target`）

当你只想调某一段时，不必全流程重跑：

```bash
# 只跑到 MHC 排名（含输入校验、免疫原性适配器、predict_mhc_ranking）
python scripts/run_all.py --run_id R001 --target mhc_ranking --mhc2_backend proxy

# 只补报告与自证
python scripts/run_all.py --run_id R001 --target report

# 只封装 SimHub 交付（可配合 --prepare_structure_inputs）
python scripts/run_all.py --run_id R001 --target simhub --structure_backend coarse

# 只跑可行性验证
python scripts/run_all.py --run_id R001 --target feasibility --feasibility_top_n 10
```

---

## 5. 分步运行（调试用）

```bash
python scripts/validate_input.py --run_id R001
python scripts/export_immuno_template.py --run_id R001
python scripts/run_immunogenicity_adapters.py --run_id R001 --backend_deepimmuno auto --backend_prime auto --backend_repitope auto
python scripts/predict_mhc_ranking.py --run_id R001 --wi_deepimmuno 1 --wi_prime 1 --wi_repitope 1 --backend_mhc1_netmhcpan auto --backend_mhc1_bigmhc auto
python scripts/select_top_peptides.py --run_id R001 --top_n 10 --min_dissimilarity 0.1
python scripts/build_multivalent_mrna.py --run_id R001 --linker AAY --poly_a_len 120 --codon_mode optimized --signal_peptide_preset mitd --tm_domain_preset cd8a_tm
python scripts/run_qc_and_report.py --run_id R001
python scripts/prepare_self_certification.py --run_id R001
python scripts/prepare_mhc_chain_sequences.py --run_id R001 --strict
python scripts/prepare_simhub_delivery.py --run_id R001 --top_k 3 --structure_backend coarse
python scripts/validate_feasibility.py --run_id R001 --top_n 10
python scripts/check_epitope_realization.py --run_id R001 --require_mhc2_real --require_mhc1_cv_real
```

说明（免疫原性真模型）：
- `run_immunogenicity_adapters.py` 支持每个工具独立后端：`auto` / `real_tsv` / `real_cmd` / `proxy`。
- `real_tsv`：读取 `results/<run_id>/tool_outputs/raw/{deepimmuno|prime|repitope}.tsv`。
- `real_cmd`：设置环境变量命令模板（示例）：
  - `IMMUNO_DEEPIMMUNO_CMD="python tools/deep_runner.py --input {input_tsv} --output {output_tsv}"`
  - `IMMUNO_PRIME_CMD="python tools/prime_runner.py --input {input_tsv} --output {output_tsv}"`
  - `IMMUNO_REPITOPE_CMD="python tools/repitope_runner.py --input {input_tsv} --output {output_tsv}"`
  命令需要产出 TSV，至少包含 `mut_peptide` 和分数字段（`score` 或对应 `immunogenicity_*`）。

说明（结构后端与 ORF 可选段）：
- `build_multivalent_mrna.py` 新增可选参数：`--signal_peptide_preset/--signal_peptide_aa`、`--tm_domain_preset/--tm_domain_aa`。
- ORF 拼接顺序为：`[signal][pep1][linker][pep2]...[TM?]`，并写入 `mrna_design.json` 字段：`signal_peptide_aa`、`tm_domain_aa`、`multivalent_core_aa`。
- `prepare_simhub_delivery.py` 支持 `--structure_backend coarse|pandora|afm`：
  - `coarse`：脚本生成粗粒度 `complex.pdb`（默认，无外部依赖）；
  - `pandora/afm`：需额外提供 `--structure_input_pdb`（外部工具结果），脚本仅做交付封装。
- GPU/队列建议：`pandora/afm` 往往是长耗时步骤，建议单独排队执行并把结果 PDB 回填，不建议在默认 `run_all.py` 里直接绑定重型结构预测。
- `prepare_mhc_chain_sequences.py` 可按 `hla_typing.json` 自动拉取并缓存 MHC-I α 链与 β2m 序列（来源：IPD-IMGT/HLA + UniProt），输出到 `results/<run_id>/structure_inputs/`。
- 若网络受限，可先跑一次联网拉取生成 `data/reference/ipd_imgt_hla_cache/`，后续离线复用缓存。

说明（MHC-I 交叉验证）：
- `predict_mhc_ranking.py` 已支持可选交叉验证列：
  - `mhc1_cv_netmhcpan_nM`（越小越好）
  - `mhc1_cv_bigmhc_score`（越大越好）
  - `mhc1_cv_source_netmhcpan`、`mhc1_cv_source_bigmhc`
- 后端参数：`--backend_mhc1_netmhcpan`、`--backend_mhc1_bigmhc`，可选 `auto/real_tsv/real_cmd/off`。
- 如需“强制实跑验收”，可增加：
  - `--require_real_mhc2`（要求 MHC-II 必须 netmhciipan）
  - `--require_real_mhc1_cv`（要求 NetMHCpan/BigMHC 至少一个 real）
- `real_tsv` 默认读取：
  - `results/<run_id>/tool_outputs/raw/mhc1_netmhcpan.tsv`
  - `results/<run_id>/tool_outputs/raw/mhc1_bigmhc.tsv`
  - 可直接参考模板：`data/examples/mhc1_cv_templates/mhc1_netmhcpan.tsv`、`data/examples/mhc1_cv_templates/mhc1_bigmhc.tsv`
- `real_tsv` / `real_cmd` 会触发 `run_all.py` 前置检查：
  - 设为 `real_tsv` 时，若对应 TSV 缺失会在执行前直接报错；
  - 设为 `real_cmd` 时，若缺少环境变量 `MHC1_NETMHCPAN_CMD` / `MHC1_BIGMHC_CMD` 会在执行前直接报错。
- `real_cmd` 通过环境变量命令模板调用（示例）：
  - `MHC1_NETMHCPAN_CMD="python tools/netmhcpan_runner.py --input {input_tsv} --output {output_tsv}"`
  - `MHC1_BIGMHC_CMD="python tools/bigmhc_runner.py --input {input_tsv} --output {output_tsv}"`
  输出 TSV 至少包含 `mut_peptide` 与对应分值列（可带 `hla_allele`）。
- 结果 CSV 新增汇总列：`mhc1_cv_source`、`mhc1_cv_tool`（便于快速确认是否真接入）。
- 每次预测会写元数据：
  - `results/<run_id>/tool_outputs/raw/mhc1_cv_netmhcpan.meta.json`
  - `results/<run_id>/tool_outputs/raw/mhc1_cv_bigmhc.meta.json`
- 跑完后可用 `scripts/check_epitope_realization.py` 一条命令验收。

快速验收示例（以 NetMHCpan TSV 为例）：

```bash
python scripts/run_all.py --run_id R001 --target mhc_ranking --mhc2_backend proxy --backend_mhc1_netmhcpan real_tsv --backend_mhc1_bigmhc off --require_real_mhc1_cv
python scripts/check_epitope_realization.py --run_id R001 --require_mhc1_cv_real
```

---

## 6. 公开数据基准

- `data/public/`：pVACtools 示例 TSV + 文献案例合并流程见 `scripts/merge_public_sources.py`、`scripts/prepare_public_dataset.py`。
- 免疫原性公开对照数据（DeepImmuno）已纳入：`data/public/immunogenicity/deepimmuno/`，可用 `scripts/import_public_deepimmuno_scores.py` 映射到 `results/<run_id>/tool_outputs/raw/deepimmuno.tsv`。
- PRIME 公开对照数据（测试输出）已纳入：`data/public/immunogenicity/prime/`，可用 `scripts/import_public_prime_scores.py` 映射到 `results/<run_id>/tool_outputs/raw/prime.tsv`。
  若拿到 PRIME 补充表（如 Table S4），可按 `data/public/immunogenicity/prime/supplementary/column_mapping.template.json` 复制成 `column_mapping.json` 后一键导入。
- 推荐验证 `run_id`：`R_public_001`。

---

## 7. 风险与边界（任务书对齐）

- MHC-II 预测准确率低于 MHC-I；当前含代理分，报告见 `SELF_CHECK.md`。
- AI / 代理免疫原性不确定性大，结论不替代湿实验（ELISpot / 四聚体等）。
- **禁止**对多肽使用 SDF + AM1-BCC 等小分子流程；SimHub 走 **AMBER14SB 纯蛋白力场**。
- 本阶段不做 LNP 详细设计。
- LNP 仅保留综述说明，见 `docs/lnp_notes.md`。

---

## 8. 相关文档

- `FINAL_CHECKLIST.md` — 交付前自检
- `RELEASE_NOTES.md` — 版本与里程碑
- `docs/TODO.md` — **后续完善清单**（NetMHCIIpan、免疫原性真工具、AF/PANDORA、信号肽/LNP 等）
- `docs/hla_typing.md` — **`hla_typing.json` 字段说明**（I 类必备、II 类可选、与上游字段对照）
- `docs/allele_naming_simple.md` — **等位基因/NetMHC 写给非生信背景**（简写、映射表、和 JSON 怎么配合）
- `data/hla_allele_map_netmhciipan.json` — **BioDriver → NetMHCIIpan** 工具有时对不上的「手工改一行」表（可空，靠脚本默认补 `HLA-`）
- `docs/netmhciipan_setup.md` — **安装 NetMHCIIpan、环境变量、与 `--mhc2_backend` 衔接**
- `docs/lnp_notes.md` — **LNP 综述级说明（与序列流水线解耦）**
- `docs/structure_backend_selection.md` — **PANDORA vs AFM 选型与落地策略（含 GPU 何时需要）**
