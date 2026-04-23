# ImmunoGen Release Notes

## v0.3.1 - 表位预测“强制实跑”验收能力

- `predict_mhc_ranking.py` 新增：
  - `--require_real_mhc2`：要求 MHC-II 必须使用 `netmhciipan`
  - `--require_real_mhc1_cv`：要求 NetMHCpan/BigMHC 至少一个为真实来源
- `run_all.py` 新增同名透传参数，支持整条流水线强制验收。
- 新增 `scripts/check_epitope_realization.py`，可对 `peptide_mhc_ranking.csv` 做后端真实性检查。
- `README.md` 与 `docs/TODO.md` 同步更新为“可强制实跑、可一键验收”的操作路径。

---

## v0.3.0 - MHC-I 交叉验证列接入（NetMHCpan / BigMHC）

- `predict_mhc_ranking.py` 新增 MHC-I 交叉验证后端参数：
  - `--backend_mhc1_netmhcpan auto|real_tsv|real_cmd|off`
  - `--backend_mhc1_bigmhc auto|real_tsv|real_cmd|off`
- 排名表新增列：`mhc1_cv_netmhcpan_nM`、`mhc1_cv_bigmhc_score`、`mhc1_cv_source_netmhcpan`、`mhc1_cv_source_bigmhc`。
- 默认 `auto` 不破坏主流程：找不到真实工具/结果时回退 `off`，保留 MHCflurry 主排序逻辑。
- `run_all.py` 与 `README.md` 已同步新增参数与输入/输出说明。

---

## v0.2.9 - 结构后端选型文档落地

- 新增 `docs/structure_backend_selection.md`，给出 PANDORA 与 ColabFold/AFM 的工程化对比与默认组合策略（PANDORA 批量 + AFM 重点复核）。
- `prepare_self_certification.py` 的 `SELF_CHECK.md` 增加结构后端策略说明，便于汇报时统一口径。

---

## v0.2.8 - MHC-I α + β2m 序列数据管线

- 新增 `scripts/prepare_mhc_chain_sequences.py`：按 `hla_typing.json` 自动选择 MHC-I 目标等位基因，并拉取/缓存 IPD-IMGT/HLA（A/B/C）蛋白序列与 UniProt β2m（P61769）。
- 输出 `results/<run_id>/structure_inputs/mhc_chain_sequences.fasta` 与 `mhc_chain_sequences.json`，用于后续 PANDORA / AFM 结构建模输入。
- `run_all.py` 新增可选参数：`--prepare_structure_inputs`、`--structure_seq_refresh_remote`、`--structure_seq_strict`。

---

## v0.2.7 - 结构后端接口 + 信号肽/TM 配置

- `build_multivalent_mrna.py` / `run_all.py` 新增信号肽与 TM 可选参数，`mrna_design.json` 记录 `signal_peptide_aa`、`tm_domain_aa`、`multivalent_core_aa`。
- `prepare_simhub_delivery.py` / `run_all.py` 新增 `--structure_backend coarse|pandora|afm` 与 `--structure_input_pdb`；`meta.json` 增加 `structure_source`、`structure_tool_version`、`replaces_coarse` 字段。
- 新增 `docs/lnp_notes.md` 与 `data/templates/lnp_formulation.json`，将 LNP 说明与序列流水线解耦。
- `run_qc_and_report.py` 的 `REPORT.md` 增补信号肽/TM 使用记录与简要免疫学说明。

---

## v0.2.6 - 免疫原性真模型可选接入（安全回退）

- `immunogenicity_adapters.py` 新增真实后端策略：`real_tsv` / `real_cmd` / `proxy`，`auto` 按 real→proxy 回退。
- 三个独立适配器与批量入口均支持后端参数；`run_all.py` 透传 `--backend_deepimmuno/--backend_prime/--backend_repitope`。
- `prepare_self_certification.py` 在 `SELF_CHECK.md` 增加 `immunogenicity_source_*` 统计，明确每轮是真模型还是代理分。

---

## v0.2.5 - 免疫原性适配器 + 预计算表合并

- 新增独立适配器模块：`scripts/immunogenicity_adapters.py`，以及 `run_deepimmuno_adapter.py` / `run_prime_adapter.py` / `run_repitope_adapter.py`。
- 新增批量入口：`scripts/run_immunogenicity_adapters.py`，输出 `results/<run_id>/tool_outputs/*.tsv`。
- `predict_mhc_ranking.py` 增加按 `run_id` 读取预计算表并 merge 的逻辑，缺失自动回退 proxy；新增来源列 `immunogenicity_source_*`。
- `run_all.py` 将适配器批量步骤纳入默认流水线（校验后先生成预计算表，再预测）。

---

## v0.2.4 - 文档：`netmhciipan_setup` 简化为 conda + 三步装 NetMHC

- **`docs/netmhciipan_setup.md`**：区分「不装也能跑」与 WSL 下 **conda + 下载 + 两个环境变量** 的最短路径。

---

## v0.2.3 - NetMHCIIpan 子进程调用层

- 新增 **`scripts/netmhciipan_runner.py`**：解析 `netMHCIIpan` 4.x 标准输出表头；`NETMHCIIPAN_BIN` / `NETMHCIIPAN_HOME` 等见 **`docs/netmhciipan_setup.md`**。
- **`predict_mhc_ranking.py` / `run_all.py`**：`--mhc2_backend auto|proxy|netmhciipan`；`peptide_mhc_ranking.csv` 增加 `mhc2_el_rank`、`mhc2_ba_nm`、`mhc2_class2_allele`、`mhc2_backend`。
- **`prepare_self_certification.py`**：按排名表中的 MHC-II 实现方式更新 **SELF_CHECK** 表格行。

---

## v0.2.2 - NetMHCIIpan 等位基因映射表 + 白话文档

- 新增 **`data/hla_allele_map_netmhciipan.json`**（`manual_overrides` 手工表，默认可空）、**`scripts/hla_allele_to_netmhciipan.py`**（常见简写补 `HLA-` 前缀）。
- 新增 **`docs/allele_naming_simple.md`**（非生物学背景可读）；**`docs/hla_typing.md`** 整篇改短、与上文互链。

---

## v0.2.1 - hla_typing.json 扩展契约（MHC-II 可选键）

- 新增 **`docs/hla_typing.md`**：必备 I 类（`HLA-A/B/C`）、可选 II 类（`HLA-DRB1`、`HLA-DQA1`、`HLA-DQB1`、`HLA-DPA1`、`HLA-DPB1`），并说明 **OptiType**（仅 I）与 **HLA-HD**（I+II）到 JSON 字段的映射。
- 新增 **`data/examples/hla_typing.class_ii.example.json`**、**`scripts/hla_typing_spec.py`**；**`validate_input.py`** 校验可选 II 类键并统计 II 类条数，未知顶键仅警告。

---

## v0.2.0 - 任务书同步（自证包 + SimHub 契约）

- **自证最小包**：`results/<run_id>/` 增加 `POSITIVE_CONTROL.md`、`SELF_CHECK.md`（与 `REPORT.md` 配套）；由 `prepare_self_certification.py` 生成，`run_all.py` 在 SimHub 交付前自动执行。
- **Simulation Hub**：按 **peptide_mhc 分支 B** 交付 `complex.pdb`（链 M/B/P）+ `hla_allele.txt` + `meta.json` + `selected_for_md.csv`；**移除 `ligand.sdf`**，禁止多肽走 SDF。
- **文档**：`README.md`、`FINAL_CHECKLIST.md` 与任务书（负责人梁心恬、上下游、风险）对齐。

---

## v0.1.0 - 从 0 到可跑通（MVP）

本版本完成了 ImmunoGen 模块从输入数据到交付结果的端到端可复现流程，支持公开数据验证与基础可行性评估。

---

## 1) 目标达成情况

- 完成 BioDriver 风格输入接入：`neoantigen_candidates.csv` + `hla_typing.json` + `meta.json`
- 完成候选肽打分、筛选、多价 mRNA 构建、质控报告、SimHub 交付包生成
- 完成公开数据验证集 `R_public_001` 跑通并给出可行性结论
- 完成 GitHub 仓库初始化与首版代码发布

---

## 2) 核心流程能力（已落地）

### 2.1 输入校验
- 脚本：`scripts/validate_input.py`
- 能力：检查输入文件是否齐全、关键列与 HLA 结构是否合规

### 2.2 肽段预测与排序
- 脚本：`scripts/predict_mhc_ranking.py`
- 能力：
  - MHC-I 亲和力预测（MHCflurry）
  - MHC-II 代理评分（可替换真实工具）
  - 免疫原性代理评分（DeepImmuno/PRIME/Repitope 接口位）
  - 综合 `rank_score` 计算

### 2.3 Top 肽筛选
- 脚本：`scripts/select_top_peptides.py`
- 能力：按 `rank_score` 筛选 Top-N，过滤与 WT 过于相似候选

### 2.4 多价 mRNA 构建
- 脚本：`scripts/build_multivalent_mrna.py`
- 能力：
  - peptide linker 串联（默认 `AAY`）
  - `basic/optimized/lineardesign` 三种密码子模式
  - 输出 `mrna_vaccine.fasta` 与 `mrna_design.json`

### 2.5 图与报告
- 脚本：`scripts/run_qc_and_report.py`
- 能力：
  - 生成亲和力热图
  - 生成 mRNA 二级结构图
  - RNAfold 真实计算（优先命令行，回退 Python ViennaRNA 绑定）
  - 输出 `qc_metrics.json` 与 `REPORT.md`

### 2.5.1 自证文档
- 脚本：`scripts/prepare_self_certification.py`
- 能力：生成 `POSITIVE_CONTROL.md`（KRAS G12D 等）、`SELF_CHECK.md`（双工具、阈值、局限）

### 2.6 Simulation Hub 交付
- 脚本：`scripts/prepare_simhub_delivery.py`
- 能力：
  - 生成 `to_simhub/<case_id>/` 标准目录
  - 输出 **`complex.pdb`**、`hla_allele.txt`、`meta.json`、`selected_for_md.csv`
  - 当前 `complex.pdb` 为粗初始构象（链 M/B/P），便于接口验收；生产请换 AlphaFold-Multimer / PANDORA

### 2.7 可行性验证
- 脚本：`scripts/validate_feasibility.py`
- 能力：
  - Top 稳定性（Jaccard）分析
  - Positive control 命中检查（KRAS G12D）
  - 输出 `feasibility_report.json` 与 `FEASIBILITY.md`

### 2.8 一键执行
- 脚本：`scripts/run_all.py`
- 能力：串行执行全流程（校验 → 预测 → 筛选 → 构建 → 报告 → **自证** → SimHub → 可行性）

---

## 3) 公开数据验证里程碑

### 3.1 数据集构建
- 使用 pVACtools 官方示例结果 + 文献经典 neoantigen（KRAS/TP53）
- 通过 `merge_public_sources.py` 与 `prepare_public_dataset.py` 构建 `R_public_001`

### 3.2 验证结果（R_public_001）
- 全流程成功跑通
- `selected_peptides.csv` 中包含 TP53 热点候选
- `FEASIBILITY.md` 结论：
  - 工程可复现：`True`
  - 方向性可行：`True`
  - KRAS G12D 命中：`True`

---

## 4) 关键产物（已可交付）

- `results/<run_id>/peptide_mhc_ranking.csv`
- `results/<run_id>/selected_peptides.csv`
- `results/<run_id>/mrna_vaccine.fasta`
- `results/<run_id>/mrna_design.json`
- `results/<run_id>/qc_metrics.json`
- `results/<run_id>/REPORT.md`
- `results/<run_id>/POSITIVE_CONTROL.md`
- `results/<run_id>/SELF_CHECK.md`
- `results/<run_id>/FEASIBILITY.md`
- `deliveries/<run_id>/to_simhub/<case_id>/complex.pdb` 等（无 SDF）

---

## 5) 已知边界与下一阶段计划

### 当前边界
- MHC-II 与免疫原性评分当前包含代理实现（用于流程验证）
- `complex.pdb` 为粗初始构象，不等同于 AlphaFold/PANDORA 精细结构
- 结论定位为“计算方向验证”，不替代湿实验结论

### 下一阶段
- 接入真实 NetMHCIIpan 结果
- 接入真实 DeepImmuno/PRIME/Repitope 推理结果
- 接入高精度 peptide-MHC 结构建模
- 增加多病例批处理与更严格 QC

---

## 6) 导师汇报一句话版本

ImmunoGen 模块已完成从候选新抗原输入到多价 mRNA 设计与 SimHub 交付包生成的端到端可复现流程，并在公开验证集上验证了工程可复现性与方向性可行性，为后续接入真实模型和高精度结构建模打下了稳定基础。
