# ImmunoGen Release Notes

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
