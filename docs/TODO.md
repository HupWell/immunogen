# ImmunoGen 任务书对照 TODO（白话版）

> 这份清单按任务书逐条对照，分成：  
> - ✅ 已完成（可直接用）  
> - 🟡 部分完成（框架有了，但真实工具/真实数据还没完全接上）  
> - ⬜ 未完成（后续迭代）

---

## 1) 当前结论（先看这个）

- ✅ 主流程已经能完整跑通：输入校验 → 排名 → Top 筛选 → mRNA 设计 → 报告/自证 → SimHub 交付。
- ✅ `R001`、`R002`、`R003`、`R_public_001` 都已产出完整结果目录。
- 🟡 多个“真实外部工具”还在“可选接入 + 自动回退”状态（不是默认强制实跑）。
- ✅ 已按任务书要求明确：peptide-MHC 分支禁止 SDF，交付 `complex.pdb`（M/B/P 链）。
- ✅ 已核对 `R_public_001`：`results` 三件套和 `to_simhub/public_case_001` 契约文件齐全，`molecule_type=peptide_mhc` 且无 `ligand.sdf`。

---

## 2) 任务书逐项对照

### 2.1 输入契约（BioDriver 输入）

- ✅ 已完成  
  - `neoantigen_candidates.csv` / `hla_typing.json` / `meta.json` 输入契约已落地并校验。
  - II 类 HLA 字段（DRB1/DQA1/DQB1/DPA1/DPB1）已纳入契约说明与示例。

### 2.2 表位预测与排序

- ✅ 已完成（可用）  
  - MHC-I 主预测：`MHCflurry`。
  - 综合排序：`rank_score = w1*HLA + w2*immunogenicity + w3*VAF + w4*dissimilarity`。
  - 已做与 WT 过于相似的过滤（`select_top_peptides.py`）。

- ✅ 已完成（可验收）  
  - 已增加“强制实跑”开关：`--require_real_mhc2`、`--require_real_mhc1_cv`。  
  - 已增加验收脚本：`scripts/check_epitope_realization.py`（可一键检查 MHC-II/MHC-I 是否真实后端）。
  - 说明：如果你的环境未接好真实工具，这些强制参数会直接报错，避免“看起来跑通但其实是回退 proxy/off”。

### 2.3 免疫原性（DeepImmuno/PRIME/Repitope）

- ✅ 已完成（工程框架）  
  - 三个工具都支持 `auto/real_tsv/real_cmd/proxy`。
  - 输出列与来源追踪列已统一。

- 🟡 部分完成（真实模型层面）  
  - 当前大多数实例仍以 `proxy` 或公开数据映射为主，不等于你本地已部署完整真模型推理链。

### 2.4 多价 ORF + mRNA 设计

- ✅ 已完成（MVP可交付）  
  - Top 肽串联、linker、UTR、polyA、密码子模式（`basic/optimized/lineardesign`）已支持。
  - 可选信号肽/TM 参数已支持并写入 `mrna_design.json`。
  - 二级结构图和基础 QC 已输出。

- ⬜ 未完成（任务书“更真实工具”层面）  
  - `LinearDesign/COOL` 真正生产级接入与严格回归验证仍待做。
  - Saluki / RNAsnp 稳定性接入还未做。

### 2.5 SimHub 交付

- ✅ 已完成（契约层）  
  - 交付目录、`complex.pdb` / `hla_allele.txt` / `meta.json` / `selected_for_md.csv` 已落地。
  - `prepare_simhub_delivery.py` 支持 `coarse/pandora/afm`，并记录结构来源元数据。
  - 链 ID 契约（M/B/P）已处理。

- 🟡 部分完成（真实结构层）  
  - 默认仍是 `coarse` 可跑通模式。
  - PANDORA/AFM 的真实批量建模与稳定评估（RMSD/H-bond/ΔG）还需你实际环境继续做。

### 2.6 自证最小包与验收

- ✅ 已完成  
  - `POSITIVE_CONTROL.md`、`SELF_CHECK.md`、`REPORT.md` 自动生成已落地。
  - IEDB 资源入口已写入 `SELF_CHECK`。

---

## 3) 目前最关键“未完全落地”项（按优先级）

### P0（优先先做）

1. 🟡 让 MHC-I 交叉验证真正“有值”  
   - 目标：`mhc1_cv_source_netmhcpan` / `mhc1_cv_source_bigmhc` 从 `off` 变成 `real_cmd` 或 `real_tsv`。
   - 现状：框架已完善（含 `run_all --target mhc_ranking`、前置校验、`mhc1_cv_source`/`mhc1_cv_tool` 汇总列、meta.json），但真实命令/真实结果文件仍待接入。

2. 🟡 让 MHC-II 在真实数据上真正跑 NetMHCIIpan  
   - 目标：`mhc2_backend` 在真实病例里出现 `netmhciipan`，不是长期 `proxy`。

3. 🟡 免疫原性从 proxy 逐步迁移到真模型  
   - 目标：`immunogenicity_source_*` 尽量由 `proxy` 转为 `real_tsv/real_cmd`。

### P1

4. 🟡 真实结构建模批量化（PANDORA/AFM）  
   - 目标：`coarse` 仅用于联调，关键交付案例换成真实 `complex.pdb`。

### P2

5. ⬜ 生产级密码子优化工具接入与评估（LinearDesign/COOL）。  
6. ⬜ mRNA 稳定性（Saluki/RNAsnp）接入 `qc_metrics.json` 和报告。  

---

## 4) 下一步执行清单（简单版）

- [ ] 在远程环境配置 `MHC1_NETMHCPAN_CMD` / `MHC1_BIGMHC_CMD`，或准备 `results/<run_id>/tool_outputs/raw/mhc1_*.tsv`，让交叉验证列不再是 `off`。  
- [ ] 在一个真实 `run_id` 上使用 `--mhc2_backend netmhciipan --require_real_mhc2 --require_real_mhc1_cv` 跑通。  
- [ ] 选择一个免疫原性工具先从 `proxy` 升级到 `real_tsv`（建议 DeepImmuno 先行）。  
- [ ] 选 Top 3-5 候选跑 PANDORA/AFM，替换 `coarse` 结构再交 SimHub。  
- [ ] 为 `R_public_001` 增加“阳性对照未命中 KRAS 时的替代阳性证据”说明（当前 `POSITIVE_CONTROL.md` 显示 KRAS 未命中）。  

建议最短路径（先打通 MHC-I）：
- [ ] `python scripts/run_all.py --run_id R001 --target mhc_ranking --mhc2_backend proxy --backend_mhc1_netmhcpan real_tsv --backend_mhc1_bigmhc off --require_real_mhc1_cv`
- [ ] `python scripts/check_epitope_realization.py --run_id R001 --require_mhc1_cv_real`

---

## 5) 说明

- 本清单以当前仓库和最近实跑日志为准，不按“是否有接口”而按“是否真的跑出真实工具结果”来区分完成度。
- 每完成一项，请同步更新：`RELEASE_NOTES.md`、`SELF_CHECK.md`、本文件。
- 本次文档核对时间：`2026-04-23`（远程服务器会话）。
