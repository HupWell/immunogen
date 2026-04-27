# ImmunoGen 阶段汇报（2026-04-24）

## 1. 本周目标与范围

本周目标是将 ImmunoGen 从“可运行”推进到“可验收的真实后端运行”，重点覆盖 `R001`、`R002`、`R003`、`R_public_001` 四个 run。  
范围聚焦 MHC-I/MHC-II 预测与免疫原性适配器真实化，不涉及结构后端与 LNP 扩展设计。

## 2. 本周已完成

- 主流程 `run_all --target mhc_ranking` 已可稳定执行并产出关键结果文件。
- MHC-I 真实后端已落地：`NetMHCpan real_cmd` 可在本机环境稳定运行。
- MHC-II 采用 `real_tsv`（`netmhciipan` 结果表）并通过严格校验。
- BigMHC 已明确为“可选增强”，默认关闭不阻塞主流程。
- PRIME 已接入真实来源（`real_tsv`），`R001` 在严格参数下可通过。
- 已增加免疫原性真实来源强制开关，防止运行时静默回退到 proxy。

## 2.1 任务书对照状态（逐条）

### 已完成

- **输入契约（任务书第 1 章）**：`neoantigen_candidates.csv`、`hla_typing.json`、`meta.json` 已接通并稳定校验。
- **表位筛选主链（任务书第 2.1 章）**：MHC-I（NetMHCpan `real_cmd`）与 MHC-II（`real_tsv`）已落地，主流程可重复运行。
- **候选排序（任务书第 2.2 章）**：综合排序公式已执行，`peptide_mhc_ranking.csv` 稳定产出。
- **多价设计与 mRNA 设计（任务书第 2.3/2.4 章）**：已产出 `selected_peptides.csv`、`mrna_vaccine.fasta`、`mrna_design.json`。
- **结果与自证交付（任务书第 3.1/3.3 章）**：`REPORT.md`、`SELF_CHECK.md`、`POSITIVE_CONTROL.md` 已有可用产物。

### 部分完成

- **免疫原性全真实化（任务书第 2.1 章）**：DeepImmuno 与 PRIME 已接入真实来源；Repitope 当前为 `real_tsv` 链路可用，但 `real_cmd` 实时推理未闭环。
- **Simulation Hub 交付链（任务书第 3.2 章）**：`to_simhub/<case_id>/` 文件结构已生成并可交付，SimHub 回传证据归档流程需继续标准化。
- **结构后端升级进度**：`R_public_001` 已切换为 `afm` 口径，`replaces_coarse=True`，已替换 coarse 占位标记。

### 未完成

- **Repitope 本机实时模型推理（任务书第 2.1 章）**：尚未完成 `real_cmd` 级闭环。
- **mRNA 稳定性高级工具（任务书第 2.5 章可选）**：Saluki/RNAsnp 未纳入主流程。
- **LNP 递送方案深化（任务书第 2.6 章可选）**：当前仅保留说明，未做工程化设计。

## 3. 当前状态分层

### 3.1 已闭环

- MHC-I：真实实时运行（NetMHCpan `real_cmd`）。
- MHC-II：真实结果表驱动（`real_tsv`）。
- PRIME：真实来源链路可用。
- 排名主链路可重复执行并稳定产出。

### 3.2 部分闭环

- Repitope 目前可走 `real_tsv` 流程链路，但“本机实时模型推理（real_cmd）”尚未完全闭环。

### 3.3 非阻塞项

- BigMHC 当前为可选增强项，不影响主流程验收。

## 4. 主要问题与处理

- 网络稳定性导致第三方仓库和依赖下载多次中断。  
  已采用本地拉取后上传服务器、压缩包导入等方式提升成功率。
- NetMHCpan 解析曾因输出列差异失败。  
  已通过包装器参数与列识别逻辑修复，`real_cmd` 运行稳定。
- Repitope 依赖链较重，环境求解和安装稳定性不足。  
  当前已保证流程可运行，下一步聚焦最小依赖闭环 `real_cmd`。

## 5. 风险与边界

- 当前“全模块实时推理”唯一缺口是 Repitope `real_cmd`。
- 其余关键链路（MHC-I/II + PRIME）已满足真实来源要求。
- 不影响阶段性流程展示，但会影响“全实时”最终验收口径。

## 6. 下周计划

### P0（必须完成）

- 完成 Repitope `real_cmd` 运行器与环境固化。
- 在 `R001/R002/R003/R_public_001` 上使用同一严格参数进行全量回归。
- 输出统一验收记录（后端来源、通过状态、失败原因）。

### P1（同步完善）

- 同步更新 `docs/TODO.md`、`SELF_CHECK.md`、`README.md` 的状态口径。
- 固化回归命令块，形成一键复检流程。

## 7. 对外汇报一句话总结

本周已完成主链路稳定运行与核心后端真实化落地，当前仅剩 Repitope 实时推理闭环，计划下周完成并给出全量严格验收结果。

