# ImmunoGen 后续完善清单（TODO）

本文档对应任务书中**尚未在主线代码中完全落地**的能力，便于排期与分工。当前 MVP 仍以 MHCflurry（MHC-I）、代理 MHC-II、代理免疫原性、粗粒度 `complex.pdb` 为主；见 `SELF_CHECK.md`、`RELEASE_NOTES.md`。

---

## 优先级说明

- **P0**：对排名/交付质量影响大、接口相对明确  
- **P1**：增强可解释性与 SimHub 可信度  
- **P2**：设计扩展与文档级能力  

建议实施顺序：**P0 MHC-II 真实工具 → P1 结构后端 → P0 免疫原性（择一先接）→ P2 信号肽配置 → P2 LNP 文档**。

---

## P0：MHC-II — NetMHCIIpan

| 状态 | 事项 |
| :---: | --- |
| [x] | 在 `hla_typing.json` 契约中明确 **II 类等位基因**（`HLA-DRB1`、`HLA-DQA1`、`HLA-DQB1`、`HLA-DPA1`、`HLA-DPB1`）及与 OptiType/HLA-HD 的字段映射 → **`docs/hla_typing.md`** + `data/examples/hla_typing.class_ii.example.json` + `scripts/hla_typing_spec.py` + `validate_input.py` |
| [x] | 维护 **等位基因名称映射表**（JSON，YAML 可用同结构）：BioDriver/契约 → NetMHCIIpan → `data/hla_allele_map_netmhciipan.json` + `scripts/hla_allele_to_netmhciipan.py`；白话说明 **`docs/allele_naming_simple.md`**，契约仍见 **`docs/hla_typing.md`** |
| [x] | 实现 **NetMHCIIpan 调用层** → `scripts/netmhciipan_runner.py`（`subprocess`，环境变量见 `docs/netmhciipan_setup.md`） |
| [x] | 批处理：同一次调用写入肽列表、`-a` 多个等位基因；**stdout 表**解析后并入 `peptide_mhc_ranking.csv`（`mhc2_el_rank` / `mhc2_ba_nm` 等；EL %Rank 优于 IC50 时以 rank 主转 `mhc2_score`） |
| [x] | `predict_mhc_ranking.py` / `run_all.py`：`--mhc2_backend auto|proxy|netmhciipan`；**auto** 在缺可执行/无 II 类时回退 `mhc2_proxy_score` |
| [x] | `prepare_self_certification.py` 根据 `mhc2_backend` 列更新 **SELF_CHECK** 中 MHC-II 说明；版本/命令以环境为准，见 `netmhciipan_setup.md` |

**相关代码：** `scripts/predict_mhc_ranking.py` 中 `mhc2_proxy_score()`。

---

## P0：免疫原性 — DeepImmuno / PRIME / Repitope

| 状态 | 事项 |
| :---: | --- |
| [x] | 统一 **输出列契约**：`immunogenicity_deepimmuno`、`immunogenicity_prime`、`immunogenicity_repitope`、综合列 `immunogenicity`；加权支持 `predict_mhc_ranking.py --wi_deepimmuno/--wi_prime/--wi_repitope`（保留旧列别名兼容历史流程） |
| [x] | 采用 **适配器 + 预计算表** 策略：新增独立脚本 `run_deepimmuno_adapter.py` / `run_prime_adapter.py` / `run_repitope_adapter.py`（批量入口 `run_immunogenicity_adapters.py`），输出 `results/<run_id>/tool_outputs/*.tsv`；`predict_mhc_ranking.py` 按 `run_id` 自动 merge，缺失回退 proxy |
| [x] | **DeepImmuno**：封装可选真模型 runner（`real_tsv` / `real_cmd` / `proxy`），命令失败自动回退（`auto`） |
| [x] | **PRIME**：同样支持可选真模型 runner（`real_tsv` / `real_cmd` / `proxy`） |
| [x] | **Repitope**：同上 |
| [x] | `predict_mhc_ranking.py`：预计算文件优先 merge，缺失回退 `deepimmuno_proxy` / `prime_proxy` / `repitope_proxy` |
| [x] | 文档与自证：`SELF_CHECK.md` 增加 `immunogenicity_source_*` 统计，标明真模型/proxy 来源 |

**相关代码：** `scripts/predict_mhc_ranking.py` 中 `deepimmuno_proxy`、`prime_proxy`、`repitope_proxy` 及 `rank_score` 组合逻辑。

---

## P1：Simulation Hub — 真实 peptide-MHC 结构（AlphaFold-Multimer / PANDORA）

| 状态 | 事项 |
| :---: | --- |
| [ ] | 从 **IPD IMGT/HLA** 或内部序列表，按分型拉取 **MHC-I α 链 + β2m** 氨基酸序列（数据管线，非单次硬编码） |
| [ ] | 评估 **PANDORA** vs **ColabFold/AlphaFold-Multimer** 二选一或组合（模板可用性、许可、批量耗时） |
| [ ] | 实现独立步骤或 CLI 标志：`prepare_simhub_delivery.py --structure_backend pandora|afm|coarse`，**默认 coarse** 保证无依赖可跑 |
| [ ] | 将工具输出 PDB **整理为契约链 ID**：Chain **M**（α）、**B**（β2m）、**P**（peptide）；写入 `deliveries/.../complex.pdb` |
| [ ] | 在 `meta.json` 中增加 `structure_source`、`structure_tool_version`、`replaces_coarse` 等字段便于 SimHub 与审计 |
| [ ] | CI/文档：说明 GPU/队列要求，避免默认绑进 `run_all.py` 长耗时步骤 |

**相关代码：** `scripts/prepare_simhub_delivery.py` 中 `write_complex_pdb()`。

---

## P2：多价 ORF — 信号肽 / 跨膜（可选）

| 状态 | 事项 |
| :---: | --- |
| [ ] | 在 `mrna_design.json` 中增加字段：`signal_peptide_aa`、`tm_domain_aa`（可选）、`multivalent_core_aa`、`linker` |
| [ ] | `build_multivalent_mrna.py`：支持 CLI 或 JSON 配置 **N 端信号肽**（如 MITD / 经典分泌肽预设名）、**C 端可选 TM** |
| [ ] | 保持氨基酸顺序文档化：`[signal][pep1][linker][pep2]...[TM?]` → 再密码子化与 UTR 拼接 |
| [ ] | `REPORT.md` 模板中简述选择信号肽的 **免疫学/表达** 依据（引用即可，不必过长） |

**相关代码：** `scripts/build_multivalent_mrna.py`。

---

## P2：LNP 与递送（第一阶段以文档为主）

| 状态 | 事项 |
| :---: | --- |
| [ ] | 在 `REPORT.md` 或独立 `docs/lnp_notes.md` 中固定 **综述级** 说明（如 SM-102、ALC-0315），明确「本仓库不生成 LNP 处方」 |
| [ ] | 若第二阶段需要结构化：新增 `lnp_formulation.json` 模板，与序列流水线 **解耦** |

---

## P2：其他工具链（任务书提及、当前为占位或简化）

| 状态 | 事项 |
| :---: | --- |
| [ ] | **MHC-I**：可选接入 NetMHCpan-4.1 / BigMHC，与 MHCflurry **交叉验证**列写入 `peptide_mhc_ranking.csv` |
| [ ] | **密码子优化**：对接 LinearDesign / COOL 等外部工具时，保留现有 `basic` / `optimized` / `lineardesign` 模式语义与回退 |
| [ ] | **IEDB**：表位相关资源以链接或导出步骤形式写入 `SELF_CHECK.md` / 操作手册（非必须自动化） |
| [ ] | **mRNA 稳定性**：可选 Saluki / RNAsnp，结果进入 `qc_metrics.json` 与图表说明 |

---

## 验收提醒（每完成一大项）

- [ ] 更新 `RELEASE_NOTES.md` 版本说明  
- [ ] 更新 `SELF_CHECK.md` / `prepare_self_certification.py` 中阈值与工具记录  
- [ ] 在干净环境中按 `README.md` 跑通至少一个 `run_id`（含 `--help` 与新参数说明）  

---

*最后更新：与任务书「尚未完全落地」章节对齐，随实现进度勾选并修订。*
