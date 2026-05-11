# ImmunoGen 任务书自检与 TODO

> **负责人：** 梁心恬（与 `README.md` 一致；若组织另有指派可改署名为「模块负责人」）  
> **文档性质：** 对照任务书做**自检**，本文**即为当前唯一 TODO**，旧版条目全部作废。  
> **本文重写日期：** 2026-05-11  
> **口径依据：** 仓库脚本、`results/`、`deliveries/` 与任务书当前版（**不要求** SimHub MD 回传闭环；**HLA-II 临床终版**上游暂未提供）。

---

## 一、执行摘要（给 Lead）

| 维度 | 状态 |
|------|------|
| **四 run 主链路** | `R001` / `R002` / `R003` / `R_public_001`：输入→排名→筛选→多价 mRNA（LinearDesign `real_cmd`）→ ViennaRNA 质控→自证→SimHub B 分支交付，已跑通。 |
| **严格验收** | `check_epitope_realization.py`（真实 MHC-I CV、MHC-II、免疫原性、真实结构）四 run 可按命令批量通过。 |
| **Positive Control** | 四 run 的 `POSITIVE_CONTROL.md` 均可满足「KRAS / `VVGADGVGK`」叙事；`R_public_001` 通过 `--ensure_positive_control_peptides` 保入并已重跑下游。 |
| **多维 MD 投递** | 工程已支持 `run_pandora_structure.py --top_k` + `prepare_simhub_delivery.py` 多子目录 `<base>_md_rXX`；**若磁盘上仅有一条 PANDORA PDB**，脚本会**提示并退回单目录**，不报错。 |
| **非阻塞项** | SimHub 回传（任务书不要求）；临床终版 HLA-II（Backlog）；可选第二套工具（Saluki 等）。 |

---

## 二、任务书章节对照（简表）

| 章节 | 结论 | 备注 |
|------|------|------|
| 输入 BioDriver | ✅ | `deliveries/<run>/to_immunogen/`；II 类终版分型见 **§五 Backlog** |
| 2.1 表位与免原 | 基本 ✅ | MHCflurry + NetMHCpan CV + NetMHCIIpan + DeepImmuno/PRIME/Repitope 真实口径；BigMHC 默认 off；**IEDB 不自动调用** |
| 2.2 排序与过滤 | ✅ | `rank_score`、WT 相似度 `select_top_peptides` |
| 2.3–2.4 多价与 mRNA | ✅ | UTR/ORF/polyA、LinearDesign；**COOL/CodonOpt 未接** |
| 2.5 二级结构 | 主线 ✅ | RNAfold/RNAeval/RNAplfold；**LinearFold/Saluki/RNAsnp 未接** |
| 2.6 LNP | 可选未深入 | 与任务书一致 |
| 3.1–3.4 报告与 SimHub 包 | ✅ | `REPORT` / 自证 / `dossier_context`；**不要求**回传写进验收 |
| 4 Top 3–5 MD case | ✅ 能力已具备 | 依赖本机生成 **K 条** PANDORA PDB；否则单目录兜底 |
| 5 独立 GitHub | 组织定 | 本仓库副本不代为打勾 |
| 6 验收与阳性对照 | ✅ | 四 run；保入参数见 **§四** |
| 7 风险声明 | 持续 | II 类局限、AI 免原性、非湿实验结论 |
| 8–9 边界与 mRNA-CAR | 不混线 / 预留 | 无联合硬 KPI |

---

## 三、已落地工程能力（不必再当「缺口」写进周会）

1. **Positive Control 保入**：`select_top_peptides.py --ensure_positive_control_peptides`；`run_all.py` 透传同名参数（**替换**当期 Top‑N 中 `rank_score` **最低**一条，**不增** `top_n`）。  
2. **多 SimHub 子 case**：`prepare_simhub_delivery.py` + `run_pandora_structure.py --top_k`；`results/<run>/meta.json` 含 `simhub_leaf_case_ids` / `simhub_delivery_dirs` 等。  
3. **`run_all.py`**：`--prepare_pandora`、`--top_k_md`、`--ensure_positive_control_peptides`。  
4. **`prepare_mhc_chain_sequences.py`**：扫描 `to_simhub/*/` 下全部 `selected_for_md.csv`。

---

## 四、可选后续（按优先级）

| 优先级 | 内容 |
|--------|------|
| **P1（运营）** | 对目标 run 跑满 `--top_k_md 5` **且** `run_pandora_structure --top_k 5`，使 `to_simhub` 真正出现 5 个子目录（算力与 Lead 配额）。 |
| **P1（材料）** | `immunogen_executive_report.html` / 老板版 md：**全量重跑或改肽后**与 `REPORT.md`、`qc_metrics.json` **手工或脚本对齐**；`mrna_candidate_qc/` 在改 `selected_peptides` 后视需要**批量重跑**。 |
| **P2** | 多 mRNA 后缀选型选定后，与默认 `mrna_vaccine.fasta` / SimHub 路径**脚本化对齐**（或扩展 `prepare_simhub_delivery`）。 |
| **P2** | 根目录 `REPORT.md`/`qc_metrics.json` 与 `mrna_design.json` **单调来源**：约定「以 design JSON 为准」或在 `run_qc` 后**强制回写**运行概况。 |
| **P3** | BigMHC 打开、第二套密码子工具、Saluki/LinearFold 等 **PoC**。 |

---

## 五、Backlog（不占当周优先级）

- **BioDriver 患者临床终版 HLA-II**：`hla_typing.json` 到位后 → 替换 `R002`/`R003`/`R_public_001` 等演示 II 类 → 全链路重跑 → 刷新报告口径。  
 **此前不阻塞**阶段验收。

---

## 六、验收与健康检查命令

```bash
# 表位与真实后端（四 run）
for r in R001 R002 R003 R_public_001; do
  python scripts/check_epitope_realization.py --run_id "$r" \
    --require_mhc1_cv_real \
    --require_mhc2_real \
    --require_real_immunogenicity \
    --require_real_structure
done
```

```bash
# SimHub 证据目录状态（可选；任务书不要求回传闭环）
python scripts/check_simhub_evidence.py --run_id R001
```

```bash
# 粗结构 / 代理残留扫描（示例）
rg "proxy|fallback_profile|coarse_initial_complex|na_for_coarse|structure_backend\"\\s*:\\s*\"coarse\"" \
  results/{R001,R002,R003,R_public_001} \
  deliveries/{R001,R002,R003,R_public_001}/to_simhub
```

---

## 七、模块边界（防串线）

- 本模块：**个性化 neoantigen → 多价 mRNA + peptide‑MHC SimHub 包**；**不做**小分子/DrugReflector 主线。  
- **mRNA‑CAR**：任务书预留；CellTherapy 与 Lead 点名后再开子项。

---

## 八、提交与仓库卫生

- 勿提交 `.autodl/`、下载缓存、未授权再分发的 `external_refs` 大二进制；`.ipynb_checkpoints/` 已列入 `.gitignore`。  
- 提交前：`git diff --cached`。

---

*以上为 2026-05-11 全文重写后的唯一基线。*
