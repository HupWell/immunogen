# ImmunoGen 任务书自检与 TODO

> **负责人：** 梁心恬（与 `README.md` 一致；若组织另有指派可改署名为「模块负责人」）  
> **文档性质：** 对照任务书做**自检**，本文**即为当前唯一 TODO**，旧版条目全部作废。  
> **本文更新日期：** 2026-05-18（在 2026-05-11 基线上增补路径 B / 转录组 zip 重构）  
> **口径依据：** 仓库脚本、`results/`、`deliveries/` 与任务书当前版（**不要求** SimHub MD 回传闭环；**HLA-II 临床终版**上游暂未提供）。

---

## 一、执行摘要（给 Lead）

| 维度 | 状态 |
|------|------|
| **当前上游侧重** | **路径 B**：`data/表达定量表格.zip` + `data/样本信息表.zip` → 一队列一 run（`T_cohort_<id>`）→ 基因 log2FC 排名；**不交付 SimHub**（见 README §1.1A 流程图）。 |
| **路径 B 工程** | ✅ 同上；P0 批量见 `results/_summary/transcriptome_batch_20260518_round2.csv`。 |
| **路径 A 运营（Top-5 MD）** | ✅ `R_public_001` 5 条 rank PDB + 5 个 `*_md_r**` SimHub 包；演示 run `R001`–`R003` 各 2–3 条（受突变肽池限制）。 |
| **探索性分组** | ✅ 支持 `death:ch1` / `meta:ch1` 等探索性排序，`provenance` 标 `exploratory` + 中文 caution；队列 103091 等已 smoke。 |
| **四 run 主链路（路径 A）** | `R001` / `R002` / `R003` / `R_public_001`：排名→多价 mRNA→SimHub peptide_mhc 包已跑通；企业服务器环境已补齐（mhcflurry pan、ViennaRNA 等）。 |
| **SimHub 交付口径** | **仅路径 A** 产生 `to_simhub/` 五件套；路径 B 经 **`b_to_a_simhub` / `run_b_to_a_simhub.py`** 过滤突变肽后走路径 A，包结构一致。 |
| **严格验收（路径 A）** | `check_epitope_realization.py` 四 run 可按命令批量通过。 |
| **非阻塞项** | SimHub 回传；临床终版 HLA-II；跨队列 run_id 合并与 Ensembl 转录本↔基因映射表。 |

---

## 二、任务书章节对照（简表）

| 章节 | 结论 | 备注 |
|------|------|------|
| 输入 BioDriver / 转录组 | ✅ 路径 B 已接 zip | `to_transcriptome/`；路径 A 仍要 `to_immunogen/` |
| 2.1 表位与免原 | 基本 ✅ | 路径 A；BigMHC 默认 off |
| 2.2 排序与过滤 | ✅ | 路径 A `rank_score`；路径 B 为基因 log2FC |
| 2.3–2.4 多价与 mRNA | ✅ | 仅路径 A |
| 2.5 二级结构 | 主线 ✅ | RNAfold/RNAeval/RNAplfold |
| 3.1–3.4 报告与 SimHub 包 | ✅ | **仅路径 A** → `to_simhub/` |
| 4 Top 3–5 MD case | ✅ 能力已具备 | 路径 A + PANDORA |
| 6 验收与阳性对照 | ✅ | 四 run 路径 A |

---

## 三、已落地工程能力（不必再当「缺口」写进周会）

### 路径 A（既有）

1. **Positive Control 保入**：`--ensure_positive_control_peptides`（`R_public_001` 已用）。  
2. **多 SimHub 子 case**：`prepare_simhub_delivery.py` + `run_pandora_structure.py --top_k`。  
3. **`run_all.py`**：`--prepare_pandora`、`--top_k_md`、`--target full` 等。  
4. **`prepare_mhc_chain_sequences.py`**：扫描 `to_simhub/*/selected_for_md.csv`。

### 路径 B（2026-05-18 新增）

5. **转录组 zip 摄入**：`scripts/ingest_transcriptome_archives.py` → `data/transcriptome/manifest/cohorts.csv`。  
6. **单队列一键**：`scripts/run_transcriptome_cohort.py --cohort_id <id>`；`run_all.py --target transcriptome_prior --cohort_id <id>`。  
7. **批量**：`scripts/run_transcriptome_batch.py`（`--dry_run` 只检查分组）。  
8. **分组配置**：`data/transcriptome/cohort_groups.yaml`（manual / skip / 探索性 death）。  
9. **契约目录**：`deliveries/T_cohort_<id>/to_transcriptome/`（`tpm_matrix.tsv`、`sample_groups.json`、`meta.json`）。  
10. **README**：§1.1A 已同步三张 Mermaid 流程图（路径 B、路径 A+SimHub、B→A→SimHub）。

---

## 四、可选后续（按优先级）

| 优先级 | 内容 |
|--------|------|
| ~~**P0（路径 B 运营）**~~ | ✅ 已跑批并汇总；后续增量队列重跑 `ingest` + `run_transcriptome_batch.py`，并在 `cohort_groups.yaml` 补 manual/skip。 |
| ~~**P1（B→A→SimHub）**~~ | ✅ `filter_neoantigen_by_gene_rank.py`、`run_b_to_a_simhub.py`、`run_all.py --target b_to_a_simhub`；见 README §1.1A B→A 命令。 |
| ~~**P1（运营）**~~ | ✅ `R_public_001` 满配 5×SimHub；`R001`–`R003` 受演示突变肽条数限制为 2–3。批量脚本 `run_path_a_simhub_top5_batch.py`；汇总 `results/_summary/path_a_simhub_top5_status.csv`。 |
| ~~**P1（材料）**~~ | ✅ 四 run 已 `sync_run_materials.py` 刷新 `REPORT.md` / `qc_metrics.json` / 自证；`R_public_001` 含 11 条 mRNA 候选 QC；老板版 `docs/ImmunoGen_管理层汇报_老板版.md` 已机器对齐。 |
| **P2** | 多 mRNA 后缀与默认 `mrna_vaccine.fasta` / SimHub 路径脚本化对齐。 |
| **P3** | BigMHC、Saluki/LinearFold 等 PoC。 |

---

## 五、Backlog（不占当周优先级）

- **BioDriver 患者临床终版 HLA-II**：到位后替换演示 II 类并重跑路径 A 四 run。  
- **路径 B 队列扩展**：新 zip 入库后重跑 `ingest`；`cohort_groups.yaml` 逐队列人工校对肿瘤/正常 vs 探索性口径。  
- **单细胞 / 10X**：独立分析，不替代路径 A 突变+HLA 输入。

---

## 六、验收与健康检查命令

### 路径 B（转录组）

```bash
# 生成 / 刷新 manifest
python scripts/ingest_transcriptome_archives.py

# 单队列（示例：肿瘤 vs 正常）
python scripts/run_transcriptome_cohort.py --cohort_id 39004

# 探索性 death 分组（示例）
python scripts/run_transcriptome_cohort.py --cohort_id 103091

# 批量（先 dry_run 看分组）
python scripts/run_transcriptome_batch.py --dry_run
python scripts/run_transcriptome_batch.py

# 一键入口
python scripts/run_all.py --run_id T_cohort_39004 --target transcriptome_prior --cohort_id 39004

# B→A→SimHub（同一 run_id 需同时具备 to_immunogen + transcriptome_prior）
python scripts/run_b_to_a_simhub.py --run_id <run_id> --top_n_genes 500
python scripts/filter_neoantigen_by_gene_rank.py --run_id <run_id> --top_n_genes 500 --dry_run

# 材料包与老板版对齐（四 run）
python scripts/sync_run_materials.py
python scripts/sync_run_materials.py --runs R_public_001 --mrna_qc
```

### 路径 A（表位与真实后端）

```bash
# Top-5 MD + SimHub 四 run 批量（需 PANDORA + biopython==1.79，见 README）
python scripts/run_path_a_simhub_top5_batch.py
python scripts/run_path_a_simhub_top5_batch.py --skip_pandora   # 仅刷新 SimHub 包

for r in R001 R002 R003 R_public_001; do
  python scripts/check_epitope_realization.py --run_id "$r" \
    --require_mhc1_cv_real \
    --require_mhc2_real \
    --require_real_immunogenicity \
    --require_real_structure
done
```

```bash
python scripts/check_simhub_evidence.py --run_id R001
```

```bash
rg "proxy|fallback_profile|coarse_initial_complex" \
  results/{R001,R002,R003,R_public_001} \
  deliveries/{R001,R002,R003,R_public_001}/to_simhub
```

---

## 七、模块边界（防串线）

- 本模块：**个性化 neoantigen → 多价 mRNA + peptide‑MHC SimHub 包**（路径 A）。  
- **路径 B**：bulk 转录组 → **基因优先级**；`T_cohort_*` run **无** `to_simhub/`。与 SimHub 相同交付必须走路径 A（或未来 B→A 衔接）。详见 `README.md` §1.1A 流程图。  
- **mRNA‑CAR**：任务书预留。

---

## 八、提交与仓库卫生

- 勿提交 `.autodl/`、下载缓存、未授权再分发的 `external_refs` 大二进制；`data/transcriptome/raw/` 体积大时可仅提交 manifest + yaml。  
- 提交前：`git diff --cached`。

---

*基线 2026-05-11；路径 B / README 流程图 / 转录组 zip 于 2026-05-18 更新。*
