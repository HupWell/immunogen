# ImmunoGen（mRNA 疫苗 / 免疫基因）

**负责人：** 梁心恬  
**模块定位：** 基于病人特异性新抗原（Neoantigen）设计多价 mRNA 疫苗。  
**一句话目标：** 从 BioDriver 的候选新抗原 + HLA 分型，输出候选 mRNA 全序列（含 UTR、密码子优化、多价串联设计），并按要求准备 Simulation Hub 交付包。

---

## 0. 一分钟导读

本模块用于病人个性化 mRNA 肿瘤疫苗设计：读取 BioDriver 的候选突变肽与 HLA 分型，综合 MHC 结合、免疫原性、VAF、与 WT 差异进行排序，选择 Top 10-20 肽段串联成多价 ORF，拼接 UTR/polyA 并生成完整 mRNA 序列，同时准备 SimHub 需要的 peptide-MHC 交付包。

第一阶段 smoke 目标：跑通 1 个病人 -> Top-10 neoantigen -> 1 条多价 mRNA -> 至少 1 个 peptide-MHC 交付包。  
交付节奏：先完成 `POSITIVE_CONTROL.md` / `SELF_CHECK.md` / `REPORT.md`，再提交 SimHub。

---

## 1. 输入输出契约

### 1.1 上游输入（BioDriver）

`deliveries/<run_id>/to_immunogen/`：

- `neoantigen_candidates.csv`（`mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf`）
- `hla_typing.json`（I 类必备：`HLA-A/B/C`；II 类可选见 `docs/hla_typing.md`）
- `meta.json`

### 1.2 本模块输出（结果产物）

`results/<run_id>/`：

- `peptide_mhc_ranking.csv`
- `selected_peptides.csv`
- `mrna_vaccine.fasta`
- `mrna_design.json`
- `qc_metrics.json`
- `figures/binding_affinity_heatmap.png`
- `figures/mrna_secondary_structure.png`
- `REPORT.md`
- `POSITIVE_CONTROL.md`
- `SELF_CHECK.md`
- `FEASIBILITY.md`
- `simhub_evidence/<case_id>/`（接收 SimHub 回传证据；含 `evidence_status.json` 与 `SIMHUB_EVIDENCE.md`）
- `meta.json`（本 run 摘要，与 SimHub 交付 `meta.json` 分工独立）

### 1.3 下游交付（Simulation Hub）

`deliveries/<run_id>/to_simhub/<case_id>/`：

- `complex.pdb`（Chain M/B/P）
- `hla_allele.txt`（可选）
- `meta.json`（`molecule_type: "peptide_mhc"`）
- `selected_for_md.csv`
- `dossier_context.json`（peptide + mRNA 最终 dossier manifest）

注意：peptide-MHC 分支禁止使用 SDF。

SimHub 回传证据检查：

```bash
python scripts/check_simhub_evidence.py --run_id R001
```

状态枚举：`not_returned`、`returned_unvalidated`、`validation_passed`、`validation_failed`。

---

## 2. 当前流程能力（已实现）

1. `validate_input.py`：输入契约检查  
2. `run_immunogenicity_adapters.py`：免疫原性适配器（auto/real_tsv/real_cmd/proxy）  
3. `predict_mhc_ranking.py`：MHC-I 主预测 + MHC-II 层 + 综合排序  
4. `select_top_peptides.py`：Top-N + WT 相似性过滤  
5. `build_multivalent_mrna.py`：多价串联 + UTR/polyA + 密码子模式  
6. `run_qc_and_report.py`：结构图、QC、报告  
7. `prepare_self_certification.py`：自证最小包  
8. `prepare_simhub_delivery.py`：SimHub 交付封装  
9. `validate_feasibility.py`：可行性验证

一键运行：

```bash
python scripts/run_all.py --run_id R001
```

真实来源常态化（推荐）：

```bash
# MHC-I：先准备 NetMHCpan real_tsv（可由真实工具输出，或由 openvax 报告桥接）
python scripts/openvax_bridge.py --run_id R001 --input "/path/to/vaccine-peptide-report.txt"

# MHC-II：支持 real_tsv（无需每次都调用 netMHCIIpan 二进制）
# 文件路径：results/R001/tool_outputs/raw/mhc2_netmhciipan.tsv
# 必需列：mut_peptide + (mhc2_el_rank|el_rank|pct_rank_el|rank)

# 运行并强制真实后端
python scripts/run_all.py --run_id R001 --target mhc_ranking --backend_mhc1_netmhcpan real_tsv --backend_mhc1_bigmhc off --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2

# 验收
python scripts/check_epitope_realization.py --run_id R001 --require_mhc1_cv_real --require_mhc2_real
```

MHC-I **`real_cmd`（本机 netMHCpan 自动跑，不写死 TSV）**：

```bash
# 依赖：已安装 NetMHCpan 4.x，并可从当前环境调用（或设 NETMHCPAN_BIN / NETMHCPAN_HOME）
export MHC1_NETMHCPAN_CMD='python tools/netmhcpan_class1_runner.py --input {input_tsv} --output {output_tsv}'

python scripts/run_all.py --run_id R001 --target mhc_ranking \
  --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off \
  --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2 \
  --backend_prime real_tsv --backend_repitope real_tsv \
  --require_real_immunogenicity_prime --require_real_immunogenicity_repitope
```

说明：若未提供 `results/<run_id>/tool_outputs/raw/prime.tsv` 或 `repitope.tsv`（或未配置对应 `real_cmd`），上述两个 `require_real_immunogenicity_*` 会阻止流程回退 proxy，并直接报错，确保 PRIME/Repitope 为真实来源。

说明：`predict_mhc_ranking` 会生成临时 `mhc1_netmhcpan_input.tsv` 并调用上述命令；包装器见 `tools/netmhcpan_class1_runner.py`。包装器**默认**为 NetMHCpan 4.2+ 追加 **`-BA`** 以得到 `Aff(nM)` 亲合力列；环境变量见 `NETMHCPAN_LEN`、`NETMHCPAN_EXTRA`、`NETMHCPAN_TIMEOUT`、`NETMHCPAN_BA`（设 `0` 可禁用 `-BA`，但表头常无 nM 会导致解析失败）。快速加载： `source scripts/env_netmhcpan.sh`（请按本机修改其中 `NETMHCPAN_HOME`）。

MHC-I **BigMHC `real_cmd`（KarchinLab [BigMHC](https://github.com/KarchinLab/bigmhc)，可选增强，默认不启用）**：

> 说明：BigMHC 不是流程必需项。当前推荐仅启用 NetMHCpan 作为 MHC-I 真实后端，BigMHC 可在网络与环境稳定后再补充为第二路交叉验证。

推荐（默认关闭 BigMHC，不阻塞主流程）：

```bash
python scripts/run_all.py --run_id R001 --target mhc_ranking \
  --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off \
  --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2
```

可选增强（启用 BigMHC 第二路）：

```bash
# 1) 克隆官方仓库（体积大，约数 GB），按官方 README 安装 PyTorch 等依赖
# 2) 设仓库根目录（其下含 src/predict.py），或显式指定 predict.py
export BIGMHC_HOME="/path/to/bigmhc"
export MHC1_BIGMHC_CMD='python tools/bigmhc_runner.py --input {input_tsv} --output {output_tsv}'

python scripts/run_all.py --run_id R001 --target mhc_ranking \
  --backend_mhc1_netmhcpan real_tsv --backend_mhc1_bigmhc real_cmd \
  --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2
```

包装器见 `tools/bigmhc_runner.py`：内部调用 `predict.py -m=el`（呈递 EL 模型），把 `BigMHC_EL` 列映射为 `mhc1_cv_bigmhc_score`。可选：`BIGMHC_MODEL`（默认 `el`）、`BIGMHC_DEVICE`（默认 `cpu`）、`BIGMHC_PREDICT_PY`（覆盖 `predict.py` 路径）。

多实例（`R002` / `R003` / `R_public_001` 等）在跑 `--backend_mhc1_netmhcpan real_tsv` 前，都需在**对应**目录放置两份表：

- `results/<run_id>/tool_outputs/raw/mhc1_netmhcpan.tsv`（I 类交叉验证；行数应覆盖本 run 的候选肽 × 病人 I 类 allele，至少肽级有值）
- `results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv`（II 类；`hla_typing.json` 里暂无 II 类时可用参考 `DRB1*15:01` 等**演示**用途，终版需换患者真实分型结果）

预检 + 可选自动跑数与验收：

```bash
python scripts/bootstrap_real_backends.py --run_id R002 --with_mhc2
python scripts/bootstrap_real_backends.py --run_id R002 --with_mhc2 --run --validate
```

I 类表亦可从外部分析经 `scripts/openvax_bridge.py` 等生成；最小列模板见 `data/examples/mhc1_cv_templates/mhc1_netmhcpan.tsv`。

---

## 3. 建议执行步骤（任务书对齐）

1. 基线验证：挑 1 个已知强免疫原性 neoantigen（如 KRAS G12D）跑通全流程  
2. 真实病例：对 Top 30-50 候选做 MHC + 免疫原性评分  
3. 排序筛选：按 `rank_score` 取 Top 10-20  
4. 多价设计：ORF 串联 + 密码子优化 + UTR 拼接 + 二级结构检查  
5. SimHub：选 Top 3-5 peptide-MHC 做 MD 稳定性验证

---

## 4. 验收标准

- 可复现：新环境按 README 可跑通
- 先自证再交付：必须有 `POSITIVE_CONTROL.md`、`SELF_CHECK.md`、`REPORT.md`
- 完整性：排名 + 多价设计 + mRNA 三件套齐全
- 可解释性：Top 候选含 HLA / 免疫原性 / VAF 证据
- Positive Control：至少 1 个已知 neoantigen 被识别为高分

---

## 5. 当前边界与风险

- MHC-II 预测精度低于 MHC-I，需要在报告注明
- 免疫原性 AI 预测不确定性较大，不替代湿实验
- 当前正式交付默认要求 PANDORA/AFM 等真实结构；`coarse` 仅保留为显式调试选项，不能作为真实交付。
- 第一阶段不做 LNP 处方设计，仅保留综述说明（`docs/lnp_notes.md`）

---

## 6. 相关文档

- `RELEASE_NOTES.md`
- `docs/hla_typing.md`
- `docs/netmhciipan_setup.md`
- `docs/structure_backend_selection.md`
- `docs/WEEKLY_REPORT_2026-04-24.md`
- `docs/TODO.md`
- `FINAL_CHECKLIST.md`
