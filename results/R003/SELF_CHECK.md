# Self Check（R003）

## 1. 双工具 / 双证据链（当前实现）

| 环节 | 工具或实现 | 用途 |
|------|----------------|------|
| MHC-I 亲和力 | **MHCflurry**（`Class1AffinityPredictor`） | 肽–HLA-I 结合强度排序 |
| MHC-II | **NetMHCIIpan**（子进程，见 `netmhciipan_runner.py`；列 `mhc2_el_rank` / `mhc2_ba_nm`） | 真实 II 类预测时 EL %Rank 等；仍须声明工具版本与等位基因写法 |
| mRNA 二级结构 / MFE | **ViennaRNA**（优先 `RNAfold` 命令，否则 Python `RNA` 绑定） | MFE 与 dot-bracket，用于质控图与 `mrna_design.json` |

- 本轮 MHC-II 行统计 netmhciipan: 18。（`peptide_mhc_ranking.csv` 中 `mhc2_backend` 列应为 netmhciipan；proxy 不再作为默认验收口径）
- 免疫原性来源统计：immunogenicity_source_deepimmuno=deepimmuno_real_cmd:18；immunogenicity_source_prime=prime_real_cmd:18；immunogenicity_source_repitope=repitope_public_dataset_knn:12，repitope_public_dataset_exact:6

## 2. 关键阈值与过滤（当前默认）

| 项目 | 默认值 | 说明 |
|------|--------|------|
| Top-N 入选肽 | `select_top_peptides.py --top_n`（默认 10） | 进入多价串联 |
| 与 WT 最小差异比例 | `--min_dissimilarity`（默认 0.1） | 过滤与野生型过相似肽 |
| 综合打分权重 | `predict_mhc_ranking.py --w1..w4` | `rank_score` 加权 |

## 3. 本轮 RNA 折叠状态

- `rnafold_status`：**fallback_profile**
- `rnafold_mfe`（若已计算）：**NA**

## 4. 已知局限（须在汇报中声明）

1. **MHC-II**：当前验收要求 **NetMHCIIpan/real_tsv** 真实来源；缺失时流程应报错，不应静默回退 proxy。
2. **免疫原性（DeepImmuno / PRIME / Repitope）**：当前验收要求 real_tsv/real_cmd 真实来源；缺失时流程应报错。
3. **NetMHCpan / BigMHC**：NetMHCpan 作为默认 MHC-I 真实交叉验证；BigMHC 为可选增强，不作为代理补位。
4. **SimHub 初始结构**：当前交付要求 PANDORA / AFM 等真实结构文件；`coarse` 仅可作为显式调试选项，不能冒充真实交付。
5. **肽–MHC 分支禁止使用 SDF + 小分子电荷模型**；本仓库 SimHub 交付已按契约仅输出 `complex.pdb` + `meta.json` + 可选 `hla_allele.txt`。

## 4.1 结构后端策略（当前建议）

- 默认交付链路要求 `pandora`/`afm` 真实 PDB，并通过 `meta.json` 标记 `replaces_coarse=true`。
- 真实结构采用“**PANDORA 批量 + AFM 重点复核**”组合策略（见 `docs/structure_backend_selection.md`）。
- 若使用 AFM，建议独立任务执行并回填 PDB，避免阻塞主流水线。

## 5. 与任务书对齐

- 先完成本目录下 **POSITIVE_CONTROL.md**、**SELF_CHECK.md** 与 **REPORT.md**，再提交 Simulation Hub 交付包（见 `deliveries/<run_id>/to_simhub/`）。

## 6. IEDB 资源（人工复核入口）

- MHC-I 结合预测工具入口：https://tools.iedb.org/mhci/
- MHC-II 结合预测工具入口：https://tools.iedb.org/mhcii/
- IEDB 主站（文献与实验数据）：https://www.iedb.org/
- 说明：本仓库当前不强制自动调用 IEDB，建议将本流程 Top 肽段导出后做交叉复核。
