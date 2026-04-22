# -*- coding: utf-8 -*-
"""
功能：生成任务书要求的「自证最小包」中的 POSITIVE_CONTROL.md 与 SELF_CHECK.md。

输入（读取已有结果）：
- results/<run_id>/selected_peptides.csv
- results/<run_id>/peptide_mhc_ranking.csv
- results/<run_id>/qc_metrics.json（可选）

输出：
- results/<run_id>/POSITIVE_CONTROL.md
- results/<run_id>/SELF_CHECK.md

说明：REPORT.md 由 run_qc_and_report.py 生成；本脚本补齐另外两份自证文档。
"""
import os
import json
import argparse
from typing import Optional

import pandas as pd


def load_json(path: str):
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_positive_control(run_id: str, out_dir: str, selected: pd.DataFrame):
    """记录已知强免疫原 neoantigen 在本流程中的复现检查。"""
    kras_hit = False
    if not selected.empty and "mut_peptide" in selected.columns:
        peptides = selected["mut_peptide"].astype(str).str.upper().tolist()
        kras_hit = any("VVGADGVGK" in p or "GADGVGK" in p for p in peptides)

    path = os.path.join(out_dir, "POSITIVE_CONTROL.md")
    content = f"""# Positive Control（{run_id}）

## 1. 对照定义

- 经典 neoantigen：**KRAS G12D** 来源 9-mer 肽段 **`VVGADGVGK`**（文献常用阳性对照之一）。
- 本流程检查：该肽是否出现在 **selected_peptides**（进入多价疫苗的 Top-N）中。

## 2. 本轮结果

- **KRAS G12D 肽段是否入选 Top 候选**：`{"是" if kras_hit else "否（请检查输入是否含该肽或提高 top_n）"}`

## 3. 如何人工复核

1. 打开 `selected_peptides.csv`，搜索 `VVGADGVGK`。
2. 若未命中，可检查 `neoantigen_candidates.csv` 是否包含该突变肽，或增大 `--top_n` 后重跑 `run_all.py`。

## 4. 说明

- 本对照用于**流程与排序方向**验证，不替代湿实验（ELISpot / 四聚体等）。
- 公开基准集建议使用 `R_public_001`（见 `data/public/` 与 `RELEASE_NOTES.md`）。
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"完成: {path}")


def write_self_check(
    run_id: str, out_dir: str, qc: dict, ranking: Optional[pd.DataFrame] = None,
):
    """双工具评分、阈值与已知局限的自证说明。"""
    mfe = qc.get("rnafold_mfe")
    rf_status = qc.get("rnafold_status", "unknown")
    mhc2_note = "（`peptide_mhc_ranking.csv` 中 `mhc2_backend` 列：proxy=代理分，netmhciipan=实跑）"
    mhc2_row = "MHC-II | **代理分**（未接 NetMHCIIpan 或无 II 类分型时） | 辅助排序的 II 类相关量纲分"
    if ranking is not None and "mhc2_backend" in ranking.columns:
        vc = ranking["mhc2_backend"].fillna("").astype(str).value_counts()
        ssum = "，".join(f"{k}: {v}" for k, v in vc.items())
        mhc2_note = f"本轮 MHC-II 行统计 {ssum}。{mhc2_note}"
        if (ranking["mhc2_backend"] == "netmhciipan").any():
            mhc2_row = "MHC-II | **NetMHCIIpan**（子进程，见 `netmhciipan_runner.py`；列 `mhc2_el_rank` / `mhc2_ba_nm`） | 真实 II 类预测时 EL %Rank 等；仍须声明工具版本与等位基因写法"
    path = os.path.join(out_dir, "SELF_CHECK.md")
    content = f"""# Self Check（{run_id}）

## 1. 双工具 / 双证据链（当前实现）

| 环节 | 工具或实现 | 用途 |
|------|----------------|------|
| MHC-I 亲和力 | **MHCflurry**（`Class1AffinityPredictor`） | 肽–HLA-I 结合强度排序 |
| {mhc2_row} |
| mRNA 二级结构 / MFE | **ViennaRNA**（优先 `RNAfold` 命令，否则 Python `RNA` 绑定） | MFE 与 dot-bracket，用于质控图与 `mrna_design.json` |

- {mhc2_note}

## 2. 关键阈值与过滤（当前默认）

| 项目 | 默认值 | 说明 |
|------|--------|------|
| Top-N 入选肽 | `select_top_peptides.py --top_n`（默认 10） | 进入多价串联 |
| 与 WT 最小差异比例 | `--min_dissimilarity`（默认 0.1） | 过滤与野生型过相似肽 |
| 综合打分权重 | `predict_mhc_ranking.py --w1..w4` | `rank_score` 加权 |

## 3. 本轮 RNA 折叠状态

- `rnafold_status`：**{rf_status}**
- `rnafold_mfe`（若已计算）：**{mfe if mfe is not None else "NA"}**

## 4. 已知局限（须在汇报中声明）

1. **MHC-II**：可通过 `--mhc2_backend` 与 `NETMHCIIPAN_BIN` 等环境变量接 **NetMHCIIpan**；未配置或无可转换的 II 类分型时行内为 **proxy**。
2. **免疫原性（DeepImmuno / PRIME / Repitope）**：当前为可复现代理分；真实模型需单独部署与校准。
3. **NetMHCpan / BigMHC**：未作为默认第二路 MHC-I；可在 `SELF_CHECK` 后续版本中补充交叉验证表。
4. **SimHub 初始结构**：当前 `complex.pdb` 为**粗粒度**多链占位，**必须**替换为 AlphaFold-Multimer / PANDORA 等生成的真实复合物。
5. **肽–MHC 分支禁止使用 SDF + 小分子电荷模型**；本仓库 SimHub 交付已按契约仅输出 `complex.pdb` + `meta.json` + 可选 `hla_allele.txt`。

## 5. 与任务书对齐

- 先完成本目录下 **POSITIVE_CONTROL.md**、**SELF_CHECK.md** 与 **REPORT.md**，再提交 Simulation Hub 交付包（见 `deliveries/<run_id>/to_simhub/`）。
"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    print(f"完成: {path}")


def main(run_id: str):
    out_dir = os.path.join("results", run_id)
    os.makedirs(out_dir, exist_ok=True)
    selected_path = os.path.join(out_dir, "selected_peptides.csv")
    if not os.path.exists(selected_path):
        raise FileNotFoundError(f"未找到: {selected_path}，请先运行至筛选步骤。")

    selected = pd.read_csv(selected_path)
    qc = load_json(os.path.join(out_dir, "qc_metrics.json"))
    rank_path = os.path.join(out_dir, "peptide_mhc_ranking.csv")
    ranking = pd.read_csv(rank_path) if os.path.exists(rank_path) else None
    write_positive_control(run_id, out_dir, selected)
    write_self_check(run_id, out_dir, qc, ranking=ranking)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
