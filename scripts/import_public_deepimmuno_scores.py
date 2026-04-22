# -*- coding: utf-8 -*-
"""
将公开 DeepImmuno 数据集映射为本项目 raw/deepimmuno.tsv（mut_peptide + score）。

用途：
- 对公开 run（如 R_public_001）快速生成“可直接被 real_tsv 读取”的 deepimmuno 输入，
  减少手工逐条填分数。

策略：
- 读取 data/public/immunogenicity/deepimmuno/ 下多个公开 CSV；
- 提取每条 peptide 的可用分数字段（优先 potential / immunogenic score / score）；
- 对同一 peptide 多来源分数取均值；
- 仅输出当前 run_id 中出现的 mut_peptide。
"""
import argparse
import os
from typing import Dict, List

import pandas as pd

from immunogenicity_adapters import ensure_tool_output_dir, load_unique_peptides


PUBLIC_DIR = os.path.join("data", "public", "immunogenicity", "deepimmuno")
PUBLIC_FILES = [
    "remove0123_sample100.csv",
    "dengue_test.csv",
    "ori_test_cells.csv",
    "sars_cov_2_result.csv",
]
SCORE_CANDIDATES = ["potential", "immunogenic score", "score", "ratio"]


def _collect_public_scores() -> Dict[str, List[float]]:
    score_map: Dict[str, List[float]] = {}
    for name in PUBLIC_FILES:
        path = os.path.join(PUBLIC_DIR, name)
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        if "peptide" not in df.columns:
            continue
        score_col = None
        for c in SCORE_CANDIDATES:
            if c in df.columns:
                score_col = c
                break
        if score_col is None:
            continue
        t = pd.DataFrame(
            {
                "mut_peptide": df["peptide"].astype(str).str.strip().str.upper(),
                "score": pd.to_numeric(df[score_col], errors="coerce"),
            }
        ).dropna(subset=["mut_peptide", "score"])
        for _, r in t.iterrows():
            pep = r["mut_peptide"]
            score_map.setdefault(pep, []).append(float(r["score"]))
    return score_map


def main(run_id: str):
    peptides = load_unique_peptides(run_id).tolist()
    score_map = _collect_public_scores()
    rows = []
    for pep in peptides:
        vals = score_map.get(pep, [])
        if not vals:
            continue
        rows.append({"mut_peptide": pep, "score": round(float(sum(vals) / len(vals)), 6)})
    out_dir = os.path.join(ensure_tool_output_dir(run_id), "raw")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "deepimmuno.tsv")
    out = pd.DataFrame(rows)
    if out.empty:
        print(f"未找到可匹配的公开 DeepImmuno 分数，未写入: {out_file}")
        return
    out.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")
    print(f"匹配肽数: {len(out)} / {len(peptides)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
