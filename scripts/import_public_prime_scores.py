# -*- coding: utf-8 -*-
"""
将 PRIME 公开数据映射为本项目 raw/prime.tsv（mut_peptide + score）。

当前默认来源：
1) data/public/immunogenicity/prime/out_compare.txt（PRIME 官方仓库测试输出）
2) data/public/immunogenicity/prime/supplementary/*.csv（可选，用户后续手工放置）
"""
import argparse
import os
from typing import Dict, List

import pandas as pd

from immunogenicity_adapters import ensure_tool_output_dir, load_unique_peptides


PRIME_DIR = os.path.join("data", "public", "immunogenicity", "prime")
OUT_COMPARE = os.path.join(PRIME_DIR, "out_compare.txt")
SUPP_DIR = os.path.join(PRIME_DIR, "supplementary")


def _collect_from_out_compare(score_map: Dict[str, List[float]]):
    if not os.path.exists(OUT_COMPARE):
        return
    # 跳过注释行，按空白分列；列2为 Score_bestAllele
    df = pd.read_csv(
        OUT_COMPARE,
        sep=r"\s+",
        comment="#",
        engine="python",
    )
    if "Peptide" not in df.columns or "Score_bestAllele" not in df.columns:
        return
    tmp = pd.DataFrame(
        {
            "mut_peptide": df["Peptide"].astype(str).str.strip().str.upper(),
            "score": pd.to_numeric(df["Score_bestAllele"], errors="coerce"),
        }
    ).dropna(subset=["mut_peptide", "score"])
    for _, r in tmp.iterrows():
        score_map.setdefault(r["mut_peptide"], []).append(float(r["score"]))


def _collect_from_supp_csv(score_map: Dict[str, List[float]]):
    """
    可选补充：读取 supplementary/*.csv（用户可自行放置 PRIME 补充表）。
    自动识别列：
    - peptide 列候选: peptide / mut_peptide
    - score 列候选: score / PRIME_score / immunogenicity / value / prediction
    """
    if not os.path.isdir(SUPP_DIR):
        return
    for name in os.listdir(SUPP_DIR):
        if not name.lower().endswith(".csv"):
            continue
        path = os.path.join(SUPP_DIR, name)
        try:
            df = pd.read_csv(path)
        except Exception:
            continue
        pep_col = None
        for c in ("peptide", "mut_peptide", "Peptide"):
            if c in df.columns:
                pep_col = c
                break
        score_col = None
        for c in ("score", "PRIME_score", "immunogenicity", "value", "prediction", "Score_bestAllele"):
            if c in df.columns:
                score_col = c
                break
        if pep_col is None or score_col is None:
            continue
        tmp = pd.DataFrame(
            {
                "mut_peptide": df[pep_col].astype(str).str.strip().str.upper(),
                "score": pd.to_numeric(df[score_col], errors="coerce"),
            }
        ).dropna(subset=["mut_peptide", "score"])
        for _, r in tmp.iterrows():
            score_map.setdefault(r["mut_peptide"], []).append(float(r["score"]))


def main(run_id: str):
    peptides = load_unique_peptides(run_id).tolist()
    score_map: Dict[str, List[float]] = {}
    _collect_from_out_compare(score_map)
    _collect_from_supp_csv(score_map)

    rows = []
    for pep in peptides:
        vals = score_map.get(pep, [])
        if vals:
            rows.append({"mut_peptide": pep, "score": round(sum(vals) / len(vals), 6)})

    out_dir = os.path.join(ensure_tool_output_dir(run_id), "raw")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, "prime.tsv")
    out = pd.DataFrame(rows)
    if out.empty:
        # 若历史存在空模板，删除以免 real_tsv 模式误读空表
        if os.path.exists(out_file):
            os.remove(out_file)
            print(f"无匹配，已删除旧文件: {out_file}")
        else:
            print(f"未找到可匹配 PRIME 公开分数，未写入: {out_file}")
        return
    out.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")
    print(f"匹配肽数: {len(out)} / {len(peptides)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
