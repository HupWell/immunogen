#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
将 external_refs/neoantigen-vaccine-pipeline（或同类报告）结果桥接为本项目 MHC-I real_tsv 输入。
"""

import argparse
import os
import pandas as pd


def _pick_col(df: pd.DataFrame, names):
    for n in names:
        if n in df.columns:
            return n
    return None


def main(run_id: str, input_path: str, peptide_col: str, affinity_col: str, allele_col: str):
    if not os.path.exists(input_path):
        raise FileNotFoundError(input_path)
    sep = "\t" if input_path.lower().endswith((".tsv", ".txt")) else ","
    src = pd.read_csv(input_path, sep=sep)
    if src.empty:
        raise ValueError("输入文件为空。")

    p_col = peptide_col or _pick_col(src, ["mut_peptide", "peptide", "Mutant peptide sequence", "vaccine_peptide"])
    a_col = affinity_col or _pick_col(src, ["mhc1_cv_netmhcpan_nM", "affinity_nM", "ic50", "IC50", "predicted_affinity"])
    h_col = allele_col or _pick_col(src, ["hla_allele", "allele", "MHC Allele"])
    if not p_col:
        raise ValueError("未识别到肽段列，请通过 --peptide_col 指定。")
    if not a_col:
        raise ValueError("未识别到亲和力列，请通过 --affinity_col 指定。")

    out = pd.DataFrame()
    out["mut_peptide"] = src[p_col].astype(str).str.strip().str.upper()
    out["mhc1_cv_netmhcpan_nM"] = pd.to_numeric(src[a_col], errors="coerce")
    out["hla_allele"] = src[h_col].astype(str).str.strip() if h_col else ""
    out = out[(out["mut_peptide"] != "") & out["mhc1_cv_netmhcpan_nM"].notna()].copy()
    out = out.drop_duplicates(subset=["mut_peptide", "hla_allele"], keep="first")
    if out.empty:
        raise ValueError("转换后为空，请检查输入列和内容。")

    raw_dir = os.path.join("results", run_id, "tool_outputs", "raw")
    os.makedirs(raw_dir, exist_ok=True)
    out_path = os.path.join(raw_dir, "mhc1_netmhcpan.tsv")
    out.to_csv(out_path, sep="\t", index=False, encoding="utf-8")
    print(f"已写入: {out_path}")
    print(f"映射列: peptide={p_col}, affinity={a_col}, allele={h_col or '(空)'}")
    print(f"行数: {len(out)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--input", required=True, help="openvax/vaxrank 报告文件路径")
    parser.add_argument("--peptide_col", default="", help="手动指定肽段列名")
    parser.add_argument("--affinity_col", default="", help="手动指定亲和力列名")
    parser.add_argument("--allele_col", default="", help="手动指定等位基因列名")
    args = parser.parse_args()
    main(args.run_id, args.input, args.peptide_col.strip(), args.affinity_col.strip(), args.allele_col.strip())
