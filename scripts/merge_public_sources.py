# -*- coding: utf-8 -*-
"""
功能：合并 pVACtools 示例候选与文献案例，生成统一公开验证集 CSV。

输入：
- data/public/pvac_example/pvac_candidates.csv
- data/public/literature/literature_cases.csv

输出：
- data/public/combined_public.csv
"""
import os
import argparse
import pandas as pd


TARGET_COLUMNS = ["gene_change", "MT_pep", "WT_pep", "vaf", "transcript", "source"]

ALIASES = {
    "gene_change": ["gene_change", "mutation", "variant", "mutation_name", "ID", "id"],
    "MT_pep": ["MT_pep", "mut_peptide", "mutant_peptide", "peptide", "tumor_peptide", "Best Peptide", "best peptide"],
    "WT_pep": ["WT_pep", "wt_peptide", "wildtype_peptide", "normal_peptide"],
    "vaf": ["vaf", "variant_vaf", "tumor_vaf", "af"],
    "transcript": ["transcript", "transcript_id", "enst_id", "tx_id"],
    "source": ["source"],
}


def detect_col(df: pd.DataFrame, names):
    lower_map = {c.lower(): c for c in df.columns}
    for n in names:
        raw = lower_map.get(n.lower())
        if raw:
            return raw
    return None


def normalize_df(df: pd.DataFrame, default_source: str) -> pd.DataFrame:
    out = pd.DataFrame()
    for col in TARGET_COLUMNS:
        raw = detect_col(df, ALIASES[col])
        if raw is None:
            if col == "WT_pep":
                out[col] = ""
            elif col == "vaf":
                out[col] = 0.2
            elif col == "transcript":
                out[col] = "NA"
            elif col == "source":
                out[col] = default_source
            else:
                out[col] = "UNKNOWN"
        else:
            out[col] = df[raw]

    out["MT_pep"] = out["MT_pep"].astype(str).str.strip().str.upper()
    out["WT_pep"] = out["WT_pep"].astype(str).str.strip().str.upper()
    out["gene_change"] = out["gene_change"].astype(str).str.strip()
    out["transcript"] = out["transcript"].astype(str).str.strip()
    out["source"] = out["source"].astype(str).str.strip()
    out["vaf"] = pd.to_numeric(out["vaf"], errors="coerce").fillna(0.2)

    missing_wt = out["WT_pep"] == ""
    out.loc[missing_wt, "WT_pep"] = out.loc[missing_wt, "MT_pep"].apply(
        lambda p: (p[:-1] + "A") if len(p) > 0 else "A"
    )
    out = out[out["MT_pep"].str.len() >= 8].copy()
    return out


def main(pvac_csv: str, literature_csv: str, output_csv: str):
    if not os.path.exists(literature_csv):
        raise FileNotFoundError(f"未找到文献案例文件: {literature_csv}")

    dfs = []
    if os.path.exists(pvac_csv):
        pvac_df = pd.read_csv(pvac_csv, sep=None, engine="python")
        dfs.append(normalize_df(pvac_df, default_source="pvac_example"))
    else:
        print(f"警告: 未找到 pVAC 示例文件，先仅使用文献案例: {pvac_csv}")

    lit_df = pd.read_csv(literature_csv, sep=None, engine="python")
    dfs.append(normalize_df(lit_df, default_source="literature"))

    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.drop_duplicates(subset=["gene_change", "MT_pep"]).reset_index(drop=True)
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    merged.to_csv(output_csv, index=False, encoding="utf-8-sig")
    print(f"完成: {output_csv}")
    print(f"合并后条目数: {len(merged)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pvac_csv", default="data/public/pvac_example/pvac_candidates.csv")
    parser.add_argument("--literature_csv", default="data/public/literature/literature_cases.csv")
    parser.add_argument("--output_csv", default="data/public/combined_public.csv")
    args = parser.parse_args()
    main(args.pvac_csv, args.literature_csv, args.output_csv)
