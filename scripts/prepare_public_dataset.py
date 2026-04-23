# -*- coding: utf-8 -*-
"""
功能：将公开数据整理为 ImmunoGen 输入格式。

输入：
- 任意 CSV（至少包含突变肽列）

输出：
- deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv
- deliveries/<run_id>/to_immunogen/hla_typing.json
- deliveries/<run_id>/to_immunogen/meta.json
"""
import os
import json
import argparse
import pandas as pd


REQUIRED_OUTPUT_COLUMNS = [
    "mutation",
    "mut_peptide",
    "wt_peptide",
    "transcript_id",
    "variant_vaf",
]

COMMON_COLUMN_ALIASES = {
    "mutation": ["mutation", "gene_change", "variant", "mutation_name", "mut"],
    "mut_peptide": ["mut_peptide", "mt_pep", "mutant_peptide", "peptide", "tumor_peptide"],
    "wt_peptide": ["wt_peptide", "wt_pep", "wildtype_peptide", "normal_peptide"],
    "transcript_id": ["transcript_id", "transcript", "enst_id", "tx_id"],
    "variant_vaf": ["variant_vaf", "vaf", "tumor_vaf", "af"],
}


def split_map_arg(mapping_text: str) -> dict:
    """
    解析列映射参数，例如：
    "mut_peptide=MT_pep,wt_peptide=WT_pep,mutation=gene_change"
    """
    result = {}
    mapping_text = (mapping_text or "").strip()
    if not mapping_text:
        return result
    for part in mapping_text.split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        result[k.strip()] = v.strip()
    return result


def normalize_hla_list(hla_text: str):
    """
    解析 HLA 参数，支持：
    "A*02:01,A*11:01;B*40:01,B*46:01;C*01:02,C*07:02"
    """
    groups = [x.strip() for x in hla_text.split(";")]
    def parse_group(i):
        if i >= len(groups) or not groups[i]:
            return []
        return [x.strip().replace("HLA-", "") for x in groups[i].split(",") if x.strip()]
    return {
        "HLA-A": parse_group(0),
        "HLA-B": parse_group(1),
        "HLA-C": parse_group(2),
    }


def auto_detect_column_map(src_df: pd.DataFrame, user_map: dict) -> dict:
    """
    自动识别常见列名；若用户显式传了 column_map，则优先用户映射。
    """
    final_map = dict(user_map)
    lower_to_raw = {c.lower(): c for c in src_df.columns}
    for target_col, aliases in COMMON_COLUMN_ALIASES.items():
        if target_col in final_map and final_map[target_col] in src_df.columns:
            continue
        for alias in aliases:
            raw = lower_to_raw.get(alias.lower())
            if raw is not None:
                final_map[target_col] = raw
                break
    return final_map


def build_output_df(src_df: pd.DataFrame, col_map: dict) -> pd.DataFrame:
    """
    按列映射将公开数据转为标准输出列。
    对缺失列填默认值，保证流程可跑。
    """
    out = pd.DataFrame()
    for col in REQUIRED_OUTPUT_COLUMNS:
        src_col = col_map.get(col, col)
        if src_col in src_df.columns:
            out[col] = src_df[src_col]
        else:
            if col == "variant_vaf":
                out[col] = 0.2
            elif col == "wt_peptide":
                out[col] = ""
            elif col == "transcript_id":
                out[col] = "NA"
            elif col == "mutation":
                out[col] = "UNKNOWN_MUTATION"
            else:
                out[col] = ""

    out["mut_peptide"] = out["mut_peptide"].astype(str).str.strip().str.upper()
    out["wt_peptide"] = out["wt_peptide"].astype(str).str.strip().str.upper()
    out["mutation"] = out["mutation"].astype(str).str.strip()
    out["transcript_id"] = out["transcript_id"].astype(str).str.strip()
    out["variant_vaf"] = pd.to_numeric(out["variant_vaf"], errors="coerce").fillna(0.2)

    # 如果 wt_peptide 缺失，则用 mut_peptide 做一个最小扰动，保证后续 dissimilarity 可计算
    missing_wt = out["wt_peptide"] == ""
    out.loc[missing_wt, "wt_peptide"] = out.loc[missing_wt, "mut_peptide"].apply(
        lambda p: (p[:-1] + "A") if len(p) > 0 else "A"
    )

    # 过滤空突变肽与过短肽段
    out = out[out["mut_peptide"].str.len() >= 8].copy()
    out = out.drop_duplicates(subset=["mutation", "mut_peptide"]).reset_index(drop=True)
    return out


def main(
    run_id: str,
    input_csv: str,
    case_id: str,
    column_map: str,
    hla: str
):
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"未找到输入 CSV: {input_csv}")

    src_df = pd.read_csv(input_csv)
    col_map = split_map_arg(column_map)
    col_map = auto_detect_column_map(src_df, col_map)
    out_df = build_output_df(src_df, col_map)
    if out_df.empty:
        raise ValueError("整理后没有有效候选肽，请检查输入列映射和数据内容。")

    out_dir = os.path.join("deliveries", run_id, "to_immunogen")
    os.makedirs(out_dir, exist_ok=True)

    neo_file = os.path.join(out_dir, "neoantigen_candidates.csv")
    hla_file = os.path.join(out_dir, "hla_typing.json")
    meta_file = os.path.join(out_dir, "meta.json")

    out_df.to_csv(neo_file, index=False, encoding="utf-8-sig")

    hla_data = normalize_hla_list(hla)
    with open(hla_file, "w", encoding="utf-8") as f:
        json.dump(hla_data, f, ensure_ascii=False, indent=2)

    meta_data = {
        "run_id": run_id,
        "case_id": case_id,
        "source": "public_dataset",
        "input_csv": input_csv,
        "rows_after_normalization": int(len(out_df)),
        "column_map": col_map,
    }
    with open(meta_file, "w", encoding="utf-8") as f:
        json.dump(meta_data, f, ensure_ascii=False, indent=2)

    print(f"完成: {neo_file}")
    print(f"完成: {hla_file}")
    print(f"完成: {meta_file}")
    print(f"候选肽条目数: {len(out_df)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R002")
    parser.add_argument("--input_csv", required=True, help="公开数据 CSV 路径")
    parser.add_argument("--case_id", default="public_case_001", help="病例 ID")
    parser.add_argument(
        "--column_map",
        default="",
        help="列映射，例如 mut_peptide=MT_pep,wt_peptide=WT_pep,mutation=gene_change,variant_vaf=vaf",
    )
    parser.add_argument(
        "--hla",
        default="A*02:01,A*11:01;B*40:01,B*46:01;C*01:02,C*07:02",
        help="HLA 输入，格式 A组;B组;C组",
    )
    args = parser.parse_args()
    main(args.run_id, args.input_csv, args.case_id, args.column_map, args.hla)
