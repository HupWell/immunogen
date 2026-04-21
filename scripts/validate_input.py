# -*- coding: utf-8 -*-
"""
功能：校验输入数据格式是否满足 ImmunoGen MVP 脚本要求。
输入：deliveries/<run_id>/to_immunogen/
输出：终端打印校验结果（通过/失败）
"""
import os
import json
import argparse
import pandas as pd


REQUIRED_CSV_COLUMNS = ["mut_peptide", "wt_peptide", "variant_vaf"]
OPTIONAL_CSV_COLUMNS = ["mutation", "transcript_id"]
REQUIRED_HLA_KEYS = ["HLA-A", "HLA-B", "HLA-C"]


def validate_csv(csv_path: str):
    """检查 neoantigen_candidates.csv 是否存在必需列与有效内容。"""
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"未找到输入文件: {csv_path}")

    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError("neoantigen_candidates.csv 为空。")

    missing = [c for c in REQUIRED_CSV_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"neoantigen_candidates.csv 缺少必需列: {missing}")

    # 检查 mut_peptide 是否至少有一条可用数据
    valid_peptides = df["mut_peptide"].astype(str).str.strip()
    if (valid_peptides == "").all():
        raise ValueError("mut_peptide 列全部为空，无法进行预测。")

    # 检查 VAF 是否可转数值（允许部分缺失）
    vaf_numeric = pd.to_numeric(df["variant_vaf"], errors="coerce")
    if vaf_numeric.isna().all():
        raise ValueError("variant_vaf 全部无法解析为数值，请检查格式。")

    return df


def validate_hla_json(hla_path: str):
    """检查 hla_typing.json 是否包含 HLA-A/B/C 且为列表。"""
    if not os.path.exists(hla_path):
        raise FileNotFoundError(f"未找到输入文件: {hla_path}")

    with open(hla_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    missing = [k for k in REQUIRED_HLA_KEYS if k not in data]
    if missing:
        raise ValueError(f"hla_typing.json 缺少必需键: {missing}")

    allele_total = 0
    for k in REQUIRED_HLA_KEYS:
        v = data.get(k)
        if not isinstance(v, list):
            raise ValueError(f"hla_typing.json 中 {k} 必须是列表。")
        allele_total += len([x for x in v if str(x).strip()])

    if allele_total == 0:
        raise ValueError("hla_typing.json 中没有有效的 HLA 等位基因。")

    return data


def validate_meta_json(meta_path: str):
    """检查 meta.json 是否存在且能正常解析。"""
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"未找到输入文件: {meta_path}")
    with open(meta_path, "r", encoding="utf-8") as f:
        return json.load(f)


def main(run_id: str):
    base = os.path.join("deliveries", run_id, "to_immunogen")
    csv_path = os.path.join(base, "neoantigen_candidates.csv")
    hla_path = os.path.join(base, "hla_typing.json")
    meta_path = os.path.join(base, "meta.json")

    df = validate_csv(csv_path)
    hla_data = validate_hla_json(hla_path)
    _ = validate_meta_json(meta_path)

    print("输入校验通过。")
    print(f"run_id: {run_id}")
    print(f"候选肽条目数: {len(df)}")
    print(f"CSV列: {list(df.columns)}")
    print(f"可选列建议: {OPTIONAL_CSV_COLUMNS}")
    print(
        "HLA计数: "
        f"A={len(hla_data.get('HLA-A', []))}, "
        f"B={len(hla_data.get('HLA-B', []))}, "
        f"C={len(hla_data.get('HLA-C', []))}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    args = parser.parse_args()
    main(args.run_id)
