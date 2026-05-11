# -*- coding: utf-8 -*-
"""
功能：从 peptide_mhc_ranking.csv 中筛选 Top-N 候选肽。
输入：results/<run_id>/peptide_mhc_ranking.csv
输出：results/<run_id>/selected_peptides.csv
"""
import os
import argparse
import pandas as pd


def hamming_dissimilarity(mut_peptide: str, wt_peptide: str) -> float:
    """
    计算突变肽与野生型肽的位点差异比例。
    - 若长度相同：按位比较不同字符占比
    - 若长度不同：返回 1.0（视为差异较大）
    """
    mut_peptide = str(mut_peptide or "").strip()
    wt_peptide = str(wt_peptide or "").strip()
    if not mut_peptide or not wt_peptide:
        return 0.0
    if len(mut_peptide) != len(wt_peptide):
        return 1.0
    diff_count = sum(1 for a, b in zip(mut_peptide, wt_peptide) if a != b)
    return diff_count / len(mut_peptide)


def _apply_positive_control_pins(
    selected: pd.DataFrame,
    filtered: pd.DataFrame,
    peptides_csv: str,
) -> pd.DataFrame:
    """
    对逗号分隔的 mut_peptide 序列做「阳性对照保入」：
    若在 filtered 池中但不在当前 Top-N 中，则用该肽替换选中集合里 rank_score 最低的一条（保持条数仍为 top_n）。
    """
    peptides = [p.strip().upper() for p in (peptides_csv or "").split(",") if p.strip()]
    if not peptides:
        return selected

    cols = list(selected.columns)
    # 转成「行列表」，便于原位替换且不破坏长度
    rows = selected.copy().reset_index(drop=True).to_dict(orient="records")
    fmp = filtered["mut_peptide"].astype(str).str.upper()

    for pep in peptides:
        present = any(str(r.get("mut_peptide", "")).strip().upper() == pep for r in rows)
        if present:
            continue
        match = filtered[fmp == pep]
        if match.empty:
            print(f"[positive_control] 跳过 {pep}：不在过滤后的 ranking 池中（可能被 min_dissimilarity 等排除）。")
            continue
        pin_ser = match.iloc[0]
        if not rows:
            break

        def _rank_key(rec: dict) -> float:
            v = pd.to_numeric(rec.get("rank_score"), errors="coerce")
            if pd.isna(v):
                return float("-inf")
            return float(v)

        worst_i = min(range(len(rows)), key=lambda i: _rank_key(rows[i]))
        base = rows[worst_i].copy()
        for col in cols:
            if col in pin_ser.index:
                base[col] = pin_ser[col]
        rows[worst_i] = base
        print(f"[positive_control] 已将 rank_score 最低条目替换为 mut_peptide={pep}。")

    out = pd.DataFrame(rows)
    if not out.empty and "rank_score" in out.columns:
        out = out.sort_values("rank_score", ascending=False, na_position="last").reset_index(drop=True)
    return out


def main(run_id: str, top_n: int, min_dissimilarity: float, ensure_positive_control_peptides: str):
    """
    主流程：
    1) 读取 peptide_mhc_ranking.csv；
    2) 以 mut_peptide 为粒度聚合（取每条肽最佳 rank_score）；
    3) 按与 WT 的差异比例过滤；
    4) 选 Top-N 输出 selected_peptides.csv。
    """
    in_file = os.path.join("results", run_id, "peptide_mhc_ranking.csv")
    out_file = os.path.join("results", run_id, "selected_peptides.csv")

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"未找到输入文件: {in_file}")

    df = pd.read_csv(in_file)
    if df.empty:
        raise ValueError("输入文件为空，无法筛选 Top 肽。")

    # 统一数值列类型，便于后续排序和导出
    df["rank_score"] = pd.to_numeric(df.get("rank_score"), errors="coerce")
    df["affinity_nM"] = pd.to_numeric(df.get("affinity_nM"), errors="coerce")
    df["variant_vaf"] = pd.to_numeric(df.get("variant_vaf"), errors="coerce")

    # 计算突变肽与 WT 差异比例（越大表示与 WT 越不相似）
    df["dissimilarity"] = df.apply(
        lambda r: hamming_dissimilarity(r.get("mut_peptide", ""), r.get("wt_peptide", "")),
        axis=1
    )

    # 聚合到“每条突变肽一行”：保留最佳 rank_score 和对应 HLA
    df = df.sort_values("rank_score", ascending=False, na_position="last")
    grouped = (
        df.groupby("mut_peptide", as_index=False)
        .first()
    )

    # 过滤与 WT 过于相似的候选肽
    filtered = grouped[grouped["dissimilarity"] >= min_dissimilarity].copy()
    selected = filtered.head(top_n).copy()
    selected = _apply_positive_control_pins(selected, filtered, ensure_positive_control_peptides)

    # 输出列顺序尽量清晰，便于后续 ORF 构建使用
    output_columns = [
        "mutation",
        "mut_peptide",
        "wt_peptide",
        "variant_vaf",
        "hla_allele",
        "affinity_nM",
        "rank_score",
        "dissimilarity",
    ]
    selected = selected[[c for c in output_columns if c in selected.columns]]
    selected.to_csv(out_file, index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")
    print(f"原始突变肽数: {grouped.shape[0]}")
    print(f"过滤后突变肽数: {filtered.shape[0]}")
    print(f"最终入选数: {selected.shape[0]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_n", type=int, default=10, help="最终选择数量，默认 10")
    parser.add_argument(
        "--min_dissimilarity",
        type=float,
        default=0.1,
        help="与 WT 最小差异比例阈值，默认 0.1",
    )
    parser.add_argument(
        "--ensure_positive_control_peptides",
        default="",
        help=(
            "逗号分隔的 mut_peptide（序列）；若在过滤池中但未进 Top-N，"
            "则替换选中集合 rank_score 最低的一条。用于任务书 Positive Control。"
        ),
    )
    args = parser.parse_args()
    main(
        args.run_id,
        args.top_n,
        args.min_dissimilarity,
        args.ensure_positive_control_peptides,
    )
