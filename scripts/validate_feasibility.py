# -*- coding: utf-8 -*-
"""
功能：自动评估 ImmunoGen 结果的基础可行性。

输入：
- results/<run_id>/peptide_mhc_ranking.csv
- results/<run_id>/selected_peptides.csv
- results/<run_id>/mrna_design.json

输出：
- results/<run_id>/feasibility_report.json
- results/<run_id>/FEASIBILITY.md
"""
import os
import json
import argparse
from typing import List, Set
import pandas as pd


def jaccard(a: Set[str], b: Set[str]) -> float:
    """计算两个集合的 Jaccard 相似度。"""
    if not a and not b:
        return 1.0
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def top_peptides_by_score(df: pd.DataFrame, top_n: int) -> List[str]:
    """按 rank_score 取 Top N 突变肽（去重后）。"""
    ranked = (
        df.sort_values("rank_score", ascending=False, na_position="last")
        .drop_duplicates(subset=["mut_peptide"])
    )
    return ranked["mut_peptide"].head(top_n).astype(str).tolist()


def evaluate_positive_controls(top_peptides: List[str]) -> dict:
    """
    Positive control 检查：
    当前示例以 KRAS G12D 常见核心片段做简易匹配。
    """
    controls = {
        "KRAS_G12D_like": ["GADGVG", "VVGADGVGK"],
    }
    result = {}
    for name, motifs in controls.items():
        hit = any(any(m in pep for m in motifs) for pep in top_peptides)
        result[name] = bool(hit)
    return result


def build_markdown(summary: dict) -> str:
    """生成人类可读的可行性报告。"""
    pc = summary["positive_control"]
    lines = [
        f"# 可行性验证报告（{summary['run_id']}）",
        "",
        "## 1. 基础统计",
        f"- 排名总行数：`{summary['ranking_rows']}`",
        f"- 去重后候选肽数：`{summary['unique_peptides']}`",
        f"- 入选肽数：`{summary['selected_count']}`",
        f"- mRNA 总长度：`{summary['mrna_length']}`",
        f"- mRNA GC%：`{summary['mrna_gc_percent']}`",
        "",
        "## 2. Top 稳定性（参数扰动）",
        f"- baseline Top{summary['top_n']} 与 variant Top{summary['top_n']} 的 Jaccard：`{summary['top_stability_jaccard']:.3f}`",
        f"- baseline Top{summary['top_n']}：`{', '.join(summary['top_baseline'])}`",
        f"- variant Top{summary['top_n']}：`{', '.join(summary['top_variant'])}`",
        "",
        "## 3. Positive Control 命中",
        f"- KRAS G12D 相关片段命中：`{pc.get('KRAS_G12D_like', False)}`",
        "",
        "## 4. 结论",
        f"- 工程可复现：`{summary['engineering_reproducible']}`",
        f"- 方向性可行：`{summary['directionally_feasible']}`",
        "- 说明：当前包含代理评分组件，结论用于流程与方向验证，不替代湿实验结论。",
        "",
    ]
    return "\n".join(lines)


def main(run_id: str, top_n: int):
    out_dir = os.path.join("results", run_id)
    ranking_file = os.path.join(out_dir, "peptide_mhc_ranking.csv")
    selected_file = os.path.join(out_dir, "selected_peptides.csv")
    design_file = os.path.join(out_dir, "mrna_design.json")

    for f in [ranking_file, selected_file, design_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"未找到输入文件: {f}")

    ranking_df = pd.read_csv(ranking_file)
    selected_df = pd.read_csv(selected_file)
    with open(design_file, "r", encoding="utf-8") as f:
        design = json.load(f)

    ranking_df["rank_score"] = pd.to_numeric(ranking_df.get("rank_score"), errors="coerce")
    ranking_df = ranking_df.dropna(subset=["mut_peptide"]).copy()

    baseline_top = top_peptides_by_score(ranking_df, top_n=top_n)

    # 简单参数扰动：将 rank_score 与 immunogenicity 做轻度混合，检查 Top 稳定性
    variant_df = ranking_df.copy()
    if "immunogenicity" in variant_df.columns:
        variant_df["immunogenicity"] = pd.to_numeric(variant_df["immunogenicity"], errors="coerce").fillna(0.0)
        variant_df["rank_score"] = variant_df["rank_score"].fillna(0.0) * 0.85 + variant_df["immunogenicity"] * 0.15
    variant_top = top_peptides_by_score(variant_df, top_n=top_n)
    stability = jaccard(set(baseline_top), set(variant_top))

    positive_control = evaluate_positive_controls(baseline_top)
    engineering_reproducible = bool(len(selected_df) > 0 and design.get("quality_metrics"))
    directionally_feasible = bool(stability >= 0.5 or any(positive_control.values()))

    summary = {
        "run_id": run_id,
        "ranking_rows": int(len(ranking_df)),
        "unique_peptides": int(ranking_df["mut_peptide"].nunique()),
        "selected_count": int(len(selected_df)),
        "mrna_length": int(design.get("quality_metrics", {}).get("full_length", 0)),
        "mrna_gc_percent": float(design.get("quality_metrics", {}).get("gc_percent", 0.0)),
        "top_n": int(top_n),
        "top_baseline": baseline_top,
        "top_variant": variant_top,
        "top_stability_jaccard": float(stability),
        "positive_control": positive_control,
        "engineering_reproducible": engineering_reproducible,
        "directionally_feasible": directionally_feasible,
    }

    json_file = os.path.join(out_dir, "feasibility_report.json")
    md_file = os.path.join(out_dir, "FEASIBILITY.md")
    with open(json_file, "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=2)
    with open(md_file, "w", encoding="utf-8") as f:
        f.write(build_markdown(summary))

    print(f"完成: {json_file}")
    print(f"完成: {md_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_n", type=int, default=10, help="用于稳定性比较的 TopN，默认 10")
    args = parser.parse_args()
    main(args.run_id, args.top_n)
