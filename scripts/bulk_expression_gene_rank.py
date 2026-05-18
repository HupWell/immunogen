# -*- coding: utf-8 -*-
"""
从 bulk 表达矩阵（基因 × 样本）在「病例 vs 对照」样本间计算均值比的对数，
输出**基因层面**优先级表，供后续 WES/Panel、功能验证或 BioDriver 选题使用。

科学边界（必读）：
- 本脚本**不**产生 mut_peptide / wt_peptide / variant_vaf，**不能**替代 neoantigen_candidates.csv；
- 无 DNA 层信息时，任何「疫苗肽段」结论都不应由本表单独支撑；
- log2FC 为转录组代理指标，受细胞组成、批次等混杂影响，需结合实验设计与质控解读。
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from typing import Dict, List, Sequence, Tuple


def _sniff_sep(first_line: str) -> str:
    """根据首行粗略判断分隔符。"""
    n_tab = first_line.count("\t")
    n_comma = first_line.count(",")
    if n_tab >= 1 and n_tab >= n_comma:
        return "\t"
    return ","


def _read_matrix(path: str) -> Tuple[str, List[str], Dict[str, List[float]]]:
    """读取矩阵：首列为基因名，其余列为数值；返回（基因列名, 样本列名列表, 基因→按样本顺序的值列表）。"""
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        first = f.readline()
        if not first:
            raise ValueError("空文件。")
        sep = _sniff_sep(first)
        header = [c.strip() for c in first.strip().split(sep)]
        if len(header) < 2:
            raise ValueError("表头至少需要「基因列 + 1 个样本列」。")
        gene_col = header[0]
        sample_cols = header[1:]
        genes: Dict[str, List[float]] = {}
        reader = csv.reader(f, delimiter=sep)
        for row in reader:
            if not row or len(row) < 2:
                continue
            g = row[0].strip()
            if not g:
                continue
            vals: List[float] = []
            for x in row[1 : len(sample_cols) + 1]:
                try:
                    vals.append(float(x))
                except ValueError:
                    vals.append(float("nan"))
            if len(vals) < len(sample_cols):
                vals.extend([float("nan")] * (len(sample_cols) - len(vals)))
            genes[g] = vals
    return gene_col, sample_cols, genes


def _mean_non_nan(xs: Sequence[float]) -> float:
    ys = [x for x in xs if not math.isnan(x)]
    if not ys:
        return float("nan")
    return sum(ys) / len(ys)


def _log2fc(case_vals: Sequence[float], ctrl_vals: Sequence[float], pseudo: float) -> Tuple[float, float, float]:
    mc = _mean_non_nan(case_vals) + pseudo
    mv = _mean_non_nan(ctrl_vals) + pseudo
    if math.isnan(mc) or math.isnan(mv) or mc <= 0 or mv <= 0:
        return float("nan"), mc - pseudo, mv - pseudo
    return math.log2(mc / mv), mc - pseudo, mv - pseudo


def main() -> None:
    parser = argparse.ArgumentParser(description="bulk 表达矩阵 → 基因优先级（非新抗原肽表）")
    parser.add_argument("--tpm_path", required=True, help="基因×样本矩阵路径（首列为基因名，其余为数值）")
    parser.add_argument(
        "--out_dir",
        required=True,
        help="输出目录（将写入 candidate_gene_ranking.csv 与 bulk_rank_provenance.json）",
    )
    parser.add_argument(
        "--gene_col_name",
        default="gene_name",
        help="输出中基因列显示名（默认 gene_name；若矩阵首列实际为 symbol 可改为 gene_symbol）",
    )
    parser.add_argument(
        "--case_columns",
        default="",
        help="病例组样本列名，逗号分隔，须与矩阵表头完全一致（如 GSM575076,GSM575077）",
    )
    parser.add_argument(
        "--control_columns",
        default="",
        help="对照组样本列名，逗号分隔，须与矩阵表头完全一致",
    )
    parser.add_argument("--pseudo", type=float, default=1.0, help="均值平滑伪计数（默认 1，避免 log(0)）")
    args = parser.parse_args()

    case_cols = [c.strip() for c in args.case_columns.split(",") if c.strip()]
    ctrl_cols = [c.strip() for c in args.control_columns.split(",") if c.strip()]
    if not case_cols or not ctrl_cols:
        print(
            "错误：必须同时提供 --case_columns 与 --control_columns（各至少 1 个样本列）。\n"
            "不设对照的「全样本排序」混杂因素过多，本脚本刻意不支持，以免误用。\n"
            "请从临床/pdata 中定义可比的两组样本，再传入对应列名。",
            file=sys.stderr,
        )
        sys.exit(2)

    file_gene_key, all_samples, genes = _read_matrix(args.tpm_path)

    def _indices(names: List[str], label: str) -> List[int]:
        """返回各样本在 vals 列表中的下标（与 sample_cols 对齐）。"""
        idxs: List[int] = []
        for n in names:
            if n not in all_samples:
                raise KeyError(f"{label} 列名不在矩阵表头中: {n!r}；可用列示例: {all_samples[:5]} ...")
            idxs.append(all_samples.index(n))
        return idxs

    try:
        case_idx = _indices(case_cols, "case")
        ctrl_idx = _indices(ctrl_cols, "control")
    except KeyError as e:
        print(f"错误: {e}", file=sys.stderr)
        sys.exit(2)

    overlap = set(case_cols) & set(ctrl_cols)
    if overlap:
        print(f"错误：病例与对照样本名不可重叠: {sorted(overlap)}", file=sys.stderr)
        sys.exit(2)

    rows_out: List[Dict[str, object]] = []
    pseudo = float(args.pseudo)
    for g, vals in genes.items():
        if len(vals) < len(all_samples):
            continue
        case_v = [vals[i] for i in case_idx]
        ctrl_v = [vals[j] for j in ctrl_idx]
        log2fc, mean_case, mean_ctrl = _log2fc(case_v, ctrl_v, pseudo)
        rows_out.append(
            {
                args.gene_col_name: g,
                "mean_tpm_case": mean_case,
                "mean_tpm_control": mean_ctrl,
                "log2_fold_change": log2fc,
            }
        )

    # 按 log2FC 降序（病例上调优先）
    rows_out.sort(
        key=lambda r: (float(r["log2_fold_change"]) if not math.isnan(float(r["log2_fold_change"])) else -1e9),
        reverse=True,
    )
    for i, r in enumerate(rows_out, start=1):
        r["rank"] = i

    os.makedirs(args.out_dir, exist_ok=True)
    out_csv = os.path.join(args.out_dir, "candidate_gene_ranking.csv")
    fieldnames = ["rank", args.gene_col_name, "mean_tpm_case", "mean_tpm_control", "log2_fold_change"]
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow({k: r[k] for k in fieldnames})

    provenance: Dict[str, object] = {
        "tool": "bulk_expression_gene_rank.py",
        "tpm_path": os.path.abspath(args.tpm_path),
        "matrix_gene_column_in_file": file_gene_key,
        "case_columns": case_cols,
        "control_columns": ctrl_cols,
        "pseudo_count": pseudo,
        "n_genes_ranked": len(rows_out),
        "caution_zh": (
            "本输出为转录组层面的基因优先级，不等价于新抗原；"
            "进入 ImmunoGen 主线仍需 neoantigen_candidates.csv 与 hla_typing.json。"
        ),
    }
    with open(os.path.join(args.out_dir, "bulk_rank_provenance.json"), "w", encoding="utf-8") as f:
        json.dump(provenance, f, ensure_ascii=False, indent=2)

    print(f"已写入: {out_csv}")
    print(f"已写入: {os.path.join(args.out_dir, 'bulk_rank_provenance.json')}")


if __name__ == "__main__":
    main()
