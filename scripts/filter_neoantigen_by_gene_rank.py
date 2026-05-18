# -*- coding: utf-8 -*-
"""
按路径 B 基因排名过滤 neoantigen_candidates.csv，供路径 A / SimHub 使用。

输入：
- results/<run_id>/<rank_subdir>/candidate_gene_ranking.csv
- deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv

输出：
- 覆盖写入 neoantigen_candidates.csv（首次备份为 neoantigen_candidates.pre_gene_filter.csv）
- results/<run_id>/<rank_subdir>/neoantigen_gene_filter_provenance.json
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Set

import pandas as pd

GENE_COLUMNS = ["gene_name", "gene_symbol", "gene", "hugo_symbol", "symbol"]
TRANSCRIPT_COLUMNS = ["transcript_id", "transcript", "enst_id", "tx_id"]
MUTATION_GENE_RE = re.compile(r"^([A-Z0-9]+)[_.]")


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def infer_gene_from_mutation(mutation: Any) -> Optional[str]:
    """从 mutation 字段解析基因符号，如 KRAS_G12D_1 -> KRAS。"""
    if mutation is None or (isinstance(mutation, float) and pd.isna(mutation)):
        return None
    text = str(mutation).strip()
    if not text or text.lower() == "nan":
        return None
    m = MUTATION_GENE_RE.match(text.upper())
    return m.group(1) if m else None


def normalize_gene_token(value: Any) -> Optional[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    text = str(value).strip().upper()
    if not text or text == "NAN":
        return None
    return text


def resolve_row_genes(row: pd.Series) -> Set[str]:
    """汇总一行候选肽可关联的基因符号集合。"""
    genes: Set[str] = set()
    for col in GENE_COLUMNS:
        if col in row.index:
            g = normalize_gene_token(row[col])
            if g:
                genes.add(g)
    g_mut = infer_gene_from_mutation(row.get("mutation"))
    if g_mut:
        genes.add(g_mut)
    return genes


def load_top_genes(
    rank_csv: str,
    top_n_genes: int,
    gene_col: str,
    min_log2fc: Optional[float],
) -> List[str]:
    rank_df = pd.read_csv(rank_csv)
    if rank_df.empty:
        raise ValueError(f"基因排名表为空: {rank_csv}")
    col = gene_col if gene_col in rank_df.columns else None
    if col is None:
        for c in ("gene_name", "gene_symbol", "gene"):
            if c in rank_df.columns:
                col = c
                break
    if col is None:
        raise ValueError(f"排名表缺少基因列（期望 {gene_col} 或 gene_name）: {rank_csv}")

    work = rank_df.copy()
    if "rank" in work.columns:
        work = work.sort_values("rank", ascending=True)
    elif "log2_fold_change" in work.columns:
        work = work.sort_values("log2_fold_change", ascending=False)

    if min_log2fc is not None and "log2_fold_change" in work.columns:
        lfc = pd.to_numeric(work["log2_fold_change"], errors="coerce")
        work = work[lfc >= float(min_log2fc)]

    if top_n_genes > 0:
        work = work.head(int(top_n_genes))

    genes: List[str] = []
    seen: Set[str] = set()
    for raw in work[col].astype(str):
        g = normalize_gene_token(raw)
        if g and g not in seen:
            seen.add(g)
            genes.append(g)
    if not genes:
        raise ValueError("Top 基因列表为空，请检查 top_n_genes / min_log2fc。")
    return genes


def filter_candidates(
    neo_df: pd.DataFrame,
    top_genes: List[str],
    match_transcript: bool,
) -> pd.DataFrame:
    allow = {normalize_gene_token(g) for g in top_genes if normalize_gene_token(g)}
    allow_tx: Set[str] = set()
    if match_transcript:
        for col in TRANSCRIPT_COLUMNS:
            if col not in neo_df.columns:
                continue
            for v in neo_df[col].astype(str):
                t = str(v).strip()
                if t and t.lower() != "nan":
                    allow_tx.add(t)

    kept: List[bool] = []
    for _, row in neo_df.iterrows():
        row_genes = resolve_row_genes(row)
        hit_gene = bool(row_genes & allow)
        hit_tx = False
        if match_transcript and allow_tx:
            for col in TRANSCRIPT_COLUMNS:
                if col not in row.index:
                    continue
                t = str(row[col]).strip()
                if t and t.lower() != "nan" and t in allow_tx:
                    hit_tx = True
                    break
        kept.append(hit_gene or hit_tx)

    return neo_df.loc[kept].copy()


def main() -> None:
    parser = argparse.ArgumentParser(description="按转录组 Top 基因过滤 neoantigen_candidates.csv")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--top_n_genes", type=int, default=500, help="取排名前 N 个基因，0 表示不截断")
    parser.add_argument(
        "--rank_subdir",
        default="transcriptome_prior",
        help="results/<run_id>/<subdir>/candidate_gene_ranking.csv",
    )
    parser.add_argument("--gene_col_name", default="gene_name", help="排名表中的基因列名")
    parser.add_argument(
        "--min_log2fc",
        type=float,
        default=None,
        help="可选：仅保留 log2_fold_change >= 该值的基因后再取 Top-N",
    )
    parser.add_argument(
        "--match_transcript",
        action="store_true",
        help="除基因符号外，允许 transcript_id 精确命中（排名表需含 transcript 列，默认关闭）",
    )
    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="只统计将保留的条数，不写文件",
    )
    parser.add_argument(
        "--allow_empty",
        action="store_true",
        help="过滤后 0 条也不报错（默认报错）",
    )
    args = parser.parse_args()

    run_id = args.run_id.strip()
    rank_dir = os.path.join("results", run_id, args.rank_subdir.strip() or "transcriptome_prior")
    rank_csv = os.path.join(rank_dir, "candidate_gene_ranking.csv")
    immuno_dir = os.path.join("deliveries", run_id, "to_immunogen")
    neo_path = os.path.join(immuno_dir, "neoantigen_candidates.csv")
    backup_path = os.path.join(immuno_dir, "neoantigen_candidates.pre_gene_filter.csv")

    if not os.path.isfile(rank_csv):
        raise FileNotFoundError(f"未找到基因排名: {rank_csv}，请先跑路径 B（transcriptome_prior）。")
    if not os.path.isfile(neo_path):
        raise FileNotFoundError(f"未找到新抗原表: {neo_path}")

    top_genes = load_top_genes(
        rank_csv, args.top_n_genes, args.gene_col_name, args.min_log2fc
    )
    neo_df = pd.read_csv(neo_path)
    n_before = len(neo_df)
    filtered = filter_candidates(neo_df, top_genes, args.match_transcript)
    n_after = len(filtered)

    sample_genes = sorted({g for _, row in neo_df.iterrows() for g in resolve_row_genes(row)})
    matched_genes = sorted(set(top_genes) & set(sample_genes))

    print(f"run_id: {run_id}")
    print(f"Top 基因数: {len(top_genes)}（取前 {args.top_n_genes}，min_log2fc={args.min_log2fc}）")
    print(f"候选肽基因可解析: {sample_genes[:20]}{'...' if len(sample_genes) > 20 else ''}")
    print(f"与 Top 池交集基因: {matched_genes}")
    print(f"过滤: {n_before} -> {n_after}")

    if n_after == 0 and not args.allow_empty:
        raise ValueError(
            "过滤后无候选肽。请确认 neoantigen 表含 gene_name/gene_symbol，"
            "或 mutation 形如 GENE_*；亦可增大 --top_n_genes。"
        )

    provenance: Dict[str, Any] = {
        "tool": "filter_neoantigen_by_gene_rank.py",
        "run_id": run_id,
        "timestamp_utc": _utc_now_iso(),
        "rank_csv": os.path.abspath(rank_csv),
        "neoantigen_csv": os.path.abspath(neo_path),
        "top_n_genes": args.top_n_genes,
        "min_log2fc": args.min_log2fc,
        "n_top_genes_loaded": len(top_genes),
        "n_candidates_before": n_before,
        "n_candidates_after": n_after,
        "matched_genes_in_neo": matched_genes,
        "match_transcript": bool(args.match_transcript),
        "caution_zh": (
            "突变肽池已按路径 B 基因排名收窄；SimHub 交付仍走路径 A 同一套脚本。"
            "若 ranking 为探索性分组，过滤结果不得单独作为疫苗终版依据。"
        ),
    }

    if args.dry_run:
        print("dry_run：未写入文件。")
        return

    os.makedirs(rank_dir, exist_ok=True)
    if not os.path.isfile(backup_path):
        shutil.copy2(neo_path, backup_path)
        print(f"已备份: {backup_path}")

    filtered.to_csv(neo_path, index=False, encoding="utf-8")
    prov_path = os.path.join(rank_dir, "neoantigen_gene_filter_provenance.json")
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(provenance, f, ensure_ascii=False, indent=2)

    print(f"已写入: {neo_path}")
    print(f"已写入: {prov_path}")


if __name__ == "__main__":
    main()
