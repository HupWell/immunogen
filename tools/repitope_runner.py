#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Repitope real_cmd 包装器（轻量实现）。

输入：--input  指向包含 mut_peptide 列的 TSV
输出：--output 指向包含 mut_peptide + score 的 TSV

实现说明：
1) 使用 external_refs/Repitope/inst 下公开数据集（Calis/MHCBN/EPIMHC/...）构建参考库；
2) 若输入肽在参考库中有完全匹配，输出该肽免疫原性标签均值（Positive=1, Negative=0）；
3) 否则在同长度参考肽中做序列近邻（difflib ratio），按相似度加权得到分数。
"""
import argparse
import csv
import os
from difflib import SequenceMatcher
from typing import Dict, List, Tuple

import pandas as pd


_DATA_FILES = [
    "Calis1.csv",
    "Calis2.csv",
    "EPIMHC.csv",
    "IMMA2.csv",
    "MHCBN.csv",
    "TANTIGEN.csv",
    "HCV.csv",
    "HIV.csv",
]


def _label_to_score(v: str) -> float:
    x = (v or "").strip().lower()
    if x in {"positive", "pos", "1", "true", "yes"}:
        return 1.0
    if x in {"negative", "neg", "0", "false", "no"}:
        return 0.0
    raise ValueError(f"未知 Immunogenicity 标签: {v}")


def _load_reference(repo_root: str) -> Tuple[Dict[str, float], Dict[int, List[Tuple[str, float]]]]:
    inst_dir = os.path.join(repo_root, "inst")
    if not os.path.isdir(inst_dir):
        raise FileNotFoundError(f"未找到 Repitope inst 目录: {inst_dir}")

    exact_bucket: Dict[str, List[float]] = {}
    by_len_bucket: Dict[int, List[Tuple[str, float]]] = {}
    for name in _DATA_FILES:
        path = os.path.join(inst_dir, name)
        if not os.path.exists(path):
            continue
        with open(path, "r", encoding="utf-8-sig", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                pep = (row.get("Peptide") or "").strip().upper()
                if not pep:
                    continue
                try:
                    sc = _label_to_score(row.get("Immunogenicity", ""))
                except Exception:
                    continue
                exact_bucket.setdefault(pep, []).append(sc)
                by_len_bucket.setdefault(len(pep), []).append((pep, sc))

    exact_mean = {k: sum(v) / max(len(v), 1) for k, v in exact_bucket.items()}
    if not exact_mean:
        raise RuntimeError("Repitope 参考库为空，无法进行 real_cmd 评分。")
    return exact_mean, by_len_bucket


def _knn_score(peptide: str, refs: List[Tuple[str, float]], top_k: int = 64) -> float:
    sims: List[Tuple[float, float]] = []
    for ref_pep, ref_score in refs:
        sim = SequenceMatcher(None, peptide, ref_pep).ratio()
        if sim <= 0.0:
            continue
        sims.append((sim, ref_score))
    if not sims:
        return 0.5
    sims.sort(key=lambda x: x[0], reverse=True)
    picked = sims[:top_k]
    w_sum = sum(w for w, _ in picked)
    if w_sum <= 0.0:
        return 0.5
    return sum(w * s for w, s in picked) / w_sum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="输入 TSV（需含 mut_peptide）")
    parser.add_argument("--output", required=True, help="输出 TSV（mut_peptide + score）")
    parser.add_argument(
        "--repo_root",
        default=os.path.join("external_refs", "Repitope"),
        help="Repitope 仓库目录（默认 external_refs/Repitope）",
    )
    args = parser.parse_args()

    src = pd.read_csv(args.input, sep="\t")
    if "mut_peptide" not in src.columns:
        raise ValueError("输入文件缺少 mut_peptide 列。")
    peptides = (
        src["mut_peptide"].astype(str).str.strip().str.upper().replace("", pd.NA).dropna().drop_duplicates()
    )
    exact_mean, by_len = _load_reference(os.path.abspath(args.repo_root))

    rows = []
    for pep in peptides.tolist():
        if pep in exact_mean:
            score = exact_mean[pep]
            source = "repitope_public_dataset_exact"
        else:
            refs = by_len.get(len(pep), [])
            if not refs:
                # 无同长度时降级到全库近邻，避免直接空结果。
                refs = [item for items in by_len.values() for item in items]
            score = _knn_score(pep, refs)
            source = "repitope_public_dataset_knn"
        rows.append(
            {
                "mut_peptide": pep,
                "score": round(float(max(0.0, min(1.0, score))), 6),
                "source": source,
                "version": "Repitope_public_dataset_v1",
            }
        )

    out = pd.DataFrame(rows).drop_duplicates(subset=["mut_peptide"], keep="first")
    out.to_csv(args.output, sep="\t", index=False, encoding="utf-8")
    print(f"Repitope real_cmd 完成，输出 {len(out)} 条。")


if __name__ == "__main__":
    main()
