# -*- coding: utf-8 -*-
"""
B → A → SimHub 一键：基因排名过滤突变肽池后，执行路径 A 至 SimHub 交付。

前置（同一 run_id）：
- results/<run_id>/transcriptome_prior/candidate_gene_ranking.csv
- deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv + hla_typing.json
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from typing import List, Optional


def _run(cmd: List[str], env: Optional[dict] = None) -> None:
    print(f"\n[执行] {' '.join(cmd)}")
    ret = subprocess.run(cmd, env=env, check=False)
    if ret.returncode != 0:
        raise RuntimeError(f"命令失败，退出码 {ret.returncode}")


def main() -> None:
    parser = argparse.ArgumentParser(description="B→A→SimHub：过滤 + run_all 路径 A")
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--top_n_genes", type=int, default=500)
    parser.add_argument("--rank_subdir", default="transcriptome_prior")
    parser.add_argument("--min_log2fc", type=float, default=None)
    parser.add_argument("--skip_filter", action="store_true", help="跳过基因过滤（已手动过滤时）")
    parser.add_argument(
        "--target",
        default="full",
        choices=["full", "mhc_ranking", "simhub", "report"],
        help="过滤后交给 run_all 的路径 A 范围",
    )
    parser.add_argument("--top_n", type=int, default=10, help="路径 A 选肽 Top-N")
    parser.add_argument("--top_k_md", type=int, default=3)
    parser.add_argument("--prepare_pandora", action="store_true")
    parser.add_argument("--pandora_python", default="")
    parser.add_argument(
        "--ensure_positive_control_peptides",
        default="",
        help="传给 run_all / select_top_peptides",
    )
    parser.add_argument(
        "--run_all_extra",
        default="",
        help="追加给 run_all.py 的参数字符串（引号内），如 '--allow_proxy_scores'",
    )
    args, unknown = parser.parse_known_args()
    if unknown:
        print(f"警告: 未识别参数将忽略: {unknown}")

    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(scripts_dir)
    python_exe = sys.executable
    run_id = args.run_id.strip()

    rank_csv = os.path.join(
        repo_root,
        "results",
        run_id,
        args.rank_subdir,
        "candidate_gene_ranking.csv",
    )
    neo_path = os.path.join(repo_root, "deliveries", run_id, "to_immunogen", "neoantigen_candidates.csv")
    if not os.path.isfile(neo_path):
        raise FileNotFoundError(f"缺少路径 A 输入: {neo_path}")
    if not args.skip_filter and not os.path.isfile(rank_csv):
        raise FileNotFoundError(f"缺少路径 B 排名: {rank_csv}")

    os.chdir(repo_root)

    if not args.skip_filter:
        filt_cmd = [
            python_exe,
            os.path.join(scripts_dir, "filter_neoantigen_by_gene_rank.py"),
            "--run_id",
            run_id,
            "--top_n_genes",
            str(args.top_n_genes),
            "--rank_subdir",
            args.rank_subdir,
        ]
        if args.min_log2fc is not None:
            filt_cmd.extend(["--min_log2fc", str(args.min_log2fc)])
        _run(filt_cmd)

    run_all_cmd = [
        python_exe,
        os.path.join(scripts_dir, "run_all.py"),
        "--run_id",
        run_id,
        "--target",
        args.target,
        "--top_n",
        str(args.top_n),
        "--top_k_md",
        str(args.top_k_md),
    ]
    if args.prepare_pandora:
        run_all_cmd.append("--prepare_pandora")
    if args.pandora_python.strip():
        run_all_cmd.extend(["--pandora_python", args.pandora_python.strip()])
    if args.ensure_positive_control_peptides.strip():
        run_all_cmd.extend(
            ["--ensure_positive_control_peptides", args.ensure_positive_control_peptides.strip()]
        )
    extra = (args.run_all_extra or "").strip()
    if extra:
        run_all_cmd.extend(extra.split())

    _run(run_all_cmd)
    print("\nB→A→SimHub 流程结束。")


if __name__ == "__main__":
    main()
