# -*- coding: utf-8 -*-
"""
四 run 路径 A：按 Top-K（默认 5）刷新 PANDORA 结构 + SimHub 交付包。

用法：
  python scripts/run_path_a_simhub_top5_batch.py
  python scripts/run_path_a_simhub_top5_batch.py --runs R_public_001 --skip_pandora
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path
from typing import List


DEFAULT_RUNS = ["R001", "R002", "R003", "R_public_001"]


def _run(cmd: List[str], env: dict) -> None:
    print(f"\n[执行] {' '.join(cmd)}")
    ret = subprocess.run(cmd, env=env, cwd=str(Path(__file__).resolve().parents[1]), check=False)
    if ret.returncode != 0:
        raise RuntimeError(f"失败: {' '.join(cmd)} (exit {ret.returncode})")


def main() -> None:
    parser = argparse.ArgumentParser(description="四 run SimHub Top-5 批量（选肽→PANDORA→SimHub）")
    parser.add_argument("--runs", default=",".join(DEFAULT_RUNS), help="逗号分隔 run_id")
    parser.add_argument("--top_k", type=int, default=5)
    parser.add_argument("--top_n_select", type=int, default=10, help="选肽池 Top-N（须 >= top_k）")
    parser.add_argument("--skip_pandora", action="store_true", help="仅刷新 SimHub 包（沿用已有 rank PDB）")
    parser.add_argument("--skip_select", action="store_true")
    parser.add_argument("--pandora_python", default="", help="含 PANDORA 的解释器，默认当前 python")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    py = (args.pandora_python or "").strip() or sys.executable
    scripts = repo / "scripts"
    runs = [x.strip() for x in args.runs.split(",") if x.strip()]

    env = os.environ.copy()
    muscle = repo / "tools" / "pandora_bin" / "muscle"
    if muscle.is_file():
        pandora_bin = repo / "tools" / "pandora_bin"
        env["PATH"] = str(pandora_bin) + os.pathsep + env.get("PATH", "")
        env["PANDORA_REAL_MUSCLE"] = str(muscle)

    for run_id in runs:
        print(f"\n{'=' * 60}\nrun_id={run_id}\n{'=' * 60}")
        if not args.skip_select:
            _run(
                [py, str(scripts / "select_top_peptides.py"), "--run_id", run_id, "--top_n", str(args.top_n_select)],
                env,
            )
        if not args.skip_pandora:
            _run(
                [
                    py,
                    str(scripts / "run_pandora_structure.py"),
                    "--run_id",
                    run_id,
                    "--top_k",
                    str(args.top_k),
                ],
                env,
            )
        _run(
            [
                py,
                str(scripts / "prepare_simhub_delivery.py"),
                "--run_id",
                run_id,
                "--top_k",
                str(args.top_k),
                "--structure_backend",
                "pandora",
            ],
            env,
        )
        _run([py, str(scripts / "check_simhub_evidence.py"), "--run_id", run_id], env)

    print("\n批量完成。")


if __name__ == "__main__":
    main()
