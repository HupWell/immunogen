#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
一键检查并驱动真实后端接入（MHC-I / MHC-II）。

用途：
1) 检查 real_tsv 文件是否存在、列是否符合要求；
2) 可选自动执行 mhc_ranking；
3) 可选自动执行真实性验收。

约定路径（与 predict_mhc_ranking / run_all 一致）：
- MHC-I：results/<run_id>/tool_outputs/raw/mhc1_netmhcpan.tsv
  必需列：mut_peptide +（mhc1_cv_netmhcpan_nM|affinity_nM|affinity|ic50|score 之一）
- MHC-II：results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv
  必需列：mut_peptide +（mhc2_el_rank|el_rank|pct_rank_el|rank 之一）

多实例说明：
- 仓库内已放齐真实 TSV 的示例包括 R001、R002、R003、R_public_001 等；新建 run 需自行生成或从
  病例工具导出上述文件。若 hla_typing.json 中暂无 II 类键名，仍可用参考等位基因（如 DRB1*15:01）
  的 mhc2 real_tsv 做演示级回填，正式交付须换成患者真实 II 类结果。
"""

import argparse
import os
import subprocess
import sys
from typing import List

import pandas as pd


def _ok(msg: str) -> None:
    print(f"[OK] {msg}")


def _warn(msg: str) -> None:
    print(f"[WARN] {msg}")


def _fail(msg: str) -> None:
    print(f"[FAIL] {msg}")


def _read_tsv(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def _check_columns(df: pd.DataFrame, required: List[str], any_of: List[List[str]], label: str) -> bool:
    ok = True
    for c in required:
        if c not in df.columns:
            _fail(f"{label} 缺少必需列: {c}")
            ok = False
    for group in any_of:
        if not any(c in df.columns for c in group):
            _fail(f"{label} 缺少候选列组之一: {group}")
            ok = False
    return ok


def main() -> None:
    parser = argparse.ArgumentParser(
        description="检查 MHC real_tsv 就绪情况，可选一键跑 mhc_ranking 与验收。",
        epilog=(
            "示例：\n"
            "  python scripts/bootstrap_real_backends.py --run_id R002 --with_mhc2\n"
            "  python scripts/bootstrap_real_backends.py --run_id R002 --with_mhc2 --run --validate\n"
            "说明：--with_mhc2 时才会检查 mhc2_netmhciipan.tsv；与 run_all 的 --mhc2_backend real_tsv 配套。"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--run_id", required=True, help="例如 R001、R002、R_public_001")
    parser.add_argument("--with_mhc2", action="store_true", help="要求同时检查 MHC-II real_tsv")
    parser.add_argument("--run", action="store_true", help="检查通过后自动执行 mhc_ranking")
    parser.add_argument("--validate", action="store_true", help="执行后自动运行 check_epitope_realization")
    args = parser.parse_args()

    run_id = args.run_id
    raw_dir = os.path.join("results", run_id, "tool_outputs", "raw")
    mhc1_path = os.path.join(raw_dir, "mhc1_netmhcpan.tsv")
    mhc2_path = os.path.join(raw_dir, "mhc2_netmhciipan.tsv")

    print(f"[INFO] run_id={run_id}")
    print(f"[INFO] raw_dir={raw_dir}")

    # 1) 检查 MHC-I real_tsv
    ready = True
    if not os.path.exists(mhc1_path):
        _fail(f"MHC-I real_tsv 文件不存在: {mhc1_path}")
        ready = False
    else:
        try:
            df1 = _read_tsv(mhc1_path)
            c_ok = _check_columns(
                df1,
                required=["mut_peptide"],
                any_of=[["mhc1_cv_netmhcpan_nM", "affinity_nM", "affinity", "ic50", "score"]],
                label="mhc1_netmhcpan.tsv",
            )
            if c_ok and not df1.empty:
                _ok(f"MHC-I real_tsv 就绪，行数={len(df1)}")
            elif c_ok and df1.empty:
                _warn("MHC-I real_tsv 存在但为空。")
            ready = ready and c_ok and (not df1.empty)
        except Exception as e:
            _fail(f"读取 MHC-I real_tsv 失败: {e}")
            ready = False

    # 2) 可选检查 MHC-II real_tsv
    mhc2_ready = True
    if args.with_mhc2:
        if not os.path.exists(mhc2_path):
            _fail(f"MHC-II real_tsv 文件不存在: {mhc2_path}")
            mhc2_ready = False
        else:
            try:
                df2 = _read_tsv(mhc2_path)
                c2_ok = _check_columns(
                    df2,
                    required=["mut_peptide"],
                    any_of=[["mhc2_el_rank", "el_rank", "pct_rank_el", "rank"]],
                    label="mhc2_netmhciipan.tsv",
                )
                if c2_ok and not df2.empty:
                    _ok(f"MHC-II real_tsv 就绪，行数={len(df2)}")
                elif c2_ok and df2.empty:
                    _warn("MHC-II real_tsv 存在但为空。")
                mhc2_ready = c2_ok and (not df2.empty)
            except Exception as e:
                _fail(f"读取 MHC-II real_tsv 失败: {e}")
                mhc2_ready = False

    all_ready = ready and (mhc2_ready if args.with_mhc2 else True)
    if not all_ready:
        _fail("真实后端输入未就绪，停止。")
        sys.exit(2)

    _ok("真实后端输入检查通过。")

    # 3) 可选自动执行 mhc_ranking
    if args.run:
        cmd = [
            sys.executable,
            os.path.join("scripts", "run_all.py"),
            "--run_id",
            run_id,
            "--target",
            "mhc_ranking",
            "--backend_mhc1_netmhcpan",
            "real_tsv",
            "--backend_mhc1_bigmhc",
            "off",
            "--require_real_mhc1_cv",
        ]
        if args.with_mhc2:
            cmd += ["--mhc2_backend", "real_tsv", "--require_real_mhc2"]
        _ok("开始执行 mhc_ranking ...")
        print("[CMD]", " ".join(cmd))
        ret = subprocess.run(cmd, check=False)
        if ret.returncode != 0:
            _fail(f"run_all 执行失败，退出码={ret.returncode}")
            sys.exit(ret.returncode)
        _ok("run_all 执行完成。")

    # 4) 可选自动验收
    if args.validate:
        vcmd = [
            sys.executable,
            os.path.join("scripts", "check_epitope_realization.py"),
            "--run_id",
            run_id,
            "--require_mhc1_cv_real",
        ]
        if args.with_mhc2:
            vcmd.append("--require_mhc2_real")
        _ok("开始执行真实性验收 ...")
        print("[CMD]", " ".join(vcmd))
        vret = subprocess.run(vcmd, check=False)
        if vret.returncode != 0:
            _fail(f"验收失败，退出码={vret.returncode}")
            sys.exit(vret.returncode)
        _ok("真实性验收通过。")


if __name__ == "__main__":
    main()
