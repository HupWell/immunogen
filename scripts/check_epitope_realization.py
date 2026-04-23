# -*- coding: utf-8 -*-
"""
功能：检查某个 run_id 的表位预测是否已使用真实后端（用于验收）。

默认检查：
- MHC-II 至少出现一行 mhc2_backend=netmhciipan
- MHC-I 交叉验证至少一个来源列出现 real_tsv/real_cmd
"""
import os
import argparse
import pandas as pd


def main(run_id: str, require_mhc2_real: bool, require_mhc1_cv_real: bool):
    path = os.path.join("results", run_id, "peptide_mhc_ranking.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError("peptide_mhc_ranking.csv 为空。")

    mhc2_ok = False
    if "mhc2_backend" in df.columns:
        mhc2_ok = (df["mhc2_backend"].fillna("").astype(str) == "netmhciipan").any()

    mhc1_ok = False
    real_tags = {"real_tsv", "real_cmd"}
    nsrc = set(df.get("mhc1_cv_source_netmhcpan", pd.Series(dtype=str)).fillna("").astype(str).tolist())
    bsrc = set(df.get("mhc1_cv_source_bigmhc", pd.Series(dtype=str)).fillna("").astype(str).tolist())
    if nsrc.intersection(real_tags) or bsrc.intersection(real_tags):
        mhc1_ok = True

    print(f"[check] run_id={run_id}")
    print(f"[check] mhc2_real={mhc2_ok}")
    print(f"[check] mhc1_cv_real={mhc1_ok}")
    print(f"[check] mhc1_cv_source_netmhcpan={sorted(x for x in nsrc if x)}")
    print(f"[check] mhc1_cv_source_bigmhc={sorted(x for x in bsrc if x)}")
    if "mhc1_cv_source" in df.columns:
        src = sorted(set(df["mhc1_cv_source"].fillna("").astype(str).tolist()) - {""})
        print(f"[check] mhc1_cv_source={src}")
    if "mhc1_cv_tool" in df.columns:
        tools = sorted(set(",".join(df["mhc1_cv_tool"].fillna("").astype(str).tolist()).split(",")) - {""})
        print(f"[check] mhc1_cv_tool={tools}")

    if require_mhc2_real and not mhc2_ok:
        raise RuntimeError("验收失败：MHC-II 未使用 netmhciipan。")
    if require_mhc1_cv_real and not mhc1_ok:
        raise RuntimeError("验收失败：MHC-I 交叉验证未使用真实来源（real_tsv/real_cmd）。")

    print("验收通过。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--require_mhc2_real", action="store_true", help="要求 MHC-II 为真实后端")
    parser.add_argument("--require_mhc1_cv_real", action="store_true", help="要求 MHC-I 交叉验证为真实后端")
    args = parser.parse_args()
    main(args.run_id, args.require_mhc2_real, args.require_mhc1_cv_real)

