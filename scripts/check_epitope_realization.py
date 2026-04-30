# -*- coding: utf-8 -*-
"""
功能：检查某个 run_id 的表位预测是否已使用真实后端（用于验收）。

默认检查：
- MHC-II 至少出现一行 mhc2_backend=netmhciipan
- MHC-I 交叉验证至少一个来源列出现 real_tsv/real_cmd
- 可选检查免疫原性三路来源不含 proxy
"""
import os
import json
import argparse
import pandas as pd


def _source_has_proxy(df: pd.DataFrame, col: str) -> bool:
    if col not in df.columns:
        return True
    return df[col].fillna("").astype(str).str.lower().str.contains("proxy").any()


def _pdb_is_real_structure(path: str) -> bool:
    if not os.path.exists(path):
        return False
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        text = f.read()

    lower = text.lower()
    coarse_markers = (
        "coarse peptide-mhc complex",
        "coarse ca trace",
        "replace with alphafold-multimer or pandora for production md",
    )
    if any(marker in lower for marker in coarse_markers):
        return False

    atom_names = []
    chains = set()
    for raw in text.splitlines():
        if raw.startswith(("ATOM", "HETATM")):
            atom_names.append(raw[12:16].strip())
            if len(raw) > 21 and raw[21].strip():
                chains.add(raw[21].strip())
    atom_count = len(atom_names)
    if not atom_count:
        return False
    ca_ratio = sum(1 for name in atom_names if name == "CA") / atom_count
    return atom_count >= 1000 and ca_ratio <= 0.4 and len(chains) >= 2


def _check_real_structure(run_id: str) -> bool:
    delivery_root = os.path.join("deliveries", run_id, "to_simhub")
    if not os.path.isdir(delivery_root):
        return False
    for case_id in os.listdir(delivery_root):
        case_dir = os.path.join(delivery_root, case_id)
        meta_path = os.path.join(case_dir, "meta.json")
        pdb_path = os.path.join(case_dir, "complex.pdb")
        if not os.path.exists(meta_path):
            continue
        with open(meta_path, "r", encoding="utf-8") as f:
            meta = json.load(f)
        if meta.get("molecule_type") == "peptide_mhc" and _pdb_is_real_structure(pdb_path):
            return True
    return False


def main(
    run_id: str,
    require_mhc2_real: bool,
    require_mhc1_cv_real: bool,
    require_real_immunogenicity: bool,
    require_real_structure: bool,
):
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
    immuno_proxy = {
        "deepimmuno": _source_has_proxy(df, "immunogenicity_source_deepimmuno"),
        "prime": _source_has_proxy(df, "immunogenicity_source_prime"),
        "repitope": _source_has_proxy(df, "immunogenicity_source_repitope"),
    }
    structure_ok = _check_real_structure(run_id)
    print(f"[check] immunogenicity_proxy={immuno_proxy}")
    print(f"[check] real_structure={structure_ok}")

    if require_mhc2_real and not mhc2_ok:
        raise RuntimeError("验收失败：MHC-II 未使用 netmhciipan。")
    if require_mhc1_cv_real and not mhc1_ok:
        raise RuntimeError("验收失败：MHC-I 交叉验证未使用真实来源（real_tsv/real_cmd）。")
    if require_real_immunogenicity and any(immuno_proxy.values()):
        raise RuntimeError("验收失败：免疫原性仍存在 proxy 或缺失来源列。")
    if require_real_structure and not structure_ok:
        raise RuntimeError("验收失败：SimHub 结构仍不是 pandora/afm 真实替换。")

    print("验收通过。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--require_mhc2_real", action="store_true", help="要求 MHC-II 为真实后端")
    parser.add_argument("--require_mhc1_cv_real", action="store_true", help="要求 MHC-I 交叉验证为真实后端")
    parser.add_argument("--require_real_immunogenicity", action="store_true", help="要求 DeepImmuno/PRIME/Repitope 均不含 proxy")
    parser.add_argument("--require_real_structure", action="store_true", help="要求 SimHub 结构为 pandora/afm 真实替换")
    args = parser.parse_args()
    main(
        args.run_id,
        args.require_mhc2_real,
        args.require_mhc1_cv_real,
        args.require_real_immunogenicity,
        args.require_real_structure,
    )

