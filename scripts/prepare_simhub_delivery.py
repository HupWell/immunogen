# -*- coding: utf-8 -*-
"""
功能：按 Simulation Hub 输入契约（分支 B：peptide_mhc）生成交付包。

输入：
- results/<run_id>/selected_peptides.csv
- deliveries/<run_id>/to_immunogen/meta.json
- deliveries/<run_id>/to_immunogen/hla_typing.json（用于 hla_allele.txt）

输出（禁止使用 SDF）：
- deliveries/<run_id>/to_simhub/<case_id>/complex.pdb   # Chain M(α) + B(β2m) + P(peptide)，粗初始构象
- deliveries/<run_id>/to_simhub/<case_id>/hla_allele.txt  # 可选，Top1 对应 HLA
- deliveries/<run_id>/to_simhub/<case_id>/meta.json
- deliveries/<run_id>/to_simhub/<case_id>/selected_for_md.csv
"""
import os
import json
import argparse
import pandas as pd
import numpy as np


def load_case_id(meta_file: str, run_id: str) -> str:
    """优先从 meta.json 读取 case_id，缺失时回退为 run_id。"""
    if not os.path.exists(meta_file):
        return run_id
    with open(meta_file, "r", encoding="utf-8") as f:
        data = json.load(f)
    case_id = str(data.get("case_id", "")).strip()
    return case_id if case_id else run_id


def residue_name(aa: str) -> str:
    """单字母转 PDB 三字母残基名。"""
    mapping = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS", "Q": "GLN", "E": "GLU",
        "G": "GLY", "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE",
        "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
    }
    return mapping.get(aa.upper(), "GLY")


def write_complex_pdb(path: str, peptide: str):
    """
    生成粗粒度 peptide-MHC 初始 complex.pdb（SimHub 契约：多链纯蛋白）。
    - Chain M：简化 MHC-I α 槽体一侧（CA 骨架）
    - Chain B：简化 β2m 轮廓（CA 骨架）
    - Chain P：候选肽 CA 骨架
    生产环境请替换为 AlphaFold-Multimer / PANDORA 输出。
    """
    atom_id = 1
    lines = [
        "HEADER    COARSE PEPTIDE-MHC COMPLEX (SimHub handoff)",
        "REMARK    Chain M = MHC class I alpha (coarse CA trace)",
        "REMARK    Chain B = beta-2-microglobulin proxy (coarse CA trace)",
        "REMARK    Chain P = neoantigen peptide (CA trace)",
        "REMARK    Replace with AlphaFold-Multimer or PANDORA for production MD",
    ]
    # Chain M：槽体一侧
    for i in range(1, 51):
        x = (i - 1) * 1.5
        y = 6.0
        z = 1.5 * np.sin(i / 3.0)
        lines.append(
            f"ATOM  {atom_id:5d}  CA  ALA M{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C"
        )
        atom_id += 1
    lines.append("TER")
    # Chain B：β2m 代理
    for i in range(1, 51):
        x = (i - 1) * 1.5
        y = -6.0
        z = 1.5 * np.cos(i / 3.0)
        lines.append(
            f"ATOM  {atom_id:5d}  CA  ALA B{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C"
        )
        atom_id += 1
    lines.append("TER")
    # Chain P：肽段
    pep = (peptide or "AAAAAAAAA")[:15]
    for i, aa in enumerate(pep, start=1):
        x = (i - 1) * 3.8
        y = 0.0
        z = 0.0
        lines.append(
            f"ATOM  {atom_id:5d}  CA  {residue_name(aa):>3s} P{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        atom_id += 1
    lines.extend(["TER", "END"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def pick_hla_allele(md_row: pd.Series, hla_json_path: str) -> str:
    """从选中行或 hla_typing.json 取一个等位基因字符串。"""
    if "hla_allele" in md_row.index and pd.notna(md_row.get("hla_allele")):
        return str(md_row["hla_allele"]).strip()
    if os.path.exists(hla_json_path):
        with open(hla_json_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        for key in ("HLA-A", "HLA-B", "HLA-C"):
            vals = data.get(key) or []
            if isinstance(vals, list) and vals:
                return str(vals[0]).strip()
    return "HLA-A*02:01"


def main(run_id: str, top_k: int):
    selected_file = os.path.join("results", run_id, "selected_peptides.csv")
    input_meta_file = os.path.join("deliveries", run_id, "to_immunogen", "meta.json")
    hla_file = os.path.join("deliveries", run_id, "to_immunogen", "hla_typing.json")

    if not os.path.exists(selected_file):
        raise FileNotFoundError(f"未找到输入文件: {selected_file}")

    case_id = load_case_id(input_meta_file, run_id)
    out_dir = os.path.join("deliveries", run_id, "to_simhub", case_id)
    os.makedirs(out_dir, exist_ok=True)

    complex_pdb = os.path.join(out_dir, "complex.pdb")
    hla_txt = os.path.join(out_dir, "hla_allele.txt")
    meta_out = os.path.join(out_dir, "meta.json")
    md_list_file = os.path.join(out_dir, "selected_for_md.csv")

    df = pd.read_csv(selected_file)
    if df.empty:
        raise ValueError("selected_peptides.csv 为空，无法准备 Simulation Hub 交付。")

    md_df = df.head(top_k).copy()
    md_df.to_csv(md_list_file, index=False, encoding="utf-8-sig")

    top_row = md_df.iloc[0]
    top_peptide = str(top_row.get("mut_peptide", "AAAAAAAAA")).strip().upper()
    write_complex_pdb(complex_pdb, top_peptide[:15] if top_peptide else "AAAAAAAAA")

    allele = pick_hla_allele(top_row, hla_file)
    with open(hla_txt, "w", encoding="utf-8") as f:
        f.write(allele + "\n")

    delivery_meta = {
        "run_id": run_id,
        "case_id": case_id,
        "molecule_type": "peptide_mhc",
        "delivery_stage": "coarse_initial_complex",
        "contract": "SimHub peptide_mhc branch B: complex.pdb + meta.json + optional hla_allele.txt; no SDF",
        "top_k_for_md": int(top_k),
        "selected_for_md_file": "selected_for_md.csv",
        "notes": [
            "complex.pdb 为粗粒度 M/B/P 三链 CA 初始构象，仅用于流程联通与接口验收。",
            "禁止使用 SDF + AM1-BCC 处理多肽；多肽须保留氨基酸残基拓扑，走 AMBER14SB 蛋白力场。",
            "生产环境请用 AlphaFold-Multimer 或 PANDORA 生成高精度 complex.pdb 后再送 MD。"
        ]
    }
    with open(meta_out, "w", encoding="utf-8") as f:
        json.dump(delivery_meta, f, ensure_ascii=False, indent=2)

    # 移除旧版文件名（若存在），避免与契约混淆
    for legacy in ("protein.pdb", "ligand.sdf"):
        legacy_path = os.path.join(out_dir, legacy)
        if os.path.exists(legacy_path):
            try:
                os.remove(legacy_path)
            except OSError:
                pass

    print(f"完成: {complex_pdb}")
    print(f"完成: {hla_txt}")
    print(f"完成: {meta_out}")
    print(f"完成: {md_list_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_k", type=int, default=3, help="用于 MD 验证的 Top 候选数量，默认 3")
    args = parser.parse_args()
    main(args.run_id, args.top_k)
