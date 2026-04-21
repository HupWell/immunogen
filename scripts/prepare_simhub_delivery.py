# -*- coding: utf-8 -*-
"""
功能：生成给 Simulation Hub 的最小交付包（MVP 版）。

输入：
- results/<run_id>/selected_peptides.csv
- deliveries/<run_id>/to_immunogen/meta.json

输出：
- deliveries/<run_id>/to_simhub/<case_id>/protein.pdb        （占位文件）
- deliveries/<run_id>/to_simhub/<case_id>/ligand.sdf         （占位文件）
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


def write_coarse_complex_pdb(path: str, peptide: str):
    """
    生成粗粒度 peptide-MHC 初始复合体（非占位空文件）：
    - 链 P：候选肽 CA 骨架
    - 链 H/I：简化 MHC 槽体两侧轮廓
    """
    atom_id = 1
    lines = [
        "HEADER    COARSE PEPTIDE-MHC INITIAL COMPLEX",
        "REMARK    coarse initial conformer for simulation handoff",
        "REMARK    replace with AlphaFold-Multimer/PANDORA for production",
    ]
    for i, aa in enumerate(peptide, start=1):
        x = (i - 1) * 3.8
        y = 0.0
        z = 0.0
        lines.append(
            f"ATOM  {atom_id:5d}  CA  {residue_name(aa):>3s} P{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C"
        )
        atom_id += 1
    for i in range(1, 51):
        x = (i - 1) * 1.5
        y = 6.0
        z = 1.5 * np.sin(i / 3.0)
        lines.append(
            f"ATOM  {atom_id:5d}  CA  ALA H{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C"
        )
        atom_id += 1
    for i in range(1, 51):
        x = (i - 1) * 1.5
        y = -6.0
        z = 1.5 * np.cos(i / 3.0)
        lines.append(
            f"ATOM  {atom_id:5d}  CA  ALA I{i:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 30.00           C"
        )
        atom_id += 1
    lines.extend(["TER", "END"])
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def write_placeholder_sdf(path: str):
    """写入最小 SDF 占位内容。"""
    content = (
        "placeholder_ligand\n"
        "  ImmunoGenMVP\n"
        "\n"
        "  0  0  0  0  0  0            999 V2000\n"
        "M  END\n"
        "$$$$\n"
    )
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)


def main(run_id: str, top_k: int):
    selected_file = os.path.join("results", run_id, "selected_peptides.csv")
    input_meta_file = os.path.join("deliveries", run_id, "to_immunogen", "meta.json")

    if not os.path.exists(selected_file):
        raise FileNotFoundError(f"未找到输入文件: {selected_file}")

    case_id = load_case_id(input_meta_file, run_id)
    out_dir = os.path.join("deliveries", run_id, "to_simhub", case_id)
    os.makedirs(out_dir, exist_ok=True)

    protein_file = os.path.join(out_dir, "protein.pdb")
    ligand_file = os.path.join(out_dir, "ligand.sdf")
    meta_file = os.path.join(out_dir, "meta.json")
    md_list_file = os.path.join(out_dir, "selected_for_md.csv")

    df = pd.read_csv(selected_file)
    if df.empty:
        raise ValueError("selected_peptides.csv 为空，无法准备 Simulation Hub 交付。")

    md_df = df.head(top_k).copy()
    md_df.to_csv(md_list_file, index=False, encoding="utf-8-sig")

    top_peptide = str(md_df.iloc[0].get("mut_peptide", "AAAAAAAAA")).strip().upper()
    write_coarse_complex_pdb(protein_file, top_peptide[:15] if top_peptide else "AAAAAAAAA")
    write_placeholder_sdf(ligand_file)

    delivery_meta = {
        "run_id": run_id,
        "case_id": case_id,
        "molecule_type": "peptide_mhc",
        "delivery_stage": "coarse_initial_complex",
        "top_k_for_md": int(top_k),
        "selected_for_md_file": "selected_for_md.csv",
        "notes": [
            "protein.pdb 当前为粗粒度初始构象，可用于流程对接，建议替换为真实 peptide-MHC 建模结果。",
            "ligand.sdf 当前为占位文件，用于占位保持接口一致。",
            "建议后续使用 AlphaFold-Multimer 或 PANDORA 构建复合物初始构象。"
        ]
    }
    with open(meta_file, "w", encoding="utf-8") as f:
        json.dump(delivery_meta, f, ensure_ascii=False, indent=2)

    print(f"完成: {protein_file}")
    print(f"完成: {ligand_file}")
    print(f"完成: {meta_file}")
    print(f"完成: {md_list_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_k", type=int, default=3, help="用于 MD 验证的 Top 候选数量，默认 3")
    args = parser.parse_args()
    main(args.run_id, args.top_k)
