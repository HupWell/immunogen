# -*- coding: utf-8 -*-
"""
功能：按 Simulation Hub 输入契约（分支 B：peptide_mhc）生成交付包。

输入：
- results/<run_id>/selected_peptides.csv
- deliveries/<run_id>/to_immunogen/meta.json
- deliveries/<run_id>/to_immunogen/hla_typing.json（用于 hla_allele.txt）

输出（禁止使用 SDF）：
- deliveries/<run_id>/to_simhub/<case_id>/complex.pdb   # Chain M(α) + B(β2m) + P(peptide)，真实结构工具输出
- deliveries/<run_id>/to_simhub/<case_id>/hla_allele.txt  # 可选，Top1 对应 HLA
- deliveries/<run_id>/to_simhub/<case_id>/meta.json
- deliveries/<run_id>/to_simhub/<case_id>/selected_for_md.csv
"""
import os
import json
import argparse
from datetime import datetime, timezone, timedelta
from typing import Dict, Tuple
import pandas as pd
import numpy as np


def normalize_chain_ids_for_contract(pdb_text: str) -> Tuple[str, Dict[str, str]]:
    """
    将输入 PDB 的前三个链统一重映射为 M/B/P。
    仅处理 ATOM/HETATM 行；其余行原样保留。
    """
    target = ["M", "B", "P"]
    seen = []
    lines_out = []
    mapping = {}
    for raw in pdb_text.splitlines():
        if raw.startswith("ATOM") or raw.startswith("HETATM"):
            line = raw
            if len(line) < 22:
                line = line.ljust(22)
            chain = line[21].strip() or "_"
            if chain not in seen:
                seen.append(chain)
            if chain in seen[:3]:
                new_chain = target[seen.index(chain)]
                mapping[chain] = new_chain
                line = f"{line[:21]}{new_chain}{line[22:]}"
            lines_out.append(line)
        else:
            lines_out.append(raw)
    if lines_out and lines_out[-1].strip() != "END":
        lines_out.append("END")
    return "\n".join(lines_out) + "\n", mapping


def assert_not_coarse_pdb(pdb_text: str, source_label: str) -> None:
    """拒绝 coarse/CA-only 结构混入生产级 SimHub 交付。"""
    lower = pdb_text.lower()
    coarse_markers = [
        "coarse peptide-mhc complex",
        "coarse ca trace",
        "replace with alphafold-multimer or pandora for production md",
    ]
    matched = [marker for marker in coarse_markers if marker in lower]
    if matched:
        raise ValueError(f"{source_label} 含 coarse 结构标记，不能作为生产 MD 的 complex.pdb：{matched}")

    atom_names = []
    chains = set()
    for raw in pdb_text.splitlines():
        if raw.startswith(("ATOM", "HETATM")):
            atom_names.append(raw[12:16].strip())
            if len(raw) > 21 and raw[21].strip():
                chains.add(raw[21].strip())
    atom_count = len(atom_names)
    ca_count = sum(1 for name in atom_names if name == "CA")
    ca_ratio = ca_count / atom_count if atom_count else 1.0

    if atom_count < 1000 or ca_ratio > 0.4 or len(chains) < 2:
        raise ValueError(
            f"{source_label} 不满足全原子 peptide-MHC 交付要求："
            f"atom_count={atom_count}, ca_ratio={ca_ratio:.3f}, chains={sorted(chains)}"
        )


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


def _value_or_not_provided(value):
    """把空值统一写成 not_provided，便于最终 dossier 显式追踪缺字段。"""
    if value is None:
        return "not_provided"
    try:
        if pd.isna(value):
            return "not_provided"
    except TypeError:
        pass
    text = str(value).strip()
    return text if text else "not_provided"


def _float_or_not_provided(value):
    if value is None:
        return "not_provided"
    try:
        if pd.isna(value):
            return "not_provided"
        return float(value)
    except (TypeError, ValueError):
        return "not_provided"


def _float_or_none(value):
    converted = _float_or_not_provided(value)
    return converted if isinstance(converted, float) else None


def _priority_from_score(score) -> str:
    """按上游综合分给 SimHub 一个稳定的优先级枚举。"""
    if score is None:
        return "medium"
    if score >= 0.7:
        return "high"
    if score >= 0.3:
        return "medium"
    return "low"


def _target_uniprot_for_hla(allele: str):
    """按 HLA-I 基因给出代表性 UniProt 条目。"""
    text = str(allele or "").upper()
    if text.startswith("HLA-A"):
        return "P01892"
    if text.startswith("HLA-B"):
        return "P01889"
    if text.startswith("HLA-C"):
        return "P04222"
    return "not_provided"


def _load_json_optional(path: str) -> Dict:
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def build_simhub_meta(
    run_id: str,
    case_id: str,
    patient_id: str,
    top_row: pd.Series,
    allele: str,
    backend: str,
    structure_input_pdb: str,
) -> Dict:
    """
    生成 SimHub SINGLE SOURCE OF TRUTH 允许的 meta.json 字段。

    注意：不要把结构调试字段写入 meta 顶层；未声明字段会导致 SimHub 拒收。
    """
    score = _float_or_none(top_row.get("rank_score"))
    target_name = f"{allele} peptide-MHC complex" if allele else "peptide-MHC complex"
    if backend == "coarse":
        structure_note = "Current complex.pdb is coarse and only suitable for connectivity testing; production MD requires AF-Multimer/PANDORA all-atom complex."
    else:
        structure_note = f"complex.pdb was prepared from {backend} structure input: {structure_input_pdb}."
    return {
        "case_id": case_id,
        "run_id": run_id,
        "patient_id": patient_id,
        "upstream_module": "ImmunoGen",
        "module_version": "ImmunoGen-v0.1.0",
        "stage_versions": {
            "immunogen": "v0.1.0",
        },
        "timestamp": datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds"),
        "molecule_type": "peptide_mhc",
        "target_name": target_name,
        "target_uniprot": _target_uniprot_for_hla(allele),
        "target_evidence_level": "medium",
        "upstream_score": score,
        "priority": _priority_from_score(score),
        "notes": [
            "Input package follows SimHub peptide_mhc contract.",
            "For peptide_mhc, parent_hit_id/scaffold_smiles/generation_source are kept for schema compatibility.",
            "No SDF is provided; peptide-MHC must use protein residue topology.",
            structure_note,
        ],
        "parent_hit_id": None,
        "scaffold_smiles": None,
        "generation_source": "de_novo",
        "carrier": None,
    }


def write_dossier_context(
    path: str,
    run_id: str,
    case_id: str,
    top_row: pd.Series,
    allele: str,
) -> None:
    """写出 SimHub 最终归档需要的 peptide + mRNA 上下文 manifest。"""
    candidate_id = f"{run_id}:{case_id}:{_value_or_not_provided(top_row.get('mut_peptide'))}"
    context = {
        "candidate_id": candidate_id,
        "therapy_line": "mRNA_vaccine",
        "peptide": {
            "mut_peptide": _value_or_not_provided(top_row.get("mut_peptide")),
            "wt_peptide": _value_or_not_provided(top_row.get("wt_peptide")),
            "hla_allele": _value_or_not_provided(allele),
            "mutation": _value_or_not_provided(top_row.get("mutation")),
            "variant_vaf": _float_or_not_provided(top_row.get("variant_vaf")),
        },
        "scores": {
            "affinity_nM": _float_or_not_provided(top_row.get("affinity_nM")),
            "immunogenicity": _float_or_not_provided(top_row.get("immunogenicity")),
            "wt_dissimilarity": _float_or_not_provided(
                top_row.get("dissimilarity", top_row.get("wt_peptide_dissimilarity"))
            ),
            "rank_score": _float_or_not_provided(top_row.get("rank_score")),
        },
        "mrna_materials": {
            "mrna_fasta": f"results/{run_id}/mrna_vaccine.fasta",
            "utr_orf_polyA_design": f"results/{run_id}/mrna_design.json",
            "modification_strategy": "N1-methylpseudouridine",
        },
        "target_rationale": (
            "Top peptide-MHC candidate selected by ImmunoGen rank_score; "
            "mRNA materials point to the multivalent vaccine design for this run."
        ),
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(context, f, ensure_ascii=False, indent=2)


def write_simhub_evidence_contract(path: str, run_id: str, case_id: str) -> None:
    """创建 SimHub 回传证据目录说明；真实回传到位后由该目录保留引用副本。"""
    os.makedirs(path, exist_ok=True)
    readme = os.path.join(path, "README.md")
    content = f"""# SimHub Evidence（{run_id}/{case_id}）

当前状态：`not_returned`

SimHub 回传后，请在本目录保留引用副本或软链接：

- `energy_report.json`
- `rmsd_profile.csv`
- `qc_flags.json`
- `summary.md`

对外交付仍以 SimHub 标准归档 `to_hospital/{run_id}/` 为准；本目录用于 ImmunoGen 自证引用。
"""
    with open(readme, "w", encoding="utf-8") as f:
        f.write(content)


def write_run_meta(
    path: str,
    run_id: str,
    case_id: str,
    selected_count: int,
    simhub_dir: str,
    evidence_dir: str,
) -> None:
    """写出 results/<run_id>/meta.json，与 SimHub meta.json 分工独立。"""
    data = {
        "run_id": run_id,
        "case_id": case_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "module": "ImmunoGen",
        "therapy_line": "mRNA_vaccine",
        "selected_peptide_count": int(selected_count),
        "outputs": {
            "ranking_csv": f"results/{run_id}/peptide_mhc_ranking.csv",
            "selected_peptides_csv": f"results/{run_id}/selected_peptides.csv",
            "mrna_fasta": f"results/{run_id}/mrna_vaccine.fasta",
            "mrna_design_json": f"results/{run_id}/mrna_design.json",
            "qc_metrics_json": f"results/{run_id}/qc_metrics.json",
            "report_md": f"results/{run_id}/REPORT.md",
            "self_check_md": f"results/{run_id}/SELF_CHECK.md",
            "simhub_delivery_dir": simhub_dir,
            "simhub_evidence_dir": evidence_dir,
        },
        "clinical_final_note": "R002/R003/R_public_001 当前仍需等待 BioDriver 患者真实 HLA-II 分型后形成临床终版。",
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)


def main(
    run_id: str,
    top_k: int,
    structure_backend: str,
    structure_input_pdb: str,
    allow_coarse_connectivity: bool,
):
    selected_file = os.path.join("results", run_id, "selected_peptides.csv")
    input_meta_file = os.path.join("deliveries", run_id, "to_immunogen", "meta.json")
    hla_file = os.path.join("deliveries", run_id, "to_immunogen", "hla_typing.json")

    if not os.path.exists(selected_file):
        raise FileNotFoundError(f"未找到输入文件: {selected_file}")

    case_id = load_case_id(input_meta_file, run_id)
    input_meta = _load_json_optional(input_meta_file)
    patient_id = str(input_meta.get("patient_id") or input_meta.get("subject_id") or run_id)
    out_dir = os.path.join("deliveries", run_id, "to_simhub", case_id)
    os.makedirs(out_dir, exist_ok=True)

    complex_pdb = os.path.join(out_dir, "complex.pdb")
    hla_txt = os.path.join(out_dir, "hla_allele.txt")
    meta_out = os.path.join(out_dir, "meta.json")
    md_list_file = os.path.join(out_dir, "selected_for_md.csv")
    dossier_context_file = os.path.join(out_dir, "dossier_context.json")

    df = pd.read_csv(selected_file)
    if df.empty:
        raise ValueError("selected_peptides.csv 为空，无法准备 Simulation Hub 交付。")

    md_df = df.head(top_k).copy()
    md_df.to_csv(md_list_file, index=False, encoding="utf-8-sig")

    top_row = md_df.iloc[0]
    top_peptide = str(top_row.get("mut_peptide", "AAAAAAAAA")).strip().upper()
    backend = structure_backend.strip().lower()
    chain_map = {}
    if backend == "coarse":
        if not allow_coarse_connectivity:
            raise ValueError("coarse 仅允许联通测试；生产交付请使用 pandora/afm，并提供 --structure_input_pdb。")
        write_complex_pdb(complex_pdb, top_peptide[:15] if top_peptide else "AAAAAAAAA")
        chain_map = {"M": "M", "B": "B", "P": "P"}
    elif backend in ("pandora", "afm"):
        if not structure_input_pdb.strip():
            raise ValueError("structure_backend 为 pandora/afm 时，必须提供 --structure_input_pdb。")
        src = structure_input_pdb.strip()
        if not os.path.exists(src):
            raise FileNotFoundError(f"未找到结构文件: {src}")
        with open(src, "r", encoding="utf-8", errors="ignore") as f:
            pdb_text = f.read()
        assert_not_coarse_pdb(pdb_text, src)
        normalized, chain_map = normalize_chain_ids_for_contract(pdb_text)
        with open(complex_pdb, "w", encoding="utf-8") as f:
            f.write(normalized)
    else:
        raise ValueError("structure_backend 仅支持 coarse/pandora/afm。")

    allele = pick_hla_allele(top_row, hla_file)
    with open(hla_txt, "w", encoding="utf-8") as f:
        f.write(allele + "\n")
    write_dossier_context(dossier_context_file, run_id, case_id, top_row, allele)

    delivery_meta = build_simhub_meta(
        run_id=run_id,
        case_id=case_id,
        patient_id=patient_id,
        top_row=top_row,
        allele=allele,
        backend=backend,
        structure_input_pdb=structure_input_pdb,
    )
    with open(meta_out, "w", encoding="utf-8") as f:
        json.dump(delivery_meta, f, ensure_ascii=False, indent=2)

    evidence_dir = os.path.join("results", run_id, "simhub_evidence", case_id)
    write_simhub_evidence_contract(evidence_dir, run_id, case_id)
    write_run_meta(
        os.path.join("results", run_id, "meta.json"),
        run_id,
        case_id,
        len(df),
        out_dir,
        evidence_dir,
    )

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
    print(f"完成: {dossier_context_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_k", type=int, default=3, help="用于 MD 验证的 Top 候选数量，默认 3")
    parser.add_argument(
        "--structure_backend",
        default="pandora",
        choices=["coarse", "pandora", "afm"],
        help="结构来源：pandora(默认)/afm；coarse 仅可显式用于联通测试",
    )
    parser.add_argument(
        "--structure_input_pdb",
        default="",
        help="structure_backend 为 pandora/afm 时，提供外部 PDB 路径。",
    )
    parser.add_argument(
        "--allow_coarse_connectivity",
        action="store_true",
        help="仅用于接口联通测试；生产 MD 交付禁止使用。",
    )
    args = parser.parse_args()
    main(
        args.run_id,
        args.top_k,
        args.structure_backend,
        args.structure_input_pdb,
        args.allow_coarse_connectivity,
    )
