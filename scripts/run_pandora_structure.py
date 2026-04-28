#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
功能：调用真实 PANDORA + MODELLER 生成 peptide-MHC-I 复合物结构。

说明：
- 只接受真实 PANDORA 模板数据库，不生成 coarse/占位结构。
- 默认读取 results/<run_id>/selected_peptides.csv 的 Top1 peptide/HLA。
- 输出 PDB 可继续交给 prepare_simhub_delivery.py --structure_backend pandora。
"""
import argparse
import json
import os
import shutil
from pathlib import Path
from typing import Dict, Iterable, Optional

import pandas as pd


def _load_case_id(run_id: str) -> str:
    meta_path = Path("deliveries") / run_id / "to_immunogen" / "meta.json"
    if not meta_path.exists():
        return run_id
    with meta_path.open("r", encoding="utf-8") as f:
        data = json.load(f)
    case_id = str(data.get("case_id", "")).strip()
    return case_id or run_id


def _load_top_target(run_id: str) -> Dict[str, str]:
    selected_path = Path("results") / run_id / "selected_peptides.csv"
    if not selected_path.exists():
        raise FileNotFoundError(f"未找到入选肽文件: {selected_path}")
    df = pd.read_csv(selected_path)
    if df.empty:
        raise ValueError(f"{selected_path} 为空，无法运行 PANDORA。")
    row = df.iloc[0]
    peptide = str(row.get("mut_peptide", "")).strip().upper()
    allele = str(row.get("hla_allele", "")).strip()
    if not peptide:
        raise ValueError("selected_peptides.csv 缺少 mut_peptide。")
    if not allele:
        raise ValueError("selected_peptides.csv 缺少 hla_allele，无法选择 PANDORA 模板。")
    return {"peptide": peptide, "allele": allele}


def _candidate_pdb_paths(pdb_root: Path, template_id: str, mhc_class: str) -> Iterable[Path]:
    if mhc_class == "I":
        yield pdb_root / "pMHCI" / f"{template_id}.pdb"
    else:
        yield pdb_root / "pMHCII" / f"{template_id}.pdb"
        yield pdb_root / "pMHCII_reversed" / f"{template_id}_reverse.pdb"


def _attach_template_paths(db, pdb_root: Path) -> None:
    """兼容新版数据库对象：为每个模板补齐当前 csb-pandora 需要的 pdb_path 字段。"""
    for template_id, template in db.MHCI_data.items():
        for path in _candidate_pdb_paths(pdb_root, template_id, "I"):
            if path.exists():
                template.pdb_path = str(path)
                break
        if not hasattr(template, "pdb_path"):
            raise FileNotFoundError(f"未找到 MHC-I 模板 PDB: {template_id}")
    for template_id, template in db.MHCII_data.items():
        for path in _candidate_pdb_paths(pdb_root, template_id, "II"):
            if path.exists():
                template.pdb_path = str(path)
                break


def _load_database(database_path: Path, pdb_root: Path):
    from PANDORA.Database import Database

    if not database_path.exists():
        raise FileNotFoundError(f"未找到 PANDORA 数据库: {database_path}")
    if not pdb_root.exists():
        raise FileNotFoundError(f"未找到 PANDORA PDB 模板目录: {pdb_root}")
    db = Database.load(str(database_path))
    _attach_template_paths(db, pdb_root)
    return db


def _match_mhci_sequence(db, allele: str) -> str:
    """从 PANDORA 官方数据库参考序列中匹配目标 HLA-I 序列。"""
    refs = getattr(db, "ref_MHCI_sequences", {}) or {}
    if allele in refs:
        return refs[allele]
    prefix = allele if allele.startswith("HLA-") else f"HLA-{allele}"
    matches = sorted(k for k in refs if k.startswith(prefix + ":"))
    if matches:
        return refs[matches[0]]
    raise ValueError(f"PANDORA 数据库中未找到 HLA-I 参考序列: {allele}")


def _write_mhci_template_without_b2m(src_pdb: Path, dst_pdb: Path) -> None:
    """旧版 PANDORA 建模 alignment 不含 B2M；建模模板需临时移除 B 链。"""
    dst_pdb.parent.mkdir(parents=True, exist_ok=True)
    with src_pdb.open("r", encoding="utf-8", errors="ignore") as src, dst_pdb.open("w", encoding="utf-8") as dst:
        for raw in src:
            if raw.startswith(("ATOM", "HETATM", "TER")) and len(raw) > 21 and raw[21] == "B":
                continue
            dst.write(raw)


def _pdb_has_chain(path: Path, chain_id: str) -> bool:
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if raw.startswith(("ATOM", "HETATM")) and len(raw) > 21 and raw[21] == chain_id:
                return True
    return False


def _split_atom_lines_by_chain(path: Path) -> Dict[str, list]:
    chains: Dict[str, list] = {}
    with path.open("r", encoding="utf-8", errors="ignore") as f:
        for raw in f:
            if raw.startswith(("ATOM", "HETATM")) and len(raw) > 21:
                chain = raw[21].strip() or "_"
                chains.setdefault(chain, []).append(raw.rstrip("\n"))
    return chains


def _write_complex_with_template_b2m(model_pdb: Path, template_pdb: Path, final_pdb: Path) -> bool:
    """
    若旧版 PANDORA 输出缺少 B2M，则从真实模板结构 graft B 链。
    返回 True 表示执行过 B2M graft。
    """
    if _pdb_has_chain(model_pdb, "B"):
        shutil.copyfile(model_pdb, final_pdb)
        return False
    model_chains = _split_atom_lines_by_chain(model_pdb)
    template_chains = _split_atom_lines_by_chain(template_pdb)
    if "B" not in template_chains:
        shutil.copyfile(model_pdb, final_pdb)
        return False
    with final_pdb.open("w", encoding="utf-8") as out:
        out.write("REMARK    PANDORA model with B2M chain grafted from selected real template\n")
        for chain in ("M", "B", "P"):
            lines = template_chains["B"] if chain == "B" else model_chains.get(chain, [])
            for line in lines:
                out.write(line + "\n")
            if lines:
                out.write("TER\n")
        out.write("END\n")
    return True


def _find_model_pdb(output_dir: Path, target_id: str) -> Path:
    candidates = []
    for path in output_dir.rglob("*.pdb"):
        name = path.name
        if name.endswith(".ini") or name.startswith("."):
            continue
        if target_id in name:
            candidates.append(path)
    if not candidates:
        for path in output_dir.rglob("*.pdb"):
            if not path.name.startswith(".") and not path.name.endswith(".ini"):
                candidates.append(path)
    if not candidates:
        raise FileNotFoundError(f"PANDORA 未产出可用 PDB，输出目录: {output_dir}")
    return max(candidates, key=lambda p: p.stat().st_size)


def main(
    run_id: str,
    database_path: str,
    pdb_root: str,
    output_root: str,
    n_loop_models: int,
    n_jobs: Optional[int],
) -> None:
    from PANDORA.PMHC import PMHC
    from PANDORA.Pandora.Pandora import Pandora

    # PANDORA 0.9 调用 MUSCLE 旧版 -in/-out 参数；项目内包装器会转成 MUSCLE 5 的真实参数。
    wrapper_dir = Path("tools") / "pandora_bin"
    if wrapper_dir.exists():
        os.environ["PATH"] = f"{wrapper_dir.resolve()}{os.pathsep}{os.environ.get('PATH', '')}"

    case_id = _load_case_id(run_id)
    top = _load_top_target(run_id)
    target_id = f"{run_id}_{case_id}_top1"
    out_dir = (Path(output_root) / run_id / case_id).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    try:
        import PANDORA

        Path(PANDORA.PANDORA_data, "outputs").mkdir(parents=True, exist_ok=True)
    except Exception:
        pass

    db = _load_database(Path(database_path), Path(pdb_root))
    mhc_sequence = _match_mhci_sequence(db, top["allele"])
    target = PMHC.Target(
        target_id,
        allele_type=[top["allele"]],
        peptide=top["peptide"],
        MHC_class="I",
        M_chain_seq=mhc_sequence,
        anchors=[],
        use_netmhcpan=False,
    )
    selector = Pandora(target, database=db, output_dir=str(out_dir))
    selector.find_template(verbose=True)
    selected_template = selector.template
    original_template_pdb = Path(selected_template.pdb_path).resolve()
    compat_template_pdb = out_dir / "template_compat" / f"{selected_template.id}.pdb"
    _write_mhci_template_without_b2m(original_template_pdb, compat_template_pdb)
    selected_template.pdb_path = str(compat_template_pdb)

    runner = Pandora(target, template=selected_template, output_dir=str(out_dir))
    runner.model(
        n_loop_models=n_loop_models,
        n_homology_models=1,
        n_jobs=n_jobs,
        loop_refinement="slow",
        pickle_out=True,
        verbose=True,
    )

    model_pdb = _find_model_pdb(out_dir, target_id)
    final_pdb = out_dir / "complex.pdb"
    b2m_grafted = _write_complex_with_template_b2m(model_pdb, original_template_pdb, final_pdb)
    meta = {
        "run_id": run_id,
        "case_id": case_id,
        "peptide": top["peptide"],
        "hla_allele": top["allele"],
        "backend": "pandora",
        "database_path": str(Path(database_path).resolve()),
        "pdb_root": str(Path(pdb_root).resolve()),
        "n_loop_models": n_loop_models,
        "template": getattr(runner.template, "id", ""),
        "template_pdb": str(original_template_pdb),
        "template_b2m_grafted": b2m_grafted,
        "raw_model_pdb": str(model_pdb),
        "final_pdb": str(final_pdb),
    }
    with (out_dir / "pandora_structure_meta.json").open("w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    print(f"完成真实 PANDORA 结构: {final_pdb}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument(
        "--database_path",
        default="external_refs/PANDORA_database_default/default/database/PANDORA_database.pkl",
        help="PANDORA 官方数据库 pkl 路径。",
    )
    parser.add_argument(
        "--pdb_root",
        default="external_refs/PANDORA_database_default/default/PDBs",
        help="PANDORA 模板 PDB 根目录，需包含 pMHCI/pMHCII。",
    )
    parser.add_argument(
        "--output_root",
        default="results/structure_models/pandora",
        help="真实结构输出根目录。",
    )
    parser.add_argument("--n_loop_models", type=int, default=5, help="MODELLER loop 模型数，默认 5。")
    parser.add_argument("--n_jobs", type=int, default=None, help="MODELLER 并行任务数，默认由 PANDORA 决定。")
    args = parser.parse_args()
    main(
        run_id=args.run_id,
        database_path=args.database_path,
        pdb_root=args.pdb_root,
        output_root=args.output_root,
        n_loop_models=args.n_loop_models,
        n_jobs=args.n_jobs,
    )
