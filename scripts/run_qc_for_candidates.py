# -*- coding: utf-8 -*-
"""
功能：对同一 run_id 下的多条 mRNA 候选逐条执行 QC，并各自产生独立报告归档。

设计原则：
1) 复用现有 `run_qc_and_report.py`（不改动其主逻辑）；
2) 候选执行时临时覆盖默认文件名 `mrna_vaccine.fasta` / `mrna_design.json`；
3) 每条候选执行完成后，将产物归档到 `results/<run_id>/mrna_candidate_qc/<candidate>/`；
4) 最后恢复 run 目录原始默认文件，避免破坏已有主线状态。

输出：
- results/<run_id>/mrna_candidate_qc/<candidate>/REPORT.md
- results/<run_id>/mrna_candidate_qc/<candidate>/qc_metrics.json
- results/<run_id>/mrna_candidate_qc/<candidate>/mrna_design.json
- results/<run_id>/mrna_candidate_qc/<candidate>/mrna_vaccine.fasta
- results/<run_id>/mrna_candidate_qc/<candidate>/figures/*
- results/<run_id>/mrna_candidate_qc/qc_summary.md
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from typing import Any, Dict, List, Tuple


def _candidate_name_from_design(path: str, design: Dict[str, Any]) -> str:
    """从设计文件推断候选名：优先 mrna_output_suffix，其次文件名后缀。"""
    suffix = str(design.get("mrna_output_suffix", "") or "").strip()
    if suffix:
        return suffix
    base = os.path.basename(path)
    if base == "mrna_design.json":
        return "default"
    prefix = "mrna_design_"
    if base.startswith(prefix) and base.endswith(".json"):
        return base[len(prefix) : -len(".json")]
    return base.replace(".json", "")


def _collect_candidates(run_dir: str) -> List[Tuple[str, str, Dict[str, Any], str]]:
    """
    收集候选：
    返回列表元素：(候选名, design_path, design_json_dict, fasta_path)
    """
    paths: List[str] = []
    default_design = os.path.join(run_dir, "mrna_design.json")
    if os.path.isfile(default_design):
        paths.append(default_design)
    paths.extend(sorted(glob.glob(os.path.join(run_dir, "mrna_design_*.json"))))

    uniq: List[str] = []
    seen = set()
    for p in paths:
        if p not in seen:
            seen.add(p)
            uniq.append(p)

    out: List[Tuple[str, str, Dict[str, Any], str]] = []
    for design_path in uniq:
        with open(design_path, "r", encoding="utf-8") as f:
            d = json.load(f)
        c_name = _candidate_name_from_design(design_path, d)
        fasta_path = str((d.get("output_files") or {}).get("fasta", "")).strip()
        if not fasta_path:
            # 兜底：按命名约定推断
            if c_name == "default":
                fasta_path = os.path.join(run_dir, "mrna_vaccine.fasta")
            else:
                fasta_path = os.path.join(run_dir, f"mrna_vaccine_{c_name}.fasta")
        if not os.path.isabs(fasta_path):
            fasta_path = os.path.join(os.getcwd(), fasta_path)
        out.append((c_name, design_path, d, fasta_path))
    return out


def _copy_if_exists(src: str, dst: str) -> None:
    if os.path.isfile(src):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)


def _to_float(v: Any) -> float | None:
    """安全转浮点，用于计算候选差异。"""
    try:
        if v is None or v == "":
            return None
        return float(v)
    except (TypeError, ValueError):
        return None


def main() -> None:
    parser = argparse.ArgumentParser(description="对多条 mRNA 候选逐条运行 QC 并归档报告")
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument(
        "--candidates",
        default="",
        help="可选：指定候选名列表（逗号分隔，如 default,v_opt,v_ld_real）；留空则自动扫描全部",
    )
    args = parser.parse_args()

    run_id = args.run_id.strip()
    run_dir = os.path.join("results", run_id)
    if not os.path.isdir(run_dir):
        raise FileNotFoundError(f"未找到 run 目录: {run_dir}")

    ranking_file = os.path.join(run_dir, "peptide_mhc_ranking.csv")
    selected_file = os.path.join(run_dir, "selected_peptides.csv")
    if not os.path.isfile(ranking_file) or not os.path.isfile(selected_file):
        raise FileNotFoundError("缺少 QC 必需输入：peptide_mhc_ranking.csv 或 selected_peptides.csv")

    candidates = _collect_candidates(run_dir)
    if not candidates:
        raise FileNotFoundError(f"{run_dir} 下未找到 mrna_design.json / mrna_design_*.json")

    requested = [x.strip() for x in args.candidates.split(",") if x.strip()]
    if requested:
        req_set = set(requested)
        candidates = [c for c in candidates if c[0] in req_set]
        if not candidates:
            raise ValueError(f"未匹配到候选名: {requested}")

    archive_root = os.path.join(run_dir, "mrna_candidate_qc")
    os.makedirs(archive_root, exist_ok=True)

    # 备份默认文件，处理完成后恢复，避免污染当前主线。
    default_fasta = os.path.join(run_dir, "mrna_vaccine.fasta")
    default_design = os.path.join(run_dir, "mrna_design.json")
    backup_fasta = os.path.join(archive_root, "_backup_mrna_vaccine.fasta")
    backup_design = os.path.join(archive_root, "_backup_mrna_design.json")
    had_default_fasta = os.path.isfile(default_fasta)
    had_default_design = os.path.isfile(default_design)
    _copy_if_exists(default_fasta, backup_fasta)
    _copy_if_exists(default_design, backup_design)

    summaries: List[Dict[str, Any]] = []
    qc_script = os.path.join("scripts", "run_qc_and_report.py")

    try:
        for c_name, design_path, design_json, fasta_path in candidates:
            if not os.path.isfile(design_path):
                print(f"[跳过] design 不存在: {design_path}")
                continue
            if not os.path.isfile(fasta_path):
                print(f"[跳过] fasta 不存在: {fasta_path}")
                continue

            print(f"\n[候选] {c_name}")
            print(f"  design: {design_path}")
            print(f"  fasta : {fasta_path}")

            # 临时覆盖默认输入，让 run_qc_and_report.py 使用该候选。
            if os.path.abspath(design_path) != os.path.abspath(default_design):
                shutil.copy2(design_path, default_design)
            if os.path.abspath(fasta_path) != os.path.abspath(default_fasta):
                shutil.copy2(fasta_path, default_fasta)

            cmd = [sys.executable, qc_script, "--run_id", run_id]
            result = subprocess.run(cmd, check=False)
            if result.returncode != 0:
                print(f"[失败] QC 执行失败: {c_name} (exit={result.returncode})")
                continue

            dst_dir = os.path.join(archive_root, c_name)
            os.makedirs(dst_dir, exist_ok=True)

            # 归档本次候选产物
            _copy_if_exists(default_fasta, os.path.join(dst_dir, "mrna_vaccine.fasta"))
            _copy_if_exists(default_design, os.path.join(dst_dir, "mrna_design.json"))
            _copy_if_exists(os.path.join(run_dir, "qc_metrics.json"), os.path.join(dst_dir, "qc_metrics.json"))
            _copy_if_exists(os.path.join(run_dir, "REPORT.md"), os.path.join(dst_dir, "REPORT.md"))

            fig_src = os.path.join(run_dir, "figures")
            fig_dst = os.path.join(dst_dir, "figures")
            if os.path.isdir(fig_dst):
                shutil.rmtree(fig_dst)
            if os.path.isdir(fig_src):
                shutil.copytree(fig_src, fig_dst)

            # 读回核心指标，用于汇总
            qm_path = os.path.join(dst_dir, "qc_metrics.json")
            qm = {}
            if os.path.isfile(qm_path):
                with open(qm_path, "r", encoding="utf-8") as f:
                    qm = json.load(f)
            summaries.append(
                {
                    "candidate": c_name,
                    "codon_mode": design_json.get("codon_mode", ""),
                    "tool": str((design_json.get("codon_optimizer") or {}).get("tool", "")),
                    "length": qm.get("sequence_length", ""),
                    "gc": qm.get("gc_percent", ""),
                    "rnafold_mfe": qm.get("rnafold_mfe", ""),
                    "stability_status": qm.get("mrna_stability_status", ""),
                    "report": os.path.join(dst_dir, "REPORT.md"),
                }
            )
            print(f"[完成] 已归档到: {dst_dir}")
    finally:
        # 恢复默认文件，尽量不影响当前工作状态
        if had_default_fasta and os.path.isfile(backup_fasta):
            shutil.copy2(backup_fasta, default_fasta)
        elif (not had_default_fasta) and os.path.isfile(default_fasta):
            os.remove(default_fasta)

        if had_default_design and os.path.isfile(backup_design):
            shutil.copy2(backup_design, default_design)
        elif (not had_default_design) and os.path.isfile(default_design):
            os.remove(default_design)

    # 写汇总
    summary_path = os.path.join(archive_root, "qc_summary.md")
    # 以 default 作为对照，补充每条候选与主线差异，便于管理层快速看变化幅度。
    base = next((x for x in summaries if x.get("candidate") == "default"), None)
    base_gc = _to_float((base or {}).get("gc", ""))
    base_mfe = _to_float((base or {}).get("rnafold_mfe", ""))

    lines = [
        f"# 候选 mRNA 批量 QC 汇总（{run_id}）",
        "",
        "| 候选 | codon_mode | 工具 | 长度(nt) | GC% | RNAfold MFE | 相对 default GC变化 | 相对 default MFE变化 | 稳定性状态 | 报告路径 |",
        "|---|---|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for s in summaries:
        gc_val = _to_float(s.get("gc", ""))
        mfe_val = _to_float(s.get("rnafold_mfe", ""))
        delta_gc = ""
        delta_mfe = ""
        if s.get("candidate") == "default":
            delta_gc = "baseline"
            delta_mfe = "baseline"
        else:
            if gc_val is not None and base_gc is not None:
                delta_gc = f"{gc_val - base_gc:+.2f}"
            if mfe_val is not None and base_mfe is not None:
                delta_mfe = f"{mfe_val - base_mfe:+.1f}"
        lines.append(
            f"| {s['candidate']} | {s['codon_mode']} | {s['tool']} | {s['length']} | {s['gc']} | {s['rnafold_mfe']} | {delta_gc} | {delta_mfe} | {s['stability_status']} | `{s['report']}` |"
        )
    if not summaries:
        lines.append("| （无成功候选） | - | - | - | - | - | - | - | - | - |")
    lines.append("")
    lines.append("说明：每条候选完整报告与图文件在 `results/<run_id>/mrna_candidate_qc/<candidate>/`。")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))

    print(f"\n汇总完成：{summary_path}")


if __name__ == "__main__":
    main()
