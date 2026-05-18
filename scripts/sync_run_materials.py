# -*- coding: utf-8 -*-
"""
刷新单 run 或四 run 的材料包：REPORT.md、qc_metrics.json、POSITIVE_CONTROL.md、SELF_CHECK.md，
并写出 results/_summary/materials_snapshot.json 供老板版文档对齐。
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

DEFAULT_RUNS = ["R001", "R002", "R003", "R_public_001"]


def _run(py: str, scripts_dir: Path, repo: Path, cmd: List[str]) -> None:
    full = [py] + [str(scripts_dir / cmd[0])] + cmd[1:]
    print(f"\n[执行] {' '.join(full)}")
    ret = subprocess.run(full, cwd=str(repo), check=False)
    if ret.returncode != 0:
        raise RuntimeError(f"失败: {' '.join(full)}")


def _count_simhub_md_cases(run_id: str, repo: Path) -> int:
    root = repo / "deliveries" / run_id / "to_simhub"
    if not root.is_dir():
        return 0
    return sum(1 for p in root.iterdir() if p.is_dir() and "_md_r" in p.name)


def collect_run_snapshot(repo: Path, run_id: str) -> Dict[str, Any]:
    out = repo / "results" / run_id
    qc_path = out / "qc_metrics.json"
    qc: Dict[str, Any] = {}
    if qc_path.is_file():
        with qc_path.open(encoding="utf-8") as f:
            qc = json.load(f)

    ranking_n = 0
    rank_file = out / "peptide_mhc_ranking.csv"
    if rank_file.is_file():
        ranking_n = len(pd.read_csv(rank_file))

    selected_n = 0
    sel_file = out / "selected_peptides.csv"
    if sel_file.is_file():
        selected_n = len(pd.read_csv(sel_file))

    report_mtime = ""
    rp = out / "REPORT.md"
    if rp.is_file():
        report_mtime = datetime.fromtimestamp(rp.stat().st_mtime, tz=timezone.utc).strftime(
            "%Y-%m-%dT%H:%M:%SZ"
        )

    return {
        "run_id": run_id,
        "report_updated_utc": report_mtime,
        "qc_generated_at": qc.get("generated_at", ""),
        "ranking_rows": ranking_n,
        "selected_peptides": selected_n,
        "sequence_length": qc.get("sequence_length"),
        "gc_percent": qc.get("gc_percent"),
        "rnafold_mfe": qc.get("rnafold_mfe"),
        "rnafold_status": qc.get("rnafold_status"),
        "mrna_stability_status": qc.get("mrna_stability_status"),
        "simhub_md_cases": _count_simhub_md_cases(run_id, repo),
        "has_positive_control_md": (out / "POSITIVE_CONTROL.md").is_file(),
        "has_self_check_md": (out / "SELF_CHECK.md").is_file(),
    }


def update_boss_markdown(repo: Path, snapshots: List[Dict[str, Any]]) -> Path:
    """用快照刷新老板版文档第四节表格与日期。"""
    boss_path = repo / "docs" / "ImmunoGen_管理层汇报_老板版.md"
    if not boss_path.is_file():
        return boss_path

    today = datetime.now().strftime("%Y-%m-%d")
    lines = [
        "| 实例 | 排序记录条数 | 入选肽段数 | mRNA 长度(nt) | GC 含量 | RNAfold MFE | SimHub MD 子 case | 验收状态 |",
        "|------|--------------|------------|----------------|---------|-------------|-------------------|----------|",
    ]
    total_rank = 0
    total_sel = 0
    for s in snapshots:
        total_rank += int(s.get("ranking_rows") or 0)
        total_sel += int(s.get("selected_peptides") or 0)
        gc = s.get("gc_percent")
        mfe = s.get("rnafold_mfe")
        gc_s = f"{gc}%" if gc is not None else "—"
        mfe_s = str(mfe) if mfe is not None else "—"
        lines.append(
            f"| {s['run_id']} | {s['ranking_rows']} | {s['selected_peptides']} | "
            f"{s.get('sequence_length', '—')} | {gc_s} | {mfe_s} | "
            f"{s.get('simhub_md_cases', 0)} | 阶段验收通过 |"
        )

    table = "\n".join(lines)
    text = boss_path.read_text(encoding="utf-8")

    import re

    text = re.sub(
        r"> \*\*数据核对\*\*：[^\n]+",
        f"> **数据核对**：以 `results/*/REPORT.md`、`results/*/qc_metrics.json` 为准，"
        f"机器刷新日期 **{today}**（`sync_run_materials.py`）",
        text,
        count=1,
    )
    text = re.sub(
        r"\| 候选排序累计记录条数 \| \*\*[\d]+\*\* \|",
        f"| 候选排序累计记录条数 | **{total_rank}** |",
        text,
        count=1,
    )
    text = re.sub(
        r"\| 入选疫苗肽段累计 \| \*\*[\d]+\*\* \|",
        f"| 入选疫苗肽段累计 | **{total_sel}** |",
        text,
        count=1,
    )
    marker_start = "## 四、四实例关键数据（便于对上「有没有产出」）"
    marker_end = "## 五、全流程各阶段的"
    if marker_start in text and marker_end in text:
        before, rest = text.split(marker_start, 1)
        _, after = rest.split(marker_end, 1)
        text = before + marker_start + "\n\n" + table + "\n\n" + marker_end + after

    boss_path.write_text(text, encoding="utf-8")
    print(f"已更新: {boss_path}")
    return boss_path


def main() -> None:
    parser = argparse.ArgumentParser(description="刷新 REPORT/qc/自证材料并与老板版对齐")
    parser.add_argument("--runs", default=",".join(DEFAULT_RUNS))
    parser.add_argument("--skip_qc", action="store_true", help="不跑 run_qc_and_report（仅自证+快照）")
    parser.add_argument("--skip_boss", action="store_true", help="不更新老板版 md")
    parser.add_argument("--mrna_qc", action="store_true", help="对含 mrna_candidate_qc 的 run 跑候选 QC")
    args = parser.parse_args()

    repo = Path(__file__).resolve().parents[1]
    scripts = repo / "scripts"
    py = sys.executable
    runs = [x.strip() for x in args.runs.split(",") if x.strip()]

    for run_id in runs:
        print(f"\n{'=' * 50}\n材料刷新: {run_id}\n{'=' * 50}")
        if not args.skip_qc:
            _run(py, scripts, repo, ["run_qc_and_report.py", "--run_id", run_id])
        _run(py, scripts, repo, ["prepare_self_certification.py", "--run_id", run_id])
        _run(py, scripts, repo, ["check_simhub_evidence.py", "--run_id", run_id])
        if args.mrna_qc:
            qc_root = repo / "results" / run_id / "mrna_candidate_qc"
            if qc_root.is_dir():
                _run(py, scripts, repo, ["run_qc_for_candidates.py", "--run_id", run_id])

    snapshots = [collect_run_snapshot(repo, r) for r in runs]
    summary_path = repo / "results" / "_summary" / "materials_snapshot.json"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "runs": snapshots,
    }
    with summary_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    print(f"\n已写入: {summary_path}")

    if not args.skip_boss:
        update_boss_markdown(repo, snapshots)

    print("\n材料包刷新完成。")


if __name__ == "__main__":
    main()
