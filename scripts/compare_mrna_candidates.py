# -*- coding: utf-8 -*-
"""
功能：汇总同一 run_id 下多条 mRNA 候选（mrna_design.json 与 mrna_design_<后缀>.json），
      生成 Markdown 对比表，便于按「序列长度、GC、密码子策略、工具」等挑选最合适的一条。

使用场景：老板要求多方案并行，技术侧用 --mrna_output_suffix 生成多条后，用本脚本一键对比，
      再选定主线复制为默认 mrna_vaccine.fasta 后跑 QC / SimHub。

用法：
  python scripts/compare_mrna_candidates.py --run_id R001
  python scripts/compare_mrna_candidates.py --run_id R001 --out results/R001/mrna_candidate_comparison.md
"""
from __future__ import annotations

import argparse
import glob
import json
import os
from typing import Any, Dict, List, Optional, Tuple


def _suffix_from_design_path(path: str) -> str:
    """从文件名推断候选后缀：mrna_design.json -> 默认；mrna_design_X.json -> X。"""
    base = os.path.basename(path)
    if base == "mrna_design.json":
        return "default"
    prefix = "mrna_design_"
    if base.startswith(prefix) and base.endswith(".json"):
        return base[len(prefix) : -len(".json")]
    return base


def _load_design(path: str) -> Tuple[str, Dict[str, Any]]:
    with open(path, "r", encoding="utf-8") as f:
        return path, json.load(f)


def _tool_label(d: Dict[str, Any]) -> str:
    co = d.get("codon_optimizer") or {}
    tool = co.get("tool") or ""
    if tool:
        return str(tool)
    return str(d.get("codon_mode", ""))


def main() -> None:
    parser = argparse.ArgumentParser(description="对比同一 run 下多条 mRNA 候选设计")
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument(
        "--out",
        default="",
        help="Markdown 输出路径，默认 results/<run_id>/mrna_candidate_comparison.md",
    )
    args = parser.parse_args()
    run_id = args.run_id.strip()
    out_dir = os.path.join("results", run_id)
    if not os.path.isdir(out_dir):
        raise FileNotFoundError(f"未找到目录: {out_dir}")

    pattern_default = os.path.join(out_dir, "mrna_design.json")
    pattern_suffix = os.path.join(out_dir, "mrna_design_*.json")
    paths: List[str] = []
    if os.path.isfile(pattern_default):
        paths.append(pattern_default)
    paths.extend(sorted(glob.glob(pattern_suffix)))
    # 去重、排序：default 优先，其余按后缀名
    seen = set()
    uniq: List[str] = []
    for p in paths:
        if p not in seen:
            seen.add(p)
            uniq.append(p)

    if not uniq:
        raise FileNotFoundError(f"{out_dir} 下未找到 mrna_design.json 或 mrna_design_*.json")

    rows: List[Dict[str, Any]] = []
    for path in uniq:
        _, data = _load_design(path)
        qm = data.get("quality_metrics") or {}
        out_files = data.get("output_files") or {}
        suffix_json = (data.get("mrna_output_suffix") or "").strip()
        suffix = suffix_json or _suffix_from_design_path(path)
        rows.append(
            {
                "候选标识": suffix,
                "design_json": path,
                "fasta": out_files.get("fasta", ""),
                "codon_mode": data.get("codon_mode", ""),
                "优化工具": _tool_label(data),
                "肽段数": data.get("selected_peptide_count", ""),
                "全长(nt)": qm.get("full_length", ""),
                "GC%": qm.get("gc_percent", ""),
            }
        )

    out_md = args.out.strip()
    if not out_md:
        out_md = os.path.join(out_dir, "mrna_candidate_comparison.md")

    lines: List[str] = [
        f"# mRNA 候选对比（{run_id}）",
        "",
        "生成后请结合业务侧重点选型：",
        "",
        "- **偏合规 / 真实工具验收**：优先 `codon_mode=real_cmd` 且优化工具为 LinearDesign 的候选。",
        "- **偏成本 / 快速迭代**：可在 `optimized` 或 `basic` 候选中先筛，再对赢家跑 `real_cmd`。",
        "- **偏稳定表达**：在候选确定后，以 `run_qc_and_report.py` 产出的 RNAfold MFE、RNAplfold 指标为准（需将选定 FASTA 置于默认 `mrna_vaccine.fasta` 再跑 QC）。",
        "",
        "| 候选标识 | codon_mode | 优化工具 | 肽段数 | 全长(nt) | GC% | FASTA | design JSON |",
        "|---|---|---|---:|---:|---:|---|---|",
    ]
    for r in rows:
        lines.append(
            "| {cid} | {cm} | {tool} | {npc} | {ln} | {gc} | `{fa}` | `{dj}` |".format(
                cid=r["候选标识"],
                cm=r["codon_mode"],
                tool=r["优化工具"],
                npc=r["肽段数"],
                ln=r["全长(nt)"],
                gc=r["GC%"],
                fa=r["fasta"] or "—",
                dj=r["design_json"],
            )
        )
    lines.extend(
        [
            "",
            "## 选定主线后的操作建议",
            "",
            f"1. 将选定候选的 FASTA 覆盖为 `{out_dir}/mrna_vaccine.fasta`，对应 JSON 覆盖为 `{out_dir}/mrna_design.json`（或自行 `cp` / 同步内容）。",
            "2. 再执行 `python scripts/run_qc_and_report.py --run_id " + run_id + "` 生成报告与稳定性指标。",
            "3. 再执行 SimHub 交付步骤（`prepare_simhub_delivery.py` 等），保证下游读取的是默认路径。",
            "",
        ]
    )

    text = "\n".join(lines)
    with open(out_md, "w", encoding="utf-8") as f:
        f.write(text)
    print(text)
    print(f"\n已写入: {out_md}")


if __name__ == "__main__":
    main()
