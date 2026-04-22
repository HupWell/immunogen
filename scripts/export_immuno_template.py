# -*- coding: utf-8 -*-
"""
导出免疫原性真模型待补全模板（全量 unique mut_peptide）。

默认生成：
results/<run_id>/tool_outputs/raw/deepimmuno.tsv

可选同时生成 prime/repitope 模板，便于后续人工或外部工具回填真实分数。
"""
import argparse
import os

import pandas as pd

from immunogenicity_adapters import ensure_tool_output_dir, load_unique_peptides


def _template_path(run_id: str, tool: str) -> str:
    return os.path.join("results", run_id, "tool_outputs", "raw", f"{tool}.tsv")


def _write_template(run_id: str, peptides: pd.Series, tool: str, overwrite: bool):
    out_dir = os.path.join(ensure_tool_output_dir(run_id), "raw")
    os.makedirs(out_dir, exist_ok=True)
    out_file = _template_path(run_id, tool)
    if os.path.exists(out_file) and not overwrite:
        print(f"跳过（已存在）: {out_file}")
        return
    df = pd.DataFrame({"mut_peptide": peptides, "score": pd.NA})
    df.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")


def main(run_id: str, include_all_tools: bool, overwrite: bool):
    peptides = load_unique_peptides(run_id)
    tools = ["deepimmuno"]
    if include_all_tools:
        tools = ["deepimmuno", "prime", "repitope"]
    for tool in tools:
        _write_template(run_id, peptides, tool, overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument(
        "--include_all_tools",
        action="store_true",
        help="若指定，则同时导出 prime.tsv 与 repitope.tsv 模板。",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="若指定，则覆盖已存在模板。",
    )
    args = parser.parse_args()
    main(args.run_id, args.include_all_tools, args.overwrite)
