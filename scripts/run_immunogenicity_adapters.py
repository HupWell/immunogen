# -*- coding: utf-8 -*-
"""
批量运行免疫原性适配器，生成预计算表：
results/<run_id>/tool_outputs/deepimmuno.tsv
results/<run_id>/tool_outputs/prime.tsv
results/<run_id>/tool_outputs/repitope.tsv
"""
import argparse
import os

from immunogenicity_adapters import (
    build_tool_df,
    ensure_tool_output_dir,
    load_unique_peptides,
    output_specs,
)


def main(run_id: str):
    peptides = load_unique_peptides(run_id)
    out_dir = ensure_tool_output_dir(run_id)
    specs = output_specs()
    for tool in ("deepimmuno", "prime", "repitope"):
        out_file = os.path.join(out_dir, specs[tool])
        df = build_tool_df(peptides, tool)
        df.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
        print(f"完成: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
