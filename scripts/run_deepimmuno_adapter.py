# -*- coding: utf-8 -*-
"""独立适配器：输出 DeepImmuno 预计算表。"""
import argparse
import os

from immunogenicity_adapters import build_tool_df, ensure_tool_output_dir, load_unique_peptides, output_specs


def main(run_id: str):
    peptides = load_unique_peptides(run_id)
    out_dir = ensure_tool_output_dir(run_id)
    out_file = os.path.join(out_dir, output_specs()["deepimmuno"])
    df = build_tool_df(peptides, "deepimmuno")
    df.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
