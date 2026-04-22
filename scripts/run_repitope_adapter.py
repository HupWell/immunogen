# -*- coding: utf-8 -*-
"""独立适配器：输出 Repitope 预计算表。"""
import argparse
import os

from immunogenicity_adapters import (
    build_tool_df_with_backend,
    ensure_tool_output_dir,
    load_unique_peptides,
    output_specs,
)


def main(run_id: str, backend: str):
    peptides = load_unique_peptides(run_id)
    out_dir = ensure_tool_output_dir(run_id)
    out_file = os.path.join(out_dir, output_specs()["repitope"])
    df, used = build_tool_df_with_backend(peptides, "repitope", run_id, backend=backend)
    df.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}（backend={backend} -> used={used}）")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument(
        "--backend",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
    )
    args = parser.parse_args()
    main(args.run_id, args.backend)
