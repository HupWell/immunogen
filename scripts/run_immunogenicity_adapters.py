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
    build_tool_df_with_backend,
    ensure_tool_output_dir,
    load_unique_peptides,
    output_specs,
)


def main(run_id: str, backend_deepimmuno: str, backend_prime: str, backend_repitope: str):
    peptides = load_unique_peptides(run_id)
    out_dir = ensure_tool_output_dir(run_id)
    specs = output_specs()
    backends = {
        "deepimmuno": backend_deepimmuno,
        "prime": backend_prime,
        "repitope": backend_repitope,
    }
    for tool in ("deepimmuno", "prime", "repitope"):
        out_file = os.path.join(out_dir, specs[tool])
        backend = backends[tool]
        try:
            df, used = build_tool_df_with_backend(peptides, tool, run_id, backend=backend)
        except Exception as e:
            raise RuntimeError(f"{tool} 适配器失败（backend={backend}）: {e}")
        df.to_csv(out_file, sep="\t", index=False, encoding="utf-8-sig")
        print(f"完成: {out_file}（{tool} backend={backend} -> used={used}）")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument(
        "--backend_deepimmuno",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="DeepImmuno 适配器模式",
    )
    parser.add_argument(
        "--backend_prime",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="PRIME 适配器模式",
    )
    parser.add_argument(
        "--backend_repitope",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="Repitope 适配器模式",
    )
    args = parser.parse_args()
    main(
        run_id=args.run_id,
        backend_deepimmuno=args.backend_deepimmuno,
        backend_prime=args.backend_prime,
        backend_repitope=args.backend_repitope,
    )
