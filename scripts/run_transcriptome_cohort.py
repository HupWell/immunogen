# -*- coding: utf-8 -*-
"""
单队列转录组流程：ingest 查找 → 解析分组 → 写 to_transcriptome → bulk 排名 → 更新 provenance。
"""
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transcriptome_cohort import (
    cohort_run_id,
    get_cohort_row,
    load_manifest,
    prepare_transcriptome_delivery,
    resolve_sample_groups,
)


def _merge_provenance(out_dir: str, resolution, cohort_id: str, run_id: str) -> None:
    prov_path = os.path.join(out_dir, "bulk_rank_provenance.json")
    if not os.path.isfile(prov_path):
        return
    with open(prov_path, encoding="utf-8") as f:
        prov = json.load(f)
    row = get_cohort_row(cohort_id)
    prov.update(
        {
            "cohort_id": str(cohort_id),
            "run_id": run_id,
            "pdata_path": row.get("pdata_path"),
            "grouping_strategy": resolution.strategy,
            "ranking_mode": resolution.ranking_mode,
            "exploratory": resolution.exploratory,
            "case_label": resolution.case_label,
            "control_label": resolution.control_label,
            "caution_zh": resolution.caution_zh,
        }
    )
    if resolution.exploratory:
        prov["caution_zh"] = (
            resolution.caution_zh
            + " 本输出为探索性基因优先级，不等价于新抗原表，不得单独作为 mRNA 疫苗设计依据。"
        )
    with open(prov_path, "w", encoding="utf-8") as f:
        json.dump(prov, f, ensure_ascii=False, indent=2)


def main() -> None:
    parser = argparse.ArgumentParser(description="运行单队列转录组基因优先级（路径 B）")
    parser.add_argument("--cohort_id", required=True, help="队列数字 ID，如 39004")
    parser.add_argument(
        "--run_id",
        default="",
        help="默认 T_cohort_<cohort_id>",
    )
    parser.add_argument(
        "--out_subdir",
        default="transcriptome_prior",
        help="results/<run_id>/ 下输出子目录",
    )
    parser.add_argument("--pseudo", type=float, default=1.0)
    parser.add_argument("--ingest", action="store_true", help="运行前先刷新 manifest（解压 zip）")
    parser.add_argument("--dry_run", action="store_true", help="只解析分组，不跑排名")
    args = parser.parse_args()

    cohort_id = str(args.cohort_id).strip()
    run_id = (args.run_id or "").strip() or cohort_run_id(cohort_id)

    if args.ingest or not load_manifest():
        from transcriptome_cohort import build_manifest, ensure_archives_extracted

        ensure_archives_extracted()
        build_manifest()

    resolution = resolve_sample_groups(cohort_id)
    print(
        f"队列 {cohort_id} | strategy={resolution.strategy} | "
        f"exploratory={resolution.exploratory} | "
        f"case={resolution.n_case} control={resolution.n_control}"
    )
    if resolution.skipped:
        print(f"跳过排名: {resolution.skip_reason}")
        prepare_transcriptome_delivery(cohort_id, resolution, run_id)
        sys.exit(0)

    prepare_transcriptome_delivery(cohort_id, resolution, run_id)

    if args.dry_run:
        print("dry_run: 已写入 deliveries，未执行 bulk_expression_gene_rank。")
        return

    out_dir = os.path.join("results", run_id, args.out_subdir)
    scripts_dir = Path(__file__).resolve().parent
    cmd = [
        sys.executable,
        str(scripts_dir / "validate_transcriptome_input.py"),
        "--run_id",
        run_id,
    ]
    subprocess.run(cmd, check=True)

    cmd = [
        sys.executable,
        str(scripts_dir / "bulk_expression_gene_rank.py"),
        "--tpm_path",
        os.path.join("deliveries", run_id, "to_transcriptome", "tpm_matrix.tsv"),
        "--out_dir",
        out_dir,
        "--case_columns",
        ",".join(resolution.case_columns),
        "--control_columns",
        ",".join(resolution.control_columns),
        "--pseudo",
        str(args.pseudo),
    ]
    print("[执行]", " ".join(cmd))
    subprocess.run(cmd, check=True)
    _merge_provenance(out_dir, resolution, cohort_id, run_id)
    print(f"完成: results/{run_id}/{args.out_subdir}/")


if __name__ == "__main__":
    main()
