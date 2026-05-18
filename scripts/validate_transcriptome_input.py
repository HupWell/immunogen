# -*- coding: utf-8 -*-
"""校验 deliveries/<run_id>/to_transcriptome/ 转录组输入契约。"""
import argparse
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transcriptome_cohort import parse_cohort_id_from_run_id


def main(run_id: str) -> None:
    base = os.path.join("deliveries", run_id, "to_transcriptome")
    meta_path = os.path.join(base, "meta.json")
    groups_path = os.path.join(base, "sample_groups.json")
    tpm_path = os.path.join(base, "tpm_matrix.tsv")

    for p in (meta_path, groups_path):
        if not os.path.isfile(p):
            raise FileNotFoundError(f"缺少: {p}")

    with open(meta_path, encoding="utf-8") as f:
        meta = json.load(f)
    with open(groups_path, encoding="utf-8") as f:
        groups = json.load(f)

    if not os.path.isfile(tpm_path) and not os.path.islink(tpm_path):
        raise FileNotFoundError(f"缺少 TPM 矩阵: {tpm_path}")

    if groups.get("skipped"):
        print(f"警告: 该 run 标记为跳过分组 — {groups.get('skip_reason', '')}")

    case_n = len(groups.get("case_columns") or [])
    ctrl_n = len(groups.get("control_columns") or [])
    if case_n == 0 or ctrl_n == 0:
        raise ValueError("case_columns 或 control_columns 为空，无法做差异排序。")

    print("转录组输入校验通过。")
    print(f"run_id: {run_id}")
    print(f"cohort_id: {meta.get('cohort_id') or parse_cohort_id_from_run_id(run_id)}")
    print(f"strategy: {groups.get('strategy')} | ranking_mode: {groups.get('ranking_mode')}")
    print(f"exploratory: {groups.get('exploratory')}")
    print(f"样本数: case={case_n}, control={ctrl_n}")
    if groups.get("exploratory"):
        print(f"注意: {groups.get('caution_zh', '')}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    args = parser.parse_args()
    main(args.run_id)
