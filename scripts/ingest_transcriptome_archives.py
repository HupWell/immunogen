# -*- coding: utf-8 -*-
"""解压表达定量/样本信息 zip，生成 data/transcriptome/manifest/cohorts.csv。"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transcriptome_cohort import (
    MANIFEST_PATH,
    PDATA_ZIP,
    TPM_ZIP,
    build_manifest,
    ensure_archives_extracted,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="摄入转录组 zip 并生成 cohort manifest")
    parser.add_argument("--tpm_zip", default=str(TPM_ZIP), help="表达定量 zip 路径")
    parser.add_argument("--pdata_zip", default=str(PDATA_ZIP), help="样本信息 zip 路径")
    parser.add_argument("--manifest", default=str(MANIFEST_PATH), help="输出 manifest CSV")
    parser.add_argument("--skip_extract", action="store_true", help="跳过解压（已解压时）")
    args = parser.parse_args()

    if not args.skip_extract:
        tpm_dir, pdata_dir = ensure_archives_extracted(
            Path(args.tpm_zip), Path(args.pdata_zip)
        )
        print(f"TPM 目录: {tpm_dir}")
        print(f"pdata 目录: {pdata_dir}")

    rows = build_manifest(manifest_path=Path(args.manifest))
    paired = sum(1 for r in rows if r["pairing_status"] == "paired")
    print(f"manifest: {args.manifest}")
    print(f"队列数: {len(rows)}，TPM+pdata 配对: {paired}")


if __name__ == "__main__":
    main()
