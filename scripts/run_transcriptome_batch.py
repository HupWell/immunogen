# -*- coding: utf-8 -*-
"""对 manifest 中可解析分组的队列批量运行转录组排名。"""
import argparse
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from transcriptome_cohort import load_manifest, resolve_sample_groups


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ingest", action="store_true", help="先刷新 manifest")
    parser.add_argument("--dry_run", action="store_true", help="只检查分组，不跑排名")
    parser.add_argument("--cohort_id", default="", help="仅跑指定队列")
    args = parser.parse_args()

    script = Path(__file__).resolve().parent / "run_transcriptome_cohort.py"
    if args.ingest:
        subprocess.run([sys.executable, str(Path(__file__).parent / "ingest_transcriptome_archives.py")], check=True)

    rows = load_manifest()
    if args.cohort_id:
        rows = [r for r in rows if r["cohort_id"] == args.cohort_id.strip()]

    summary = []
    for row in rows:
        cid = row["cohort_id"]
        if row["pairing_status"] != "paired":
            summary.append((cid, "skip_manifest", row["pairing_status"]))
            continue
        res = resolve_sample_groups(cid)
        if res.skipped:
            summary.append((cid, "skip_groups", res.skip_reason))
            continue
        cmd = [sys.executable, str(script), "--cohort_id", cid]
        if args.dry_run:
            cmd.append("--dry_run")
        print(f"\n=== cohort {cid} ({res.strategy}, exploratory={res.exploratory}) ===")
        ret = subprocess.run(cmd).returncode
        summary.append((cid, "ok" if ret == 0 else "fail", res.strategy))

    print("\n--- 批量汇总 ---")
    for cid, status, detail in summary:
        print(f"{cid}\t{status}\t{detail}")


if __name__ == "__main__":
    main()
