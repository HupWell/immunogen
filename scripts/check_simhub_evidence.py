# -*- coding: utf-8 -*-
"""
功能：检查并归档 SimHub 回传证据状态。

输入：
- results/<run_id>/simhub_evidence/<case_id>/

期望回传：
- energy_report.json
- rmsd_profile.csv
- qc_flags.json
- summary.md
- 轨迹文件（*.xtc / *.dcd / *.nc / *.mdcrd / *.trr）
- 日志文件（*.log，或 logs/*.log）

输出：
- evidence_status.json
- SIMHUB_EVIDENCE.md
"""
import argparse
import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional


TRAJECTORY_PATTERNS = ("*.xtc", "*.dcd", "*.nc", "*.mdcrd", "*.trr")
LOG_PATTERNS = ("*.log", "logs/*.log")


def _load_json(path: Path) -> Dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def _case_id_for_run(run_id: str) -> str:
    run_meta = _load_json(Path("results") / run_id / "meta.json")
    if run_meta.get("case_id"):
        return str(run_meta["case_id"])
    simhub_root = Path("deliveries") / run_id / "to_simhub"
    case_dirs = sorted(p for p in simhub_root.glob("*") if p.is_dir())
    if case_dirs:
        return case_dirs[0].name
    return run_id


def _glob_any(root: Path, patterns) -> List[str]:
    matches: List[str] = []
    for pattern in patterns:
        matches.extend(str(p) for p in sorted(root.glob(pattern)) if p.is_file())
    return matches


def _qc_flags_failed(path: Path) -> Optional[bool]:
    if not path.exists():
        return None
    data = _load_json(path)
    if not data:
        return None
    text = json.dumps(data, ensure_ascii=False).lower()
    failure_tokens = ("fail", "failed", "error", "unstable", "reject")
    pass_tokens = ("pass", "passed", "ok", "stable")
    if any(token in text for token in failure_tokens):
        return True
    if any(token in text for token in pass_tokens):
        return False
    return None


def inspect_evidence(run_id: str, case_id: Optional[str] = None) -> Dict:
    case = case_id or _case_id_for_run(run_id)
    evidence_dir = Path("results") / run_id / "simhub_evidence" / case
    evidence_dir.mkdir(parents=True, exist_ok=True)

    files = {
        "energy_report": evidence_dir / "energy_report.json",
        "rmsd_profile": evidence_dir / "rmsd_profile.csv",
        "qc_flags": evidence_dir / "qc_flags.json",
        "summary": evidence_dir / "summary.md",
    }
    trajectories = _glob_any(evidence_dir, TRAJECTORY_PATTERNS)
    logs = _glob_any(evidence_dir, LOG_PATTERNS)
    returned_payload = [
        str(path)
        for path in files.values()
        if path.exists() and path.stat().st_size > 0
    ] + trajectories + logs

    missing = [
        name
        for name, path in files.items()
        if not path.exists() or path.stat().st_size == 0
    ]
    if not trajectories:
        missing.append("trajectory")
    if not logs:
        missing.append("log")

    qc_failed = _qc_flags_failed(files["qc_flags"])
    if not returned_payload:
        status = "not_returned"
        label = "未返回"
    elif missing:
        status = "validation_failed"
        label = "已验收失败"
    elif qc_failed is True:
        status = "validation_failed"
        label = "已验收失败"
        missing.append("qc_flags_failed")
    elif qc_failed is None:
        status = "returned_unvalidated"
        label = "已返回未验收"
    else:
        status = "validation_passed"
        label = "已验收通过"

    return {
        "run_id": run_id,
        "case_id": case,
        "checked_at": datetime.now().isoformat(timespec="seconds"),
        "status": status,
        "status_label": label,
        "evidence_dir": str(evidence_dir),
        "required_files": {name: str(path) for name, path in files.items()},
        "trajectory_files": trajectories,
        "log_files": logs,
        "missing_items": missing,
        "qc_flags_failed": qc_failed,
    }


def write_outputs(status: Dict) -> None:
    evidence_dir = Path(status["evidence_dir"])
    evidence_dir.mkdir(parents=True, exist_ok=True)

    status_path = evidence_dir / "evidence_status.json"
    with status_path.open("w", encoding="utf-8") as f:
        json.dump(status, f, ensure_ascii=False, indent=2)

    report_path = evidence_dir / "SIMHUB_EVIDENCE.md"
    missing = status.get("missing_items") or []
    content = f"""# SimHub Evidence Status（{status['run_id']}/{status['case_id']}）

- 状态：`{status['status']}`（{status['status_label']}）
- 检查时间：`{status['checked_at']}`
- 证据目录：`{status['evidence_dir']}`

## 必需证据

- energy report：`{status['required_files']['energy_report']}`
- RMSD profile：`{status['required_files']['rmsd_profile']}`
- QC flags：`{status['required_files']['qc_flags']}`
- summary：`{status['required_files']['summary']}`
- trajectory files：`{', '.join(status['trajectory_files']) if status['trajectory_files'] else 'not_provided'}`
- log files：`{', '.join(status['log_files']) if status['log_files'] else 'not_provided'}`

## 缺失或失败项

{chr(10).join(f'- `{item}`' for item in missing) if missing else '- 无'}

## 状态定义

- `not_returned`：SimHub 尚未回传证据。
- `returned_unvalidated`：证据已返回，但 QC flags 未给出明确通过/失败信号。
- `validation_passed`：证据完整，且 QC flags 未提示失败。
- `validation_failed`：证据缺失或 QC flags 提示失败。
"""
    with report_path.open("w", encoding="utf-8") as f:
        f.write(content)

    readme_path = evidence_dir / "README.md"
    with readme_path.open("w", encoding="utf-8") as f:
        f.write(content)


def main(run_id: str, case_id: Optional[str]) -> None:
    status = inspect_evidence(run_id, case_id=case_id)
    write_outputs(status)
    print(f"{run_id}/{status['case_id']}: {status['status']} ({status['status_label']})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True)
    parser.add_argument("--case_id", default="")
    args = parser.parse_args()
    main(args.run_id, args.case_id.strip() or None)
