# -*- coding: utf-8 -*-
"""
转录组队列公共工具：解压归档、manifest、样本分组解析（含探索性 death/meta）。
"""
from __future__ import annotations

import csv
import json
import os
import re
import zipfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    import yaml
except ImportError:
    yaml = None  # type: ignore

REPO_ROOT = Path(__file__).resolve().parents[1]
TRANSCRIPTOME_ROOT = REPO_ROOT / "data" / "transcriptome"
MANIFEST_PATH = TRANSCRIPTOME_ROOT / "manifest" / "cohorts.csv"
COHORT_GROUPS_YAML = TRANSCRIPTOME_ROOT / "cohort_groups.yaml"
TPM_ZIP = REPO_ROOT / "data" / "表达定量表格.zip"
PDATA_ZIP = REPO_ROOT / "data" / "样本信息表.zip"

RUN_ID_PREFIX = "T_cohort_"


def cohort_run_id(cohort_id: str) -> str:
    """一队列一个 run_id。"""
    cid = str(cohort_id).strip()
    if cid.startswith(RUN_ID_PREFIX):
        return cid
    return f"{RUN_ID_PREFIX}{cid}"


def parse_cohort_id_from_run_id(run_id: str) -> Optional[str]:
    rid = (run_id or "").strip()
    if rid.startswith(RUN_ID_PREFIX):
        return rid[len(RUN_ID_PREFIX) :]
    return None


def ensure_archives_extracted(
    tpm_zip: Path = TPM_ZIP,
    pdata_zip: Path = PDATA_ZIP,
    out_root: Path = TRANSCRIPTOME_ROOT,
) -> Tuple[Path, Path]:
    """解压两个 zip 到 data/transcriptome/raw/tpm 与 raw/pdata。"""
    tpm_dir = out_root / "raw" / "tpm"
    pdata_dir = out_root / "raw" / "pdata"
    tpm_dir.mkdir(parents=True, exist_ok=True)
    pdata_dir.mkdir(parents=True, exist_ok=True)

    if tpm_zip.is_file():
        with zipfile.ZipFile(tpm_zip, "r") as zf:
            zf.extractall(tpm_dir)
    if pdata_zip.is_file():
        with zipfile.ZipFile(pdata_zip, "r") as zf:
            zf.extractall(pdata_dir)

    return tpm_dir, pdata_dir


def _cohort_id_from_tpm(name: str) -> Optional[str]:
    m = re.match(r"TPM_(\d+)_cleaning\.txt$", name, re.I)
    return m.group(1) if m else None


def _cohort_id_from_pdata(name: str) -> Optional[str]:
    m = re.match(r"pdata_(\d+)(?:_[^./]+)?\.csv$", name, re.I)
    return m.group(1) if m else None


def pick_pdata_for_cohort(pdata_root: Path, cohort_id: str) -> Optional[Path]:
    """优先主表 pdata_{id}.csv，否则取同 id 最大文件变体。"""
    cid = str(cohort_id)
    all_csv = list(pdata_root.rglob(f"pdata_{cid}.csv"))
    if all_csv:
        return max(all_csv, key=lambda p: p.stat().st_size)
    candidates = list(pdata_root.rglob(f"pdata_{cid}_*.csv"))
    if not candidates:
        candidates = [p for p in pdata_root.rglob("pdata_*.csv") if _cohort_id_from_pdata(p.name) == cid]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_size)


def build_manifest(
    tpm_dir: Optional[Path] = None,
    pdata_dir: Optional[Path] = None,
    manifest_path: Path = MANIFEST_PATH,
) -> List[Dict[str, str]]:
    """扫描 TPM / pdata 配对并写入 manifest（递归搜索解压子目录）。"""
    if tpm_dir is None or pdata_dir is None:
        tpm_dir, pdata_dir = ensure_archives_extracted()

    tpm_by_id: Dict[str, Path] = {}
    for p in tpm_dir.rglob("TPM_*_cleaning.txt"):
        cid = _cohort_id_from_tpm(p.name)
        if cid:
            tpm_by_id[cid] = p

    pdata_by_id: Dict[str, Path] = {}
    seen_ids: set = set()
    for p in sorted(pdata_dir.rglob("pdata_*.csv")):
        cid = _cohort_id_from_pdata(p.name)
        if not cid or cid in seen_ids:
            continue
        chosen = pick_pdata_for_cohort(pdata_dir, cid)
        if chosen:
            pdata_by_id[cid] = chosen
            seen_ids.add(cid)

    all_ids = sorted(set(tpm_by_id) | set(pdata_by_id))
    rows: List[Dict[str, str]] = []
    for cid in all_ids:
        tpm_p = tpm_by_id.get(cid)
        pdata_p = pdata_by_id.get(cid)
        if tpm_p and pdata_p:
            status = "paired"
        elif tpm_p:
            status = "tpm_only"
        elif pdata_p:
            status = "pdata_only"
        else:
            status = "unknown"
        rows.append(
            {
                "cohort_id": cid,
                "run_id": cohort_run_id(cid),
                "tpm_path": str(tpm_p.resolve()) if tpm_p else "",
                "pdata_path": str(pdata_p.resolve()) if pdata_p else "",
                "pairing_status": status,
            }
        )

    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    with open(manifest_path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["cohort_id", "run_id", "tpm_path", "pdata_path", "pairing_status"],
        )
        w.writeheader()
        w.writerows(rows)
    return rows


def load_manifest(manifest_path: Path = MANIFEST_PATH) -> List[Dict[str, str]]:
    if not manifest_path.is_file():
        return build_manifest()
    with open(manifest_path, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def get_cohort_row(cohort_id: str, manifest_path: Path = MANIFEST_PATH) -> Dict[str, str]:
    cid = str(cohort_id).strip()
    for row in load_manifest(manifest_path):
        if row.get("cohort_id") == cid:
            return row
    raise KeyError(f"manifest 中未找到 cohort_id={cid!r}，请先运行 ingest_transcriptome_archives.py")


def load_cohort_groups_yaml(path: Path = COHORT_GROUPS_YAML) -> Dict[str, Any]:
    if not path.is_file():
        return {"default_strategy_order": ["tumor_normal", "exploratory_death", "exploratory_meta"], "cohorts": {}}
    if yaml is None:
        raise RuntimeError("需要 PyYAML：pip install pyyaml")
    with open(path, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def read_pdata(pdata_path: Path) -> Tuple[List[str], List[Dict[str, str]]]:
    """返回 (字段名列表, 行字典列表)。样本 ID 列优先 geo_accession。"""
    with open(pdata_path, encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames or []
        rows = [dict(r) for r in reader]
    return list(fields), rows


def _sample_id_from_row(row: Dict[str, str]) -> str:
    for key in ("geo_accession", "GEO_accession", ""):
        if key and row.get(key):
            return str(row[key]).strip().strip('"')
    first = row.get("") or row.get(None)
    if first:
        return str(first).strip().strip('"')
    return ""


def _read_tpm_header(tpm_path: Path) -> List[str]:
    with open(tpm_path, encoding="utf-8", errors="replace") as f:
        line = f.readline().strip()
    sep = "\t" if line.count("\t") >= line.count(",") else ","
    return [c.strip() for c in line.split(sep)][1:]


def _filter_gsm_in_tpm(gsms: List[str], tpm_samples: List[str]) -> List[str]:
    tpm_set = set(tpm_samples)
    return [g for g in gsms if g in tpm_set]


@dataclass
class GroupResolution:
    cohort_id: str
    strategy: str
    ranking_mode: str  # tumor_normal | clinical_proxy
    case_columns: List[str]
    control_columns: List[str]
    case_label: str
    control_label: str
    exploratory: bool
    caution_zh: str
    n_case: int
    n_control: int
    skipped: bool = False
    skip_reason: str = ""

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _strategy_from_yaml(cohort_id: str, cfg: Dict[str, Any]) -> Optional[GroupResolution]:
    cohorts = cfg.get("cohorts") or {}
    spec = cohorts.get(str(cohort_id)) or cohorts.get(int(cohort_id))  # type: ignore
    if not spec:
        return None
    if spec.get("strategy") == "skip":
        return GroupResolution(
            cohort_id=str(cohort_id),
            strategy="skip",
            ranking_mode="",
            case_columns=[],
            control_columns=[],
            case_label="",
            control_label="",
            exploratory=False,
            caution_zh=spec.get("note_zh", "配置为跳过。"),
            n_case=0,
            n_control=0,
            skipped=True,
            skip_reason=spec.get("note_zh", "cohort_groups.yaml 指定 skip"),
        )
    if spec.get("strategy") == "manual":
        case_col = spec["case_column"]
        ctrl_col = spec.get("control_column", case_col)
        case_vals = spec["case_values"]
        ctrl_vals = spec["control_values"]
        return None  # 需要 pdata 行，在外层用 _resolve_manual 处理
    return None


def _resolve_manual(
    cohort_id: str,
    rows: List[Dict[str, str]],
    tpm_samples: List[str],
    spec: Dict[str, Any],
) -> GroupResolution:
    case_col = spec["case_column"]
    ctrl_col = spec.get("control_column", case_col)
    case_vals = {str(v).strip() for v in spec["case_values"]}
    ctrl_vals = {str(v).strip() for v in spec["control_values"]}
    exploratory = bool(spec.get("exploratory", False))
    case_gsms: List[str] = []
    ctrl_gsms: List[str] = []
    for row in rows:
        sid = _sample_id_from_row(row)
        if not sid:
            continue
        cv = (row.get(case_col) or "").strip()
        rv = (row.get(ctrl_col) or "").strip()
        if cv in case_vals:
            case_gsms.append(sid)
        if rv in ctrl_vals:
            ctrl_gsms.append(sid)
    case_gsms = _filter_gsm_in_tpm(case_gsms, tpm_samples)
    ctrl_gsms = _filter_gsm_in_tpm(ctrl_gsms, tpm_samples)
    mode = "clinical_proxy" if exploratory else "tumor_normal"
    caution = spec.get(
        "caution_zh",
        "手动分组；探索性排序不得单独作为疫苗依据。" if exploratory else "手动肿瘤/对照分组。",
    )
    return GroupResolution(
        cohort_id=str(cohort_id),
        strategy="manual",
        ranking_mode=mode,
        case_columns=case_gsms,
        control_columns=ctrl_gsms,
        case_label=f"{case_col} in {sorted(case_vals)[:3]}",
        control_label=f"{ctrl_col} in {sorted(ctrl_vals)[:3]}",
        exploratory=exploratory,
        caution_zh=caution,
        n_case=len(case_gsms),
        n_control=len(ctrl_gsms),
        skipped=len(case_gsms) == 0 or len(ctrl_gsms) == 0,
        skip_reason="" if (case_gsms and ctrl_gsms) else "手动规则未匹配到足够样本",
    )


def _tumor_normal_label_from_text(text: str) -> Optional[str]:
    """从 source/title/characteristics 文本判断 case/control；无法判断则 None。"""
    if not text or not str(text).strip():
        return None
    t = str(text).strip().lower()
    if any(x in t for x in ("reference", "stratagene", "human universal", "human reference")):
        return None
    if "normal-like" in t or "normal breast-like" in t or "normal like" in t:
        return None
    ctrl_hints = (
        "normal breast tissue",
        "breast tissue, normal",
        "tissue: normal breast",
        "tissue type: normal breast",
        "sample: normal breast",
        "cancer-adjacent normal",
    )
    for h in ctrl_hints:
        if h in t:
            return "control"
    if t == "human breast":
        return "control"
    if re.search(r"\bnormal\b", t) and "tumor" not in t and "cancer" not in t and "carcinoma" not in t:
        if t.startswith("normal-") or "normal breast" in t:
            return "control"
    case_hints = (
        "breast tumor",
        "breast tissue, cancer",
        "invasive breast",
        "ductal carcinoma",
        "metastasis",
        "primary tumor",
        "primary breast tumors",
        "frozen tissue of primary breast",
        "tissue: tumor",
        "sample: tumor",
        "tissue type: invasive",
        "tissue type: ductal",
    )
    for h in case_hints:
        if h in t:
            return "case"
    if "tumor" in t or ("cancer" in t and "adjacent" not in t):
        return "case"
    if "adenocarcinoma" in t:
        return "case"
    return None


def _label_row_tumor_normal(row: Dict[str, str]) -> Optional[str]:
    """单样本肿瘤/正常标签；双通道阵列在 ch1 为 reference 时用 ch2 判定。"""
    src1 = (row.get("source_name_ch1") or "").strip().lower()
    src2 = (row.get("source_name_ch2") or "").strip()
    if src1 and any(x in src1 for x in ("reference", "stratagene", "universal")):
        lab = _tumor_normal_label_from_text(src2)
        if lab:
            return lab
        if src2 and not re.search(r"normal", src2, re.I):
            return "case"
    for col in ("source_name_ch1", "source_name_ch2", "characteristics_ch1", "characteristics_ch2", "title"):
        val = (row.get(col) or "").strip()
        lab = _tumor_normal_label_from_text(val)
        if lab:
            return lab
    return None


def _resolve_tumor_normal(
    cohort_id: str,
    rows: List[Dict[str, str]],
    tpm_samples: List[str],
) -> Optional[GroupResolution]:
    case_gsms: List[str] = []
    ctrl_gsms: List[str] = []
    for row in rows:
        sid = _sample_id_from_row(row)
        if not sid:
            continue
        label = _label_row_tumor_normal(row)
        if label == "case":
            case_gsms.append(sid)
        elif label == "control":
            ctrl_gsms.append(sid)
    case_gsms = _filter_gsm_in_tpm(case_gsms, tpm_samples)
    ctrl_gsms = _filter_gsm_in_tpm(ctrl_gsms, tpm_samples)
    if len(case_gsms) < 2 or len(ctrl_gsms) < 2:
        return None
    return GroupResolution(
        cohort_id=str(cohort_id),
        strategy="tumor_normal",
        ranking_mode="tumor_normal",
        case_columns=case_gsms,
        control_columns=ctrl_gsms,
        case_label="tumor (title/source 启发式)",
        control_label="normal (title/source 启发式)",
        exploratory=False,
        caution_zh="基于 title/source_name 的肿瘤 vs 正常组织分组；非个体配对时需结合队列设计解读。",
        n_case=len(case_gsms),
        n_control=len(ctrl_gsms),
    )


def _resolve_binary_column(
    cohort_id: str,
    rows: List[Dict[str, str]],
    tpm_samples: List[str],
    col: str,
    strategy: str,
    case_value: str,
    control_value: str,
) -> Optional[GroupResolution]:
    case_gsms: List[str] = []
    ctrl_gsms: List[str] = []
    for row in rows:
        sid = _sample_id_from_row(row)
        if not sid:
            continue
        v = (row.get(col) or "").strip()
        if v == case_value:
            case_gsms.append(sid)
        elif v == control_value:
            ctrl_gsms.append(sid)
    case_gsms = _filter_gsm_in_tpm(case_gsms, tpm_samples)
    ctrl_gsms = _filter_gsm_in_tpm(ctrl_gsms, tpm_samples)
    if len(case_gsms) < 3 or len(ctrl_gsms) < 3:
        return None
    label = col.replace(":ch1", "")
    return GroupResolution(
        cohort_id=str(cohort_id),
        strategy=strategy,
        ranking_mode="clinical_proxy",
        case_columns=case_gsms,
        control_columns=ctrl_gsms,
        case_label=f"{label}={case_value}",
        control_label=f"{label}={control_value}",
        exploratory=True,
        caution_zh=(
            f"探索性排序：以 {label} 二分（{case_value} vs {control_value}）代理「病例/对照」，"
            "反映临床结局或转移等混杂，不等价于肿瘤组织 vs 正常组织，"
            "不得单独作为新抗原或疫苗设计依据。"
        ),
        n_case=len(case_gsms),
        n_control=len(ctrl_gsms),
    )


def resolve_sample_groups(
    cohort_id: str,
    pdata_path: Optional[Path] = None,
    tpm_path: Optional[Path] = None,
    yaml_path: Path = COHORT_GROUPS_YAML,
) -> GroupResolution:
    """解析单个队列的 case/control GSM 列表。"""
    row = get_cohort_row(cohort_id)
    if pdata_path is None:
        pdata_path = Path(row["pdata_path"])
    if tpm_path is None:
        tpm_path = Path(row["tpm_path"])
    if not pdata_path.is_file() or not tpm_path.is_file():
        return GroupResolution(
            cohort_id=str(cohort_id),
            strategy="unpaired",
            ranking_mode="",
            case_columns=[],
            control_columns=[],
            case_label="",
            control_label="",
            exploratory=False,
            caution_zh="缺少 TPM 或 pdata 文件。",
            n_case=0,
            n_control=0,
            skipped=True,
            skip_reason=row.get("pairing_status", "not_paired"),
        )

    _, rows = read_pdata(pdata_path)
    tpm_samples = _read_tpm_header(tpm_path)
    cfg = load_cohort_groups_yaml(yaml_path)
    cohorts_cfg = cfg.get("cohorts") or {}
    spec = cohorts_cfg.get(str(cohort_id))
    if spec:
        if spec.get("strategy") == "skip":
            return GroupResolution(
                cohort_id=str(cohort_id),
                strategy="skip",
                ranking_mode="",
                case_columns=[],
                control_columns=[],
                case_label="",
                control_label="",
                exploratory=False,
                caution_zh=spec.get("note_zh", ""),
                n_case=0,
                n_control=0,
                skipped=True,
                skip_reason=spec.get("note_zh", "skip"),
            )
        if spec.get("strategy") == "manual":
            return _resolve_manual(cohort_id, rows, tpm_samples, spec)

    order = cfg.get("default_strategy_order") or [
        "tumor_normal",
        "exploratory_death",
        "exploratory_meta",
    ]
    for strat in order:
        if strat == "tumor_normal":
            res = _resolve_tumor_normal(cohort_id, rows, tpm_samples)
            if res:
                return res
        elif strat == "exploratory_death":
            res = _resolve_binary_column(
                cohort_id, rows, tpm_samples, "death:ch1", "exploratory_death", "1", "0"
            )
            if res:
                return res
        elif strat == "exploratory_meta":
            res = _resolve_binary_column(
                cohort_id, rows, tpm_samples, "meta:ch1", "exploratory_meta", "1", "0"
            )
            if res:
                return res

    return GroupResolution(
        cohort_id=str(cohort_id),
        strategy="none",
        ranking_mode="",
        case_columns=[],
        control_columns=[],
        case_label="",
        control_label="",
        exploratory=False,
        caution_zh="未能自动解析分组；请在 cohort_groups.yaml 中补充 manual 规则。",
        n_case=0,
        n_control=0,
        skipped=True,
        skip_reason="no_strategy_matched",
    )


def delivery_transcriptome_dir(run_id: str) -> Path:
    return REPO_ROOT / "deliveries" / run_id / "to_transcriptome"


def prepare_transcriptome_delivery(
    cohort_id: str,
    resolution: GroupResolution,
    run_id: Optional[str] = None,
) -> Path:
    """写入 deliveries/<run_id>/to_transcriptome/ 契约文件。"""
    rid = run_id or cohort_run_id(cohort_id)
    row = get_cohort_row(cohort_id)
    out = delivery_transcriptome_dir(rid)
    out.mkdir(parents=True, exist_ok=True)

    tpm_src = Path(row["tpm_path"])
    tpm_link = out / "tpm_matrix.tsv"
    if tpm_link.exists() or tpm_link.is_symlink():
        tpm_link.unlink()
    try:
        os.symlink(os.path.abspath(tpm_src), tpm_link)
    except OSError:
        import shutil
        shutil.copy2(tpm_src, tpm_link)

    groups = {
        "cohort_id": str(cohort_id),
        "run_id": rid,
        "strategy": resolution.strategy,
        "ranking_mode": resolution.ranking_mode,
        "exploratory": resolution.exploratory,
        "case_columns": resolution.case_columns,
        "control_columns": resolution.control_columns,
        "case_label": resolution.case_label,
        "control_label": resolution.control_label,
        "caution_zh": resolution.caution_zh,
        "skipped": resolution.skipped,
        "skip_reason": resolution.skip_reason,
    }
    with open(out / "sample_groups.json", "w", encoding="utf-8") as f:
        json.dump(groups, f, ensure_ascii=False, indent=2)

    meta = {
        "run_id": rid,
        "cohort_id": str(cohort_id),
        "input_type": "bulk_transcriptome",
        "tpm_path": str(tpm_link.resolve()),
        "pdata_path": row["pdata_path"],
        "pairing_status": row["pairing_status"],
        "ranking_mode": resolution.ranking_mode,
        "exploratory": resolution.exploratory,
    }
    with open(out / "meta.json", "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    return out
