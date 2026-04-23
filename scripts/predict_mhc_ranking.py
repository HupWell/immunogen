# -*- coding: utf-8 -*-
"""
脚本名称：predict_mhc_ranking.py
主要功能：
1) 读取 BioDriver 提供的候选突变肽和 HLA 分型；
2) 使用 MHCflurry 预测每条突变肽与病人 HLA-I 的结合亲和力；
3) 可选使用 NetMHCIIpan（`--mhc2_backend`）做 MHC-II 层评分，否则使用代理分；
4) 生成排名文件 peptide_mhc_ranking.csv。

输入文件：
- deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv
- deliveries/<run_id>/to_immunogen/hla_typing.json

输出文件：
- results/<run_id>/peptide_mhc_ranking.csv
"""

import os
import json
import argparse
import subprocess
from typing import Any, Dict, Optional, Tuple

import pandas as pd
import numpy as np

from netmhciipan_runner import (
    build_peptide_mhc2_lookup,
    resolve_netmhciipan_bin,
    collect_netmhcii_alleles,
)
from immunogenicity_adapters import (
    deepimmuno_proxy,
    prime_proxy,
    repitope_proxy,
)


def flatten_hla(hla_json: dict):
    """
    把 HLA-A/HLA-B/HLA-C 三类等位基因合并成一个列表。
    例如：
    {
      "HLA-A": ["A*02:01"],
      "HLA-B": ["B*40:01"],
      "HLA-C": ["C*07:02"]
    }
    -> ["A*02:01", "B*40:01", "C*07:02"]
    """
    alleles = []
    for k in ["HLA-A", "HLA-B", "HLA-C"]:
        values = hla_json.get(k, [])
        if isinstance(values, list):
            alleles.extend(values)
    return alleles


def normalize_allele(a: str):
    """
    把 HLA 格式标准化为 MHCflurry 可识别格式。
    - 输入可能是 A*02:01 / HLA-A*02:01
    - 统一输出 HLA-A*02:01
    若格式不符合 A*/B*/C*，返回 None。
    """
    a = a.strip().replace("HLA-", "")
    # 统一为 MHCflurry 可识别格式，例如 HLA-A*02:01
    if not a.startswith(("A*", "B*", "C*")):
        return None
    return "HLA-" + a


def minmax_series(s: pd.Series, reverse: bool = False) -> pd.Series:
    """对分值做 0-1 归一化；reverse=True 表示原值越小越好。"""
    s = pd.to_numeric(s, errors="coerce")
    if s.notna().sum() == 0:
        return pd.Series([0.0] * len(s), index=s.index)
    if reverse:
        s = -s
    vmin, vmax = s.min(), s.max()
    if pd.isna(vmin) or pd.isna(vmax) or np.isclose(vmax, vmin):
        return pd.Series([0.5] * len(s), index=s.index)
    return (s - vmin) / (vmax - vmin)


def hamming_dissimilarity(mut_peptide: str, wt_peptide: str) -> float:
    """计算突变肽与 WT 的位点差异比例。"""
    mut_peptide = str(mut_peptide or "").strip().upper()
    wt_peptide = str(wt_peptide or "").strip().upper()
    if not mut_peptide or not wt_peptide:
        return 0.0
    if len(mut_peptide) != len(wt_peptide):
        return 1.0
    diff = sum(1 for a, b in zip(mut_peptide, wt_peptide) if a != b)
    return diff / len(mut_peptide)


def mhc2_score_from_el_rank(pct_rank_el: float) -> float:
    """%Rank EL 越小越好，转成 0~1 分（越大表示 II 类呈递预测越好，便于与 proxy 同量纲）。"""
    v = 1.0 - float(pct_rank_el) / 100.0
    return round(max(0.0, min(1.0, v)), 4)


def mhc2_proxy_score(peptide: str, hla_allele: str) -> float:
    """
    MHC-II 代理评分（0-1）：
    - 偏好长度 13-18
    - 带电与疏水残基均衡
    说明：当前为可复现代理模型，后续可替换为 NetMHCIIpan 实测值。
    """
    pep = peptide.upper()
    length = len(pep)
    length_score = max(0.0, 1.0 - abs(length - 15) / 10.0)
    hydrophobic = sum(aa in "AILMFWVY" for aa in pep) / max(length, 1)
    charged = sum(aa in "KRDEH" for aa in pep) / max(length, 1)
    balance_score = 1.0 - min(1.0, abs(hydrophobic - charged))
    allele_bias = 0.55 if "DRB1" in hla_allele else 0.5
    return round(max(0.0, min(1.0, (length_score * 0.5 + balance_score * 0.4 + allele_bias * 0.1))), 4)


def load_precomputed_immunogenicity(run_id: str) -> Dict[str, Dict[str, float]]:
    """
    读取 results/<run_id>/tool_outputs/*.tsv 作为预计算表。
    缺失文件不报错，返回可用部分；主流程会对缺失值回退 proxy。
    """
    base = os.path.join("results", run_id, "tool_outputs")
    files = {
        "immunogenicity_deepimmuno": "deepimmuno.tsv",
        "immunogenicity_prime": "prime.tsv",
        "immunogenicity_repitope": "repitope.tsv",
    }
    merged = {}
    for col, name in files.items():
        path = os.path.join(base, name)
        if not os.path.exists(path):
            continue
        try:
            tdf = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if "mut_peptide" not in tdf.columns or col not in tdf.columns:
            continue
        for _, r in tdf.iterrows():
            pep = str(r.get("mut_peptide", "")).strip().upper()
            if not pep:
                continue
            if pep not in merged:
                merged[pep] = {}
            v = pd.to_numeric(r.get(col), errors="coerce")
            if pd.notna(v):
                merged[pep][col] = float(v)
    return merged


def _cmd_input_pairs(df: pd.DataFrame, alleles: list) -> pd.DataFrame:
    """生成交叉验证工具可复用的输入对：mut_peptide × hla_allele。"""
    peps = []
    for _, r in df.iterrows():
        p = str(r.get("mut_peptide", "")).strip().upper()
        if 8 <= len(p) <= 14:
            peps.append(p)
    peps = list(dict.fromkeys(peps))
    rows = []
    for p in peps:
        for a in alleles:
            rows.append({"mut_peptide": p, "hla_allele": a})
    return pd.DataFrame(rows)


def _tool_raw_tsv(run_id: str, tool: str) -> str:
    base = os.path.join("results", run_id, "tool_outputs", "raw")
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, f"{tool}.tsv")


def _tool_meta_json(run_id: str, tool: str) -> str:
    base = os.path.join("results", run_id, "tool_outputs", "raw")
    os.makedirs(base, exist_ok=True)
    return os.path.join(base, f"mhc1_cv_{tool}.meta.json")


def _write_mhc1_cv_meta(
    run_id: str,
    tool: str,
    backend_req: str,
    source_used: str,
    rows: int,
    input_rows: int,
) -> None:
    """
    记录交叉验证来源元数据，便于后续定位“为什么是 off/real_*”。
    """
    payload = {
        "run_id": run_id,
        "tool": tool,
        "backend_requested": backend_req,
        "source_used": source_used,
        "rows": int(rows),
        "input_rows": int(input_rows),
        "raw_tsv_path": _tool_raw_tsv(run_id, tool),
    }
    with open(_tool_meta_json(run_id, tool), "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)


def _summarize_mhc1_cv_source(src_netmhcpan: str, src_bigmhc: str) -> str:
    real_tags = {"real_tsv", "real_cmd"}
    tags = [x for x in (src_netmhcpan, src_bigmhc) if x in real_tags]
    if not tags:
        return "off"
    if "real_cmd" in tags:
        return "real_cmd"
    return "real_tsv"


def _summarize_mhc1_cv_tool(cv_n: Any, cv_b: Any, src_netmhcpan: str, src_bigmhc: str) -> str:
    tools = []
    if pd.notna(pd.to_numeric(cv_n, errors="coerce")):
        tools.append("netmhcpan")
    if pd.notna(pd.to_numeric(cv_b, errors="coerce")):
        tools.append("bigmhc")
    if not tools:
        if src_netmhcpan in {"real_tsv", "real_cmd"}:
            tools.append("netmhcpan")
        if src_bigmhc in {"real_tsv", "real_cmd"}:
            tools.append("bigmhc")
    return ",".join(tools)


def _normalize_mhc1_cv_table(tool: str, raw_df: pd.DataFrame) -> pd.DataFrame:
    """
    统一交叉验证表格式：
    - 必须有 mut_peptide
    - hla_allele 可选（缺失时按 peptide 级匹配）
    """
    if "mut_peptide" not in raw_df.columns:
        raise ValueError(f"{tool} 结果缺少 mut_peptide 列。")
    out = pd.DataFrame()
    out["mut_peptide"] = raw_df["mut_peptide"].astype(str).str.strip().str.upper()
    if "hla_allele" in raw_df.columns:
        out["hla_allele"] = raw_df["hla_allele"].astype(str).str.strip()
    else:
        out["hla_allele"] = ""

    if tool == "mhc1_netmhcpan":
        cands = ["mhc1_cv_netmhcpan_nM", "affinity_nM", "affinity", "ic50", "score"]
        metric_col = next((c for c in cands if c in raw_df.columns), None)
        if metric_col is None:
            raise ValueError("mhc1_netmhcpan 结果缺少亲和力列（mhc1_cv_netmhcpan_nM/affinity_nM/affinity/ic50/score）。")
        out["mhc1_cv_netmhcpan_nM"] = pd.to_numeric(raw_df[metric_col], errors="coerce")
    elif tool == "mhc1_bigmhc":
        cands = ["mhc1_cv_bigmhc_score", "bigmhc_score", "presentation_score", "score"]
        metric_col = next((c for c in cands if c in raw_df.columns), None)
        if metric_col is None:
            raise ValueError("mhc1_bigmhc 结果缺少分值列（bigmhc_score/presentation_score/score）。")
        out["mhc1_cv_bigmhc_score"] = pd.to_numeric(raw_df[metric_col], errors="coerce")
    else:
        raise ValueError(f"未知交叉验证工具: {tool}")
    out = out[out["mut_peptide"] != ""].copy()
    return out


def _run_mhc1_cv_cmd(tool: str, run_id: str, input_pairs: pd.DataFrame) -> pd.DataFrame:
    """按环境变量命令模板调用外部工具并读取结果 TSV。"""
    env_key = "MHC1_NETMHCPAN_CMD" if tool == "mhc1_netmhcpan" else "MHC1_BIGMHC_CMD"
    cmd_tpl = os.environ.get(env_key, "").strip()
    if not cmd_tpl:
        raise RuntimeError(f"未设置环境变量 {env_key}")
    tmp_in = _tool_raw_tsv(run_id, f"{tool}_input")
    tmp_out = _tool_raw_tsv(run_id, tool)
    input_pairs.to_csv(tmp_in, sep="\t", index=False, encoding="utf-8")
    cmd = cmd_tpl.format(input_tsv=tmp_in, output_tsv=tmp_out, run_id=run_id)
    ret = subprocess.run(cmd, shell=True, check=False, capture_output=True, text=True)
    if ret.returncode != 0:
        stderr = (ret.stderr or "").strip()
        raise RuntimeError(f"{tool} real_cmd 执行失败: {stderr}")
    if not os.path.exists(tmp_out):
        raise FileNotFoundError(tmp_out)
    rdf = pd.read_csv(tmp_out, sep="\t")
    return _normalize_mhc1_cv_table(tool, rdf)


def _load_mhc1_cv_with_backend(tool: str, run_id: str, backend: str, input_pairs: pd.DataFrame) -> Tuple[pd.DataFrame, str]:
    """
    返回：标准化表 + 实际来源（real_tsv/real_cmd/proxy）
    注：此处 proxy 表示“未接入真实工具，仅占位”，不会影响主排序。
    """
    b = backend.strip().lower()
    if b not in ("auto", "real_tsv", "real_cmd", "off"):
        raise ValueError(f"{tool} backend 非法: {backend}")
    tsv_path = _tool_raw_tsv(run_id, tool)
    if b == "off":
        return pd.DataFrame(), "off"
    if b in ("auto", "real_tsv"):
        if os.path.exists(tsv_path):
            try:
                df = pd.read_csv(tsv_path, sep="\t")
                return _normalize_mhc1_cv_table(tool, df), "real_tsv"
            except Exception as e:
                if b == "real_tsv":
                    raise RuntimeError(f"{tool} real_tsv 读取失败: {tsv_path}; 原因: {e}") from e
                print(f"{tool} auto: real_tsv 不可用，回退。原因: {e}")
        elif b == "real_tsv":
            raise FileNotFoundError(tsv_path)
    if b in ("auto", "real_cmd"):
        try:
            return _run_mhc1_cv_cmd(tool, run_id, input_pairs), "real_cmd"
        except Exception as e:
            if b == "real_cmd":
                raise
            print(f"{tool} auto: real_cmd 不可用，回退。原因: {e}")
    return pd.DataFrame(), "off"


def _build_cv_lookup(cv_df: pd.DataFrame, metric_col: str) -> Tuple[Dict[Tuple[str, str], float], Dict[str, float]]:
    """构建 (peptide, allele) 与 peptide 级别的双层查表。"""
    pair_lookup: Dict[Tuple[str, str], float] = {}
    pep_lookup: Dict[str, float] = {}
    if cv_df.empty or metric_col not in cv_df.columns:
        return pair_lookup, pep_lookup
    for _, r in cv_df.iterrows():
        pep = str(r.get("mut_peptide", "")).strip().upper()
        alle = str(r.get("hla_allele", "")).strip()
        val = pd.to_numeric(r.get(metric_col), errors="coerce")
        if not pep or pd.isna(val):
            continue
        if alle:
            pair_lookup[(pep, alle)] = float(val)
        # peptide 级别记录取最优值：NetMHCpan nM 越低越好，BigMHC score 越高越好
        if pep not in pep_lookup:
            pep_lookup[pep] = float(val)
        else:
            if metric_col.endswith("_nM"):
                pep_lookup[pep] = min(pep_lookup[pep], float(val))
            else:
                pep_lookup[pep] = max(pep_lookup[pep], float(val))
    return pair_lookup, pep_lookup


def _prepare_mhc2_lookup(
    run_id: str,
    hla_json: dict,
    candidate_peptides: list,
    mhc2_backend: str,
) -> Tuple[Optional[Dict[str, Any]], str]:
    """
    按 backend 预计算 (peptide->NetMHCIIpan 行, 或 None)；返回 (lookup 或 None, 实际模式说明)。
    auto：能调用真实工具且 hla 含可转换的 II 类时尝试，否则 None（下游用 proxy）。
    """
    def _load_mhc2_real_tsv(_run_id: str) -> Optional[Dict[str, Any]]:
        """
        读取预计算的 NetMHCIIpan 结果，便于在无二进制环境下也保持真实来源常态化。
        约定路径：results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv
        """
        tsv_path = os.path.join("results", _run_id, "tool_outputs", "raw", "mhc2_netmhciipan.tsv")
        if not os.path.exists(tsv_path):
            return None
        df_tsv = pd.read_csv(tsv_path, sep="\t")
        if "mut_peptide" not in df_tsv.columns:
            raise ValueError(f"mhc2 real_tsv 缺少 mut_peptide 列: {tsv_path}")
        rank_col = None
        for c in ("mhc2_el_rank", "el_rank", "pct_rank_el", "rank"):
            if c in df_tsv.columns:
                rank_col = c
                break
        if rank_col is None:
            raise ValueError(f"mhc2 real_tsv 缺少 EL 排名列（mhc2_el_rank/el_rank/pct_rank_el/rank）: {tsv_path}")
        ba_col = "mhc2_ba_nm" if "mhc2_ba_nm" in df_tsv.columns else ("ba_nm" if "ba_nm" in df_tsv.columns else None)
        allele_col = None
        for c in ("mhc2_class2_allele", "mhc2_allele", "allele", "hla_allele"):
            if c in df_tsv.columns:
                allele_col = c
                break
        lookup: Dict[str, Any] = {}
        for _, row in df_tsv.iterrows():
            pep = str(row.get("mut_peptide", "")).strip().upper()
            if not pep:
                continue
            if pep not in {str(x).strip().upper() for x in candidate_peptides}:
                continue
            el_rank = pd.to_numeric(row.get(rank_col), errors="coerce")
            if pd.isna(el_rank):
                continue
            ba_nm = pd.to_numeric(row.get(ba_col), errors="coerce") if ba_col else np.nan
            lookup[pep] = {
                "mhc2_source": "netmhciipan_tsv",
                "mhc2_el_rank": float(el_rank),
                "mhc2_ba_nm": (None if pd.isna(ba_nm) else float(ba_nm)),
                "mhc2_allele": (str(row.get(allele_col, "")).strip() if allele_col else ""),
            }
        return lookup or None

    b = mhc2_backend.lower().strip()
    if b == "proxy":
        return None, "proxy"
    if b not in ("auto", "real_tsv", "netmhciipan"):
        raise ValueError("mhc2_backend 须为 auto / proxy / real_tsv / netmhciipan")

    if b == "netmhciipan":
        if not resolve_netmhciipan_bin():
            raise FileNotFoundError("mhc2_backend=netmhciipan 但未找到可执行文件，请设置 NETMHCIIPAN_BIN 或 PATH。")
        if not collect_netmhcii_alleles(hla_json):
            raise ValueError("hla_typing 中无可用 II 类（自动转换失败），请补充 HLA-DRB1*… 等或改用 proxy。")
        lu, st = build_peptide_mhc2_lookup(hla_json, candidate_peptides)
        if st != "ok" or not lu:
            raise RuntimeError(f"NetMHCIIpan 未得到可用结果: {st}")
        return lu, "netmhciipan"

    if b == "real_tsv":
        lu = _load_mhc2_real_tsv(run_id)
        if not lu:
            raise FileNotFoundError(
                f"mhc2_backend=real_tsv 但未找到或未解析到可用结果: results/{run_id}/tool_outputs/raw/mhc2_netmhciipan.tsv"
            )
        return lu, "netmhciipan"

    # auto
    if resolve_netmhciipan_bin() and collect_netmhcii_alleles(hla_json):
        lu, st = build_peptide_mhc2_lookup(hla_json, candidate_peptides)
        if st == "ok" and lu:
            return lu, "netmhciipan"
        print(f"NetMHCIIpan（auto）未使用，将回退 proxy。原因: {st}")
    # auto 下优先尝试真实 TSV 回填，避免长期 proxy
    lu = _load_mhc2_real_tsv(run_id)
    if lu:
        return lu, "netmhciipan"
    return None, "proxy"


def main(
    run_id: str,
    w1: float,
    w2: float,
    w3: float,
    w4: float,
    mhc2_backend: str = "auto",
    wi_deepimmuno: float = 1.0,
    wi_prime: float = 1.0,
    wi_repitope: float = 1.0,
    backend_mhc1_netmhcpan: str = "auto",
    backend_mhc1_bigmhc: str = "auto",
    require_real_mhc2: bool = False,
    require_real_mhc1_cv: bool = False,
):
    """
    主流程：
    1) 定位输入输出路径；
    2) 读取候选肽和 HLA；
    3) 对每条突变肽做 MHC-I 亲和力预测；
    4) 可选：NetMHCIIpan 对 II 类（见 --mhc2_backend）；
    5) 按亲和力生成初版 rank_score 并输出 CSV。
    """
    try:
        from mhcflurry import Class1AffinityPredictor
    except ModuleNotFoundError as e:
        raise ModuleNotFoundError(
            "未安装 mhcflurry，请先安装并执行 mhcflurry-downloads fetch 后再运行预测。"
        ) from e

    # 根据 run_id 拼接输入目录与输出目录
    base = os.path.join("deliveries", run_id, "to_immunogen")
    out_dir = os.path.join("results", run_id)
    os.makedirs(out_dir, exist_ok=True)

    csv_path = os.path.join(base, "neoantigen_candidates.csv")
    hla_path = os.path.join(base, "hla_typing.json")

    # 输入文件存在性检查，便于第一时间发现路径或命名问题
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"未找到输入文件: {csv_path}")
    if not os.path.exists(hla_path):
        raise FileNotFoundError(f"未找到输入文件: {hla_path}")

    # 读取候选肽表和 HLA 分型
    df = pd.read_csv(csv_path)
    with open(hla_path, "r", encoding="utf-8") as f:
        hla_json = json.load(f)

    # 标准化 HLA 格式并过滤非法条目
    alleles = [normalize_allele(x) for x in flatten_hla(hla_json)]
    alleles = [x for x in alleles if x is not None]
    if len(alleles) == 0:
        raise ValueError("没有可用的 HLA 等位基因，请检查 hla_typing.json 格式。")

    # 供 NetMHCIIpan：去重后的、长度 9–14 且会进入 I 类预测的肽（与下循环条件一致，略去 8 mer）
    cands_mhc2 = []
    for _, r in df.iterrows():
        p0 = str(r.get("mut_peptide", "")).strip()
        if 8 <= len(p0) <= 14 and len(p0) >= 9:
            cands_mhc2.append(p0)
    cands_mhc2 = list(dict.fromkeys(cands_mhc2))
    mhc2_lookup, mhc2_mode = _prepare_mhc2_lookup(run_id, hla_json, cands_mhc2, mhc2_backend)
    print(f"MHC-II 层: 模式={mhc2_mode}，查表条目数={len(mhc2_lookup) if mhc2_lookup else 0}")
    precomputed_immuno = load_precomputed_immunogenicity(run_id)
    print(f"免疫原性预计算表: 命中肽数={len(precomputed_immuno)}（缺失时自动回退 proxy）")
    cv_input = _cmd_input_pairs(df, alleles)
    cv_netmhcpan_df, cv_netmhcpan_src = _load_mhc1_cv_with_backend(
        tool="mhc1_netmhcpan",
        run_id=run_id,
        backend=backend_mhc1_netmhcpan,
        input_pairs=cv_input,
    )
    cv_bigmhc_df, cv_bigmhc_src = _load_mhc1_cv_with_backend(
        tool="mhc1_bigmhc",
        run_id=run_id,
        backend=backend_mhc1_bigmhc,
        input_pairs=cv_input,
    )
    n_pair, n_pep = _build_cv_lookup(cv_netmhcpan_df, "mhc1_cv_netmhcpan_nM")
    b_pair, b_pep = _build_cv_lookup(cv_bigmhc_df, "mhc1_cv_bigmhc_score")
    _write_mhc1_cv_meta(
        run_id=run_id,
        tool="netmhcpan",
        backend_req=backend_mhc1_netmhcpan,
        source_used=cv_netmhcpan_src,
        rows=len(cv_netmhcpan_df),
        input_rows=len(cv_input),
    )
    _write_mhc1_cv_meta(
        run_id=run_id,
        tool="bigmhc",
        backend_req=backend_mhc1_bigmhc,
        source_used=cv_bigmhc_src,
        rows=len(cv_bigmhc_df),
        input_rows=len(cv_input),
    )
    print(f"MHC-I 交叉验证 NetMHCpan: source={cv_netmhcpan_src}, rows={len(cv_netmhcpan_df)}")
    print(f"MHC-I 交叉验证 BigMHC: source={cv_bigmhc_src}, rows={len(cv_bigmhc_df)}")
    if require_real_mhc2 and mhc2_mode != "netmhciipan":
        raise RuntimeError(
            "已启用 --require_real_mhc2，但当前 MHC-II 不是 netmhciipan（可能回退到了 proxy）。"
        )
    if require_real_mhc1_cv and (cv_netmhcpan_src not in ("real_tsv", "real_cmd")) and (cv_bigmhc_src not in ("real_tsv", "real_cmd")):
        raise RuntimeError(
            "已启用 --require_real_mhc1_cv，但 NetMHCpan/BigMHC 均未使用真实结果（当前 source=off）。"
        )

    # 加载 MHCflurry MHC-I 亲和力模型
    predictor = Class1AffinityPredictor.load()
    rows = []

    # 对每条突变肽与病人全部 HLA 做预测
    for _, row in df.iterrows():
        pep = str(row.get("mut_peptide", "")).strip()
        # MHC-I 常见肽长范围：8~14，超出范围先跳过
        if len(pep) < 8 or len(pep) > 14:
            continue

        # 同一条肽对多个 HLA 同时预测
        pred = predictor.predict_to_dataframe(
            peptides=[pep] * len(alleles),
            alleles=alleles,
            include_percentile_ranks=False,
            include_confidence_intervals=False
        )
        # 逐行收集结果，保留原始输入中的关键字段
        for _, pr in pred.iterrows():
            key = (pep.upper(), str(pr["allele"]))
            cv_n = n_pair.get(key, n_pep.get(pep.upper()))
            cv_b = b_pair.get(key, b_pep.get(pep.upper()))
            mhc1_cv_source = _summarize_mhc1_cv_source(cv_netmhcpan_src, cv_bigmhc_src)
            mhc1_cv_tool = _summarize_mhc1_cv_tool(cv_n, cv_b, cv_netmhcpan_src, cv_bigmhc_src)
            lu = mhc2_lookup.get(pep) if mhc2_lookup else None
            if lu and lu.get("mhc2_source") in ("netmhciipan", "netmhciipan_tsv") and lu.get("mhc2_el_rank") is not None:
                m2s = mhc2_score_from_el_rank(lu["mhc2_el_rank"])
                m2src = "netmhciipan"
                m2el = lu.get("mhc2_el_rank")
                m2ba = lu.get("mhc2_ba_nm")
                m2a2 = lu.get("mhc2_allele", "")
            else:
                m2s = mhc2_proxy_score(pep, str(pr["allele"]))
                m2src = "proxy"
                m2el, m2ba, m2a2 = None, None, ""
            pre = precomputed_immuno.get(pep, {})
            immunogenicity_deepimmuno = pre.get("immunogenicity_deepimmuno", deepimmuno_proxy(pep))
            immunogenicity_prime = pre.get("immunogenicity_prime", prime_proxy(pep))
            immunogenicity_repitope = pre.get("immunogenicity_repitope", repitope_proxy(pep))
            src_deepimmuno = "precomputed" if "immunogenicity_deepimmuno" in pre else "proxy"
            src_prime = "precomputed" if "immunogenicity_prime" in pre else "proxy"
            src_repitope = "precomputed" if "immunogenicity_repitope" in pre else "proxy"
            rows.append({
                "mutation": row.get("mutation", ""),
                "mut_peptide": pep,
                "wt_peptide": row.get("wt_peptide", ""),
                "variant_vaf": row.get("variant_vaf", ""),
                "hla_allele": pr["allele"],
                "affinity_nM": pr.get("affinity", None),
                "mhc1_cv_netmhcpan_nM": cv_n,
                "mhc1_cv_bigmhc_score": cv_b,
                "mhc1_cv_source_netmhcpan": cv_netmhcpan_src,
                "mhc1_cv_source_bigmhc": cv_bigmhc_src,
                "mhc1_cv_source": mhc1_cv_source,
                "mhc1_cv_tool": mhc1_cv_tool,
                "mhc2_score": m2s,
                "mhc2_el_rank": m2el,
                "mhc2_ba_nm": m2ba,
                "mhc2_class2_allele": m2a2,
                "mhc2_backend": m2src,
                # 新契约列：免疫原性分项
                "immunogenicity_deepimmuno": immunogenicity_deepimmuno,
                "immunogenicity_prime": immunogenicity_prime,
                "immunogenicity_repitope": immunogenicity_repitope,
                "immunogenicity_source_deepimmuno": src_deepimmuno,
                "immunogenicity_source_prime": src_prime,
                "immunogenicity_source_repitope": src_repitope,
                # 兼容旧列名（后续逐步淘汰）
                "deepimmuno_score": immunogenicity_deepimmuno,
                "prime_score": immunogenicity_prime,
                "repitope_score": immunogenicity_repitope,
                "wt_peptide_dissimilarity": hamming_dissimilarity(pep, row.get("wt_peptide", "")),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("没有得到预测结果，请检查肽长度、HLA 格式或输入内容。")

    out["affinity_nM"] = pd.to_numeric(out["affinity_nM"], errors="coerce")
    out["mhc1_cv_netmhcpan_nM"] = pd.to_numeric(out["mhc1_cv_netmhcpan_nM"], errors="coerce")
    out["mhc1_cv_bigmhc_score"] = pd.to_numeric(out["mhc1_cv_bigmhc_score"], errors="coerce")
    if "mhc2_el_rank" in out.columns:
        out["mhc2_el_rank"] = pd.to_numeric(out["mhc2_el_rank"], errors="coerce")
    if "mhc2_ba_nm" in out.columns:
        out["mhc2_ba_nm"] = pd.to_numeric(out["mhc2_ba_nm"], errors="coerce")
    out["variant_vaf"] = pd.to_numeric(out["variant_vaf"], errors="coerce").fillna(0.0)
    immuno_weight_sum = wi_deepimmuno + wi_prime + wi_repitope
    if np.isclose(immuno_weight_sum, 0.0):
        raise ValueError("免疫原性子权重和不能为 0，请调整 --wi_* 参数。")
    out["immunogenicity"] = (
        wi_deepimmuno * out["immunogenicity_deepimmuno"]
        + wi_prime * out["immunogenicity_prime"]
        + wi_repitope * out["immunogenicity_repitope"]
    ) / immuno_weight_sum

    # 标准化综合打分：rank_score = w1*HLA_affinity + w2*immunogenicity + w3*VAF + w4*dissimilarity
    out["hla_affinity_norm"] = minmax_series(out["affinity_nM"], reverse=True)
    out["immunogenicity_norm"] = minmax_series(out["immunogenicity"])
    out["vaf_norm"] = minmax_series(out["variant_vaf"])
    out["dissimilarity_norm"] = minmax_series(out["wt_peptide_dissimilarity"])
    out["rank_score"] = (
        w1 * out["hla_affinity_norm"]
        + w2 * out["immunogenicity_norm"]
        + w3 * out["vaf_norm"]
        + w4 * out["dissimilarity_norm"]
    )
    out = out.sort_values("rank_score", ascending=False)

    out_file = os.path.join(out_dir, "peptide_mhc_ranking.csv")
    # utf-8-sig 便于 Windows 下用 Excel 直接打开不乱码
    out.to_csv(out_file, index=False, encoding="utf-8-sig")
    print(f"完成: {out_file}")
    print(
        "后端摘要: "
        f"mhc2_mode={mhc2_mode}, "
        f"mhc1_netmhcpan_source={cv_netmhcpan_src}, "
        f"mhc1_bigmhc_source={cv_bigmhc_src}"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--w1", type=float, default=0.45, help="HLA_affinity 权重")
    parser.add_argument("--w2", type=float, default=0.25, help="immunogenicity 权重")
    parser.add_argument("--w3", type=float, default=0.15, help="VAF 权重")
    parser.add_argument("--w4", type=float, default=0.15, help="wt_dissimilarity 权重")
    parser.add_argument("--wi_deepimmuno", type=float, default=1.0, help="immunogenicity_deepimmuno 子权重")
    parser.add_argument("--wi_prime", type=float, default=1.0, help="immunogenicity_prime 子权重")
    parser.add_argument("--wi_repitope", type=float, default=1.0, help="immunogenicity_repitope 子权重")
    parser.add_argument(
        "--mhc2_backend",
        default="auto",
        choices=["auto", "proxy", "real_tsv", "netmhciipan"],
        help="MHC-II 层: auto=优先 netmhciipan，其次 real_tsv，最后 proxy；real_tsv=读取 results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv",
    )
    parser.add_argument(
        "--backend_mhc1_netmhcpan",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 NetMHCpan 后端：auto/real_tsv/real_cmd/off",
    )
    parser.add_argument(
        "--backend_mhc1_bigmhc",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 BigMHC 后端：auto/real_tsv/real_cmd/off",
    )
    parser.add_argument(
        "--require_real_mhc2",
        action="store_true",
        help="要求 MHC-II 必须使用 netmhciipan；若回退 proxy 则直接报错退出。",
    )
    parser.add_argument(
        "--require_real_mhc1_cv",
        action="store_true",
        help="要求 NetMHCpan/BigMHC 至少一个使用真实结果（real_tsv/real_cmd）；否则报错。",
    )
    args = parser.parse_args()
    main(
        args.run_id,
        args.w1,
        args.w2,
        args.w3,
        args.w4,
        mhc2_backend=args.mhc2_backend,
        wi_deepimmuno=args.wi_deepimmuno,
        wi_prime=args.wi_prime,
        wi_repitope=args.wi_repitope,
        backend_mhc1_netmhcpan=args.backend_mhc1_netmhcpan,
        backend_mhc1_bigmhc=args.backend_mhc1_bigmhc,
        require_real_mhc2=args.require_real_mhc2,
        require_real_mhc1_cv=args.require_real_mhc1_cv,
    )
