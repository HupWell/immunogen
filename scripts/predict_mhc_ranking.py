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
from typing import Any, Dict, Optional, Tuple

import pandas as pd
from mhcflurry import Class1AffinityPredictor
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


def _prepare_mhc2_lookup(
    hla_json: dict,
    candidate_peptides: list,
    mhc2_backend: str,
) -> Tuple[Optional[Dict[str, Any]], str]:
    """
    按 backend 预计算 (peptide->NetMHCIIpan 行, 或 None)；返回 (lookup 或 None, 实际模式说明)。
    auto：能调用真实工具且 hla 含可转换的 II 类时尝试，否则 None（下游用 proxy）。
    """
    b = mhc2_backend.lower().strip()
    if b == "proxy":
        return None, "proxy"
    if b not in ("auto", "netmhciipan"):
        raise ValueError("mhc2_backend 须为 auto / proxy / netmhciipan")

    if b == "netmhciipan":
        if not resolve_netmhciipan_bin():
            raise FileNotFoundError("mhc2_backend=netmhciipan 但未找到可执行文件，请设置 NETMHCIIPAN_BIN 或 PATH。")
        if not collect_netmhcii_alleles(hla_json):
            raise ValueError("hla_typing 中无可用 II 类（自动转换失败），请补充 HLA-DRB1*… 等或改用 proxy。")
        lu, st = build_peptide_mhc2_lookup(hla_json, candidate_peptides)
        if st != "ok" or not lu:
            raise RuntimeError(f"NetMHCIIpan 未得到可用结果: {st}")
        return lu, "netmhciipan"

    # auto
    if resolve_netmhciipan_bin() and collect_netmhcii_alleles(hla_json):
        lu, st = build_peptide_mhc2_lookup(hla_json, candidate_peptides)
        if st == "ok" and lu:
            return lu, "netmhciipan"
        print(f"NetMHCIIpan（auto）未使用，将回退 proxy。原因: {st}")
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
):
    """
    主流程：
    1) 定位输入输出路径；
    2) 读取候选肽和 HLA；
    3) 对每条突变肽做 MHC-I 亲和力预测；
    4) 可选：NetMHCIIpan 对 II 类（见 --mhc2_backend）；
    5) 按亲和力生成初版 rank_score 并输出 CSV。
    """
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
    mhc2_lookup, mhc2_mode = _prepare_mhc2_lookup(hla_json, cands_mhc2, mhc2_backend)
    print(f"MHC-II 层: 模式={mhc2_mode}，查表条目数={len(mhc2_lookup) if mhc2_lookup else 0}")
    precomputed_immuno = load_precomputed_immunogenicity(run_id)
    print(f"免疫原性预计算表: 命中肽数={len(precomputed_immuno)}（缺失时自动回退 proxy）")

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
            lu = mhc2_lookup.get(pep) if mhc2_lookup else None
            if lu and lu.get("mhc2_source") == "netmhciipan" and lu.get("mhc2_el_rank") is not None:
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
        choices=["auto", "proxy", "netmhciipan"],
        help="MHC-II 层: auto=有工具与 II 等位则 NetMHCIIpan 否则代理; proxy=仅代理; netmhciipan=强制（失败则报错）",
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
    )
