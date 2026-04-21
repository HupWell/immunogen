# -*- coding: utf-8 -*-
"""
脚本名称：predict_mhc_ranking.py
主要功能：
1) 读取 BioDriver 提供的候选突变肽和 HLA 分型；
2) 使用 MHCflurry 预测每条突变肽与病人 HLA 的结合亲和力；
3) 生成最基础的排名文件 peptide_mhc_ranking.csv。

输入文件：
- deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv
- deliveries/<run_id>/to_immunogen/hla_typing.json

输出文件：
- results/<run_id>/peptide_mhc_ranking.csv
"""

import os
import json
import argparse
import pandas as pd
from mhcflurry import Class1AffinityPredictor
import numpy as np


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


def deepimmuno_proxy(peptide: str) -> float:
    """DeepImmuno 代理分：偏好多样性与带电残基比例。"""
    pep = peptide.upper()
    uniq = len(set(pep)) / max(len(pep), 1)
    charged = sum(aa in "KRDEH" for aa in pep) / max(len(pep), 1)
    return round(max(0.0, min(1.0, uniq * 0.7 + charged * 0.3)), 4)


def prime_proxy(peptide: str) -> float:
    """PRIME 代理分：偏好锚定位点多样性与芳香族残基。"""
    pep = peptide.upper()
    anchor = pep[1] + pep[-1] if len(pep) >= 2 else pep
    anchor_div = len(set(anchor)) / max(len(anchor), 1)
    aromatic = sum(aa in "FWY" for aa in pep) / max(len(pep), 1)
    return round(max(0.0, min(1.0, anchor_div * 0.6 + aromatic * 0.4)), 4)


def repitope_proxy(peptide: str) -> float:
    """Repitope 代理分：偏好中等疏水与脯氨酸/甘氨酸比例。"""
    pep = peptide.upper()
    hydrophobic = sum(aa in "AILMFWVY" for aa in pep) / max(len(pep), 1)
    pg_ratio = sum(aa in "PG" for aa in pep) / max(len(pep), 1)
    score = (1.0 - abs(hydrophobic - 0.45)) * 0.7 + pg_ratio * 0.3
    return round(max(0.0, min(1.0, score)), 4)


def main(run_id: str, w1: float, w2: float, w3: float, w4: float):
    """
    主流程：
    1) 定位输入输出路径；
    2) 读取候选肽和 HLA；
    3) 对每条突变肽做 MHC-I 亲和力预测；
    4) 按亲和力生成初版 rank_score 并输出 CSV。
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
            rows.append({
                "mutation": row.get("mutation", ""),
                "mut_peptide": pep,
                "wt_peptide": row.get("wt_peptide", ""),
                "variant_vaf": row.get("variant_vaf", ""),
                "hla_allele": pr["allele"],
                "affinity_nM": pr.get("affinity", None),
                "mhc2_score": mhc2_proxy_score(pep, str(pr["allele"])),
                "deepimmuno_score": deepimmuno_proxy(pep),
                "prime_score": prime_proxy(pep),
                "repitope_score": repitope_proxy(pep),
                "wt_peptide_dissimilarity": hamming_dissimilarity(pep, row.get("wt_peptide", "")),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        raise ValueError("没有得到预测结果，请检查肽长度、HLA 格式或输入内容。")

    out["affinity_nM"] = pd.to_numeric(out["affinity_nM"], errors="coerce")
    out["variant_vaf"] = pd.to_numeric(out["variant_vaf"], errors="coerce").fillna(0.0)
    out["immunogenicity"] = (
        out["deepimmuno_score"] + out["prime_score"] + out["repitope_score"]
    ) / 3.0

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
    args = parser.parse_args()
    main(args.run_id, args.w1, args.w2, args.w3, args.w4)
