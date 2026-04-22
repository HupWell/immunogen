# -*- coding: utf-8 -*-
"""
免疫原性适配器公共逻辑（当前为可复现代理实现）。

说明：
- 本模块统一维护三个工具的“占位可运行”评分函数，便于后续替换为真实模型推理。
- 各独立脚本会输出 results/<run_id>/tool_outputs/*.tsv，主流程按 run_id 合并。
"""
import os
from typing import Dict

import pandas as pd


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


def load_unique_peptides(run_id: str) -> pd.Series:
    """读取 run_id 对应输入，返回去重后的 mut_peptide 序列。"""
    path = os.path.join("deliveries", run_id, "to_immunogen", "neoantigen_candidates.csv")
    if not os.path.exists(path):
        raise FileNotFoundError(f"未找到输入文件: {path}")
    df = pd.read_csv(path)
    if "mut_peptide" not in df.columns:
        raise ValueError("neoantigen_candidates.csv 缺少 mut_peptide 列。")
    peps = (
        df["mut_peptide"]
        .astype(str)
        .str.strip()
        .str.upper()
        .replace("", pd.NA)
        .dropna()
        .drop_duplicates()
    )
    if peps.empty:
        raise ValueError("mut_peptide 无有效值。")
    return peps


def ensure_tool_output_dir(run_id: str) -> str:
    out_dir = os.path.join("results", run_id, "tool_outputs")
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def build_tool_df(peptides: pd.Series, tool_name: str) -> pd.DataFrame:
    """按工具名生成标准输出 DataFrame。"""
    if tool_name == "deepimmuno":
        fn = deepimmuno_proxy
        col = "immunogenicity_deepimmuno"
    elif tool_name == "prime":
        fn = prime_proxy
        col = "immunogenicity_prime"
    elif tool_name == "repitope":
        fn = repitope_proxy
        col = "immunogenicity_repitope"
    else:
        raise ValueError(f"未知工具名: {tool_name}")

    out = pd.DataFrame({"mut_peptide": peptides})
    out[col] = out["mut_peptide"].map(fn)
    out["source"] = f"{tool_name}_proxy"
    out["version"] = "proxy_v1"
    return out


def output_specs() -> Dict[str, str]:
    """返回工具名到文件名的映射。"""
    return {
        "deepimmuno": "deepimmuno.tsv",
        "prime": "prime.tsv",
        "repitope": "repitope.tsv",
    }
