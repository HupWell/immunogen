# -*- coding: utf-8 -*-
"""
免疫原性适配器公共逻辑（可选真模型 + 自动回退 proxy）。

说明：
- 本模块支持三种来源：
  1) proxy：内置可复现代理分；
  2) real_tsv：读取 results/<run_id>/tool_outputs/raw/<tool>.tsv；
  3) real_cmd：调用外部命令生成 TSV（通过环境变量配置）。
- 各适配器输出统一文件：results/<run_id>/tool_outputs/*.tsv
  列至少包含：mut_peptide + immunogenicity_* + source + version。
"""
import os
import shlex
import subprocess
import tempfile
from typing import Dict, Tuple

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


def output_column(tool_name: str) -> str:
    if tool_name == "deepimmuno":
        return "immunogenicity_deepimmuno"
    if tool_name == "prime":
        return "immunogenicity_prime"
    if tool_name == "repitope":
        return "immunogenicity_repitope"
    raise ValueError(f"未知工具名: {tool_name}")


def proxy_builder(tool_name: str):
    if tool_name == "deepimmuno":
        return deepimmuno_proxy
    if tool_name == "prime":
        return prime_proxy
    if tool_name == "repitope":
        return repitope_proxy
    raise ValueError(f"未知工具名: {tool_name}")


def _normalize_real_table(df: pd.DataFrame, tool_name: str) -> pd.DataFrame:
    """
    真实工具 TSV 归一化为统一列：
    - mut_peptide
    - immunogenicity_<tool>
    - source
    - version
    """
    col = output_column(tool_name)
    if "mut_peptide" not in df.columns:
        raise ValueError("真实工具 TSV 缺少 mut_peptide 列。")
    # 支持常见分数字段别名
    candidate_score_cols = [col, "score", "immunogenicity", "value", "pred", "prediction"]
    score_col = None
    for c in candidate_score_cols:
        if c in df.columns:
            score_col = c
            break
    if score_col is None:
        raise ValueError(f"真实工具 TSV 缺少分数列（可用: {candidate_score_cols}）。")
    out = pd.DataFrame()
    out["mut_peptide"] = (
        df["mut_peptide"].astype(str).str.strip().str.upper().replace("", pd.NA).dropna()
    )
    out[col] = pd.to_numeric(df[score_col], errors="coerce")
    out["source"] = df.get("source", f"{tool_name}_real")
    out["version"] = df.get("version", "unknown")
    out = out.dropna(subset=["mut_peptide", col]).drop_duplicates(subset=["mut_peptide"], keep="first")
    if out.empty:
        raise ValueError("真实工具 TSV 归一化后为空。")
    return out


def _raw_tsv_path(run_id: str, tool_name: str) -> str:
    return os.path.join("results", run_id, "tool_outputs", "raw", f"{tool_name}.tsv")


def _cmd_env_key(tool_name: str) -> str:
    return f"IMMUNO_{tool_name.upper()}_CMD"


def _run_real_cmd(peptides: pd.Series, run_id: str, tool_name: str) -> pd.DataFrame:
    """
    调用外部真模型命令。
    环境变量示例（Windows/POSIX 都可）：
    IMMUNO_DEEPIMMUNO_CMD="python tools/deep_runner.py --input {input_tsv} --output {output_tsv}"
    """
    cmd_tmpl = (os.environ.get(_cmd_env_key(tool_name)) or "").strip()
    if not cmd_tmpl:
        raise RuntimeError(f"未设置 {_cmd_env_key(tool_name)}。")
    with tempfile.TemporaryDirectory(prefix=f"immuno_{tool_name}_") as td:
        in_tsv = os.path.join(td, "input.tsv")
        out_tsv = os.path.join(td, "output.tsv")
        pd.DataFrame({"mut_peptide": peptides}).to_csv(in_tsv, sep="\t", index=False, encoding="utf-8")
        cmd = cmd_tmpl.format(input_tsv=in_tsv, output_tsv=out_tsv, run_id=run_id)
        proc = subprocess.run(
            shlex.split(cmd, posix=False),
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            msg = (proc.stdout or "") + "\n" + (proc.stderr or "")
            raise RuntimeError(f"命令执行失败 code={proc.returncode}: {msg[:1500]}")
        if not os.path.exists(out_tsv):
            raise RuntimeError("命令执行后未产出 output.tsv。")
        rdf = pd.read_csv(out_tsv, sep="\t")
    return _normalize_real_table(rdf, tool_name)


def _try_real_tsv(run_id: str, tool_name: str) -> pd.DataFrame:
    path = _raw_tsv_path(run_id, tool_name)
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    rdf = pd.read_csv(path, sep="\t")
    return _normalize_real_table(rdf, tool_name)


def build_tool_df_with_backend(
    peptides: pd.Series,
    tool_name: str,
    run_id: str,
    backend: str = "auto",
) -> Tuple[pd.DataFrame, str]:
    """
    生成工具输出并返回 (df, backend_used)。
    backend:
    - auto：优先 real_tsv -> real_cmd -> proxy
    - real_tsv：仅读 raw/<tool>.tsv，失败抛错
    - real_cmd：仅调用命令，失败抛错
    - proxy：仅代理
    """
    b = backend.lower().strip()
    if b not in ("auto", "real_tsv", "real_cmd", "proxy"):
        raise ValueError("backend 必须是 auto/real_tsv/real_cmd/proxy")

    if b == "proxy":
        return build_tool_df(peptides, tool_name), "proxy"
    if b == "real_tsv":
        return _try_real_tsv(run_id, tool_name), "real_tsv"
    if b == "real_cmd":
        return _run_real_cmd(peptides, run_id, tool_name), "real_cmd"

    # auto
    try:
        return _try_real_tsv(run_id, tool_name), "real_tsv"
    except Exception:
        pass
    try:
        return _run_real_cmd(peptides, run_id, tool_name), "real_cmd"
    except Exception:
        pass
    return build_tool_df(peptides, tool_name), "proxy"
