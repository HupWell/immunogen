#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DeepImmuno real_cmd 包装器。

输入：--input  指向包含 mut_peptide 列的 TSV
输出：--output 指向包含 mut_peptide + score 的 TSV

说明：
1) 调用官方仓库 external_refs/DeepImmuno-main/deepimmuno-cnn.py（multiple 模式）；
2) DeepImmuno-CNN 原生支持 9/10 mer；8mer 会走扩展近似（A+pep / pep+A 后取均值）；
3) 为兼容本仓库 real_cmd 协议，输出统一列名 score（后续会被适配器归一化）。
"""
import argparse
import os
import subprocess
import tempfile

import pandas as pd


def _normalize_hla(hla: str) -> str:
    """将常见 HLA 写法转换为 DeepImmuno 更稳的格式（例如 HLA-A*0201）。"""
    x = (hla or "").strip().upper().replace("HLA-", "")
    if "*" in x and ":" in x:
        gene, rest = x.split("*", 1)
        a, b = rest.split(":", 1)
        return f"HLA-{gene}*{a}{b}"
    if x and not x.startswith("HLA-"):
        return f"HLA-{x}"
    return hla


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="输入 TSV（需含 mut_peptide）")
    parser.add_argument("--output", required=True, help="输出 TSV（mut_peptide + score）")
    args = parser.parse_args()

    repo_root = os.path.abspath(
        os.environ.get("DEEPIMMUNO_REPO", os.path.join("external_refs", "DeepImmuno-main"))
    )
    script_path = os.path.join(repo_root, "deepimmuno-cnn.py")
    if not os.path.exists(script_path):
        raise FileNotFoundError(f"未找到 DeepImmuno 脚本: {script_path}")

    src = pd.read_csv(args.input, sep="\t")
    if "mut_peptide" not in src.columns:
        raise ValueError("输入文件缺少 mut_peptide 列。")
    src["mut_peptide"] = src["mut_peptide"].astype(str).str.strip().str.upper()

    # DeepImmuno-CNN 官方仅支持 9/10 mer。
    # 为避免 8mer 全量回退 proxy，这里对 8mer 做可追溯的“最小扩展近似”：
    # - 生成 A+pep 与 pep+A 两个 9mer
    # - 分别跑 DeepImmuno，再对同一原始 8mer 取均值
    # 该策略会在 source/version 中明确标注，便于下游审计。
    base = src.drop_duplicates(subset=["mut_peptide"]).copy()
    records = []
    for pep in base["mut_peptide"].tolist():
        n = len(pep)
        if n in (9, 10):
            records.append({"query_peptide": pep, "origin_peptide": pep, "mode": "native"})
        elif n == 8:
            records.append({"query_peptide": f"A{pep}", "origin_peptide": pep, "mode": "expand8_leftA"})
            records.append({"query_peptide": f"{pep}A", "origin_peptide": pep, "mode": "expand8_rightA"})
    runnable = pd.DataFrame(records)
    default_hla = _normalize_hla(os.environ.get("DEEPIMMUNO_DEFAULT_HLA", "HLA-A*0201"))
    if runnable.empty:
        out = pd.DataFrame(columns=["mut_peptide", "score", "source", "version"])
        out.to_csv(args.output, sep="\t", index=False, encoding="utf-8")
        print("DeepImmuno: 无可运行肽（仅支持 9/10 mer），已输出空结果。")
        return

    with tempfile.TemporaryDirectory(prefix="deepimmuno_cmd_") as td:
        in_csv = os.path.join(td, "deepimmuno_input.csv")
        out_dir = os.path.join(td, "out")
        os.makedirs(out_dir, exist_ok=True)
        pd.DataFrame({"peptide": runnable["query_peptide"], "HLA": default_hla}).to_csv(
            in_csv, index=False, header=False
        )

        cmd = [
            "python",
            script_path,
            "--mode",
            "multiple",
            "--intdir",
            in_csv,
            "--outdir",
            out_dir,
        ]
        # 强制 CPU 跑 DeepImmuno，避免服务器 CUDA/cuDNN 版本不匹配导致真实推理失败。
        run_env = dict(os.environ)
        run_env["CUDA_VISIBLE_DEVICES"] = "-1"
        proc = subprocess.run(
            cmd,
            cwd=repo_root,
            env=run_env,
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0:
            msg = (proc.stdout or "") + "\n" + (proc.stderr or "")
            raise RuntimeError(f"DeepImmuno 运行失败: {msg[:2000]}")

        pred_file = os.path.join(out_dir, "deepimmuno-cnn-result.txt")
        if not os.path.exists(pred_file):
            raise FileNotFoundError(f"未找到输出文件: {pred_file}")

        pred = pd.read_csv(pred_file, sep="\t")
        if "peptide" not in pred.columns or "immunogenicity" not in pred.columns:
            raise ValueError("DeepImmuno 输出缺少 peptide/immunogenicity 列。")

        pred["query_peptide"] = pred["peptide"].astype(str).str.strip().str.upper()
        merged = runnable.merge(pred[["query_peptide", "immunogenicity"]], on="query_peptide", how="left")
        merged["score"] = pd.to_numeric(merged["immunogenicity"], errors="coerce")
        out = (
            merged.groupby("origin_peptide", as_index=False)["score"].mean()
            .rename(columns={"origin_peptide": "mut_peptide"})
        )
        out["source"] = "deepimmuno_real_cmd"
        if (merged["mode"] != "native").any():
            out["version"] = "DeepImmuno-CNN(8mer_expand_A_mean)"
        else:
            out["version"] = "DeepImmuno-CNN"
        out = out.dropna(subset=["mut_peptide", "score"]).drop_duplicates(subset=["mut_peptide"], keep="first")
        out.to_csv(args.output, sep="\t", index=False, encoding="utf-8")
        print(f"DeepImmuno real_cmd 完成，输出 {len(out)} 条。")


if __name__ == "__main__":
    main()

