#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
BigMHC real_cmd 包装器（KarchinLab BigMHC：src/predict.py）。

与 predict_mhc_ranking 的 MHC1_BIGMHC_CMD 协议一致：
  --input  含列 mut_peptide, hla_allele 的 TSV
  --output  含列 mut_peptide, hla_allele, mhc1_cv_bigmhc_score 的 TSV

依赖：克隆并安装
  https://github.com/KarchinLab/bigmhc
仓库较大（约数 GB），需本机/容器内可运行：
  python <bigmhc>/src/predict.py -i=... -m=el -a=0 -p=1 -c=1 -d=cpu

环境变量（建议）：
- BIGMHC_HOME   BigMHC 仓库根目录（其下需存在 src/predict.py），默认 external_refs/bigmhc
- BIGMHC_MODEL  默认 el（BigMHC EL 呈递模型；亦可 im，依论文场景）
- BIGMHC_DEVICE 默认 cpu（或 all / 0 等，传给 predict.py -d=）
- BIGMHC_MAXBAT  可选，调小以防 OOM（对应 predict.py -b=）
- BIGMHC_PREDICT_PY  若设置则直接使用该 predict.py 绝对路径，覆盖 BIGMHC_HOME
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Tuple

import pandas as pd


def _find_predict_script() -> str:
    p = (os.environ.get("BIGMHC_PREDICT_PY") or "").strip()
    if p and os.path.isfile(p):
        return p
    home = (os.environ.get("BIGMHC_HOME") or "").strip() or os.path.join(
        "external_refs", "bigmhc"
    )
    cand = os.path.join(os.path.abspath(home), "src", "predict.py")
    if os.path.isfile(cand):
        return cand
    raise FileNotFoundError(
        f"未找到 BigMHC predict.py。请克隆 KarchinLab/bigmhc 并设置 BIGMHC_HOME，"
        f"或设置 BIGMHC_PREDICT_PY。期望路径: {cand}"
    )


def mhc_key(s: str) -> str:
    x = re.sub(r"^HLA-", "", (s or "").strip(), flags=re.I)
    m = re.match(r"^([ABC])\*([\d:]+)$", x, re.I)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    m2 = re.match(r"^([ABC])([\d:]+)$", x, re.I)
    if m2:
        return f"{m2.group(1).upper()}{m2.group(2)}"
    return x


def _pick_score_col(df: pd.DataFrame) -> str:
    for c in df.columns:
        u = c.upper()
        if "BIGMHC_EL" in u or c == "BigMHC_EL":
            return c
    for c in df.columns:
        u = c.upper()
        if "PRESENT" in u:
            return c
    for c in df.columns:
        if c.lower() in ("score", "el", "prd") or c.endswith("_EL"):
            return c
    raise ValueError(
        f"BigMHC 输出中未找到分数字段（预期 BigMHC_EL 等），实际列: {list(df.columns)}"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="将 mut_peptide×HLA 表转为 BigMHC EL 分并写本仓库所需 TSV。",
    )
    parser.add_argument("--input", required=True, help="含 mut_peptide, hla_allele")
    parser.add_argument("--output", required=True, help="输出 mhc1_cv_bigmhc_score")
    args = parser.parse_args()

    src = pd.read_csv(args.input, sep="\t")
    if "mut_peptide" not in src.columns or "hla_allele" not in src.columns:
        raise ValueError("input 需含列 mut_peptide, hla_allele")
    src = src.copy()
    src["mut_peptide"] = src["mut_peptide"].astype(str).str.strip().str.upper()
    src["hla_allele"] = src["hla_allele"].astype(str).str.strip()
    if src.empty:
        pd.DataFrame(
            columns=["mut_peptide", "hla_allele", "mhc1_cv_bigmhc_score", "source", "version"]
        ).to_csv(args.output, sep="\t", index=False, encoding="utf-8")
        return

    predict_py = _find_predict_script()
    pdir = os.path.dirname(os.path.dirname(predict_py))
    model = (os.environ.get("BIGMHC_MODEL") or "el").strip()
    dev = (os.environ.get("BIGMHC_DEVICE") or "cpu").strip()
    maxbat = (os.environ.get("BIGMHC_MAXBAT") or "").strip()

    # BigMHC 需要 CSV；官方示例列顺序：allele, peptide + 可带表头
    with tempfile.TemporaryDirectory(prefix="bigmhc_cmd_") as td:
        in_csv = os.path.join(td, "bigmhc_in.csv")
        out_csv = os.path.join(td, "bigmhc_out.csv")
        # 与 README 中 example1 一致：首列 mhc、次列 pep，一行表头
        work = pd.DataFrame(
            {
                "mhc": src["hla_allele"],
                "pep": src["mut_peptide"],
            }
        )
        work.to_csv(in_csv, index=False, encoding="utf-8")
        # -a=0 mhc, -p=1 pep, -c=1 跳一行表头（见 KarchinLab README: columns zero-indexed）
        cmd: List[str] = [
            shutil.which("python3") or shutil.which("python") or "python",
            predict_py,
            f"-i={in_csv}",
            f"-m={model}",
            "-a=0",
            "-p=1",
            "-c=1",
            f"-d={dev}",
            f"-o={out_csv}",
        ]
        if maxbat:
            cmd.append(f"-b={maxbat}")

        wdir = pdir if os.path.isdir(pdir) else None
        run_env = dict(os.environ)
        if dev.lower() in ("cpu", "CPU"):
            run_env.setdefault("CUDA_VISIBLE_DEVICES", "-1")

        proc = subprocess.run(
            cmd, cwd=wdir, env=run_env, capture_output=True, text=True, check=False
        )
        if proc.returncode != 0:
            msg = (proc.stdout or "") + "\n" + (proc.stderr or "")
            raise RuntimeError(
                f"BigMHC predict 失败(退出码 {proc.returncode})。"
                f"请确认已按官方 README 安装依赖、模型文件齐全。\n{msg[:4000]}"
            )

        if not os.path.isfile(out_csv):
            # 回退：部分版本可能仍写 in_csv + .prd
            alt = in_csv + ".prd"
            if os.path.isfile(alt):
                out_csv = alt
            else:
                raise FileNotFoundError(f"未找到 BigMHC 输出: {out_csv}")

        out_df = pd.read_csv(out_csv)
        if "pep" in out_df.columns:
            pcol = "pep"
        elif "Peptide" in out_df.columns:
            pcol = "Peptide"
        else:
            pcol = next(
                (
                    c
                    for c in out_df.columns
                    if c.lower() in ("peptide", "mut_peptide", "p")
                ),
                None,
            )
        if "mhc" in out_df.columns:
            mcol = "mhc"
        elif "MHC" in out_df.columns:
            mcol = "MHC"
        else:
            mcol = next(
                (
                    c
                    for c in out_df.columns
                    if "mhc" in c.lower() or c.lower() in ("allele", "hla")
                ),
                None,
            )
        if pcol is None or mcol is None:
            raise ValueError(f"无法识别肽/等位基因列: {list(out_df.columns)}")

        score_c = _pick_score_col(out_df)
        out_df["mut_peptide"] = out_df[pcol].astype(str).str.strip().str.upper()
        out_df["hla_allele_in"] = out_df[mcol].astype(str).str.strip()
        out_df["mhc1_cv_bigmhc_score"] = pd.to_numeric(out_df[score_c], errors="coerce")

    # 与原始输入行对齐（顺序可能因 BigMHC 重排，故按键合并）
    look: Dict[Tuple[str, str], float] = {}
    for _, row in out_df.iterrows():
        pp = str(row["mut_peptide"]).strip().upper()
        mm = str(row["hla_allele_in"]).strip()
        val = row["mhc1_cv_bigmhc_score"]
        if pd.isna(val):
            continue
        look[(pp, mhc_key(mm))] = float(val)

    rows_out = []
    for _, r in src.iterrows():
        p = r["mut_peptide"]
        h = r["hla_allele"]
        s = look.get((p, mhc_key(h)))
        if s is None or (isinstance(s, float) and pd.isna(s)):
            raise RuntimeError(
                f"未在 BigMHC 结果中匹配: mut_peptide={p}, hla_allele={h} "
                f"(可检查 BigMHC 是否支持该等位/肽长)。"
            )
        rows_out.append(
            {
                "mut_peptide": p,
                "hla_allele": h,
                "mhc1_cv_bigmhc_score": float(s),
                "source": "bigmhc_real_cmd",
                "version": f"KarchinLab-BigMHC-{model}",
            }
        )

    pd.DataFrame(rows_out).to_csv(args.output, sep="\t", index=False, encoding="utf-8")
    print(f"BigMHC real_cmd 完成，{len(rows_out)} 行，输出 {args.output}")


if __name__ == "__main__":
    main()
