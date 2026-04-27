#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NetMHCpan 4.x（MHC-I）real_cmd 包装器，与 predict_mhc_ranking 的 MHC1_NETMHCPAN_CMD 协议对接。

入参：--input  含列 mut_peptide, hla_allele 的 TSV（与主流程写出的 mhc1_netmhcpan_input 一致）
出参：--output  含列 mut_peptide, hla_allele, affinity_nM（及可选 mhc1_cv_netmhcpan_nM 同义列）

环境变量（可选）：
- NETMHCPAN_BIN      可执行文件路径；未设则 PATH 中查找 netMHCpan / netMHCpan.pl
- NETMHCPAN_HOME     子进程工作目录（与官方安装包一致，部分环境需要）
- NETMHCPAN_LEN      默认 8,9,10,11（I 类常见肽长；可改为 9 等）
- NETMHCPAN_EXTRA    附加参数，空格分割（如 -t -500）
- NETMHCPAN_TIMEOUT  秒，默认 1200
- NETMHCPAN_BA       默认 1：自动追加 -BA，输出 Aff(nM) 亲合力；设为 0 则不加（此时表头常无 nM，解析会失败）
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import tempfile
from typing import Dict, List, Optional, Tuple

import pandas as pd


def resolve_netmhcpan_bin() -> str:
    """定位 NetMHCpan 可执行文件。"""
    env = (os.environ.get("NETMHCPAN_BIN") or "").strip()
    if env and os.path.isfile(env):
        return env
    if env:
        w = shutil.which(env)
        if w:
            return w
    for name in ("netMHCpan", "netMHCpan.pl"):
        p = shutil.which(name)
        if p:
            return p
    raise FileNotFoundError(
        "未找到 netMHCpan。请设置 NETMHCPAN_BIN 或将可执行文件加入 PATH。"
    )


def to_netmhcpan_class1_allele(allele: str) -> str:
    """
    将仓库内常用写法转为 NetMHCpan 4.x 常见的 -a 形式（HLA-A02:01）。
    入参示例：HLA-A*02:01、A*02:01
    """
    a = (allele or "").strip()
    a = re.sub(r"^HLA-", "", a, flags=re.I)
    if re.match(r"^[ABC]\*\d+:\d+", a, re.I):
        gene, rest = a[0], a[2:]
        return f"HLA-{gene}{rest}"
    if re.match(r"^[ABC]\d+:\d+", a, re.I):
        return f"HLA-{a[0]}{a[1:]}"
    if a.upper().startswith("HLA-"):
        return a if a[0].isupper() else a
    return f"HLA-{a}" if not a.upper().startswith("HLA-") else a


def _find_header_line(lines: List[str]) -> Optional[Tuple[int, List[str]]]:
    for i, line in enumerate(lines):
        if "Peptide" in line and (
            "HLA" in line or "Allele" in line or "MHC" in line or "allele" in line
        ):
            parts = line.split()
            if "Peptide" in parts:
                return i, parts
    return None


def parse_netmhcpan_stdout(stdout: str) -> pd.DataFrame:
    """
    解析 netMHCpan 4.x 标准输出表格（空格分列）。
    取 BA 亲合力 nM（列名含 nM / Aff / IC50 等）。
    """
    lines = stdout.splitlines()
    found = _find_header_line(lines)
    if not found:
        raise ValueError(
            "未解析到含 Peptide 的表头，请确认为 netMHCpan 4.x 且输出完整。"
        )
    hidx, headers = found
    hmap = {name: j for j, name in enumerate(headers)}

    def idx(*names) -> int:
        for n in names:
            if n in hmap:
                return hmap[n]
        return -1

    i_pep = idx("Peptide")
    i_allele = idx("Allele", "MHC", "HLA")
    # 4.2c 带 -BA 时常为 Aff(nM)；无 -BA 时往往只有 Score_EL / %Rank_EL
    i_nm = -1
    for j, h in enumerate(headers):
        if h in ("nM", "IC50", "IC50(nM)", "BA_IC50", "Aff(nM)", "Affinity(nM)"):
            i_nm = j
            break
    if i_nm < 0:
        for j, h in enumerate(headers):
            hl = h.lower()
            if "aff" in hl and "nm" in hl.replace("(", "").replace(")", ""):
                i_nm = j
                break
    if i_nm < 0:
        for j, h in enumerate(headers):
            if "nM" in h and "core" not in h.lower() and "rank" not in h.lower():
                i_nm = j
                break
    if i_nm < 0:
        i_nm = idx("IC50", "NetMHC", "Exp", "Affinity", "BA")

    if i_pep < 0 or i_allele < 0:
        raise ValueError(f"表头缺 Peptide/Allele，实际列: {headers}")
    if i_nm < 0:
        raise ValueError(
            f"表头中未找到亲合力 nM 列，实际列: {headers}。"
            "NetMHCpan 4.2 需加 -BA 才有 Aff(nM)；请使用默认 NETMHCPAN_BA=1 或在 NETMHCPAN_EXTRA 中加入 -BA。"
        )

    rows = []
    for line in lines[hidx + 1 :]:
        s = line.strip()
        if not s or s.startswith("---") or s.startswith("#"):
            continue
        toks = line.split()
        if not toks:
            continue
        try:
            int(toks[0])
        except ValueError:
            continue
        if i_pep >= len(toks) or i_allele >= len(toks) or i_nm >= len(toks):
            continue
        pep = toks[i_pep]
        alle = toks[i_allele]
        try:
            t_nm = toks[i_nm]
            if t_nm in ("NA", "N/A", "nan"):
                continue
            val = float(t_nm)
        except (ValueError, IndexError):
            continue
        rows.append({"peptide": pep, "mhc_allele": alle, "nM": val})
    if not rows:
        raise ValueError("表头后无有效 nM 数据行。请检查肽长是否与 -l 参数匹配。")
    return pd.DataFrame(rows)


def mhc_match_key(allele: str) -> str:
    """
    将不同写法的 I 类 allele 压成同一 key，用于对齐工具输出与输入。
    例：HLA-A*02:01、HLA-A02:01、A*02:01 -> A02:01
    """
    x = (allele or "").strip()
    x = re.sub(r"^HLA-", "", x, flags=re.I)
    m = re.match(r"^([ABC])\*([\d:]+)$", x, re.I)
    if m:
        return f"{m.group(1).upper()}{m.group(2)}"
    m2 = re.match(r"^([ABC])([\d:]+)$", x, re.I)
    if m2:
        return f"{m2.group(1).upper()}{m2.group(2)}"
    return x


def run_netMHCpan_batch(
    peptides: List[str],
    alleles_netm: List[str],
) -> str:
    """一次调用：所有肽、所有 I 类 allele（逗号分隔）。"""
    exe = resolve_netmhcpan_bin()
    wdir = (os.environ.get("NETMHCPAN_HOME") or "").strip() or None
    to_sec = int((os.environ.get("NETMHCPAN_TIMEOUT") or "1200").strip())
    lens = (os.environ.get("NETMHCPAN_LEN") or "8,9,10,11").strip()
    extra = [x for x in (os.environ.get("NETMHCPAN_EXTRA") or "").split() if x]
    a_str = ",".join(alleles_netm)
    env = {**os.environ, "LANG": "C", "LC_ALL": "C"}
    # 4.2 默认无 BA 列；必须 -BA 才有 Aff(nM)（与 netMHCIIpan 包装器一致）
    use_ba = (os.environ.get("NETMHCPAN_BA", "1") or "1").strip().lower() not in (
        "0",
        "false",
        "no",
        "n",
    )
    extra_upper = {x.upper() for x in extra}
    ba_in_extra = "-BA" in extra_upper

    with tempfile.TemporaryDirectory(prefix="nmp1_") as td:
        pfile = os.path.join(td, "in.pep")
        with open(pfile, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(peptides) + "\n")
        cmd = [exe, "-f", pfile, "-a", a_str, "-l", lens, "-inptype", "1"]
        if use_ba and not ba_in_extra:
            cmd.append("-BA")
        cmd = cmd + extra
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=wdir,
            env=env,
            timeout=to_sec,
        )
    out = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if proc.returncode != 0:
        raise RuntimeError(
            f"netMHCpan 退出码 {proc.returncode}。输出前 2000 字符：\n{out[:2000]}"
        )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="从 mut_peptide×HLA 表调用 netMHCpan，写 affinity_nM TSV。",
    )
    parser.add_argument("--input", required=True, help="含 mut_peptide, hla_allele 的 TSV")
    parser.add_argument(
        "--output", required=True, help="输出 TSV：mut_peptide, hla_allele, affinity_nM"
    )
    args = parser.parse_args()

    src = pd.read_csv(args.input, sep="\t")
    if "mut_peptide" not in src.columns or "hla_allele" not in src.columns:
        raise ValueError("input 需包含列 mut_peptide, hla_allele")
    src["mut_peptide"] = src["mut_peptide"].astype(str).str.strip().str.upper()
    src["hla_allele"] = src["hla_allele"].astype(str).str.strip()

    # 去重后的 (pep, 原始 hla) 与 netMHCpan allele 形式
    rows_u = src.drop_duplicates(subset=["mut_peptide", "hla_allele"])
    if rows_u.empty:
        pd.DataFrame(
            columns=["mut_peptide", "hla_allele", "affinity_nM", "mhc1_cv_netmhcpan_nM"]
        ).to_csv(args.output, sep="\t", index=False, encoding="utf-8")
        print("无输入行，已写空 TSV。")
        return

    peps = list(dict.fromkeys(rows_u["mut_peptide"].tolist()))
    alleles_raw = list(dict.fromkeys(rows_u["hla_allele"].tolist()))
    alleles_nm = [to_netmhcpan_class1_allele(x) for x in alleles_raw]
    # 同一 netM 形式去重
    a_unique: List[str] = []
    for a in alleles_nm:
        if a not in a_unique:
            a_unique.append(a)

    if len(a_unique) > 50:
        raise ValueError("一次 -a 中 allele 数量过多，请分批或联系维护者扩展脚本。")

    raw_stdout = run_netMHCpan_batch(peps, a_unique)
    pred = parse_netmhcpan_stdout(raw_stdout)
    # (peptide, allele_key) -> nM（后写覆盖先写，同 key 应一致）
    lookup: Dict[Tuple[str, str], float] = {}
    for _, pr in pred.iterrows():
        p1 = str(pr["peptide"]).strip().upper()
        ak = mhc_match_key(str(pr["mhc_allele"]))
        lookup[(p1, ak)] = float(pr["nM"])

    out_rows = []
    for _, r in rows_u.iterrows():
        p = r["mut_peptide"]
        h_in = r["hla_allele"]
        h_key = mhc_match_key(h_in)
        m = lookup.get((p, h_key))
        if m is None:
            h_alt = mhc_match_key(to_netmhcpan_class1_allele(h_in))
            m = lookup.get((p, h_alt))
        if m is None:
            raise RuntimeError(
                f"未在 netMHCpan 结果中匹配: peptide={p}, hla_allele={h_in} "
                f"(key={h_key})。请检查 NETMHCPAN_LEN 是否含该肽长，或 HLA 写法是否与 netMHCpan 一致。"
            )
        out_rows.append(
            {
                "mut_peptide": p,
                "hla_allele": h_in,
                "affinity_nM": m,
                "mhc1_cv_netmhcpan_nM": m,
            }
        )

    out = pd.DataFrame(out_rows)
    out.to_csv(args.output, sep="\t", index=False, encoding="utf-8")
    print(f"netMHCpan class1 完成，{len(out)} 行，输出 {args.output}")


if __name__ == "__main__":
    main()
