# -*- coding: utf-8 -*-
"""
NetMHCIIpan 调用层：subprocess 运行可执行文件、解析标准输出表格。

- 等位基因：自动将 **HLA-DRB1*15:01** 转为 **DRB1_1501**（与 openvax mhctools 中 DRB 单链习惯一致）；
- DQ/DP 双链需写成工具认可的一条名字（可用手工映射 + `data/hla_allele_map_netmhciipan.json`），未自动配对。

环境变量（可选）：
- NETMHCIIPAN_BIN   可执行文件；未设则在 PATH 中找 netMHCIIpan / netMHCIIpan.pl
- NETMHCIIPAN_HOME  工作目录（多数字符包需从安装目录起算路径）
- NETMHCIIPAN_BA=0  不追加 -BA（默认追加 -BA 以便解析 nM 列）
- NETMHCIIPAN_EXTRA  附加参数，空格分隔
- NETMHCIIPAN_TIMEOUT  秒，默认 600
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile
from typing import List, Optional, Tuple

import pandas as pd

from hla_allele_to_netmhciipan import load_mapping, to_netmhciipan
from hla_typing_spec import flatten_hla_class_ii

_DRB_RE = re.compile(
    r"^(?:HLA-)?(DRB1)\*(\d{1,3}):(\d{1,3})(?::\d+)*$",
    re.IGNORECASE,
)

# 行内已含 DQA1…-DQB1… 等
_COMBO = re.compile(
    r"^(HLA-)?(DPA1|DQA1)\*.*-(DPB1|DQB1)",
    re.IGNORECASE,
)


def resolve_netmhciipan_bin() -> Optional[str]:
    env = (os.environ.get("NETMHCIIPAN_BIN") or "").strip()
    if env and os.path.isfile(env):
        return env
    if env:
        w = shutil.which(env)
        if w:
            return w
    for name in ("netMHCIIpan", "netMHCIIpan.pl"):
        p = shutil.which(name)
        if p:
            return p
    return None


def ii_allele_to_netmhcii_token(allele: str) -> Optional[str]:
    """
    将 hla_typing 中一条 II 类转为 -a 中的单个 token；无法识别则返回 None。
    """
    _map = load_mapping()
    s0 = to_netmhciipan(allele, _map)
    s = s0.strip()
    if re.match(r"^DRB1_\d{3,4}$", s, re.I):
        return s.upper()

    m1 = _DRB_RE.match(s)
    if m1:
        fam, code = m1.group(2), m1.group(3)
        return f"DRB1_{int(fam):02d}{int(code):02d}"

    if _COMBO.search(s) or (s.count("-") >= 1 and "DQA1" in s.upper() and "DQB1" in s.upper()):
        if not s.upper().startswith("HLA-"):
            return "HLA-" + s
        return s
    if "DPA1" in s.upper() and "DPB1" in s.upper() and "-" in s:
        if not s.upper().startswith("HLA-"):
            return "HLA-" + s
        return s
    return None


def collect_netmhcii_alleles(hla_json: dict) -> List[str]:
    seen = set()
    out: List[str] = []
    for a in flatten_hla_class_ii(hla_json):
        tok = ii_allele_to_netmhcii_token(a)
        if tok and tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def _find_header_line(lines: List[str]) -> Optional[Tuple[int, List[str]]]:
    for i, line in enumerate(lines):
        if "Peptide" in line and ("MHC" in line or "allele" in line.lower()):
            parts = line.split()
            if "Peptide" in parts:
                return i, parts
    return None


def parse_netmhciipan_stdout(stdout: str) -> pd.DataFrame:
    """
    解析 NetMHCIIpan 4.x 标准输出。列名来自表头，兼容列顺序差异。
    """
    lines = stdout.splitlines()
    found = _find_header_line(lines)
    if not found:
        raise ValueError("未找到含 Peptide / MHC 的表头，请确认可执行文件为 NetMHCIIpan 4.x 且已输出预测表。")
    hidx, headers = found
    hmap = {name: j for j, name in enumerate(headers)}

    def idx(*names) -> int:
        for n in names:
            if n in hmap:
                return hmap[n]
        return -1

    i_pep = idx("Peptide")
    i_mhc = idx("MHC")
    i_rkel = idx("%Rank_EL", "EL_Rank", "%Rank")
    i_nm = idx("Affinity(nM)", "nM")

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
        if i_pep < 0 or i_mhc < 0 or i_pep >= len(toks) or i_mhc >= len(toks):
            continue
        peptide, mhc = toks[i_pep], toks[i_mhc]
        r_el = None
        if i_rkel >= 0 and i_rkel < len(toks):
            try:
                t = toks[i_rkel]
                if t not in ("NA", "N/A"):
                    r_el = float(t)
            except ValueError:
                pass
        ba_nm = None
        if i_nm >= 0 and i_nm < len(toks):
            try:
                t = toks[i_nm]
                if t not in ("NA", "N/A"):
                    ba_nm = float(t)
            except ValueError:
                pass
        rows.append(
            {
                "peptide": peptide,
                "mhc": mhc,
                "pct_rank_el": r_el,
                "ba_nm": ba_nm,
            }
        )
    if not rows:
        raise ValueError("表头后无有效数据行。")
    return pd.DataFrame(rows)


def run_netmhciipan(
    peptide_list: List[str],
    allele_tokens: List[str],
) -> str:
    """
    单次子进程：PEPTIDE 模式（-inptype 1），返回合并后的 stdout+stderr 文本供解析。
    """
    exe = resolve_netmhciipan_bin()
    if not exe:
        raise FileNotFoundError(
            "未找到 NetMHCIIpan。请设置 NETMHCIIPAN_BIN 或将可执行文件放入 PATH。"
        )
    if not allele_tokens:
        raise ValueError("等位基因为空。")
    if len(allele_tokens) > 20:
        raise ValueError("一次最多 20 个 MHC 等位基因（工具限制）。")

    peps = [p.strip().upper() for p in peptide_list if p and len(p.strip()) >= 9]
    if not peps:
        raise ValueError("无长度≥9 的肽段，NetMHCIIpan PEPTIDE 模式无法计算。")

    a_str = ",".join(allele_tokens)
    extra: List[str] = []
    if (os.environ.get("NETMHCIIPAN_BA", "1").strip() not in ("0", "false", "no", "N")):
        extra.append("-BA")
    extra += (os.environ.get("NETMHCIIPAN_EXTRA") or "").split()

    wdir = (os.environ.get("NETMHCIIPAN_HOME") or "").strip() or None
    to_sec = int(os.environ.get("NETMHCIIPAN_TIMEOUT", "600"))
    env = {**os.environ, "LANG": "C", "LC_ALL": "C"}

    with tempfile.TemporaryDirectory(prefix="nm2_") as td:
        pfile = os.path.join(td, "in.pep")
        with open(pfile, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(peps) + "\n")
        cmd = [exe, "-f", pfile, "-a", a_str, "-inptype", "1"] + extra
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
            f"NetMHCIIpan 退出码 {proc.returncode}。输出前 2000 字符：\n{out[:2000]}"
        )
    return out


def build_peptide_mhc2_lookup(
    hla_json: dict,
    candidate_peptides: List[str],
) -> Tuple[dict, str]:
    """
    对去重后的候选肽运行 NetMHCIIpan，返回 (peptide -> 字段字典, 状态说明)。
    字典项：mhc2_el_rank（越小越好，取所有等位基因中最低 rank）、mhc2_ba_nm、mhc2_allele、source=netmhciipan
    失败时返回 ({}, 错误原因) 由上层回退到 proxy。
    """
    alleles = collect_netmhcii_alleles(hla_json)
    if not alleles:
        return {}, "no_class2_alleles"
    upep = list(dict.fromkeys([p for p in candidate_peptides if p and len(p) >= 9]))
    if not upep:
        return {}, "no_peptide_len_ge_9"
    try:
        raw = run_netmhciipan(upep, alleles)
    except (FileNotFoundError, ValueError, RuntimeError) as e:
        return {}, str(e)

    try:
        df = parse_netmhciipan_stdout(raw)
    except Exception as e:
        return {}, f"parse: {e}"

    lookup: dict = {}
    for pep in upep:
        sub = df[df["peptide"] == pep]
        if sub.empty:
            sub = df[df["peptide"].str.upper() == pep.upper()]
        if sub.empty:
            continue
        sub2 = sub.dropna(subset=["pct_rank_el"])
        if not sub2.empty:
            best = sub2.loc[sub2["pct_rank_el"].idxmin()]
        else:
            best = sub.iloc[0]
        pr = best["pct_rank_el"]
        bnm = best.get("ba_nm")
        lookup[pep] = {
            "mhc2_el_rank": None if pr is None or pd.isna(pr) else float(pr),
            "mhc2_ba_nm": None if bnm is None or pd.isna(bnm) else float(bnm),
            "mhc2_allele": str(best.get("mhc", "")),
            "mhc2_source": "netmhciipan",
        }
    if not lookup:
        return {}, "no_overlap_peptide"
    return lookup, "ok"
