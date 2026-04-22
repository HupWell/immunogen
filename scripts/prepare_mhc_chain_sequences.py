# -*- coding: utf-8 -*-
"""
功能：按 hla_typing.json 自动准备结构建模所需的 MHC-I α 链与 β2m 氨基酸序列。

数据来源（默认）：
- IPD-IMGT/HLA（A/B/C 蛋白 FASTA）
- UniProt 人源 β2m（P61769）

输出：
- results/<run_id>/structure_inputs/mhc_chain_sequences.fasta
- results/<run_id>/structure_inputs/mhc_chain_sequences.json
"""
import os
import re
import json
import argparse
from urllib import request
from urllib.error import URLError, HTTPError


IPD_FASTA_URLS = {
    "A": "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/A_prot.fasta",
    "B": "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/B_prot.fasta",
    "C": "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/C_prot.fasta",
}
B2M_FASTA_URL = "https://rest.uniprot.org/uniprotkb/P61769.fasta"


def normalize_hla_allele(raw: str) -> str:
    s = str(raw or "").strip().upper()
    if not s:
        return ""
    if s.startswith("HLA-"):
        return s
    if re.match(r"^[ABC]\*", s):
        return f"HLA-{s}"
    return s


def parse_fasta(text: str):
    seqs = []
    header = None
    chunks = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                seqs.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        seqs.append((header, "".join(chunks).upper()))
    return seqs


def fetch_with_cache(url: str, cache_path: str, refresh: bool):
    if (not refresh) and os.path.exists(cache_path):
        with open(cache_path, "r", encoding="utf-8") as f:
            return f.read(), "cache"
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    try:
        with request.urlopen(url, timeout=40) as resp:
            text = resp.read().decode("utf-8", errors="ignore")
        with open(cache_path, "w", encoding="utf-8") as f:
            f.write(text)
        return text, "remote"
    except (URLError, HTTPError, TimeoutError):
        if os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8") as f:
                return f.read(), "cache_fallback"
        raise


def extract_gene_and_allele(header: str):
    """
    尝试从 IPD FASTA 头提取类似 A*02:01:01:01 的字段。
    """
    m = re.search(r"\b([ABC]\*\d{2,3}:\d{2,3}(?::\d{2,3}){0,2}[A-Z]?)\b", header.upper())
    if not m:
        return "", ""
    token = m.group(1)
    gene = token.split("*", 1)[0]
    return gene, token


def build_alpha_index(ipd_text_by_gene: dict):
    """
    索引键：
    - HLA-A*02:01:01:01（全字段）
    - HLA-A*02:01（两字段）
    """
    index = {}
    for gene, fasta_text in ipd_text_by_gene.items():
        for header, seq in parse_fasta(fasta_text):
            g, allele_full = extract_gene_and_allele(header)
            if g != gene or not allele_full or not seq:
                continue
            parts = allele_full.split(":")
            two_field = ":".join(parts[:2]) if len(parts) >= 2 else allele_full
            k_full = f"HLA-{allele_full}"
            k_two = f"HLA-{two_field}"
            if k_full not in index:
                index[k_full] = seq
            if k_two not in index:
                index[k_two] = seq
    return index


def pick_target_alpha_allele(hla_json: dict):
    """
    当前策略：优先 HLA-A，其次 HLA-B，再 HLA-C；各键取第一个可用等位基因。
    """
    for key in ("HLA-A", "HLA-B", "HLA-C"):
        vals = hla_json.get(key) or []
        if isinstance(vals, list):
            for v in vals:
                norm = normalize_hla_allele(v)
                if norm:
                    return norm
    return ""


def load_json(path: str):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_fasta(path: str, records: list):
    with open(path, "w", encoding="utf-8") as f:
        for seq_id, seq in records:
            f.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i : i + 80] + "\n")


def main(run_id: str, refresh_remote: bool, strict: bool):
    hla_path = os.path.join("deliveries", run_id, "to_immunogen", "hla_typing.json")
    if not os.path.exists(hla_path):
        raise FileNotFoundError(f"未找到输入文件: {hla_path}")
    hla_json = load_json(hla_path)
    target_allele = pick_target_alpha_allele(hla_json)
    if not target_allele:
        raise ValueError("未在 hla_typing.json 中找到可用的 HLA-A/B/C 等位基因。")

    cache_dir = os.path.join("data", "reference", "ipd_imgt_hla_cache")
    ipd_text_by_gene = {}
    source_map = {}
    try:
        for gene, url in IPD_FASTA_URLS.items():
            cache_path = os.path.join(cache_dir, f"{gene}_prot.fasta")
            txt, src = fetch_with_cache(url, cache_path, refresh=refresh_remote)
            ipd_text_by_gene[gene] = txt
            source_map[f"ipd_{gene}"] = src
    except Exception as e:
        if strict:
            raise RuntimeError(f"拉取 IPD-IMGT/HLA 数据失败: {e}")
        print(f"[WARN] 拉取 IPD-IMGT/HLA 数据失败，跳过本步: {e}")
        return

    b2m_cache = os.path.join("data", "reference", "b2m_human_P61769.fasta")
    try:
        b2m_text, b2m_src = fetch_with_cache(B2M_FASTA_URL, b2m_cache, refresh=refresh_remote)
    except Exception as e:
        if strict:
            raise RuntimeError(f"拉取 β2m 序列失败: {e}")
        print(f"[WARN] 拉取 β2m 序列失败，跳过本步: {e}")
        return
    source_map["b2m"] = b2m_src
    b2m_parsed = parse_fasta(b2m_text)
    if not b2m_parsed:
        raise RuntimeError("β2m FASTA 解析失败。")
    b2m_seq = b2m_parsed[0][1]

    alpha_index = build_alpha_index(ipd_text_by_gene)
    alpha_seq = alpha_index.get(target_allele)
    if (not alpha_seq) and "*" in target_allele:
        # 二次容错：尝试两字段匹配
        gene, rest = target_allele.split("*", 1)
        two_field = ":".join(rest.split(":")[:2])
        alpha_seq = alpha_index.get(f"{gene}*{two_field}")
    if not alpha_seq:
        msg = f"未在 IPD 数据中匹配到 {target_allele} 的 α 链序列。"
        if strict:
            raise RuntimeError(msg)
        print(f"[WARN] {msg}")
        return

    out_dir = os.path.join("results", run_id, "structure_inputs")
    os.makedirs(out_dir, exist_ok=True)
    out_fasta = os.path.join(out_dir, "mhc_chain_sequences.fasta")
    out_json = os.path.join(out_dir, "mhc_chain_sequences.json")

    write_fasta(
        out_fasta,
        [
            (f"MHC_ALPHA|{target_allele}|source=IPD-IMGT-HLA", alpha_seq),
            ("B2M|P61769|source=UniProt", b2m_seq),
        ],
    )
    meta = {
        "run_id": run_id,
        "target_alpha_allele": target_allele,
        "output_fasta": out_fasta,
        "alpha_length": len(alpha_seq),
        "b2m_length": len(b2m_seq),
        "data_source": {
            "alpha": "IPD-IMGT/HLA GitHub FASTA",
            "b2m": "UniProt P61769 FASTA",
            "fetch_mode": source_map,
        },
    }
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

    print(f"完成: {out_fasta}")
    print(f"完成: {out_json}")
    print(f"alpha: {target_allele} (len={len(alpha_seq)})")
    print(f"b2m: P61769 (len={len(b2m_seq)})")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--refresh_remote", action="store_true", help="强制刷新远程 FASTA 并更新缓存")
    parser.add_argument("--strict", action="store_true", help="匹配不到序列时直接报错退出")
    args = parser.parse_args()
    main(run_id=args.run_id, refresh_remote=args.refresh_remote, strict=args.strict)

