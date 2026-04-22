# -*- coding: utf-8 -*-
"""
功能：将入选肽段串联并生成多价 mRNA 序列。
输入：results/<run_id>/selected_peptides.csv
输出：results/<run_id>/mrna_vaccine.fasta, results/<run_id>/mrna_design.json
"""
import os
import json
import argparse
import subprocess
import pandas as pd

SIGNAL_PEPTIDE_PRESETS = {
    # 常用分泌信号肽示例（教学/流程用途）
    "igkv": "METDTLLLWVLLLWVPGSTG",
    # 任务书提到 MITD，此处给出流程占位序列（生产请以实验方案为准）
    "mitd": "MASSLAFLLLLLLQVSAA",
}

TM_DOMAIN_PRESETS = {
    # 常见跨膜锚定占位（流程/接口用途）
    "cd8a_tm": "LFWLFFPVVLTVD",
}


CODON_TABLE = {
    "A": "GCT", "R": "CGT", "N": "AAT", "D": "GAT", "C": "TGT",
    "Q": "CAA", "E": "GAA", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTT", "P": "CCT",
    "S": "TCT", "T": "ACT", "W": "TGG", "Y": "TAT", "V": "GTT",
    "*": "TAA",
}

PREFERRED_CODONS = {
    "A": ["GCC", "GCT", "GCA", "GCG"],
    "R": ["CGC", "CGT", "AGA", "AGG", "CGA", "CGG"],
    "N": ["AAC", "AAT"],
    "D": ["GAC", "GAT"],
    "C": ["TGC", "TGT"],
    "Q": ["CAG", "CAA"],
    "E": ["GAG", "GAA"],
    "G": ["GGC", "GGT", "GGA", "GGG"],
    "H": ["CAC", "CAT"],
    "I": ["ATC", "ATT", "ATA"],
    "L": ["CTG", "CTC", "CTT", "CTA", "TTA", "TTG"],
    "K": ["AAG", "AAA"],
    "M": ["ATG"],
    "F": ["TTC", "TTT"],
    "P": ["CCC", "CCT", "CCA", "CCG"],
    "S": ["AGC", "TCC", "TCT", "AGT", "TCA", "TCG"],
    "T": ["ACC", "ACT", "ACA", "ACG"],
    "W": ["TGG"],
    "Y": ["TAC", "TAT"],
    "V": ["GTG", "GTC", "GTT", "GTA"],
    "*": ["TAA", "TAG", "TGA"],
}


def aa_to_dna(aa_seq: str) -> str:
    """
    将氨基酸序列按简单密码子表转换为 DNA 序列。
    这是 MVP 版本的密码子优化，占位可跑通，后续可替换为 LinearDesign 等工具。
    """
    dna = []
    for aa in aa_seq:
        if aa not in CODON_TABLE:
            raise ValueError(f"发现不支持的氨基酸字符: {aa}")
        dna.append(CODON_TABLE[aa])
    return "".join(dna)


def aa_to_dna_optimized(aa_seq: str, gc_target: float = 0.55) -> str:
    """
    简化高级密码子优化：
    - 在常用密码子集合中选择更接近目标 GC 的密码子
    - 保留原有 ORF 氨基酸不变
    """
    seq = []
    for aa in aa_seq:
        choices = PREFERRED_CODONS.get(aa)
        if not choices:
            raise ValueError(f"发现不支持的氨基酸字符: {aa}")
        # 选择 GC 含量最接近目标的密码子
        scored = []
        for codon in choices:
            gc = (codon.count("G") + codon.count("C")) / 3.0
            scored.append((abs(gc - gc_target), -gc, codon))
        scored.sort(key=lambda x: (x[0], x[1]))
        seq.append(scored[0][2])
    return "".join(seq)


def run_lineardesign_if_available(lineardesign_bin: str, aa_seq: str) -> str:
    """
    如果提供了 LinearDesign 可执行文件，则调用其输出优化序列。
    返回 DNA/RNA 字母序列；失败时抛异常。
    """
    if not lineardesign_bin:
        raise ValueError("未提供 LinearDesign 路径。")
    if not os.path.exists(lineardesign_bin):
        raise FileNotFoundError(f"LinearDesign 可执行文件不存在: {lineardesign_bin}")

    cmd = [lineardesign_bin, aa_seq]
    result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"LinearDesign 执行失败: {result.stderr.strip()}")
    seq = result.stdout.strip().splitlines()[-1].strip().upper()
    if not seq:
        raise RuntimeError("LinearDesign 输出为空。")
    return seq.replace("U", "T")


def write_fasta(path: str, seq_id: str, sequence: str, line_width: int = 80):
    """将序列写入 FASTA 文件。"""
    with open(path, "w", encoding="utf-8") as f:
        f.write(f">{seq_id}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + "\n")


def gc_content(seq: str) -> float:
    """计算 GC 含量百分比。"""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return (gc / len(seq)) * 100


def resolve_optional_aa(
    signal_peptide_preset: str,
    signal_peptide_aa: str,
    tm_domain_preset: str,
    tm_domain_aa: str,
):
    """解析信号肽与 TM 来源（预设优先，其次手工 AA）。"""
    sig_src = "none"
    tm_src = "none"
    sig = ""
    tm = ""

    if signal_peptide_preset:
        key = signal_peptide_preset.strip().lower()
        if key not in SIGNAL_PEPTIDE_PRESETS:
            raise ValueError(f"未知 signal_peptide_preset: {signal_peptide_preset}")
        sig = SIGNAL_PEPTIDE_PRESETS[key]
        sig_src = f"preset:{key}"
    elif signal_peptide_aa.strip():
        sig = signal_peptide_aa.strip().upper()
        sig_src = "manual"

    if tm_domain_preset:
        key = tm_domain_preset.strip().lower()
        if key not in TM_DOMAIN_PRESETS:
            raise ValueError(f"未知 tm_domain_preset: {tm_domain_preset}")
        tm = TM_DOMAIN_PRESETS[key]
        tm_src = f"preset:{key}"
    elif tm_domain_aa.strip():
        tm = tm_domain_aa.strip().upper()
        tm_src = "manual"

    return sig, sig_src, tm, tm_src


def main(
    run_id: str,
    linker: str,
    poly_a_len: int,
    codon_mode: str,
    lineardesign_bin: str,
    signal_peptide_preset: str,
    signal_peptide_aa: str,
    tm_domain_preset: str,
    tm_domain_aa: str,
):
    """
    主流程：
    1) 读取 selected_peptides.csv；
    2) 按 linker 串联肽段得到多价 ORF（氨基酸）；
    3) 转换为 DNA 编码序列；
    4) 拼接 5'UTR + ORF + 3'UTR + polyA；
    5) 输出 FASTA 和设计说明 JSON。
    """
    in_file = os.path.join("results", run_id, "selected_peptides.csv")
    out_dir = os.path.join("results", run_id)
    fasta_file = os.path.join(out_dir, "mrna_vaccine.fasta")
    design_file = os.path.join(out_dir, "mrna_design.json")

    if not os.path.exists(in_file):
        raise FileNotFoundError(f"未找到输入文件: {in_file}")

    df = pd.read_csv(in_file)
    if df.empty or "mut_peptide" not in df.columns:
        raise ValueError("selected_peptides.csv 缺少有效的 mut_peptide 数据。")

    # 按文件顺序取候选肽，去除空值
    peptides = [
        str(x).strip().upper()
        for x in df["mut_peptide"].tolist()
        if str(x).strip()
    ]
    if not peptides:
        raise ValueError("没有可用于构建 ORF 的肽段。")

    # 构建多价核心 ORF：Pep1 + linker + Pep2 + ...
    multivalent_core_aa = linker.join(peptides)
    signal_aa, signal_src, tm_aa, tm_src = resolve_optional_aa(
        signal_peptide_preset, signal_peptide_aa, tm_domain_preset, tm_domain_aa
    )
    # 最终 ORF：N端信号肽 + 多价核心 + C端可选TM
    orf_aa = f"{signal_aa}{multivalent_core_aa}{tm_aa}"
    if codon_mode == "basic":
        orf_dna = aa_to_dna(orf_aa)
    elif codon_mode == "optimized":
        orf_dna = aa_to_dna_optimized(orf_aa)
    elif codon_mode == "lineardesign":
        orf_dna = run_lineardesign_if_available(lineardesign_bin, orf_aa)
    else:
        raise ValueError(f"不支持的 codon_mode: {codon_mode}")

    # MVP 版本 UTR（后续可替换为更成熟的人源优化 UTR）
    utr5 = "GCCACC"
    utr3 = "CCTG"
    poly_a = "A" * poly_a_len

    # 这里输出 DNA 字母形式，作为 mRNA 设计占位结果（流程可跑通）
    full_seq = utr5 + orf_dna + utr3 + poly_a

    write_fasta(
        path=fasta_file,
        seq_id=f"{run_id}_multivalent_mrna_design",
        sequence=full_seq,
    )

    design = {
        "run_id": run_id,
        "design_stage": "MVP",
        "note": "当前版本支持 basic/optimized/lineardesign 三种密码子模式。",
        "input_file": in_file,
        "selected_peptide_count": len(peptides),
        "selected_peptides": peptides,
        "linker": linker,
        "signal_peptide_aa": signal_aa,
        "signal_peptide_source": signal_src,
        "tm_domain_aa": tm_aa,
        "tm_domain_source": tm_src,
        "multivalent_core_aa": multivalent_core_aa,
        "codon_mode": codon_mode,
        "segments": {
            "utr5": {"sequence": utr5, "length": len(utr5)},
            "orf_aa": {"sequence": orf_aa, "length": len(orf_aa)},
            "orf_dna": {"length": len(orf_dna)},
            "utr3": {"sequence": utr3, "length": len(utr3)},
            "poly_a": {"length": len(poly_a)},
        },
        "quality_metrics": {
            "full_length": len(full_seq),
            "gc_percent": round(gc_content(full_seq), 2),
        },
        "modification_strategy": "N1-methylpseudouridine（文档说明，未在序列字母中直接编码）",
        "output_files": {
            "fasta": fasta_file,
            "design_json": design_file,
        },
    }

    with open(design_file, "w", encoding="utf-8") as f:
        json.dump(design, f, ensure_ascii=False, indent=2)

    print(f"完成: {fasta_file}")
    print(f"完成: {design_file}")
    print(f"肽段数: {len(peptides)}")
    print(f"总长度: {len(full_seq)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--linker", default="AAY", help="肽段连接子，默认 AAY")
    parser.add_argument("--poly_a_len", type=int, default=120, help="polyA 长度，默认 120")
    parser.add_argument(
        "--codon_mode",
        default="optimized",
        choices=["basic", "optimized", "lineardesign"],
        help="密码子模式：basic/optimized/lineardesign，默认 optimized",
    )
    parser.add_argument(
        "--lineardesign_bin",
        default="",
        help="LinearDesign 可执行文件路径（codon_mode=lineardesign 时需要）",
    )
    parser.add_argument(
        "--signal_peptide_preset",
        default="",
        choices=["", "igkv", "mitd"],
        help="N端信号肽预设（可选）：igkv/mitd",
    )
    parser.add_argument(
        "--signal_peptide_aa",
        default="",
        help="N端信号肽手工氨基酸序列（与 preset 二选一）",
    )
    parser.add_argument(
        "--tm_domain_preset",
        default="",
        choices=["", "cd8a_tm"],
        help="C端TM预设（可选）：cd8a_tm",
    )
    parser.add_argument(
        "--tm_domain_aa",
        default="",
        help="C端TM手工氨基酸序列（与 preset 二选一）",
    )
    args = parser.parse_args()
    main(
        args.run_id,
        args.linker,
        args.poly_a_len,
        args.codon_mode,
        args.lineardesign_bin,
        args.signal_peptide_preset,
        args.signal_peptide_aa,
        args.tm_domain_preset,
        args.tm_domain_aa,
    )
