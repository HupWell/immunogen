# -*- coding: utf-8 -*-
"""
功能：将入选肽段串联并生成多价 mRNA 序列。
输入：results/<run_id>/selected_peptides.csv
输出：results/<run_id>/mrna_vaccine.fasta, results/<run_id>/mrna_design.json
"""
import os
import json
import argparse
import re
import shlex
import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional
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

DNA_TRANSLATION_TABLE = {
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    "AAT": "N", "AAC": "N",
    "GAT": "D", "GAC": "D",
    "TGT": "C", "TGC": "C",
    "CAA": "Q", "CAG": "Q",
    "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    "CAT": "H", "CAC": "H",
    "ATT": "I", "ATC": "I", "ATA": "I",
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "AAA": "K", "AAG": "K",
    "ATG": "M",
    "TTT": "F", "TTC": "F",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "TGG": "W",
    "TAT": "Y", "TAC": "Y",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TAA": "*", "TAG": "*", "TGA": "*",
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


def _read_fasta_sequence(path: Path) -> str:
    """读取 FASTA 中的首条序列。"""
    seq = []
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    return "".join(seq).strip().upper()


def _normalize_coding_sequence(seq: str) -> str:
    """把真实工具输出规整为 DNA 字母序列，并拒绝非编码字符。"""
    normalized = re.sub(r"\s+", "", seq).upper().replace("U", "T")
    if not normalized:
        raise RuntimeError("真实密码子优化工具没有输出序列。")
    invalid = sorted(set(normalized) - set("ACGT"))
    if invalid:
        raise RuntimeError(f"真实密码子优化工具输出包含非 DNA/RNA 字符: {''.join(invalid)}")
    if len(normalized) % 3 != 0:
        raise RuntimeError(f"真实密码子优化工具输出长度不是 3 的倍数: {len(normalized)}")
    return normalized


def _translate_dna(seq: str) -> str:
    """翻译 DNA 编码序列，用于确认真实工具没有改变蛋白 ORF。"""
    aa = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        if codon not in DNA_TRANSLATION_TABLE:
            raise RuntimeError(f"无法翻译密码子: {codon}")
        aa.append(DNA_TRANSLATION_TABLE[codon])
    return "".join(aa).rstrip("*")


def _validate_optimized_orf(orf_dna: str, expected_aa: str) -> None:
    """验收优化后 ORF 与输入氨基酸序列一致。"""
    translated = _translate_dna(orf_dna)
    if translated != expected_aa:
        raise RuntimeError(
            "真实密码子优化工具输出翻译结果与输入 ORF 不一致: "
            f"expected_len={len(expected_aa)}, translated_len={len(translated)}"
        )


def _parse_lineardesign_stdout(stdout: str) -> str:
    """解析 LinearDesign 标准输出中的 mRNA sequence 行。"""
    for line in stdout.splitlines():
        if line.lower().startswith("mrna sequence:"):
            return line.split(":", 1)[1].strip()
    lines = [x.strip() for x in stdout.splitlines() if x.strip()]
    return lines[-1] if lines else ""


def _format_cmd(template: str, values: Dict[str, str]) -> List[str]:
    """展开 real_cmd 模板，保留命令参数边界。"""
    if not template.strip():
        raise ValueError("codon_mode=real_cmd 时必须提供 --codon_real_cmd。")
    return shlex.split(template.format(**values))


def _run_version_cmd(version_cmd: str) -> str:
    """记录真实工具版本；版本命令缺失时不影响主流程。"""
    if not version_cmd.strip():
        return ""
    result = subprocess.run(
        shlex.split(version_cmd),
        capture_output=True,
        text=True,
        check=False,
    )
    output = (result.stdout or result.stderr or "").strip()
    return output.splitlines()[0] if output else ""


def run_real_codon_optimizer(
    run_id: str,
    aa_seq: str,
    codon_real_tool: str,
    codon_real_cmd: str,
    codon_real_version_cmd: str,
    out_dir: str,
) -> Dict[str, object]:
    """
    调用真实命令行密码子优化工具；失败时直接报错，不回退到内部简化优化。

    LinearDesign 可通过 stdin 输出 mRNA sequence；VaxPress 等工具可通过
    {input_fasta}/{output_dir}/{output_fasta} 占位符写出 FASTA。
    """
    work_dir = Path(out_dir) / "codon_optimization"
    work_dir.mkdir(parents=True, exist_ok=True)
    input_fasta = work_dir / "input_protein.fasta"
    output_fasta = work_dir / "optimized_orf.fasta"
    stdout_log = work_dir / "real_cmd.stdout.txt"
    stderr_log = work_dir / "real_cmd.stderr.txt"
    tool_output_dir = work_dir / "tool_output"
    for stale_file in (output_fasta, stdout_log, stderr_log):
        if stale_file.exists():
            stale_file.unlink()
    if tool_output_dir.exists():
        shutil.rmtree(tool_output_dir)
    tool_output_dir.mkdir(parents=True, exist_ok=True)

    write_fasta(str(input_fasta), f"{run_id}_orf_protein", aa_seq)
    values = {
        "run_id": run_id,
        "aa_seq": aa_seq,
        "input_fasta": str(input_fasta),
        "output_fasta": str(output_fasta),
        "output_dir": str(tool_output_dir),
    }
    cmd = _format_cmd(codon_real_cmd, values)
    result = subprocess.run(
        cmd,
        input=f"{aa_seq}\n",
        capture_output=True,
        text=True,
        check=False,
    )
    stdout_log.write_text(result.stdout or "", encoding="utf-8")
    stderr_log.write_text(result.stderr or "", encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            "真实密码子优化命令执行失败: "
            f"returncode={result.returncode}, stderr={stderr_log}"
        )

    output_candidates = [
        output_fasta,
        tool_output_dir / "best-sequence.fasta",
        tool_output_dir / "best_sequence.fasta",
    ]
    raw_seq: Optional[str] = None
    selected_output = ""
    for candidate in output_candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            raw_seq = _read_fasta_sequence(candidate)
            selected_output = str(candidate)
            break
    if raw_seq is None:
        raw_seq = _parse_lineardesign_stdout(result.stdout or "")
        selected_output = str(stdout_log)

    orf_dna = _normalize_coding_sequence(raw_seq)
    _validate_optimized_orf(orf_dna, aa_seq)
    if not output_fasta.exists():
        write_fasta(str(output_fasta), f"{run_id}_real_cmd_optimized_orf", orf_dna)

    return {
        "backend": "real_cmd",
        "tool": codon_real_tool.strip() or "external_codon_optimizer",
        "tool_version": _run_version_cmd(codon_real_version_cmd),
        "command": cmd,
        "returncode": result.returncode,
        "input_fasta": str(input_fasta),
        "selected_output": selected_output,
        "normalized_output_fasta": str(output_fasta),
        "stdout_log": str(stdout_log),
        "stderr_log": str(stderr_log),
        "output_dir": str(tool_output_dir),
    }


def run_lineardesign_if_available(lineardesign_bin: str, aa_seq: str) -> str:
    """
    如果提供了 LinearDesign 可执行文件，则调用其输出优化序列。
    返回 DNA/RNA 字母序列；失败时抛异常。
    """
    if not lineardesign_bin:
        raise ValueError("未提供 LinearDesign 路径。")
    if not os.path.exists(lineardesign_bin):
        raise FileNotFoundError(f"LinearDesign 可执行文件不存在: {lineardesign_bin}")

    cmd = [lineardesign_bin]
    result = subprocess.run(cmd, input=f"{aa_seq}\n", capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"LinearDesign 执行失败: {result.stderr.strip()}")
    return _normalize_coding_sequence(_parse_lineardesign_stdout(result.stdout))


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
    codon_real_tool: str,
    codon_real_cmd: str,
    codon_real_version_cmd: str,
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
    codon_optimizer = {
        "backend": codon_mode,
        "tool": "internal",
        "tool_version": "",
        "command": [],
        "input_fasta": "",
        "selected_output": "",
        "normalized_output_fasta": "",
        "stdout_log": "",
        "stderr_log": "",
        "output_dir": "",
    }
    if codon_mode == "basic":
        orf_dna = aa_to_dna(orf_aa)
    elif codon_mode == "optimized":
        orf_dna = aa_to_dna_optimized(orf_aa)
    elif codon_mode == "lineardesign":
        orf_dna = run_lineardesign_if_available(lineardesign_bin, orf_aa)
        _validate_optimized_orf(orf_dna, orf_aa)
        codon_optimizer.update({
            "backend": "lineardesign_legacy",
            "tool": "LinearDesign",
            "command": [lineardesign_bin],
        })
    elif codon_mode == "real_cmd":
        real_result = run_real_codon_optimizer(
            run_id=run_id,
            aa_seq=orf_aa,
            codon_real_tool=codon_real_tool,
            codon_real_cmd=codon_real_cmd,
            codon_real_version_cmd=codon_real_version_cmd,
            out_dir=out_dir,
        )
        orf_dna = _read_fasta_sequence(Path(str(real_result["normalized_output_fasta"])))
        codon_optimizer.update(real_result)
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
        "note": "当前版本支持 basic/optimized/lineardesign/real_cmd；real_cmd 为真实命令行工具接入且失败不回退。",
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
        "codon_optimizer": codon_optimizer,
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
        choices=["basic", "optimized", "lineardesign", "real_cmd"],
        help="密码子模式：basic/optimized/lineardesign/real_cmd，默认 optimized",
    )
    parser.add_argument(
        "--lineardesign_bin",
        default="",
        help="LinearDesign 可执行文件路径（codon_mode=lineardesign 时需要）",
    )
    parser.add_argument(
        "--codon_real_tool",
        default="",
        help="真实密码子优化工具名称，用于记录验收元数据",
    )
    parser.add_argument(
        "--codon_real_cmd",
        default="",
        help=(
            "真实密码子优化命令模板（codon_mode=real_cmd 时需要），"
            "可使用 {input_fasta}/{output_fasta}/{output_dir}/{aa_seq}/{run_id}"
        ),
    )
    parser.add_argument(
        "--codon_real_version_cmd",
        default="",
        help="真实密码子优化工具版本命令，用于记录验收元数据",
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
        args.codon_real_tool,
        args.codon_real_cmd,
        args.codon_real_version_cmd,
        args.signal_peptide_preset,
        args.signal_peptide_aa,
        args.tm_domain_preset,
        args.tm_domain_aa,
    )
