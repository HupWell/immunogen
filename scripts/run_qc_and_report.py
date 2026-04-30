# -*- coding: utf-8 -*-
"""
功能：执行基础质控并生成报告文件。
输入：results/<run_id>/mrna_vaccine.fasta 等中间结果
输出：results/<run_id>/figures/*, results/<run_id>/REPORT.md
"""
import os
import json
import argparse
import subprocess
import tempfile
from datetime import datetime
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def read_fasta_sequence(path: str) -> str:
    """读取 FASTA 文件中的主序列（忽略标题行）。"""
    seq_lines = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()


def gc_content(seq: str) -> float:
    """计算序列 GC 含量百分比。"""
    if not seq:
        return 0.0
    gc = sum(1 for b in seq if b in ("G", "C"))
    return (gc / len(seq)) * 100


def build_qc_json(run_id: str, out_dir: str, sequence: str, design: dict) -> str:
    """生成基础质控 JSON，便于后续流程读取。"""
    stability = design.get("mrna_stability", {})
    qc = {
        "run_id": run_id,
        "generated_at": datetime.now().isoformat(timespec="seconds"),
        "sequence_length": len(sequence),
        "gc_percent": round(gc_content(sequence), 2),
        "segment_lengths": {
            "utr5": design.get("segments", {}).get("utr5", {}).get("length"),
            "orf_aa": design.get("segments", {}).get("orf_aa", {}).get("length"),
            "orf_dna": design.get("segments", {}).get("orf_dna", {}).get("length"),
            "utr3": design.get("segments", {}).get("utr3", {}).get("length"),
            "poly_a": design.get("segments", {}).get("poly_a", {}).get("length"),
        },
        "rnafold_status": design.get("rnafold_status", "not_run"),
        "rnafold_note": design.get("rnafold_note", ""),
        "rnafold_mfe": design.get("rnafold_mfe"),
        "mrna_stability_status": stability.get("status", "not_run"),
        "mrna_stability_tools": stability.get("tools", []),
        "mrna_stability_metrics": stability.get("metrics", {}),
        "mrna_stability_files": stability.get("files", {}),
    }
    out_file = os.path.join(out_dir, "qc_metrics.json")
    with open(out_file, "w", encoding="utf-8") as f:
        json.dump(qc, f, ensure_ascii=False, indent=2)
    return out_file


def draw_binding_affinity_heatmap(out_dir: str, ranking_df: pd.DataFrame):
    """绘制真实亲和力热图。"""
    fig_dir = os.path.join(out_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    heatmap_file = os.path.join(fig_dir, "binding_affinity_heatmap.png")
    data = ranking_df.copy()
    data["affinity_nM"] = pd.to_numeric(data.get("affinity_nM"), errors="coerce")
    if data["affinity_nM"].isna().all():
        data["affinity_nM"] = pd.to_numeric(data.get("rank_score"), errors="coerce")
    pivot = data.pivot_table(
        index="mut_peptide",
        columns="hla_allele",
        values="affinity_nM",
        aggfunc="min"
    )
    plt.figure(figsize=(10, max(3, 0.7 * max(1, len(pivot.index)))))
    sns.heatmap(pivot, cmap="viridis", linewidths=0.5, annot=True, fmt=".2g")
    plt.title("Peptide-HLA Binding Affinity Heatmap")
    plt.xlabel("HLA Allele")
    plt.ylabel("Mutant Peptide")
    plt.tight_layout()
    plt.savefig(heatmap_file, dpi=180)
    plt.close()
    return heatmap_file


def run_rnafold(sequence: str):
    """优先调用 RNAfold 命令；不可用时回退到 ViennaRNA Python 绑定。"""
    cmd = ["RNAfold", "--noPS"]
    try:
        result = subprocess.run(cmd, input=sequence + "\n", capture_output=True, text=True, check=False)
    except FileNotFoundError:
        result = None
    if result is not None and result.returncode == 0:
        lines = [x.strip() for x in result.stdout.splitlines() if x.strip()]
        if len(lines) >= 2:
            structure_line = lines[1]
            parts = structure_line.rsplit("(", 1)
            structure = parts[0].strip()
            mfe = None
            if len(parts) == 2 and ")" in parts[1]:
                try:
                    mfe = float(parts[1].replace(")", "").strip())
                except ValueError:
                    mfe = None
            return structure, mfe
    # 回退：使用 Python 绑定 RNA.fold_compound
    try:
        import RNA  # type: ignore
        fc = RNA.fold_compound(sequence)
        structure, mfe = fc.mfe()
        return structure, float(mfe)
    except Exception:
        return None, None


def _run_version(cmd):
    """读取外部工具版本，便于下游复检。"""
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except FileNotFoundError:
        return ""
    output = (result.stdout or result.stderr or "").strip()
    return output.splitlines()[0] if output else ""


def _parse_rnaeval_energy(stdout: str):
    for line in stdout.splitlines():
        if "(" not in line or ")" not in line:
            continue
        try:
            return float(line.rsplit("(", 1)[1].split(")", 1)[0].strip())
        except (ValueError, IndexError):
            continue
    return None


def _summarize_lunp(path: str):
    """汇总 RNAplfold 非配对概率，作为局部可及性真实指标。"""
    values_l1 = []
    values_l10 = []
    values_all = []
    if not os.path.exists(path):
        return {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2 or not parts[0].isdigit():
                continue
            nums = []
            for item in parts[1:]:
                if item == "NA":
                    continue
                try:
                    nums.append(float(item))
                except ValueError:
                    continue
            if nums:
                values_l1.append(nums[0])
                values_all.extend(nums)
                if len(nums) >= 10:
                    values_l10.append(nums[9])

    def avg(values):
        return round(sum(values) / len(values), 6) if values else None

    return {
        "mean_unpaired_l1": avg(values_l1),
        "mean_unpaired_l10": avg(values_l10),
        "mean_unpaired_all_windows": avg(values_all),
        "positions_with_l1": len(values_l1),
    }


def run_mrna_stability_tools(out_dir: str, sequence: str, structure: str, mfe):
    """
    运行真实 mRNA 稳定性工具链。

    当前采用 ViennaRNA 增强口径：RNAfold 给出全局 MFE，RNAeval 复核该结构能量，
    RNAplfold 输出局部非配对概率，作为局部可及性/结构稳定性指标。
    """
    stability_dir = os.path.join(out_dir, "mrna_stability")
    os.makedirs(stability_dir, exist_ok=True)
    input_fasta = os.path.join(stability_dir, "input_mrna.fasta")
    input_rna = os.path.join(stability_dir, "input_mrna.rna")
    rnafold_stdout = os.path.join(stability_dir, "rnafold.stdout.txt")
    rnaeval_stdout = os.path.join(stability_dir, "rnaeval.stdout.txt")
    rnaplfold_stdout = os.path.join(stability_dir, "rnaplfold.stdout.txt")
    rnaplfold_stderr = os.path.join(stability_dir, "rnaplfold.stderr.txt")

    rna_seq = sequence.replace("T", "U")
    with open(input_fasta, "w", encoding="utf-8") as f:
        f.write(">mrna_design\n")
        for i in range(0, len(rna_seq), 80):
            f.write(rna_seq[i:i + 80] + "\n")
    with open(input_rna, "w", encoding="utf-8") as f:
        f.write(rna_seq + "\n")

    tools = [
        {"name": "RNAfold", "version": _run_version(["RNAfold", "--version"])},
        {"name": "RNAeval", "version": _run_version(["RNAeval", "--version"])},
        {"name": "RNAplfold", "version": _run_version(["RNAplfold", "--version"])},
    ]
    commands = {}
    metrics = {
        "rnafold_mfe": mfe,
        "rnaeval_energy": None,
    }

    rnafold_cmd = ["RNAfold", "--noPS"]
    commands["rnafold"] = rnafold_cmd
    try:
        rnafold_result = subprocess.run(
            rnafold_cmd, input=rna_seq + "\n", capture_output=True, text=True, check=False
        )
        with open(rnafold_stdout, "w", encoding="utf-8") as f:
            f.write(rnafold_result.stdout or "")
        if rnafold_result.returncode != 0:
            raise RuntimeError(rnafold_result.stderr.strip())
    except FileNotFoundError as exc:
        raise RuntimeError("未找到 RNAfold，无法完成真实 mRNA 稳定性评估。") from exc

    if structure:
        rnaeval_cmd = ["RNAeval"]
        commands["rnaeval"] = rnaeval_cmd
        try:
            rnaeval_result = subprocess.run(
                rnaeval_cmd,
                input=f"{rna_seq}\n{structure}\n",
                capture_output=True,
                text=True,
                check=False,
            )
            with open(rnaeval_stdout, "w", encoding="utf-8") as f:
                f.write(rnaeval_result.stdout or "")
            if rnaeval_result.returncode != 0:
                raise RuntimeError(rnaeval_result.stderr.strip())
            metrics["rnaeval_energy"] = _parse_rnaeval_energy(rnaeval_result.stdout or "")
        except FileNotFoundError as exc:
            raise RuntimeError("未找到 RNAeval，无法完成真实 mRNA 稳定性复核。") from exc

    with tempfile.TemporaryDirectory() as tmpdir:
        window = min(240, max(80, len(rna_seq)))
        span = min(160, window)
        rnaplfold_cmd = ["RNAplfold", "-u", "30", "-W", str(window), "-L", str(span)]
        commands["rnaplfold"] = rnaplfold_cmd
        try:
            rnaplfold_result = subprocess.run(
                rnaplfold_cmd,
                input=rna_seq + "\n",
                capture_output=True,
                text=True,
                check=False,
                cwd=tmpdir,
            )
        except FileNotFoundError as exc:
            raise RuntimeError("未找到 RNAplfold，无法完成真实局部可及性评估。") from exc
        with open(rnaplfold_stdout, "w", encoding="utf-8") as f:
            f.write(rnaplfold_result.stdout or "")
        with open(rnaplfold_stderr, "w", encoding="utf-8") as f:
            f.write(rnaplfold_result.stderr or "")
        if rnaplfold_result.returncode != 0:
            raise RuntimeError(f"RNAplfold 执行失败: {rnaplfold_stderr}")
        lunp_src = os.path.join(tmpdir, "plfold_lunp")
        lunp_out = os.path.join(stability_dir, "rnaplfold_lunp.tsv")
        if not os.path.exists(lunp_src):
            raise RuntimeError("RNAplfold 未生成 plfold_lunp 输出。")
        with open(lunp_src, "r", encoding="utf-8") as src, open(lunp_out, "w", encoding="utf-8") as dst:
            dst.write(src.read())
        metrics.update(_summarize_lunp(lunp_out))

    files = {
        "input_fasta": input_fasta,
        "input_rna": input_rna,
        "rnafold_stdout": rnafold_stdout,
        "rnaeval_stdout": rnaeval_stdout,
        "rnaplfold_stdout": rnaplfold_stdout,
        "rnaplfold_stderr": rnaplfold_stderr,
        "rnaplfold_lunp": os.path.join(stability_dir, "rnaplfold_lunp.tsv"),
    }
    return {
        "status": "ok",
        "note": "ViennaRNA 真实稳定性工具链计算成功：RNAfold + RNAeval + RNAplfold。",
        "tools": tools,
        "commands": commands,
        "metrics": metrics,
        "files": files,
    }


def draw_secondary_structure_figure(out_dir: str, sequence: str, structure: str, mfe):
    """绘制二级结构轮廓图（真实 RNAfold 或回退 profile）。"""
    fig_dir = os.path.join(out_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    structure_file = os.path.join(fig_dir, "mrna_secondary_structure.png")
    if not structure:
        structure = "." * len(sequence)
    y = []
    for c in structure:
        if c == "(":
            y.append(1)
        elif c == ")":
            y.append(-1)
        else:
            y.append(0)
    plt.figure(figsize=(12, 3))
    plt.plot(y, linewidth=1.0)
    if mfe is None:
        title = "mRNA Secondary Structure (dot-bracket profile, MFE=NA)"
    else:
        title = f"mRNA Secondary Structure (dot-bracket profile, MFE={mfe:.2f} kcal/mol)"
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Pairing state")
    plt.tight_layout()
    plt.savefig(structure_file, dpi=180)
    plt.close()
    return structure_file


def write_report(
    run_id: str,
    out_dir: str,
    ranking_df: pd.DataFrame,
    selected_df: pd.DataFrame,
    design: dict,
    qc_file: str,
    heatmap_file: str,
    structure_file: str
):
    """生成 REPORT.md，汇总本轮运行输入、过程和结果。"""
    report_file = os.path.join(out_dir, "REPORT.md")
    full_length = design.get("quality_metrics", {}).get("full_length", "NA")
    gc_percent = design.get("quality_metrics", {}).get("gc_percent", "NA")
    selected_count = design.get("selected_peptide_count", len(selected_df))
    stability = design.get("mrna_stability", {})
    stability_metrics = stability.get("metrics", {})
    stability_tools = ", ".join(
        f"{x.get('name')} {x.get('version', '')}".strip()
        for x in stability.get("tools", [])
    ) or "NA"
    stability_files = stability.get("files", {})
    # 避免依赖 tabulate，使用 CSV 文本块展示前 10 条
    top_preview = selected_df.head(10).to_csv(index=False)

    content = f"""# ImmunoGen 报告（{run_id}）

## 1. 运行概况
- 运行编号：`{run_id}`
- 生成时间：`{datetime.now().isoformat(timespec='seconds')}`
- 排名结果条目数：`{len(ranking_df)}`
- 入选肽段数量：`{selected_count}`

## 2. mRNA 设计摘要
- 总长度：`{full_length}`
- GC 含量：`{gc_percent}%`
- linker：`{design.get("linker", "NA")}`
- signal_peptide：`{design.get("signal_peptide_source", "none")}` / `len={len(str(design.get("signal_peptide_aa", "")))}`
- tm_domain：`{design.get("tm_domain_source", "none")}` / `len={len(str(design.get("tm_domain_aa", "")))}`
- 修饰策略说明：`{design.get("modification_strategy", "NA")}`

## 3. Top 候选肽（最多展示前 10 条）
```csv
{top_preview}
```

## 4. 文件产物清单
- `peptide_mhc_ranking.csv`
- `selected_peptides.csv`
- `mrna_vaccine.fasta`
- `mrna_design.json`
- `qc_metrics.json`
- `POSITIVE_CONTROL.md`、`SELF_CHECK.md`：完整执行 `run_all.py` 时由 `prepare_self_certification.py` 写入；若仅单独运行 `run_qc_and_report.py` 则可能尚未生成。
- `figures/{os.path.basename(heatmap_file)}`
- `figures/{os.path.basename(structure_file)}`

## 5. 说明与后续
- 当前为 MVP 版本，已经可复现跑通“预测 -> 筛选 -> 构建 -> 报告”主流程。
- `affinity_nM` 若为空，通常由示例数据或等位基因覆盖范围导致，建议换真实病例数据复核。
- 当前已输出真实热图；RNAfold 若可用会写入 MFE，否则使用回退 profile 图。
- 信号肽/TM 可选项为工程接口，具体序列应由实验方案确定（本报告仅记录参数与长度）。
- 基础质控 JSON：`{os.path.basename(qc_file)}`

## 5.1 mRNA 稳定性真实工具链
- 状态：`{stability.get("status", "not_run")}`
- 工具：`{stability_tools}`
- RNAfold MFE：`{stability_metrics.get("rnafold_mfe", "NA")}`
- RNAeval 复核能量：`{stability_metrics.get("rnaeval_energy", "NA")}`
- RNAplfold 平均单碱基非配对概率：`{stability_metrics.get("mean_unpaired_l1", "NA")}`
- RNAplfold 平均 10nt 非配对概率：`{stability_metrics.get("mean_unpaired_l10", "NA")}`
- 稳定性输入：`{stability_files.get("input_fasta", "NA")}`
- RNAplfold 输出：`{stability_files.get("rnaplfold_lunp", "NA")}`

## 6. 信号肽与表达定位（简述）
- 如启用 N 端信号肽，常见目的是促进分泌/加工路径，提高抗原递呈机会；最终方案需结合实验体系验证。
- 如启用 C 端 TM，常用于增强膜定位或改变加工路径；是否采用以实验设计与安全性评估为准。
- 本流程只记录“是否启用 + 使用何种序列来源（preset/manual）”，不对免疫学效果做结论性承诺。
"""
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(content)
    return report_file


def main(run_id: str):
    """
    主流程：
    1) 读取 ranking/selected/mRNA/design 文件；
    2) 计算基础 QC 指标并写入 qc_metrics.json；
    3) 生成图文件占位；
    4) 生成 REPORT.md。
    """
    out_dir = os.path.join("results", run_id)
    ranking_file = os.path.join(out_dir, "peptide_mhc_ranking.csv")
    selected_file = os.path.join(out_dir, "selected_peptides.csv")
    fasta_file = os.path.join(out_dir, "mrna_vaccine.fasta")
    design_file = os.path.join(out_dir, "mrna_design.json")

    for f in [ranking_file, selected_file, fasta_file, design_file]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"未找到必需输入文件: {f}")

    ranking_df = pd.read_csv(ranking_file)
    selected_df = pd.read_csv(selected_file)
    with open(design_file, "r", encoding="utf-8") as f:
        design = json.load(f)
    sequence = read_fasta_sequence(fasta_file)

    structure, mfe = run_rnafold(sequence)
    if structure is None:
        design["rnafold_status"] = "fallback_profile"
        design["rnafold_note"] = "未检测到 RNAfold，使用回退结构轮廓图。"
        design["rnafold_mfe"] = None
    else:
        design["rnafold_status"] = "ok"
        design["rnafold_note"] = "RNAfold 计算成功。"
        design["rnafold_mfe"] = mfe
        design["mrna_stability"] = run_mrna_stability_tools(out_dir, sequence, structure, mfe)
    with open(design_file, "w", encoding="utf-8") as f:
        json.dump(design, f, ensure_ascii=False, indent=2)

    heatmap_file = draw_binding_affinity_heatmap(out_dir, ranking_df)
    structure_file = draw_secondary_structure_figure(out_dir, sequence, structure, mfe)

    qc_file = build_qc_json(run_id, out_dir, sequence, design)
    report_file = write_report(
        run_id, out_dir, ranking_df, selected_df, design, qc_file, heatmap_file, structure_file
    )

    print(f"完成: {qc_file}")
    print(f"完成: {report_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    args = parser.parse_args()
    main(args.run_id)
