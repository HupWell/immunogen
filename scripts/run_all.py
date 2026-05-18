# -*- coding: utf-8 -*-
"""
功能：一键串行执行 ImmunoGen MVP 全流程。
执行顺序：
0)（可选，`--target transcriptome_prior`）bulk 基因优先级侧车，见 README §1.1A  
0b)（可选，`--target b_to_a_simhub`）filter_neoantigen_by_gene_rank.py → 路径 A 全流程  
1) validate_input.py
2) run_immunogenicity_adapters.py
3) predict_mhc_ranking.py
4) select_top_peptides.py
5) build_multivalent_mrna.py
6) run_qc_and_report.py
7) prepare_self_certification.py（POSITIVE_CONTROL.md / SELF_CHECK.md）
8)（可选，`--prepare_structure_inputs`）prepare_mhc_chain_sequences.py  
9)（可选，`--prepare_pandora`）run_pandora_structure.py --top_k 与 `--top_k_md` 对齐（可用 `--pandora_python` 指定含 csb-pandora 的解释器）  
10) prepare_simhub_delivery.py（支持多 `<case>_md_rXX` SimHub 子目录）  
11) validate_feasibility.py  
12)（目标含 full/report/simhub 时）check_simhub_evidence.py
"""
import os
import sys
import argparse
import subprocess
from typing import Dict, Optional


def run_step(command: list, env: Optional[Dict[str, str]] = None):
    """执行单个步骤并打印命令，失败时立即退出。"""
    print(f"\n[执行] {' '.join(command)}")
    result = subprocess.run(command, check=False, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"步骤执行失败，退出码: {result.returncode}")


def _pandora_step_env(pandora_exe: str) -> Optional[Dict[str, str]]:
    """
    为 PANDORA 子进程补齐 PATH 与 MUSCLE：tools/pandora_bin/muscle 依赖真实 muscle，
    且 shebang 的 python3 应优先落在含 csb-pandora 的 conda env/bin。
    """
    bin_dir = os.path.dirname(os.path.abspath(pandora_exe))
    if not bin_dir or not os.path.isdir(bin_dir):
        return None
    muscle = os.path.join(bin_dir, "muscle")
    if not os.path.isfile(muscle):
        return None
    env = os.environ.copy()
    env["PATH"] = bin_dir + os.pathsep + env.get("PATH", "")
    env["PANDORA_REAL_MUSCLE"] = muscle
    return env


def _preflight_mhc1_backend(run_id: str, backend_netmhcpan: str, backend_bigmhc: str):
    """
    在执行 predict_mhc_ranking 前做 MHC-I 交叉验证后端检查，尽早给出可读错误。
    """
    raw_dir = os.path.join("results", run_id, "tool_outputs", "raw")
    checks = [
        ("mhc1_netmhcpan", backend_netmhcpan, "MHC1_NETMHCPAN_CMD"),
        ("mhc1_bigmhc", backend_bigmhc, "MHC1_BIGMHC_CMD"),
    ]
    for tool, backend, env_key in checks:
        b = (backend or "").strip().lower()
        if b == "real_tsv":
            tsv_path = os.path.join(raw_dir, f"{tool}.tsv")
            if not os.path.exists(tsv_path):
                raise FileNotFoundError(
                    f"{tool} 后端设为 real_tsv，但未找到文件: {tsv_path}"
                )
        elif b == "real_cmd":
            if not os.environ.get(env_key, "").strip():
                raise RuntimeError(
                    f"{tool} 后端设为 real_cmd，但未设置环境变量 {env_key}"
                )


def _preflight_mhc2_backend(run_id: str, mhc2_backend: str):
    """
    在执行 predict_mhc_ranking 前做 MHC-II 后端检查（仅 real_tsv 需要本地文件）。
    """
    b = (mhc2_backend or "").strip().lower()
    if b != "real_tsv":
        return
    tsv_path = os.path.join("results", run_id, "tool_outputs", "raw", "mhc2_netmhciipan.tsv")
    if not os.path.exists(tsv_path):
        raise FileNotFoundError(
            f"mhc2_backend 设为 real_tsv，但未找到文件: {tsv_path}"
        )


def _preflight_immunogenicity_backend(run_id: str, backend: str, tool: str):
    """
    在执行 run_immunogenicity_adapters 前做免疫原性后端检查。
    """
    b = (backend or "").strip().lower()
    if b == "real_tsv":
        tsv_path = os.path.join("results", run_id, "tool_outputs", "raw", f"{tool}.tsv")
        if not os.path.exists(tsv_path):
            raise FileNotFoundError(
                f"{tool} 后端设为 real_tsv，但未找到文件: {tsv_path}"
            )
    elif b == "real_cmd":
        env_key = f"IMMUNO_{tool.upper()}_CMD"
        if not os.environ.get(env_key, "").strip():
            raise RuntimeError(
                f"{tool} 后端设为 real_cmd，但未设置环境变量 {env_key}"
            )


def main(
    run_id: str,
    top_n: int,
    min_dissimilarity: float,
    linker: str,
    poly_a_len: int,
    codon_mode: str,
    codon_real_tool: str,
    codon_real_cmd: str,
    codon_real_version_cmd: str,
    signal_peptide_preset: str,
    signal_peptide_aa: str,
    tm_domain_preset: str,
    tm_domain_aa: str,
    top_k_md: int,
    prepare_pandora: bool,
    prepare_structure_inputs: bool,
    structure_seq_refresh_remote: bool,
    structure_seq_strict: bool,
    structure_backend: str,
    structure_input_pdb: str,
    feasibility_top_n: int,
    mhc2_backend: str,
    wi_deepimmuno: float,
    wi_prime: float,
    wi_repitope: float,
    backend_mhc1_netmhcpan: str,
    backend_mhc1_bigmhc: str,
    require_real_mhc2: bool,
    require_real_mhc1_cv: bool,
    backend_deepimmuno: str,
    backend_prime: str,
    backend_repitope: str,
    require_real_immunogenicity_deepimmuno: bool,
    require_real_immunogenicity_prime: bool,
    require_real_immunogenicity_repitope: bool,
    allow_proxy_scores: bool,
    target: str,
    mrna_output_suffix: str,
    ensure_positive_control_peptides: str,
    pandora_python: str,
    bulk_tpm_path: str,
    bulk_case_columns: str,
    bulk_control_columns: str,
    bulk_gene_col_name: str,
    bulk_pseudo: float,
    bulk_out_subdir: str,
    cohort_id: str = "",
    b_to_a_top_n_genes: int = 500,
    b_to_a_rank_subdir: str = "transcriptome_prior",
    b_to_a_min_log2fc: Optional[float] = None,
    skip_b_to_a_filter: bool = False,
):
    """按固定顺序执行全流程。"""
    python_exe = sys.executable
    scripts_dir = os.path.dirname(os.path.abspath(__file__))
    # PANDORA 常安装在独立 conda 环境；其余步骤仍用当前解释器（如含 TensorFlow 的排名环境）。
    pandora_exe = (pandora_python or "").strip() or python_exe

    out_dir_bulk = os.path.join("results", run_id, (bulk_out_subdir or "transcriptome_prior").strip() or "transcriptome_prior")
    step0_bulk_rank = [
        python_exe,
        os.path.join(scripts_dir, "bulk_expression_gene_rank.py"),
        "--tpm_path",
        (bulk_tpm_path or "").strip(),
        "--out_dir",
        out_dir_bulk,
        "--case_columns",
        (bulk_case_columns or "").strip(),
        "--control_columns",
        (bulk_control_columns or "").strip(),
    ]
    gcn = (bulk_gene_col_name or "gene_name").strip() or "gene_name"
    if gcn != "gene_name":
        step0_bulk_rank.extend(["--gene_col_name", gcn])
    if float(bulk_pseudo) != 1.0:
        step0_bulk_rank.extend(["--pseudo", str(float(bulk_pseudo))])

    rank_sub = (b_to_a_rank_subdir or "transcriptome_prior").strip() or "transcriptome_prior"
    step0_gene_filter = [
        python_exe,
        os.path.join(scripts_dir, "filter_neoantigen_by_gene_rank.py"),
        "--run_id",
        run_id,
        "--top_n_genes",
        str(int(b_to_a_top_n_genes)),
        "--rank_subdir",
        rank_sub,
    ]
    if b_to_a_min_log2fc is not None:
        step0_gene_filter.extend(["--min_log2fc", str(float(b_to_a_min_log2fc))])

    step1 = [
        python_exe,
        os.path.join(scripts_dir, "validate_input.py"),
        "--run_id",
        run_id,
    ]
    step2 = [
        python_exe,
        os.path.join(scripts_dir, "run_immunogenicity_adapters.py"),
        "--run_id",
        run_id,
        "--backend_deepimmuno",
        backend_deepimmuno,
        "--backend_prime",
        backend_prime,
        "--backend_repitope",
        backend_repitope,
    ]
    step3 = [
        python_exe,
        os.path.join(scripts_dir, "predict_mhc_ranking.py"),
        "--run_id",
        run_id,
        "--mhc2_backend",
        mhc2_backend,
        "--wi_deepimmuno",
        str(wi_deepimmuno),
        "--wi_prime",
        str(wi_prime),
        "--wi_repitope",
        str(wi_repitope),
        "--backend_mhc1_netmhcpan",
        backend_mhc1_netmhcpan,
        "--backend_mhc1_bigmhc",
        backend_mhc1_bigmhc,
    ]
    if require_real_mhc2:
        step3.append("--require_real_mhc2")
    if require_real_mhc1_cv:
        step3.append("--require_real_mhc1_cv")
    if require_real_immunogenicity_deepimmuno:
        step3.append("--require_real_immunogenicity_deepimmuno")
    if require_real_immunogenicity_prime:
        step3.append("--require_real_immunogenicity_prime")
    if require_real_immunogenicity_repitope:
        step3.append("--require_real_immunogenicity_repitope")
    if allow_proxy_scores:
        step3.append("--allow_proxy_scores")
    step4 = [
        python_exe,
        os.path.join(scripts_dir, "select_top_peptides.py"),
        "--run_id",
        run_id,
        "--top_n",
        str(top_n),
        "--min_dissimilarity",
        str(min_dissimilarity),
    ]
    if (ensure_positive_control_peptides or "").strip():
        step4.extend(
            [
                "--ensure_positive_control_peptides",
                ensure_positive_control_peptides.strip(),
            ]
        )
    step5 = [
        python_exe,
        os.path.join(scripts_dir, "build_multivalent_mrna.py"),
        "--run_id",
        run_id,
        "--linker",
        linker,
        "--poly_a_len",
        str(poly_a_len),
        "--codon_mode",
        codon_mode,
        "--codon_real_tool",
        codon_real_tool,
        "--codon_real_cmd",
        codon_real_cmd,
        "--codon_real_version_cmd",
        codon_real_version_cmd,
        "--signal_peptide_preset",
        signal_peptide_preset,
        "--signal_peptide_aa",
        signal_peptide_aa,
        "--tm_domain_preset",
        tm_domain_preset,
        "--tm_domain_aa",
        tm_domain_aa,
    ]
    if (mrna_output_suffix or "").strip():
        step5.extend(["--mrna_output_suffix", mrna_output_suffix.strip()])
    step6 = [
        python_exe,
        os.path.join(scripts_dir, "run_qc_and_report.py"),
        "--run_id",
        run_id,
    ]
    step7 = [
        python_exe,
        os.path.join(scripts_dir, "prepare_self_certification.py"),
        "--run_id",
        run_id,
    ]
    step8 = [
        python_exe,
        os.path.join(scripts_dir, "prepare_simhub_delivery.py"),
        "--run_id",
        run_id,
        "--top_k",
        str(top_k_md),
        "--structure_backend",
        structure_backend,
        "--structure_input_pdb",
        structure_input_pdb,
    ]
    step8_1 = [
        python_exe,
        os.path.join(scripts_dir, "prepare_mhc_chain_sequences.py"),
        "--run_id",
        run_id,
    ]
    if structure_seq_refresh_remote:
        step8_1.append("--refresh_remote")
    if structure_seq_strict:
        step8_1.append("--strict")
    step8_pandora = [
        pandora_exe,
        os.path.join(scripts_dir, "run_pandora_structure.py"),
        "--run_id",
        run_id,
        "--top_k",
        str(top_k_md),
    ]
    step9 = [
        python_exe,
        os.path.join(scripts_dir, "validate_feasibility.py"),
        "--run_id",
        run_id,
        "--top_n",
        str(feasibility_top_n),
    ]
    step10 = [
        python_exe,
        os.path.join(scripts_dir, "check_simhub_evidence.py"),
        "--run_id",
        run_id,
    ]

    if target != "transcriptome_prior" and require_real_mhc1_cv and backend_mhc1_netmhcpan == "off" and backend_mhc1_bigmhc == "off":
        raise RuntimeError("已启用 --require_real_mhc1_cv，但两个 MHC-I 交叉验证后端均为 off。")

    path_a_full = [step1, step2, step3, step4, step5, step6, step8, step10, step7, step9]
    pipelines = {
        "full": path_a_full,
        "mhc_ranking": [step1, step2, step3],
        "report": [step6, step10, step7],
        "simhub": [step8, step10, step7],
        "feasibility": [step9],
        "transcriptome_prior": [step0_bulk_rank],
        "b_to_a_simhub": path_a_full,
    }
    selected = list(pipelines[target])

    if target == "b_to_a_simhub":
        rank_csv = os.path.join("results", run_id, rank_sub, "candidate_gene_ranking.csv")
        neo_csv = os.path.join("deliveries", run_id, "to_immunogen", "neoantigen_candidates.csv")
        if not skip_b_to_a_filter:
            if not os.path.isfile(rank_csv):
                raise FileNotFoundError(
                    f"target=b_to_a_simhub 需要路径 B 排名: {rank_csv}；"
                    "可先 --target transcriptome_prior 或提供 --cohort_id。"
                )
            if not os.path.isfile(neo_csv):
                raise FileNotFoundError(f"target=b_to_a_simhub 需要: {neo_csv}")
            selected = [step0_gene_filter] + selected
        print(f"开始执行 B→A→SimHub（run_id={run_id}）...", flush=True)

    if target == "transcriptome_prior":
        cid = (cohort_id or "").strip()
        if cid:
            # 一队列一个 run_id：由 run_transcriptome_cohort 写契约并排名
            cohort_script = os.path.join(scripts_dir, "run_transcriptome_cohort.py")
            cohort_cmd = [
                python_exe,
                cohort_script,
                "--cohort_id",
                cid,
                "--run_id",
                run_id,
                "--out_subdir",
                (bulk_out_subdir or "transcriptome_prior").strip() or "transcriptome_prior",
                "--pseudo",
                str(float(bulk_pseudo)),
            ]
            print(f"开始执行 ImmunoGen 流程（target={target}, cohort_id={cid}）...", flush=True)
            run_step(cohort_cmd)
            print("\n所选步骤执行完成。")
            return
        if not (bulk_tpm_path or "").strip():
            raise RuntimeError(
                "target=transcriptome_prior 时必须提供 --bulk_tpm_path，或改用 --cohort_id。"
            )
        if not (bulk_case_columns or "").strip() or not (bulk_control_columns or "").strip():
            raise RuntimeError(
                "target=transcriptome_prior 时必须同时提供 --bulk_case_columns 与 "
                "--bulk_control_columns（未使用 --cohort_id 时）。"
            )

    print(f"开始执行 ImmunoGen 流程（target={target}）...", flush=True)
    if target == "transcriptome_prior":
        for step in selected:
            run_step(step)
        print("\n所选步骤执行完成。")
        return

    if step2 in selected:
        _preflight_immunogenicity_backend(run_id, backend_deepimmuno, "deepimmuno")
        _preflight_immunogenicity_backend(run_id, backend_prime, "prime")
        _preflight_immunogenicity_backend(run_id, backend_repitope, "repitope")
    if step3 in selected:
        _preflight_mhc1_backend(run_id, backend_mhc1_netmhcpan, backend_mhc1_bigmhc)
        _preflight_mhc2_backend(run_id, mhc2_backend)
    if prepare_pandora and step8 in selected:
        probe = subprocess.run(
            [pandora_exe, "-c", "import PANDORA"],
            capture_output=True,
            text=True,
            check=False,
        )
        if probe.returncode != 0:
            raise RuntimeError(
                "已启用 --prepare_pandora，但当前 PANDORA 解释器无法 import PANDORA。"
                f" 解释器: {pandora_exe}；stderr: {(probe.stderr or '').strip()}。"
                " 请安装 csb-pandora 或使用 --pandora_python 指向已安装该包的环境（例如 pandora_src 的 python）。"
            )
    for step in selected:
        if step is step8 and prepare_structure_inputs:
            run_step(step8_1)
        if step is step8 and prepare_pandora:
            run_step(step8_pandora, env=_pandora_step_env(pandora_exe))
        run_step(step)
    print("\n所选步骤执行完成。")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--run_id", required=True, help="例如 R001")
    parser.add_argument("--top_n", type=int, default=10, help="Top 候选肽数量，默认 10")
    parser.add_argument(
        "--min_dissimilarity",
        type=float,
        default=0.1,
        help="与 WT 的最小差异比例阈值，默认 0.1",
    )
    parser.add_argument("--linker", default="AAY", help="肽段连接子，默认 AAY")
    parser.add_argument("--poly_a_len", type=int, default=120, help="polyA 长度，默认 120")
    parser.add_argument(
        "--codon_mode",
        default="optimized",
        choices=["basic", "optimized", "lineardesign", "real_cmd"],
        help="密码子模式，默认 optimized",
    )
    parser.add_argument(
        "--codon_real_tool",
        default="",
        help="codon_mode=real_cmd 时记录的真实密码子优化工具名称",
    )
    parser.add_argument(
        "--codon_real_cmd",
        default="",
        help="codon_mode=real_cmd 时的真实密码子优化命令模板",
    )
    parser.add_argument(
        "--codon_real_version_cmd",
        default="",
        help="真实密码子优化工具版本命令，用于记录验收元数据",
    )
    parser.add_argument("--signal_peptide_preset", default="", choices=["", "igkv", "mitd"], help="N端信号肽预设")
    parser.add_argument("--signal_peptide_aa", default="", help="N端信号肽手工AA")
    parser.add_argument("--tm_domain_preset", default="", choices=["", "cd8a_tm"], help="C端TM预设")
    parser.add_argument("--tm_domain_aa", default="", help="C端TM手工AA")
    parser.add_argument(
        "--top_k_md",
        type=int,
        default=3,
        help="SimHub 多维 MD 候选数（与 selected_for_md 对齐）；与 run_pandora_structure --top_k 一致，建议 3–5。",
    )
    parser.add_argument(
        "--prepare_pandora",
        action="store_true",
        help="在 prepare_simhub_delivery 前调用 run_pandora_structure.py（--top_k = top_k_md）；需 PANDORA 可用，可用 --pandora_python 指定解释器。",
    )
    parser.add_argument(
        "--pandora_python",
        default="",
        help=(
            "仅用于 run_pandora_structure：已安装 csb-pandora（及 MODELLER 等）的 python 可执行文件路径；"
            "留空则与当前流程使用同一解释器。"
        ),
    )
    parser.add_argument(
        "--prepare_structure_inputs",
        action="store_true",
        help="生成结构建模输入序列（MHC-I alpha + b2m），输出到 results/<run_id>/structure_inputs/",
    )
    parser.add_argument("--structure_seq_refresh_remote", action="store_true", help="刷新 IPD/UniProt 远程序列缓存")
    parser.add_argument("--structure_seq_strict", action="store_true", help="拉取/匹配失败时中断流程")
    parser.add_argument(
        "--structure_backend",
        default="pandora",
        choices=["coarse", "pandora", "afm"],
        help="SimHub 结构来源：pandora(默认)/afm/coarse；真实交付需提供 --structure_input_pdb",
    )
    parser.add_argument(
        "--structure_input_pdb",
        default="",
        help=(
            "structure_backend 为 afm 时必须提供；为 pandora 且留空时，"
            "prepare_simhub_delivery 从 results/structure_models/pandora/.../rank_XX/complex.pdb 自动解析。"
        ),
    )
    parser.add_argument("--feasibility_top_n", type=int, default=10, help="可行性验证 TopN，默认 10")
    parser.add_argument(
        "--mhc2_backend",
        default="real_tsv",
        choices=["auto", "proxy", "real_tsv", "netmhciipan"],
        help="MHC-II：默认 real_tsv，读取 results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv；auto 也不会静默回退 proxy。",
    )
    parser.add_argument("--wi_deepimmuno", type=float, default=1.0, help="immunogenicity_deepimmuno 子权重")
    parser.add_argument("--wi_prime", type=float, default=1.0, help="immunogenicity_prime 子权重")
    parser.add_argument("--wi_repitope", type=float, default=1.0, help="immunogenicity_repitope 子权重")
    parser.add_argument(
        "--backend_mhc1_netmhcpan",
        default="real_tsv",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 NetMHCpan 后端",
    )
    parser.add_argument(
        "--backend_mhc1_bigmhc",
        default="off",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 BigMHC 后端（默认 off；NetMHCpan real_tsv 作为必需真实交叉验证）",
    )
    parser.add_argument("--require_real_mhc2", action="store_true", default=True, help="要求 MHC-II 必须 real（默认开启）")
    parser.add_argument("--require_real_mhc1_cv", action="store_true", default=True, help="要求 MHC-I 交叉验证至少一个 real（默认开启）")
    parser.add_argument(
        "--backend_deepimmuno",
        default="real_tsv",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="DeepImmuno 适配器后端",
    )
    parser.add_argument(
        "--backend_prime",
        default="real_tsv",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="PRIME 适配器后端",
    )
    parser.add_argument(
        "--backend_repitope",
        default="real_tsv",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="Repitope 适配器后端",
    )
    parser.add_argument(
        "--require_real_immunogenicity_deepimmuno",
        action="store_true",
        default=True,
        help="要求 DeepImmuno 必须使用真实来源（默认开启）。",
    )
    parser.add_argument(
        "--require_real_immunogenicity_prime",
        action="store_true",
        default=True,
        help="要求 PRIME 必须使用真实来源（real_tsv/real_cmd），若回退 proxy 则报错。",
    )
    parser.add_argument(
        "--require_real_immunogenicity_repitope",
        action="store_true",
        default=True,
        help="要求 Repitope 必须使用真实来源（real_tsv/real_cmd），若回退 proxy 则报错。",
    )
    parser.add_argument(
        "--allow_proxy_scores",
        action="store_true",
        help="显式允许 MHC-II 或免疫原性缺失时使用代理分；默认禁止。",
    )
    parser.add_argument(
        "--target",
        default="full",
        choices=[
            "full",
            "mhc_ranking",
            "report",
            "simhub",
            "feasibility",
            "transcriptome_prior",
            "b_to_a_simhub",
        ],
        help=(
            "执行范围：full=路径 A 全流程；b_to_a_simhub=基因池过滤后路径 A 全流程；"
            "transcriptome_prior=仅 bulk 基因侧车（见 README §1.1A）。"
        ),
    )
    parser.add_argument(
        "--b_to_a_top_n_genes",
        type=int,
        default=500,
        help="target=b_to_a_simhub：取路径 B Top 基因数过滤突变肽，默认 500。",
    )
    parser.add_argument(
        "--b_to_a_rank_subdir",
        default="transcriptome_prior",
        help="target=b_to_a_simhub：基因排名子目录（位于 results/<run_id>/ 下）。",
    )
    parser.add_argument(
        "--b_to_a_min_log2fc",
        type=float,
        default=None,
        help="target=b_to_a_simhub：过滤前仅保留 log2FC >= 该值的基因。",
    )
    parser.add_argument(
        "--skip_b_to_a_filter",
        action="store_true",
        help="target=b_to_a_simhub：跳过基因过滤（neoantigen 已手动收窄时）。",
    )
    parser.add_argument(
        "--bulk_tpm_path",
        default="",
        help="仅 target=transcriptome_prior：基因×样本表达矩阵路径（首列基因名，其余为数值）。",
    )
    parser.add_argument(
        "--bulk_case_columns",
        default="",
        help="仅 target=transcriptome_prior：病例组样本列名，逗号分隔，须与矩阵表头一致。",
    )
    parser.add_argument(
        "--bulk_control_columns",
        default="",
        help="仅 target=transcriptome_prior：对照组样本列名，逗号分隔。",
    )
    parser.add_argument(
        "--bulk_gene_col_name",
        default="gene_name",
        help="仅 target=transcriptome_prior：输出 CSV 中基因列显示名（默认 gene_name）。",
    )
    parser.add_argument(
        "--bulk_pseudo",
        type=float,
        default=1.0,
        help="仅 target=transcriptome_prior：log2FC 均值平滑伪计数，默认 1。",
    )
    parser.add_argument(
        "--bulk_out_subdir",
        default="transcriptome_prior",
        help="仅 target=transcriptome_prior：输出子目录名，位于 results/<run_id>/<subdir>/，默认 transcriptome_prior。",
    )
    parser.add_argument(
        "--cohort_id",
        default="",
        help=(
            "仅 target=transcriptome_prior：队列数字 ID（如 39004），"
            "自动从 data/表达定量表格.zip 与 样本信息表.zip 解析分组并运行；"
            "建议 run_id=T_cohort_<cohort_id>。"
        ),
    )
    parser.add_argument(
        "--mrna_output_suffix",
        default="",
        help=(
            "传给 build_multivalent_mrna：可选输出后缀，同一 run 生成 mrna_vaccine_<后缀>.fasta 等多条候选；"
            "留空则保持默认 mrna_vaccine.fasta（与 QC/SimHub 一致）。"
        ),
    )
    parser.add_argument(
        "--ensure_positive_control_peptides",
        default="",
        help=(
            "传给 select_top_peptides：逗号分隔 mut_peptide；若在过滤池中但未进 Top-N，则用其替换选中集中 rank_score 最低的一条。"
        ),
    )
    args = parser.parse_args()

    main(
        run_id=args.run_id,
        top_n=args.top_n,
        min_dissimilarity=args.min_dissimilarity,
        linker=args.linker,
        poly_a_len=args.poly_a_len,
        codon_mode=args.codon_mode,
        codon_real_tool=args.codon_real_tool,
        codon_real_cmd=args.codon_real_cmd,
        codon_real_version_cmd=args.codon_real_version_cmd,
        signal_peptide_preset=args.signal_peptide_preset,
        signal_peptide_aa=args.signal_peptide_aa,
        tm_domain_preset=args.tm_domain_preset,
        tm_domain_aa=args.tm_domain_aa,
        top_k_md=args.top_k_md,
        prepare_pandora=args.prepare_pandora,
        prepare_structure_inputs=args.prepare_structure_inputs,
        structure_seq_refresh_remote=args.structure_seq_refresh_remote,
        structure_seq_strict=args.structure_seq_strict,
        structure_backend=args.structure_backend,
        structure_input_pdb=args.structure_input_pdb,
        feasibility_top_n=args.feasibility_top_n,
        mhc2_backend=args.mhc2_backend,
        wi_deepimmuno=args.wi_deepimmuno,
        wi_prime=args.wi_prime,
        wi_repitope=args.wi_repitope,
        backend_mhc1_netmhcpan=args.backend_mhc1_netmhcpan,
        backend_mhc1_bigmhc=args.backend_mhc1_bigmhc,
        require_real_mhc2=args.require_real_mhc2,
        require_real_mhc1_cv=args.require_real_mhc1_cv,
        backend_deepimmuno=args.backend_deepimmuno,
        backend_prime=args.backend_prime,
        backend_repitope=args.backend_repitope,
        require_real_immunogenicity_deepimmuno=args.require_real_immunogenicity_deepimmuno,
        require_real_immunogenicity_prime=args.require_real_immunogenicity_prime,
        require_real_immunogenicity_repitope=args.require_real_immunogenicity_repitope,
        allow_proxy_scores=args.allow_proxy_scores,
        target=args.target,
        mrna_output_suffix=args.mrna_output_suffix,
        ensure_positive_control_peptides=args.ensure_positive_control_peptides,
        pandora_python=args.pandora_python,
        bulk_tpm_path=args.bulk_tpm_path,
        bulk_case_columns=args.bulk_case_columns,
        bulk_control_columns=args.bulk_control_columns,
        bulk_gene_col_name=args.bulk_gene_col_name,
        bulk_pseudo=args.bulk_pseudo,
        bulk_out_subdir=args.bulk_out_subdir,
        cohort_id=args.cohort_id,
        b_to_a_top_n_genes=args.b_to_a_top_n_genes,
        b_to_a_rank_subdir=args.b_to_a_rank_subdir,
        b_to_a_min_log2fc=args.b_to_a_min_log2fc,
        skip_b_to_a_filter=args.skip_b_to_a_filter,
    )
