# -*- coding: utf-8 -*-
"""
功能：一键串行执行 ImmunoGen MVP 全流程。
执行顺序：
1) validate_input.py
2) run_immunogenicity_adapters.py
3) predict_mhc_ranking.py
4) select_top_peptides.py
5) build_multivalent_mrna.py
6) run_qc_and_report.py
7) prepare_self_certification.py（POSITIVE_CONTROL.md / SELF_CHECK.md）
8) prepare_simhub_delivery.py
9) validate_feasibility.py
"""
import os
import sys
import argparse
import subprocess


def run_step(command: list):
    """执行单个步骤并打印命令，失败时立即退出。"""
    print(f"\n[执行] {' '.join(command)}")
    result = subprocess.run(command, check=False)
    if result.returncode != 0:
        raise RuntimeError(f"步骤执行失败，退出码: {result.returncode}")


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


def main(
    run_id: str,
    top_n: int,
    min_dissimilarity: float,
    linker: str,
    poly_a_len: int,
    codon_mode: str,
    signal_peptide_preset: str,
    signal_peptide_aa: str,
    tm_domain_preset: str,
    tm_domain_aa: str,
    top_k_md: int,
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
    target: str,
):
    """按固定顺序执行全流程。"""
    python_exe = sys.executable
    scripts_dir = os.path.dirname(os.path.abspath(__file__))

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
        "--signal_peptide_preset",
        signal_peptide_preset,
        "--signal_peptide_aa",
        signal_peptide_aa,
        "--tm_domain_preset",
        tm_domain_preset,
        "--tm_domain_aa",
        tm_domain_aa,
    ]
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
    step9 = [
        python_exe,
        os.path.join(scripts_dir, "validate_feasibility.py"),
        "--run_id",
        run_id,
        "--top_n",
        str(feasibility_top_n),
    ]

    if require_real_mhc1_cv and backend_mhc1_netmhcpan == "off" and backend_mhc1_bigmhc == "off":
        raise RuntimeError("已启用 --require_real_mhc1_cv，但两个 MHC-I 交叉验证后端均为 off。")

    pipelines = {
        "full": [step1, step2, step3, step4, step5, step6, step7, step8, step9],
        "mhc_ranking": [step1, step2, step3],
        "report": [step6, step7],
        "simhub": [step8],
        "feasibility": [step9],
    }
    selected = pipelines[target]

    print(f"开始执行 ImmunoGen 流程（target={target}）...")
    if step3 in selected:
        _preflight_mhc1_backend(run_id, backend_mhc1_netmhcpan, backend_mhc1_bigmhc)
        _preflight_mhc2_backend(run_id, mhc2_backend)
    for step in selected:
        if step is step8 and prepare_structure_inputs:
            run_step(step8_1)
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
        choices=["basic", "optimized", "lineardesign"],
        help="密码子模式，默认 optimized",
    )
    parser.add_argument("--signal_peptide_preset", default="", choices=["", "igkv", "mitd"], help="N端信号肽预设")
    parser.add_argument("--signal_peptide_aa", default="", help="N端信号肽手工AA")
    parser.add_argument("--tm_domain_preset", default="", choices=["", "cd8a_tm"], help="C端TM预设")
    parser.add_argument("--tm_domain_aa", default="", help="C端TM手工AA")
    parser.add_argument("--top_k_md", type=int, default=3, help="提交 SimHub 的 Top K，默认 3")
    parser.add_argument(
        "--prepare_structure_inputs",
        action="store_true",
        help="生成结构建模输入序列（MHC-I alpha + b2m），输出到 results/<run_id>/structure_inputs/",
    )
    parser.add_argument("--structure_seq_refresh_remote", action="store_true", help="刷新 IPD/UniProt 远程序列缓存")
    parser.add_argument("--structure_seq_strict", action="store_true", help="拉取/匹配失败时中断流程")
    parser.add_argument(
        "--structure_backend",
        default="coarse",
        choices=["coarse", "pandora", "afm"],
        help="SimHub 结构来源：coarse(默认)/pandora/afm",
    )
    parser.add_argument("--structure_input_pdb", default="", help="若 structure_backend 非 coarse，提供外部 PDB 路径")
    parser.add_argument("--feasibility_top_n", type=int, default=10, help="可行性验证 TopN，默认 10")
    parser.add_argument(
        "--mhc2_backend",
        default="auto",
        choices=["auto", "proxy", "real_tsv", "netmhciipan"],
        help="MHC-II：auto=优先 netmhciipan，其次 real_tsv，最后 proxy；real_tsv 读取 results/<run_id>/tool_outputs/raw/mhc2_netmhciipan.tsv。",
    )
    parser.add_argument("--wi_deepimmuno", type=float, default=1.0, help="immunogenicity_deepimmuno 子权重")
    parser.add_argument("--wi_prime", type=float, default=1.0, help="immunogenicity_prime 子权重")
    parser.add_argument("--wi_repitope", type=float, default=1.0, help="immunogenicity_repitope 子权重")
    parser.add_argument(
        "--backend_mhc1_netmhcpan",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 NetMHCpan 后端",
    )
    parser.add_argument(
        "--backend_mhc1_bigmhc",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "off"],
        help="MHC-I 交叉验证 BigMHC 后端",
    )
    parser.add_argument("--require_real_mhc2", action="store_true", help="要求 MHC-II 必须 real（netmhciipan）")
    parser.add_argument("--require_real_mhc1_cv", action="store_true", help="要求 MHC-I 交叉验证至少一个 real")
    parser.add_argument(
        "--backend_deepimmuno",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="DeepImmuno 适配器后端",
    )
    parser.add_argument(
        "--backend_prime",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="PRIME 适配器后端",
    )
    parser.add_argument(
        "--backend_repitope",
        default="auto",
        choices=["auto", "real_tsv", "real_cmd", "proxy"],
        help="Repitope 适配器后端",
    )
    parser.add_argument(
        "--target",
        default="full",
        choices=["full", "mhc_ranking", "report", "simhub", "feasibility"],
        help="执行范围：full=全流程；mhc_ranking=到表位排名；report=仅报告/自证；simhub=仅交付封装；feasibility=仅可行性。",
    )
    args = parser.parse_args()

    main(
        run_id=args.run_id,
        top_n=args.top_n,
        min_dissimilarity=args.min_dissimilarity,
        linker=args.linker,
        poly_a_len=args.poly_a_len,
        codon_mode=args.codon_mode,
        signal_peptide_preset=args.signal_peptide_preset,
        signal_peptide_aa=args.signal_peptide_aa,
        tm_domain_preset=args.tm_domain_preset,
        tm_domain_aa=args.tm_domain_aa,
        top_k_md=args.top_k_md,
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
        target=args.target,
    )
