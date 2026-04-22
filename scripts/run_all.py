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


def main(
    run_id: str,
    top_n: int,
    min_dissimilarity: float,
    linker: str,
    poly_a_len: int,
    codon_mode: str,
    top_k_md: int,
    feasibility_top_n: int,
    mhc2_backend: str,
    wi_deepimmuno: float,
    wi_prime: float,
    wi_repitope: float,
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
    ]
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
    ]
    step9 = [
        python_exe,
        os.path.join(scripts_dir, "validate_feasibility.py"),
        "--run_id",
        run_id,
        "--top_n",
        str(feasibility_top_n),
    ]

    print("开始执行 ImmunoGen MVP 全流程...")
    run_step(step1)
    run_step(step2)
    run_step(step3)
    run_step(step4)
    run_step(step5)
    run_step(step6)
    run_step(step7)
    run_step(step8)
    run_step(step9)
    print("\n全部步骤执行完成。")


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
    parser.add_argument("--top_k_md", type=int, default=3, help="提交 SimHub 的 Top K，默认 3")
    parser.add_argument("--feasibility_top_n", type=int, default=10, help="可行性验证 TopN，默认 10")
    parser.add_argument(
        "--mhc2_backend",
        default="auto",
        choices=["auto", "proxy", "netmhciipan"],
        help="MHC-II：auto=有 NetMHCIIpan 与 II 分型则实跑否则代理；见 predict_mhc_ranking 说明。",
    )
    parser.add_argument("--wi_deepimmuno", type=float, default=1.0, help="immunogenicity_deepimmuno 子权重")
    parser.add_argument("--wi_prime", type=float, default=1.0, help="immunogenicity_prime 子权重")
    parser.add_argument("--wi_repitope", type=float, default=1.0, help="immunogenicity_repitope 子权重")
    args = parser.parse_args()

    main(
        run_id=args.run_id,
        top_n=args.top_n,
        min_dissimilarity=args.min_dissimilarity,
        linker=args.linker,
        poly_a_len=args.poly_a_len,
        codon_mode=args.codon_mode,
        top_k_md=args.top_k_md,
        feasibility_top_n=args.feasibility_top_n,
        mhc2_backend=args.mhc2_backend,
        wi_deepimmuno=args.wi_deepimmuno,
        wi_prime=args.wi_prime,
        wi_repitope=args.wi_repitope,
    )
