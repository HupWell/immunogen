#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LinearDesign 命令包装器：固定工作目录，避免官方二进制找不到相对路径依赖。"""
import argparse
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lineardesign_dir", default="external_refs/LinearDesign")
    parser.add_argument("--lambda_value", default="0.3")
    parser.add_argument("--version", action="store_true")
    args = parser.parse_args()

    root = Path(args.lineardesign_dir).resolve()
    binary = root / "lineardesign"
    if args.version:
        print("LinearDesign official main; source=LinearDesignSoftware/LinearDesign; wrapper=tools/lineardesign_runner.py")
        return 0
    if not binary.exists():
        raise FileNotFoundError(f"未找到 LinearDesign 可执行文件: {binary}")

    protein = sys.stdin.read()
    result = subprocess.run(
        [str(binary), "--lambda", str(args.lambda_value)],
        input=protein,
        capture_output=False,
        text=True,
        check=False,
        cwd=str(root),
    )
    return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
