#!/usr/bin/env bash
# DeepImmuno / Repitope 免疫原性真实后端环境变量示例。
# 用法：在仓库根目录执行 source scripts/env_immunogenicity.sh

if [[ -n "${BASH_VERSION:-}" ]] || [[ -n "${ZSH_VERSION:-}" ]]; then
  _IMMUNO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
else
  _IMMUNO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

export IMMUNO_DEEPIMMUNO_CMD="python \"${_IMMUNO_ROOT}/tools/deepimmuno_runner.py\" --input {input_tsv} --output {output_tsv}"
export IMMUNO_REPITOPE_CMD="python \"${_IMMUNO_ROOT}/tools/repitope_runner.py\" --input {input_tsv} --output {output_tsv}"

echo "[env_immunogenicity] IMMUNO_DEEPIMMUNO_CMD 已设置（仓库: ${_IMMUNO_ROOT}）"
echo "[env_immunogenicity] IMMUNO_REPITOPE_CMD 已设置（仓库: ${_IMMUNO_ROOT}）"

unset _IMMUNO_ROOT
