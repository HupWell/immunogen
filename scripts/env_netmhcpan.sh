#!/usr/bin/env bash
# NetMHCpan 4.2+ 与本仓库 MHC-I real_cmd 对接示例。
# 用法：在仓库根目录执行  source scripts/env_netmhcpan.sh
# 请按本机实际路径修改 NETMHCPAN_HOME。

if [[ -n "${BASH_VERSION:-}" ]] || [[ -n "${ZSH_VERSION:-}" ]]; then
  _MHC_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")/.." && pwd)"
else
  _MHC_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
fi

export NETMHCPAN_HOME="${NETMHCPAN_HOME:-/root/netMHCpan-4.2}"
export NETMHCPAN_BIN="${NETMHCPAN_BIN:-$NETMHCPAN_HOME/netMHCpan}"
export MHC1_NETMHCPAN_CMD="python \"${_MHC_ROOT}/tools/netmhcpan_class1_runner.py\" --input {input_tsv} --output {output_tsv}"

echo "[env_netmhcpan] NETMHCPAN_HOME=$NETMHCPAN_HOME"
echo "[env_netmhcpan] NETMHCPAN_BIN=$NETMHCPAN_BIN"
echo "[env_netmhcpan] MHC1_NETMHCPAN_CMD 已设置（含仓库: ${_MHC_ROOT}）"

unset _MHC_ROOT
