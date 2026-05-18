#!/usr/bin/env bash
# MHCflurry 2.x 默认需要 models_class1_pan（体积大）。
# 若已通过「mhcflurry-downloads fetch … models_class1」或本机已有 models_class1/models/manifest.csv，
# 可设置 MHCFLURRY_DEFAULT_CLASS1_MODELS，使 Class1AffinityPredictor.load() 直接加载该目录。
# 用法：在仓库根目录执行  source scripts/env_mhcflurry.sh

_DEFAULT_MODELS="${HOME}/.local/share/mhcflurry/4/2.2.0/models_class1/models"

if [[ -f "${_DEFAULT_MODELS}/manifest.csv" ]]; then
  export MHCFLURRY_DEFAULT_CLASS1_MODELS="${_DEFAULT_MODELS}"
  echo "[env_mhcflurry] MHCFLURRY_DEFAULT_CLASS1_MODELS=${MHCFLURRY_DEFAULT_CLASS1_MODELS}"
else
  echo "[env_mhcflurry] 未找到 ${_DEFAULT_MODELS}/manifest.csv，跳过设置（请运行 mhcflurry-downloads fetch models_class1）"
fi

unset _DEFAULT_MODELS
