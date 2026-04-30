# P3-4 ColabFold GPU 环境增强验证

## 结论

- GPU 版 JAX 已在 `colabfold` 环境中修复并验证。
- `colabfold_batch` 已真正进入 GPU 计算阶段，并完成 `R_public_001/public_case_001` 的一次 AFM 复核。
- AlphaFold-multimer v3 参数已在 `external_refs/colabfold_data` 缓存中可用。
- AFM 已产出 PDB、PAE、pLDDT 和 score JSON。
- 本轮 AFM 置信度较低：mean pLDDT `31.79`，pTM `0.18`，ipTM `0.11`；因此不建议替换当前 PANDORA 交付结构。
- 详细复核记录见：`docs/P3_4_AFM_REVIEW_R_PUBLIC_001_2026-04-30.md`。

## GPU / JAX 验证

执行命令：

```bash
conda run -n colabfold python -c "import jax; print(jax.__version__); print(jax.devices()); print(jax.default_backend())"
```

结果：

```text
0.6.2
[CudaDevice(id=0)]
gpu
```

说明：此前 `jax 0.6.2` 只能识别 CPU；已通过安装匹配版本的 CUDA 12 JAX 组件修复：

```bash
conda run -n colabfold python -m pip install -U "jax[cuda12]==0.6.2"
```

## 重点候选输入

重点复核对象：

- run_id：`R_public_001`
- case_id：`public_case_001`
- peptide：`AAAAVIPTV`
- HLA：`HLA-C*01:02`

已重新生成 MHC-I alpha / B2M 输入：

```bash
python scripts/prepare_mhc_chain_sequences.py --run_id R_public_001 --strict
```

输出：

- `results/R_public_001/structure_inputs/mhc_chain_sequences.fasta`
- `results/R_public_001/structure_inputs/mhc_chain_sequences.json`
- `results/R_public_001/structure_inputs/public_case_001_colabfold_multimer.fasta`

ColabFold multimer 输入格式：

```text
MHC_ALPHA:B2M:PEPTIDE
```

## ColabFold 运行状态

尝试运行：

```bash
conda run -n colabfold colabfold_batch \
  --msa-mode single_sequence \
  --model-type alphafold2_multimer_v3 \
  --num-models 1 \
  --num-recycle 1 \
  --num-seeds 1 \
  --num-ensemble 1 \
  --rank multimer \
  --overwrite-existing-results \
  results/R_public_001/structure_inputs/public_case_001_colabfold_multimer.fasta \
  results/structure_models/afm/R_public_001/public_case_001
```

现象：

- `colabfold_batch` 启动并写入 `log.txt`
- 日志记录 `Running on GPU`
- 完成 1 个 model、1 次 recycle 的 AFM multimer 预测
- 输出 unrelaxed PDB、score JSON、PAE 图、pLDDT 图

参数源检查：

```text
https://storage.googleapis.com/alphafold/alphafold_params_colab_2022-12-06.tar
status=200
Content-Length=4099624960
```

## 与 PANDORA 对比状态

PANDORA 结构已存在：

- `results/structure_models/pandora/R_public_001/public_case_001/complex.pdb`

AFM / ColabFold PDB 已生成：

- `results/structure_models/afm/R_public_001/public_case_001/R_public_001_public_case_001_top1_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb`
- `results/structure_models/afm/R_public_001/public_case_001/complex_afm_chains_MBP.pdb`
- `results/structure_models/afm/R_public_001/public_case_001/afm_vs_pandora_comparison.json`

当前结论：

- AFM receptor 链长度与 PANDORA 交付结构不同，不直接给出 full-chain RMSD。
- peptide 链 CA Kabsch RMSD 为 `8.697 Å`。
- 因 AFM pLDDT / pTM / ipTM 偏低，当前仅作为低置信二级复核证据，不替换 PANDORA 交付结构。

## 下一步

1. 完成 AlphaFold-multimer v3 参数下载到大容量目录，例如：

   ```bash
   conda run -n colabfold python -c "from pathlib import Path; from colabfold.download import download_alphafold_params; download_alphafold_params('alphafold2_multimer_v3', Path('external_refs/colabfold_data'))"
   ```

2. 使用显式数据目录重跑 ColabFold：

   ```bash
   conda run -n colabfold colabfold_batch \
     --data external_refs/colabfold_data \
     --msa-mode single_sequence \
     --model-type alphafold2_multimer_v3 \
     --num-models 1 \
     --num-recycle 1 \
     --num-seeds 1 \
     --num-ensemble 1 \
     --rank multimer \
     --overwrite-existing-results \
     results/R_public_001/structure_inputs/public_case_001_colabfold_multimer.fasta \
     results/structure_models/afm/R_public_001/public_case_001
   ```

3. 生成 AFM PDB 后，再与 PANDORA PDB 做结构差异复核并更新本报告。
