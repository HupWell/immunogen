# P3-4 AFM / ColabFold GPU 复核记录（R_public_001）

## 结论

- 已真正使用 GPU 跑完一次 ColabFold / AlphaFold-Multimer v3 复核。
- 复核对象：`R_public_001/public_case_001`，peptide 为 `AAAAVIPTV`。
- AFM 已生成 PDB、PAE、pLDDT 图和 score JSON。
- 但本次 AFM 置信度较低：mean pLDDT `31.79`，pTM `0.18`，ipTM `0.11`。
- 因此当前不建议用该 AFM 结果替换 SimHub 交付用 PANDORA 结构；AFM 结果仅作为低置信二级复核证据归档。

## GPU 运行证据

日志文件：

- `results/structure_models/afm/R_public_001/public_case_001/log.txt`

关键日志：

```text
Running on GPU
rank_001_alphafold2_multimer_v3_model_1_seed_000 pLDDT=31.8 pTM=0.18 ipTM=0.107
Done
```

## 输出文件

- AFM 原始 PDB：`results/structure_models/afm/R_public_001/public_case_001/R_public_001_public_case_001_top1_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb`
- AFM 链 ID 归一化 PDB：`results/structure_models/afm/R_public_001/public_case_001/complex_afm_chains_MBP.pdb`
- 对比摘要 JSON：`results/structure_models/afm/R_public_001/public_case_001/afm_vs_pandora_comparison.json`
- PAE 图：`results/structure_models/afm/R_public_001/public_case_001/R_public_001_public_case_001_top1_pae.png`
- pLDDT 图：`results/structure_models/afm/R_public_001/public_case_001/R_public_001_public_case_001_top1_plddt.png`

## 与 PANDORA 对比

PANDORA 当前交付结构：

- `results/structure_models/pandora/R_public_001/public_case_001/complex.pdb`

链统计：

- PANDORA：M 链 `274` residues，B 链 `100` residues，P 链 `9` residues。
- AFM：A 链 `366` residues，B 链 `119` residues，C 链 `9` residues。

比较口径：

- AFM 输入使用完整 HLA-C alpha / B2M 序列。
- PANDORA 交付结构为模板建模/截断后的 pMHC 结构。
- 两者 receptor 链长度不同，因此不直接给出 full-chain RMSD。
- peptide 链按 CA 原子 Kabsch 对齐后 RMSD 为 `8.697 Å`。

## 交付建议

- SimHub 当前生产交付仍使用 PANDORA 全原子 pMHC：`deliveries/R_public_001/to_simhub/public_case_001/complex.pdb`。
- AFM 结果已证明 GPU 复核链路可运行，但由于本轮置信度低，不作为替换结构。
- 如果后续要提升 AFM 可信度，建议增加 MSA、提高 recycle 数、运行多个 seed/model 后再筛选。
