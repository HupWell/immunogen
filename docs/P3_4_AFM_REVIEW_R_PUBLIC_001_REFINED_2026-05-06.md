# P3-4 AFM 复核质量提升记录（R_public_001）

## 结论（相对 2026-04-30 本轮）

- 本轮把 AFM 复核参数升级为：`msa-mode=mmseqs2_uniref`、`recycle=3`、`num-models=2`（并使用更多 seed）。
- 复核置信度显著提升：从先前平均 `pLDDT≈31.79 / pTM≈0.18` 提升到 `pLDDT≈83+`、`pTM≈0.78-0.79`、`ipTM≈0.85-0.86`（以 log 为准）。
- 生产交付（SimHub）仍建议继续使用 PANDORA 的结构作为主结构；本轮 AFM 结果更适合作为“更强证据归档/二级复核”，以降低未来替换结构的风险。

## AFM 运行参数

对象：
- `run_id=R_public_001`
- `case_id=public_case_001`
- peptide：`AAAAVIPTV`
- HLA：`HLA-C*01:02`

关键输出目录：
- `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/`

关键运行参数（colabfold_batch）：
- `--msa-mode mmseqs2_uniref`
- `--num-recycle 3`
- `--num-models 2`
- `--num-seeds 3`
- `--rank multimer`
- `--model-type alphafold2_multimer_v3`

## 关键分数（来自 log）

（seed/model 均在同一目录的 `log.txt` 中可核对；下面列 recycle=3 的代表性结果，并补充 seed=000 的中间迭代用于说明稳定性。）

- `alphafold2_multimer_v3_model_1_seed_000`
  - recycle=1：`pLDDT=83.9 pTM=0.791 ipTM=0.862`
  - recycle=3：`pLDDT=83.9 pTM=0.787 ipTM=0.854`
- `alphafold2_multimer_v3_model_2_seed_000`
  - recycle=1：`pLDDT=83.1 pTM=0.781 ipTM=0.851`
  - recycle=3：`pLDDT=83.9 pTM=0.784 ipTM=0.852`
- `alphafold2_multimer_v3_model_2_seed_002`
  - recycle=3：`pLDDT=83.8 pTM=0.784 ipTM=0.852`
- `alphafold2_multimer_v3_model_1_seed_001`
  - recycle=3：`pLDDT=83.6 pTM=0.787 ipTM=0.855`
- `alphafold2_multimer_v3_model_2_seed_001`
  - recycle=3：`pLDDT=83.4 pTM=0.782 ipTM=0.852`
- `alphafold2_multimer_v3_model_1_seed_002`
  - recycle=3：`pLDDT=83.6 pTM=0.787 ipTM=0.854`

## 输出文件（已生成）

- score JSON / AFM PDB（按 `multimer` 指标的 rank_001~rank_006 结果）：
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_001_alphafold2_multimer_v3_model_1_seed_002.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_002.pdb`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_002_alphafold2_multimer_v3_model_1_seed_001.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_002_alphafold2_multimer_v3_model_1_seed_001.pdb`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_003_alphafold2_multimer_v3_model_1_seed_000.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_003_alphafold2_multimer_v3_model_1_seed_000.pdb`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_004_alphafold2_multimer_v3_model_2_seed_002.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_004_alphafold2_multimer_v3_model_2_seed_002.pdb`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_005_alphafold2_multimer_v3_model_2_seed_000.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_005_alphafold2_multimer_v3_model_2_seed_000.pdb`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_scores_rank_006_alphafold2_multimer_v3_model_2_seed_001.json`
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/R_public_001_public_case_001_top1_unrelaxed_rank_006_alphafold2_multimer_v3_model_2_seed_001.pdb`
- log：
  - `results/structure_models/afm/R_public_001/public_case_001/afm_refined_msarecycle3_models2_seeds3/log.txt`

## 与 PANDORA 的交付建议

- 本轮 AFM 置信度提升明显，可以作为“更强的二级结构证据”归档。
- 是否替换 SimHub 生产交付结构，仍需要做一次与 PANDORA 的严格结构对齐/对比确认（尤其是 peptide 链与 receptor 关键链段）。

