# P2 PANDORA 真实结构生成与验收记录（2026-04-28）

## 1. 数据库与环境

- PANDORA 环境：`pandora_src`
- MODELLER：已通过真实 license 校验
- PANDORA 官方数据库：Zenodo `PANDORA Database`，文件 `default.tar.gz`
- 下载校验：`md5=8e75734b7ed6394742e3c7d2d674396b`
- 模板库规模：MHC-I 864 个模板，MHC-II 136 个模板

## 2. 执行命令

```bash
for r in R001 R002 R003 R_public_001; do
  conda run -n pandora_src python scripts/run_pandora_structure.py \
    --run_id "$r" --n_loop_models 5 --n_jobs 5
done
```

随后使用真实 PDB 刷新 SimHub 交付：

```bash
python scripts/prepare_simhub_delivery.py \
  --run_id <run_id> \
  --structure_backend pandora \
  --structure_input_pdb results/structure_models/pandora/<run_id>/<case_id>/complex.pdb
```

## 3. 结构结果

| Run | Case | Peptide | HLA-I | PANDORA template | MODELLER models | Best molpdf |
|---|---|---|---|---|---:|---:|
| `R001` | `demo_case_001` | `HMTEVVRHC` | `HLA-A*02:01` | `6VR5` | 5 | -4.8826 |
| `R002` | `public_case_002` | `ARHGGWTTK` | `HLA-A*02:01` | `1LDP` | 5 | 5.0982 |
| `R003` | `public_case_003` | `ARHGGWTTK` | `HLA-A*02:01` | `1LDP` | 5 | 5.0982 |
| `R_public_001` | `public_case_001` | `AAAAVIPTV` | `HLA-C*01:02` | `6GH4` | 5 | 35.4779 |

## 4. 交付验收

已执行：

```bash
for r in R001 R002 R003 R_public_001; do
  python scripts/check_epitope_realization.py --run_id "$r" \
    --require_mhc1_cv_real \
    --require_mhc2_real \
    --require_real_immunogenicity \
    --require_real_structure
done
```

四个 run 均通过：

- `mhc2_real=True`
- `mhc1_cv_real=True`
- `immunogenicity_proxy=False`
- `real_structure=True`

## 5. PDB 链检查

| Run | Chain M residues | Chain B residues | Chain P residues |
|---|---:|---:|---:|
| `R001` | 274 | 98 | 9 |
| `R002` | 272 | 99 | 9 |
| `R003` | 272 | 99 | 9 |
| `R_public_001` | 274 | 100 | 9 |

说明：当前 `csb-pandora` 旧版建模流程只对 M 链与 peptide 建模；B2M 链来自 PANDORA 选中的真实模板结构，并在 `pandora_structure_meta.json` 中记录 `template_b2m_grafted=true`。
