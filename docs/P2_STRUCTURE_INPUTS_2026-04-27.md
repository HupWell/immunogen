# P2 结构输入准备归档（2026-04-27）

## 1. 执行目的

为 `R001`、`R002`、`R003`、`R_public_001` 准备 peptide-MHC 高精度结构建模输入，完成 P2-1 的结构输入准备。

## 2. 执行命令

```bash
for r in R001 R002 R003 R_public_001; do
  python scripts/prepare_mhc_chain_sequences.py --run_id "$r" --strict
done
```

## 3. 输出文件

每个 run 均已生成：

- `results/<run_id>/structure_inputs/mhc_chain_sequences.fasta`
- `results/<run_id>/structure_inputs/mhc_chain_sequences.json`

## 4. Top peptide-MHC 组合

| Run | Case | Top peptide | Top HLA-I | alpha input | alpha length | B2M length |
|---|---|---|---|---|---:|---:|
| `R001` | `demo_case_001` | `HMTEVVRHC` | `HLA-A*02:01` | `HLA-A*02:01` | 365 | 119 |
| `R002` | `public_case_002` | `ARHGGWTTK` | `HLA-A*02:01` | `HLA-A*02:01` | 365 | 119 |
| `R003` | `public_case_003` | `ARHGGWTTK` | `HLA-A*02:01` | `HLA-A*02:01` | 365 | 119 |
| `R_public_001` | `public_case_001` | `AAAAVIPTV` | `HLA-C*01:02` | `HLA-C*01:02` | 366 | 119 |

## 5. 修正记录

- `scripts/prepare_mhc_chain_sequences.py` 已从“默认取 hla_typing.json 中第一个 HLA-A/B/C”等位基因，修正为优先读取：
  1. `deliveries/<run_id>/to_simhub/<case_id>/selected_for_md.csv`
  2. `results/<run_id>/selected_peptides.csv`
  3. `deliveries/<run_id>/to_immunogen/hla_typing.json`
- 这样可保证结构输入 alpha 链与 Top peptide-HLA 组合一致。
- 本次修正后，`R_public_001` 的 alpha input 已从 `HLA-A*02:01` 更正为 `HLA-C*01:02`。

## 6. 当前结构交付状态

四个 run 当前 `deliveries/<run_id>/to_simhub/<case_id>/meta.json` 仍为：

- `structure_backend`: `coarse`
- `delivery_stage`: `coarse_initial_complex`
- `replaces_coarse`: `false`

因此 P2-1 已完成，但 P2-2/P2-3 仍需继续：使用 PANDORA/AFM 生成真实 `complex.pdb`，替换 coarse 本体并做结构质量复核。
