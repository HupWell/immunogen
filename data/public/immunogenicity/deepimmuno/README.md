# DeepImmuno Public Data (Imported)

此目录存放从 DeepImmuno 官方仓库公开下载的数据文件，用于本项目公开基准的免疫原性打分对照。

来源（GitHub）：
- https://github.com/frankligy/DeepImmuno
- `reproduce/data/remove0123_sample100.csv`
- `reproduce/data/dengue_test.csv`
- `reproduce/data/ori_test_cells.csv`
- `reproduce/data/sars_cov_2_result.csv`

使用方式：
1. 运行 `python scripts/import_public_deepimmuno_scores.py --run_id R_public_001`
2. 生成 `results/R_public_001/tool_outputs/raw/deepimmuno.tsv`
3. 再运行适配器：
   - `python scripts/run_immunogenicity_adapters.py --run_id R_public_001 --backend_deepimmuno real_tsv --backend_prime proxy --backend_repitope proxy`

注意：
- 本项目只做工程映射与流程验证，不代表这些公开数据可直接替代真实病人实验结论。
- 分数字段在不同公开文件中的定义略有差异，本项目默认按 `potential` / `immunogenic score` / `score` / `ratio` 的优先级读取并做同肽均值。
