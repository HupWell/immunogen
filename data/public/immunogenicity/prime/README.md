# PRIME Public Data (Imported)

此目录用于存放 PRIME 公开数据，供本项目导入为 `raw/prime.tsv`。

已包含：
- `out_compare.txt`（来自 https://github.com/GfellerLab/PRIME 的测试输出）
- `test_input.txt`（对应测试输入）

导入命令：
- `python scripts/import_public_prime_scores.py --run_id R_public_001`

导入逻辑：
1. 先读取 `out_compare.txt` 的 `Peptide` 与 `Score_bestAllele`；
2. 若存在 `supplementary/column_mapping.json`，按映射显式读取指定 CSV 列；
3. 若存在 `supplementary/*.csv`（你后续可手工放 PRIME 补充表），自动尝试识别并追加；
3. 对同一 peptide 的多来源分数取均值，输出到：
   - `results/<run_id>/tool_outputs/raw/prime.tsv`

注意：
- `out_compare.txt` 是 PRIME 测试集，和公开基准 `R_public_001` 重叠可能较低；
- 匹配不足时，流程会由 `auto` 回退 proxy，不会中断。
- 映射模板：`supplementary/column_mapping.template.json`。将其复制为 `column_mapping.json` 并按真实列名修改即可。
