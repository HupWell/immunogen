# ImmunoGen 最终交付检查单

> 使用方式：每完成一项把 `[ ]` 改为 `[x]`。  
> 默认示例：`run_id=R001`（可替换）。

## 1. 环境与复现

- [ ] `conda activate immunogen` 正常
- [ ] `python scripts/run_all.py --run_id R001` 一键跑通
- [ ] 全链路无关键报错（校验 -> 预测 -> 筛选 -> 构建 -> 报告 -> 自证 -> SimHub -> 可行性）

## 2. 输入完整性

- [ ] `deliveries/R001/to_immunogen/neoantigen_candidates.csv` 存在
- [ ] `deliveries/R001/to_immunogen/hla_typing.json` 存在
- [ ] `deliveries/R001/to_immunogen/meta.json` 存在
- [ ] `neoantigen_candidates.csv` 包含必需列

## 3. 结果产物（results）

- [ ] `results/R001/peptide_mhc_ranking.csv` 存在且非空
- [ ] `results/R001/selected_peptides.csv` 存在且非空
- [ ] `results/R001/mrna_vaccine.fasta` 存在且非空
- [ ] `results/R001/mrna_design.json` 存在
- [ ] `results/R001/REPORT.md` 存在
- [ ] `results/R001/POSITIVE_CONTROL.md` 存在
- [ ] `results/R001/SELF_CHECK.md` 存在

## 4. 质控与图

- [ ] `results/R001/figures/binding_affinity_heatmap.png` 可打开
- [ ] `results/R001/figures/mrna_secondary_structure.png` 可打开
- [ ] `results/R001/qc_metrics.json` 中 `rnafold_status=ok` 或有明确原因说明

## 5. SimHub 交付（分支 B）

- [ ] `deliveries/R001/to_simhub/<case_id>/complex.pdb` 存在
- [ ] `deliveries/R001/to_simhub/<case_id>/meta.json` 中 `molecule_type=peptide_mhc`
- [ ] `deliveries/R001/to_simhub/<case_id>/selected_for_md.csv` 存在
- [ ] 交付目录不包含 `ligand.sdf`
- [ ] 报告中注明当前是否为 coarse 初始结构

## 6. 真实性验收（关键）

- [ ] 若声明真实 MHC-I 交叉验证：`mhc1_cv_source_*` 至少一个为 `real_tsv/real_cmd`
- [ ] 若声明真实 MHC-II：`mhc2_backend=netmhciipan`
- [ ] 执行 `python scripts/check_epitope_realization.py --run_id R001` 并通过目标验收项

## 7. 解释性与边界声明

- [ ] Top 候选有 HLA / 免疫原性 / VAF 三类证据
- [ ] `SELF_CHECK.md` 写明 proxy 与真实工具边界
- [ ] 报告明确“计算验证不替代湿实验”

## 8. 提交前最后确认

- [ ] 文档（README/TODO/SELF_CHECK/REPORT）已同步
- [ ] 不包含环境噪声文件（如 `.autodl/*.db`）
- [ ] `external_refs` 的管理方式（子模块/普通目录）已明确
