# ImmunoGen 最终交付检查单

> 使用方法：每完成一项就把 `[ ]` 改成 `[x]`。

## 1. 环境与复现

- [ ] `conda activate immunogen` 可正常激活环境
- [ ] `python scripts/run_all.py --run_id R_public_001` 可一键跑通
- [ ] 无关键报错（输入校验、预测、构建、报告、交付、可行性都成功）

## 2. 输入数据完整性

- [ ] `deliveries/R_public_001/to_immunogen/neoantigen_candidates.csv` 存在
- [ ] `deliveries/R_public_001/to_immunogen/hla_typing.json` 存在
- [ ] `deliveries/R_public_001/to_immunogen/meta.json` 存在
- [ ] `neoantigen_candidates.csv` 包含关键列：`mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf`

## 3. ImmunoGen 主结果（3.1）

- [ ] `results/R_public_001/peptide_mhc_ranking.csv` 存在且非空
- [ ] `results/R_public_001/selected_peptides.csv` 存在且非空
- [ ] `results/R_public_001/mrna_vaccine.fasta` 存在且非空
- [ ] `results/R_public_001/mrna_design.json` 存在且含 `rnafold_status`
- [ ] `results/R_public_001/REPORT.md` 存在
- [ ] `results/R_public_001/FEASIBILITY.md` 存在

## 4. 图与质控

- [ ] `results/R_public_001/figures/binding_affinity_heatmap.png` 可正常打开
- [ ] `results/R_public_001/figures/mrna_secondary_structure.png` 可正常打开
- [ ] `results/R_public_001/qc_metrics.json` 中 `rnafold_status` 为 `ok`（或说明回退原因）
- [ ] `results/R_public_001/qc_metrics.json` 中有合理的 `gc_percent` 与 `sequence_length`

## 5. Simulation Hub 交付（3.2）

- [ ] `deliveries/R_public_001/to_simhub/public_case_001/protein.pdb` 存在
- [ ] `deliveries/R_public_001/to_simhub/public_case_001/ligand.sdf` 存在
- [ ] `deliveries/R_public_001/to_simhub/public_case_001/meta.json` 存在
- [ ] `deliveries/R_public_001/to_simhub/public_case_001/selected_for_md.csv` 存在
- [ ] 明确说明 `protein.pdb` 当前为粗初始构象（非最终高精度结构）

## 6. 可行性与对照验证

- [ ] `results/R_public_001/FEASIBILITY.md` 中 `engineering_reproducible` 为 `True`
- [ ] `results/R_public_001/FEASIBILITY.md` 中 `directionally_feasible` 为 `True`
- [ ] `results/R_public_001/FEASIBILITY.md` 中 `KRAS_G12D_like` 命中为 `True`
- [ ] `results/R_public_001/selected_peptides.csv` 中可见 TP53 热点候选（如 `TP53_R175H` / `TP53_R248Q`）

## 7. 报告边界说明（必须）

- [ ] 报告中写清数据来源：`pVACtools 示例 + 文献补充`
- [ ] 报告中写清 MHC-II/免疫原性当前为代理评分（若尚未替换真实模型）
- [ ] 报告中写清当前结果用于流程与方向验证，不替代湿实验结论
- [ ] 报告中写清后续计划：真实 MHC-II、真实免疫原性模型、真实 peptide-MHC 精细结构建模

## 8. 最终交付物打包前确认

- [ ] 脚本目录 `scripts/` 保持可运行状态（无临时调试残留）
- [ ] `README.md` 命令与当前代码一致
- [ ] 所有路径使用 UTF-8 正常显示，无中文乱码
- [ ] 交付前再次执行一次 `python scripts/run_all.py --run_id R_public_001` 复核

---

## 一句话结论（交付时可复制）

本次 ImmunoGen 模块已完成从候选新抗原输入到多价 mRNA 序列输出及 Simulation Hub 交付包生成的端到端流程验证；在公开验证集上具备可复现性和方向性可行性，后续将接入真实 MHC-II/免疫原性模型与高精度结构建模以提升生物学可信度。
