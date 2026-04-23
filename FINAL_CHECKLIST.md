# ImmunoGen 最终交付检查单（对齐更新任务书）

> 使用方法：每完成一项就把 `[ ]` 改成 `[x]`。  
> 示例 `run_id`：`R_public_001`；请按需替换。

## 1. 环境与复现

- [ ] `conda activate immunogen` 可正常激活环境（需在当前远程会话手动确认）
- [ ] `python scripts/run_all.py --run_id R_public_001` 可一键跑通（建议本次会话复跑）
- [ ] 无关键报错（校验 → 预测 → 筛选 → 构建 → 报告 → **自证** → SimHub → 可行性）（需复跑确认）

## 2. 输入数据完整性

- [x] `deliveries/R_public_001/to_immunogen/neoantigen_candidates.csv` 存在
- [x] `deliveries/R_public_001/to_immunogen/hla_typing.json` 存在
- [x] `deliveries/R_public_001/to_immunogen/meta.json` 存在（含 `case_id` 更佳）
- [x] `neoantigen_candidates.csv` 含列：`mutation, mut_peptide, wt_peptide, transcript_id, variant_vaf`

## 3. ImmunoGen 主结果 + 自证最小包（results）

- [x] `results/R_public_001/peptide_mhc_ranking.csv` 存在且非空
- [x] `results/R_public_001/selected_peptides.csv` 存在且非空
- [x] `results/R_public_001/mrna_vaccine.fasta` 存在且非空
- [x] `results/R_public_001/mrna_design.json` 存在且含 `rnafold_status`
- [x] **`results/R_public_001/REPORT.md` 存在**
- [x] **`results/R_public_001/POSITIVE_CONTROL.md` 存在**
- [x] **`results/R_public_001/SELF_CHECK.md` 存在**
- [x] `results/R_public_001/FEASIBILITY.md` 存在（可选但推荐）

## 4. 图与质控

- [x] `results/R_public_001/figures/binding_affinity_heatmap.png` 可正常打开
- [x] `results/R_public_001/figures/mrna_secondary_structure.png` 可正常打开（标题含 MFE）
- [x] `results/R_public_001/qc_metrics.json` 中 `rnafold_status` 为 `ok`（或 `SELF_CHECK` 中已说明原因）

## 5. Simulation Hub 交付（契约分支 B：无 SDF）

- [x] **`deliveries/R_public_001/to_simhub/<case_id>/complex.pdb` 存在**
- [x] **`deliveries/R_public_001/to_simhub/<case_id>/hla_allele.txt` 存在**（可选但本脚本会生成）
- [x] **`deliveries/R_public_001/to_simhub/<case_id>/meta.json` 中 `molecule_type` 为 `peptide_mhc`**
- [x] `deliveries/R_public_001/to_simhub/<case_id>/selected_for_md.csv` 存在
- [x] **确认未包含 `ligand.sdf`**（多肽禁止走 SDF）
- [x] 汇报中说明：`complex.pdb` 当前为粗初始构象，送 MD 前需换 AlphaFold-Multimer / PANDORA 高精度结构

## 6. 可行性与 Positive Control

- [x] `FEASIBILITY.md` 中工程/方向性结论符合预期（若使用公开基准）
- [x] `POSITIVE_CONTROL.md` 中 KRAS G12D 对照说明与人工复核步骤可读
- [x] `selected_peptides.csv` 中可见文献热点（如 TP53）或等价阳性证据

## 7. 边界说明（汇报必提）

- [x] `SELF_CHECK.md` 已写明：MHC-II / 免疫原性代理、与 NetMHCpan 等差异
- [x] 已写明：计算验证 ≠ 湿实验结论

## 8. 提交前再跑一遍

- [ ] `python scripts/run_all.py --run_id R_public_001` 复核通过（建议本次远程环境执行后再勾选）

---

## 一句话结论（交付时可复制）

ImmunoGen 已完成从 BioDriver 风格输入到多价 mRNA 与 SimHub（peptide_mhc、无 SDF）交付目录的端到端流程；`results/` 下具备 **POSITIVE_CONTROL.md、SELF_CHECK.md、REPORT.md** 自证最小包，公开基准可复现，后续将替换代理分为真实工具链并升级 `complex.pdb` 为高精度初始构象。
