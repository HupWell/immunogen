# ImmunoGen 任务书对照 TODO（2026-04-23）

> 说明：按“是否已在仓库真实落地”划分，不按“是否仅有接口”划分。

---

## 一、已完成

### 1) 输入契约
- [x] BioDriver 输入三件套校验（`neoantigen_candidates.csv` / `hla_typing.json` / `meta.json`）
- [x] HLA-I 必填与 HLA-II 可选字段说明（见 `docs/hla_typing.md`）

### 2) 预测与排序主流程
- [x] MHC-I 主流程预测（MHCflurry）
- [x] 综合打分 `rank_score = w1*HLA + w2*immunogenicity + w3*VAF + w4*dissimilarity`
- [x] 与 WT 过相似过滤（`select_top_peptides.py`）

### 3) 工程化交付链路
- [x] 多价 ORF + mRNA 输出（FASTA + 设计 JSON）
- [x] QC 图与报告输出
- [x] 自证最小包自动生成（`POSITIVE_CONTROL.md` / `SELF_CHECK.md` / `REPORT.md`）
- [x] SimHub 分支 B 交付目录（`complex.pdb` / `meta.json` / `selected_for_md.csv`）
- [x] `run_all.py` 一键串联执行

### 4) 验收保护
- [x] 强制实跑开关：`--require_real_mhc2`、`--require_real_mhc1_cv`
- [x] 验收脚本：`scripts/check_epitope_realization.py`
- [x] MHC-II 新增 `real_tsv` 后端并接入 `auto` 回退链（`predict_mhc_ranking.py` / `run_all.py`）
- [x] MHC-I 新增 `openvax_bridge.py`，支持把外部报告桥接为 `mhc1_netmhcpan.tsv`

---

## 二、部分完成（有框架，真实能力未全面落地）

### P0（优先级最高）
- [ ] MHC-I 交叉验证真实来源常态化（当前很多 run 仍可能是 `off`）
  - 目标：`mhc1_cv_source_netmhcpan` / `mhc1_cv_source_bigmhc` 出现 `real_tsv/real_cmd`
  - 现状：代码层已支持 `real_tsv/real_cmd` + 报告桥接，主要缺口在真实结果文件的持续供给与数据运维。
- [ ] MHC-II 真实 NetMHCIIpan 常态化（当前常见 `proxy`）
  - 目标：`mhc2_backend=netmhciipan` 在真实病例稳定出现
  - 现状：代码层已支持 `netmhciipan` 与 `real_tsv`，主要缺口在 II 类分型与真实结果覆盖率。
- [ ] 免疫原性真实模型常态化
  - 目标：`immunogenicity_source_*` 从 `proxy/precomputed` 逐步提升到 `real_tsv/real_cmd`

### P1
- [ ] 真实结构后端批量化（PANDORA/AFM）
  - 目标：关键交付 case 的 `complex.pdb` 不再使用 coarse 占位

---

## 三、未完成（任务书明确但尚未接入）

### P2
- [ ] 生产级密码子优化工具（LinearDesign/COOL）真实接入与回归验证
- [ ] mRNA 稳定性工具（Saluki / RNAsnp）接入 `qc_metrics.json` 与报告
- [ ] SimHub 回传证据自动归档流程标准化（`simhub_evidence/<case_id>/` 全字段闭环）

---

## 四、当前实际缺口（你现在最该补的）

1. **文档恢复与持续维护**
   - [x] `README.md` 已恢复
   - [x] 本 TODO 已恢复
   - [ ] `FINAL_CHECKLIST.md` 需随每轮 run 更新勾选状态

2. **真实工具接入优先顺序**
   - [ ] 先打通 MHC-I 交叉验证 real_tsv（最容易落地）
   - [ ] 再打通 NetMHCIIpan real
   - [ ] 最后替换免疫原性 proxy 为真实推理

3. **提交质量风险**
   - [ ] 清理 `.autodl/` 这类环境噪声文件，避免污染主分支
   - [ ] 明确 `external_refs/neoantigen-vaccine-pipeline` 是子模块还是普通目录

---

## 五、建议你马上执行的 3 条命令

```bash
python scripts/run_all.py --run_id R001 --target mhc_ranking --backend_mhc1_netmhcpan real_tsv --backend_mhc1_bigmhc off --require_real_mhc1_cv
python scripts/check_epitope_realization.py --run_id R001 --require_mhc1_cv_real
python scripts/run_all.py --run_id R001 --target full
```

---

最后更新：2026-04-23
