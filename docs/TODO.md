# ImmunoGen TODO（当前状态版）

> 负责人：梁心恬  
> 文档目的：只保留当前有效状态、剩余阻塞和下一步动作。  
> 最近一次核对：2026-04-30  
> 当前结论：`R001/R002/R003/R_public_001` 已完成阶段版严格重跑；机器可读结果未发现 `proxy` / `fallback` / `coarse` 残留。

---

## 一、当前可交付结果

以下内容可以作为**阶段验收版/下游技术检查版**交付：

- [x] `deliveries/<run_id>/to_simhub/<case_id>/`：SimHub 交付包，包含 PANDORA 真实结构 `complex.pdb`
- [x] `results/<run_id>/mrna_vaccine.fasta`：LinearDesign `real_cmd` 密码子优化后的 mRNA 设计
- [x] `results/<run_id>/mrna_design.json`：mRNA 设计元数据，含工具、命令、输入输出路径
- [x] `results/<run_id>/qc_metrics.json`：RNAfold / RNAeval / RNAplfold 稳定性指标
- [x] `results/<run_id>/REPORT.md`、`SELF_CHECK.md`、`POSITIVE_CONTROL.md`：报告与自证材料
- [x] `deliveries/<run_id>/to_simhub/<case_id>/dossier_context.json`：peptide + mRNA dossier manifest
- [x] `results/<run_id>/meta.json`：本 run 自追溯摘要

适用实例：

- [x] `R001`
- [x] `R002`
- [x] `R003`
- [x] `R_public_001`

注意：`R002/R003/R_public_001` 的 HLA-II 仍是公开队列演示分型，不是 BioDriver 患者临床终版分型。

---

## 二、已完成

### 1) 主流程

- [x] 输入契约已落地：`deliveries/<run_id>/to_immunogen/neoantigen_candidates.csv`、`hla_typing.json`、`meta.json`
- [x] 主流程已打通：输入校验 -> 免疫原性适配 -> MHC 排名 -> Top 肽筛选 -> mRNA 构建 -> QC/报告 -> SimHub 交付
- [x] 自证材料自动生成：`POSITIVE_CONTROL.md`、`SELF_CHECK.md`、`REPORT.md`
- [x] SimHub 分支 B 交付格式稳定：`complex.pdb` + `hla_allele.txt` + `meta.json` + `selected_for_md.csv`

### 2) 真实后端与防回退

- [x] MHC-I 交叉验证：NetMHCpan `real_tsv`
- [x] MHC-II：NetMHCIIpan `real_tsv`，运行后结果列为 `netmhciipan`
- [x] DeepImmuno：真实来源，当前结果列为 `deepimmuno_real_cmd`
- [x] PRIME：真实来源，当前结果列为 `prime_real_cmd`
- [x] Repitope：公开数据集真实来源，当前结果列为 `repitope_public_dataset_exact/knn`
- [x] BigMHC：显式 `off`，不作为代理补位
- [x] 默认禁止 proxy / fallback 静默回退

### 3) 结构真实后端

- [x] PANDORA 官方模板数据库已接入：Zenodo `default.tar.gz`
- [x] `scripts/run_pandora_structure.py` 已落地：真实 PANDORA + MODELLER 生成 peptide-MHC-I 结构
- [x] `R001/R002/R003/R_public_001` 均已生成真实 PANDORA `complex.pdb`
- [x] SimHub `meta.json` 已改为 SINGLE SOURCE OF TRUTH 白名单字段，不再写入未声明的结构调试字段
- [x] 已完成结构复核：Chain M / B / P、peptide 对齐、无 coarse 标记、非 CA-only
- [x] `prepare_simhub_delivery.py` 已增加生产交付硬校验：coarse 标记、CA-only 比例过高、原子数过少或链数不足均会拒绝

### 4) P3-1 生产级密码子优化

- [x] 选定官方 LinearDesign 本机命令行接入
- [x] `build_multivalent_mrna.py --codon_mode real_cmd` 已支持真实命令适配器
- [x] 已记录工具名、版本命令、真实命令、输入 FASTA、输出 FASTA、stdout/stderr 日志路径
- [x] 四个实例已完成 LinearDesign 真实验收：`codon_mode=real_cmd`、`codon_optimizer.tool=LinearDesign`

### 5) P3-2 mRNA 稳定性真实工具链

- [x] 选定当前本机可复检口径：ViennaRNA `RNAfold` + `RNAeval` + `RNAplfold`
- [x] 已写入 `qc_metrics.json` / `mrna_design.json` / `REPORT.md`
- [x] 已输出可复检字段：工具版本、命令、输入 RNA/FASTA、关键分数、日志路径
- [x] 四个实例已完成真实稳定性验收
  - 指标：`rnafold_mfe`、`rnaeval_energy`、`mean_unpaired_l1`、`mean_unpaired_l10`
  - 输出目录：`results/<run_id>/mrna_stability/`

### 6) 最近一次全量严格复核

- [x] 已重新运行四个实例完整链路
- [x] 已补齐任务书 3.4 要求的 `dossier_context.json`
- [x] 已补齐 `results/<run_id>/meta.json` 与 `simhub_evidence/<case_id>/README.md` 回传占位契约
- [x] 已按 SimHub SINGLE SOURCE OF TRUTH 修复 peptide-MHC `meta.json`：补齐共享必填字段，并移除未声明字段，避免 `E_META_INCOMPLETE`
- [x] 已执行严格验收：
  ```bash
  for r in R001 R002 R003 R_public_001; do
    python scripts/check_epitope_realization.py --run_id "$r" \
      --require_mhc1_cv_real \
      --require_mhc2_real \
      --require_real_immunogenicity \
      --require_real_structure
  done
  ```
- [x] 已字段级扫描关键交付文件，未发现实际交付包中的 coarse 文件头或联通测试说明

---

## 三、仍未完成

### 1) 临床终版数据闭环【上游阻塞】

- [ ] 等待 BioDriver 提供患者真实 HLA-II 分型
- [ ] 将 `R002/R003/R_public_001` 的公开队列演示分型替换为 BioDriver 患者真实 II 类分型
- [ ] 替换后自动重跑：MHC-II -> 排名 -> 筛选 -> mRNA -> QC -> SimHub -> 自证
- [ ] 形成“BioDriver 真实 II 类分型 -> 自动重跑 -> 自动更新验收”的临床终版闭环

### 2) P3-3 SimHub 回传证据自动归档

- [x] 定义归档契约：`results/<run_id>/simhub_evidence/<case_id>/`
- [x] 当前未回传状态已用 `README.md` 标记为 `not_returned`
- [x] 增加回传文件完整性检查：轨迹、能量、结构稳定性摘要、日志
- [x] 在 `SELF_CHECK.md` 和独立 `SIMHUB_EVIDENCE.md` 中自动引用 SimHub 回传证据
- [x] 标记下游检查状态：`not_returned` / `returned_unvalidated` / `validation_passed` / `validation_failed`

### 3) P3-4 ColabFold GPU 环境增强【非主链路阻塞】

- [x] 修复 GPU 版 JAX 配置，使 `colabfold_batch` 能识别 GPU
- [x] 用 AFM / ColabFold 对重点候选 `R_public_001/public_case_001` 做结构复核
- [x] 将 AFM 复核 PDB 与 PANDORA PDB 的差异记录进验证报告：`docs/P3_4_AFM_REVIEW_R_PUBLIC_001_2026-04-30.md`
- [x] 已记录当前验证状态：`docs/P3_4_COLABFOLD_GPU_VALIDATION_2026-04-29.md`

### 4) 可选增强

- [ ] Saluki / RNAsnp 后续可作为额外稳定性工具接入；当前 ViennaRNA 增强口径已满足本机真实稳定性验收
- [ ] BigMHC 可作为额外 MHC-I 交叉验证增强；当前 NetMHCpan `real_tsv` 已满足必需真实交叉验证
- [ ] 整理 Git 工作区中的本机噪声和外部工具目录，避免误提交 `.autodl/`、`external_refs/`、下载缓存等

---

## 四、下一步优先级

1. **等待 SimHub 回传 MD 结果**
   输入包已交付，后续需要接收轨迹、能量、RMSD、QC flags 和 summary。

2. **等待 BioDriver HLA-II 临床终版数据**  
   数据到位后做全量重跑和临床终版自证。

3. **可选：提升 AFM 复核质量**  
   当前 AFM 已跑通但置信度低；如需替换结构，建议增加 MSA、recycle、model/seed 后重新筛选。

---

## 五、防回退检查命令

严格验收：

```bash
for r in R001 R002 R003 R_public_001; do
  python scripts/check_epitope_realization.py --run_id "$r" \
    --require_mhc1_cv_real \
    --require_mhc2_real \
    --require_real_immunogenicity \
    --require_real_structure
done
```

关键结果字段扫描：

```bash
rg "proxy|fallback_profile|coarse_initial_complex|na_for_coarse|structure_backend\"\\s*:\\s*\"coarse\"" \
  results/{R001,R002,R003,R_public_001} \
  deliveries/{R001,R002,R003,R_public_001}/to_simhub
rg "COARSE PEPTIDE-MHC COMPLEX|coarse CA trace|Replace with AlphaFold-Multimer|仅用于流程联通|生产环境请用 AlphaFold-Multimer" \
  deliveries/{R001,R002,R003,R_public_001}/to_simhub
```

提交前注意：

- [ ] 不提交 `.autodl/*`
- [ ] 不提交未授权再分发的 `external_refs/LinearDesign/`
- [ ] 不提交下载缓存和大型外部数据库
- [ ] 若提交代码，先检查 `git diff --cached`
