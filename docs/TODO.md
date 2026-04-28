# ImmunoGen 任务书对照 TODO（重写版）

> 负责人：梁心恬  
> 文档目的：写清楚“已经做完什么、还没做什么、下一步先做什么”。  
> 最近一次核对：2026-04-27

---

## 一、已经完成（Done）

### 1) 主流程与交付链路
- [x] 输入契约已落地：`deliveries/<run_id>/to_immunogen/` 的 `neoantigen_candidates.csv`、`hla_typing.json`、`meta.json`
- [x] 一键主流程已打通：预测 -> 排序 -> 筛选 -> 多价 mRNA -> QC/报告 -> SimHub 交付
- [x] 自证材料自动生成：`POSITIVE_CONTROL.md`、`SELF_CHECK.md`、`REPORT.md`
- [x] SimHub 分支 B 交付格式稳定：`complex.pdb` + `hla_allele.txt` + `meta.json` + `selected_for_md.csv`（无 SDF）

### 2) P0 验收项（已全部完成）
- [x] `R_public_001` 具备 MHC-I 真实交叉验证输入，`mhc1_cv_source_netmhcpan` 达到 `real_cmd`
  - 路径已落实：`MHC1_NETMHCPAN_CMD` + `tools/netmhcpan_class1_runner.py`
  - 复核命令：`source scripts/env_netmhcpan.sh && python scripts/run_all.py --run_id R_public_001 --target mhc_ranking --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2`
- [x] `R_public_001` 具备 MHC-II 真实输入（当前为 `real_tsv` 路径，运行后后端为 `netmhciipan`）
- [x] `R_public_001` 严格真实性验收已通过
  - 验证命令：`python scripts/check_epitope_realization.py --run_id R_public_001 --require_mhc1_cv_real --require_mhc2_real`
  - 结果要点：`mhc1_cv_source_netmhcpan=['real_cmd']`、`mhc2_real=True`

### 3) 免疫原性真实后端
- [x] DeepImmuno `real_cmd` 已工程化接入：`tools/deepimmuno_runner.py` + `IMMUNO_DEEPIMMUNO_CMD`
- [x] Repitope `real_cmd` 已工程化接入：`tools/repitope_runner.py` + `IMMUNO_REPITOPE_CMD`
- [x] 免疫原性一键环境脚本已提供：`scripts/env_immunogenicity.sh`
- [x] `R001` / `R002` / `R003` / `R_public_001` 已通过免疫原性严格来源校验（`--require_real_immunogenicity_prime --require_real_immunogenicity_repitope`）
- [x] `R_public_001` 历史 `raw/repitope.tsv` 缺分数字段问题已修复并重建

### 4) 多实例实跑状态
- [x] `R001` 全流程可跑通并完成严格验收
- [x] `R002` 全流程可跑通并完成严格验收
- [x] `R003` 全流程可跑通并完成严格验收
- [x] `R_public_001` 已达成当前阶段严格口径（MHC-I real_cmd + MHC-II real + 免疫原性真实来源）

---

## 二、还没完成（Not Done）

### 1) 数据真实性与临床终版缺口
- [ ] `R_public_001`、`R002`、`R003` 的 MHC-II 当前已完成公开文献口径演示替代，但尚未替换为 BioDriver 上游患者真实 II 类分型
- [ ] 尚未形成“BioDriver 真实 II 类分型 -> 自动重跑 -> 自动更新验收”的临床终版闭环

### 2) 结构交付本体缺口
- [x] `R001/R002/R003/R_public_001` 已使用 PANDORA 真实生成 `complex.pdb`，并完成 M/B/P 三链质量复核

### 3) 任务书增强项未工程化
- [ ] 生产级密码子优化真实接入（LinearDesign / COOL）
- [ ] mRNA 稳定性工具接入（Saluki / RNAsnp）
- [ ] SimHub 回传证据自动归档闭环（`simhub_evidence/<case_id>/` 自动化）

---

## 三、下一步执行清单（按优先级）

### P1（已完成演示替代，临床终版等待上游）

> **2026-04-27 更新**：BioDriver 上游数据暂时无法获取；P1 替代方案 A 已完成，可用于流程展示。临床终版等待 BioDriver 后续提供真实 II 类分型。

#### P1-1: 向 BioDriver 请求 HLA-II 分型数据【阻塞】
- [x] 已创建数据请求文档：`docs/biodriver_hla_ii_request.md`
- [ ] ~~等待 BioDriver 提供三个 run 的真实 HLA-II 分型~~（上游暂时无法提供，任务挂起）
- [ ] ~~收到数据后验证格式并更新 `hla_typing.json`~~

#### P1-替代方案 A：使用公开队列真实 II 类分型（演示级）
若需在当前阶段展示"真实 II 类分型"流程（非终版临床数据），可从公开队列获取真实分型并绑定到 run：

- [x] 方案 A1：选用 TCGA-SKCM 公开研究口径（Dhall et al., 2020）的 HLA-II 类型范围，完成三个 run 的演示绑定
  - 来源：Dhall et al. 2020（PMC7113398）公开说明（II 类包含 DRB1/DQB1/DPB1）
  - 已执行：
    1. 更新 `deliveries/R_public_001/to_immunogen/hla_typing.json`
    2. 更新 `deliveries/R002/to_immunogen/hla_typing.json`
    3. 更新 `deliveries/R003/to_immunogen/hla_typing.json`
    4. 重跑严格链路并归档日志到 `results/<run>/validation_logs/`
- [x] 方案 A2：明确声明"演示数据替代"口径
  - 已在三个 run 的 `SELF_CHECK.md` 中注明："公开队列数据，非患者真实临床终版分型"

#### P1-2: 重建 MHC-II 预测并重新验收【已完成（演示口径）】
- [x] 使用替代方案 A 的 II 类分型完成 `mhc2_netmhciipan.tsv` 重建并重跑
- [x] 三个 run 重新执行严格链路并归档验收结果
  - 归档文档：`docs/P1_ALTERNATIVE_A_VALIDATION_2026-04-27.md`
  - 执行命令（预备）：
    ```bash
    for r in R_public_001 R002 R003; do
      source scripts/env_netmhcpan.sh
      python scripts/run_all.py --run_id $r --target full \
        --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off \
        --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2 \
        --require_real_immunogenicity_prime --require_real_immunogenicity_repitope
    done
    ```
  - 验收命令（预备）：
    ```bash
    for r in R_public_001 R002 R003; do
      python scripts/check_epitope_realization.py --run_id $r \
        --require_mhc1_cv_real --require_mhc2_real
    done
    ```

#### P1-3: 更新自证文档【已完成（演示口径）】
- [x] 更新 `R_public_001/SELF_CHECK.md`（添加 II 类来源声明）
- [x] 更新 `R002/SELF_CHECK.md`
- [x] 更新 `R003/SELF_CHECK.md`
- [x] 更新本 TODO.md（记录 P1 替代方案 A 执行结果）

---

**当前结论**：
- P1 替代方案 A 已完成（可用于流程展示与验收归档）
- 临床终版仍需 BioDriver 真实 II 类分型，届时需再执行一次全量重跑

### P2（当前优先推进）：结构本体替换与质量复核

#### P2-1: 结构输入准备
- [x] 为 `R001/R002/R003/R_public_001` 逐例确认 Top peptide-MHC 组合（来自 `selected_for_md.csv`）
- [x] 生成/复核结构建模输入序列（MHC-I alpha、B2M、peptide）
- [x] 明确每个 case 的建模后端：优先 PANDORA 批量，AFM 重点复核
  - 归档文档：`docs/P2_STRUCTURE_INPUTS_2026-04-27.md`
  - 修正记录：`prepare_mhc_chain_sequences.py` 已改为优先使用 Top HLA-I 等位基因生成 alpha 链输入

#### P2-2: 高精度结构本体替换
- [x] 安装结构建模环境（只安装真实工具，不生成伪造结构）
  - 归档文档：`docs/P2_STRUCTURE_ENV_SETUP_2026-04-28.md`
  - 已可用：`pandora_src` 环境、`MUSCLE 5.2`、系统 `blastp/makeblastdb 2.9.0+`
  - 已可用：`colabfold` 环境与 `colabfold_batch` 命令入口
  - 已配置：MODELLER 真实 license key（不记录明文）
  - 已修复：`csb-pandora 0.9` 的 Biopython 兼容问题与 IMGT HTTPS 下载问题
  - 已完成：下载 Zenodo 官方 PANDORA 数据库 `default.tar.gz`，并接入 864 个 MHC-I / 136 个 MHC-II 模板
  - 待增强：ColabFold GPU 版 JAX（当前 `jax 0.6.2` 仅识别 CPU）
- [x] 用 PANDORA 生成 `R001/R002/R003/R_public_001` 真实 peptide-MHC 复合物
- [x] 替换 `deliveries/<run_id>/to_simhub/<case_id>/complex.pdb`
- [x] 更新 `meta.json`：`structure_backend=pandora`、`structure_tool_version`、`replaces_coarse=true`

#### P2-3: 结构质量复核
- [x] 检查 `complex.pdb` 是否包含 Chain M / B / P
- [x] 检查 peptide 是否与 `selected_for_md.csv` 对齐
- [x] 检查是否仍有 `coarse_initial_complex` / `na_for_coarse` / `replaces_coarse=false`
- [x] 归档结构复核记录：`docs/P2_PANDORA_STRUCTURE_VALIDATION_2026-04-28.md`

#### P2-4: SimHub 交付刷新
- [x] 重新运行 `prepare_simhub_delivery.py`（非 coarse 后端）
- [x] 重新检查 `deliveries/<run_id>/to_simhub/<case_id>/` 文件完整性
- [x] 更新 `SELF_CHECK.md` 中结构后端说明

### P3（第二阶段增强）
- [ ] 接入 LinearDesign/COOL 并纳入主流程参数
- [ ] 接入 Saluki/RNAsnp 并产出可复检指标
- [ ] 打通 SimHub 回传证据自动归档

---

## 四、执行约束（防回退）

- [ ] 每次 `run_all` 后同步更新本文件（仅保留当前有效状态，不累计历史噪声）
- [ ] 每次提交前至少执行一次：`python scripts/check_epitope_realization.py --run_id <run_id> --require_mhc1_cv_real --require_mhc2_real`
- [ ] 提交时剔除噪声变更：`.autodl/*`、`external_refs` 子仓库状态
