# BioDriver HLA Class II 分型数据请求文档

> 请求日期：2026-04-27  
> 请求方：ImmunoGen（梁心恬）  
> 用途：完成 R_public_001 / R002 / R003 的 MHC-II 真实分型替换

---

## 1. 背景说明

当前 ImmunoGen 模块中三个关键 run（R_public_001、R002、R003）的 MHC-II 预测仍使用参考等位基因（如 `DRB1*15:01`）进行回填演示，尚未替换为患者真实 II 类分型。

为达成"真实 II 类分型 -> 自动重跑 -> 自动更新验收"的标准化闭环，需要 BioDriver 上游提供这三个 run 对应患者的 HLA Class II 分型数据。

---

## 2. 请求数据清单

### 2.1 目标 Run 列表

| Run ID | Case ID | 当前状态 | 需补充字段 |
|--------|---------|----------|-----------|
| R_public_001 | public_case_001 | 仅含 I 类分型 | HLA-DRB1、HLA-DQB1、HLA-DPB1（至少 DRB1） |
| R002 | public_case_002 | 仅含 I 类分型 | HLA-DRB1、HLA-DQB1、HLA-DPB1（至少 DRB1） |
| R003 | public_case_003 | 仅含 I 类分型 | HLA-DRB1、HLA-DQB1、HLA-DPB1（至少 DRB1） |

### 2.2 数据格式要求

**文件位置**：`deliveries/<run_id>/to_immunogen/hla_typing.json`

**格式模板**（现有 I 类分型基础上增加 II 类）：

```json
{
  "HLA-A": ["A*02:01", "A*11:01"],
  "HLA-B": ["B*40:01", "B*46:01"],
  "HLA-C": ["C*01:02", "C*07:02"],
  "HLA-DRB1": ["DRB1*15:01", "DRB1*04:05"],
  "HLA-DQB1": ["DQB1*06:02", "DQB1*03:01"],
  "HLA-DPB1": ["DPB1*04:01", "DPB1*02:01"]
}
```

**分辨率要求**：
- 至少 2-field（如 `DRB1*15:01`）
- 优选 4-field（如 `DRB1*15:01:01:01`）

**来源优先级**（按可靠性排序）：
1. **临床 HLA 检测报告**（NGS/SSO/SSP）
2. **高质量 WES/WGS/RNA-seq 推断**（如 arcasHLA、HLA-HD、xHLA 等工具，需附 QC 指标）
3. **公开队列官方注释**（如 TCGA 样本的已知分型）

---

## 3. 下游使用流程

收到数据后，ImmunoGen 将执行以下步骤：

1. **更新 `hla_typing.json`**：在现有 I 类分型基础上增加 II 类分型字段
2. **重建 MHC-II 预测表**：使用真实 II 类分型重新运行 NetMHCIIpan，生成新的 `mhc2_netmhciipan.tsv`
3. **重跑严格链路**：
   ```bash
   python scripts/run_all.py --run_id <run_id> --target full \
     --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off \
     --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2 \
     --require_real_immunogenicity_prime --require_real_immunogenicity_repitope
   ```
4. **真实性验收**：
   ```bash
   python scripts/check_epitope_realization.py --run_id <run_id> \
     --require_mhc1_cv_real --require_mhc2_real
   ```
5. **更新交付文档**：归档验收结果，更新 `SELF_CHECK.md`

---

## 4. 验收标准

- BioDriver 提供的 HLA-II 分型需与患者临床记录或测序数据一致
- 若来源为工具推断，需明确工具名称（如 HLA-HD v1.4.0）、版本及 QC 指标
- ImmunoGen 将验证 II 类分型格式是否符合 NetMHCIIpan 输入要求

---

## 5. 时间计划（建议）

| 阶段 | 预计时间 | 交付物 |
|------|----------|--------|
| BioDriver 提供 HLA-II 分型 | T+3 工作日 | 3 个更新的 `hla_typing.json` |
| ImmunoGen 重建 MHC-II 表并重跑 | T+5 工作日 | 新的 `mhc2_netmhciipan.tsv` + 验收报告 |
| 更新 SELF_CHECK.md 并归档 | T+6 工作日 | 完整的 P1 完成报告 |

---

## 6. 联系方式

如有问题，请联系 ImmunoGen 负责人：梁心恬

