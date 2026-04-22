# hla_typing.json 输入契约（ImmunoGen / BioDriver）

本文档与仓库内 `deliveries/<run_id>/to_immunogen/hla_typing.json` 对齐，供上游 BioDriver 写出稳定字段、下游 NetMHCIIpan 等工具对接。**MHC-I 为必需；MHC-II 为可选**，便于仅完成 WES I 类分型的历史管线继续运行。

---

## 1. 顶层结构

- 根节点为 **JSON 对象**。
- 所有 HLA 键对应的值均为 **JSON 字符串数组**；单个等位基因占数组中一项。
- 等位基因字符串建议使用 **HLA 官方命名**（多数字段如 `*02:01`，长度基因可用 `*01:01` 等；与上游工具版本一致即可）。

### 1.1 必备键（MHC-I）

| 键 | 说明 |
| --- | --- |
| `HLA-A` | A 座两条（或一条，依生物学与上游策略），最常见形式 `A*02:01` |
| `HLA-B` | B 座 |
| `HLA-C` | C 座 |

### 1.2 可选键（MHC-II，供 NetMHCIIpan / 辅助 T 细胞表位等）

未提供时，ImmunoGen 仍只做 **MHC-I 亲和力**（MHCflurry）与 **MHC-II 代理分**；提供后便于后续脚本直接消费真实 II 类分型，无需再改字段名。

| 键 | 说明 |
| --- | --- |
| `HLA-DRB1` | DRβ 链，常见如 `DRB1*15:01` |
| `HLA-DQA1` | DQα |
| `HLA-DQB1` | DQβ（与 DQA1 配对信息由生物学家或上游表维护；JSON 中仅平铺等位基因列表） |
| `HLA-DPA1` | DPα |
| `HLA-DPB1` | DPβ |

**约定：** 可省略某基因键，或显式给空数组 `[]` 表示「本 patient 经上游判定该座缺失/未 call」；若键存在，必须为数组（不可为 `null` 或字符串）。

### 1.3 与示例文件

- 仅 I 类最小示例：各 `deliveries/R*/to_immunogen/hla_typing.json` 现有内容。
- **含 II 类完整示例**（供对照）：`data/examples/hla_typing.class_ii.example.json`

---

## 2. 上游工具 → 本契约字段映射

### 2.1 OptiType（经典仅 MHC-I）

- OptiType 标准输出为 **A/B/C 三座的推断分型**，**不包含** MHC-II。
- 常见结果列为等位基因字符串（可能带 `HLA-` 前缀，如 `HLA-B*15:25` 或 `B*15:25`）。
- **映射到本 JSON：**
  - OptiType A → `HLA-A` 列表（一条或按输出拆成两条）；
  - 同理 B → `HLA-B`，C → `HLA-C`。
- **II 类**：需由 BioDriver 使用 **HLA-HD**、HLA*LA、seq2HLA 等**单独**流程写入本文档第二节可选键，**不可**从 OptiType 直接得到。

### 2.2 HLA-HD（可含 I + II）

- HLA-HD 报告常见包含 **A/B/C/DRB1/DQA1/DQB1/DPA1/DPB1** 等（具体以所用版本与 `report`/`result` 实际格式为准）。
- **映射到本 JSON：**
  - I 类基因名 `A` / `B` / `C` 对应 → `HLA-A` / `HLA-B` / `HLA-C`（值为去除多余空格后的等位基因名，与 IMGT/HLA 一致）；
  - `DRB1` → `HLA-DRB1`；
  - `DQA1` → `HLA-DQA1`；`DQB1` → `HLA-DQB1`；
  - `DPA1` → `HLA-DPA1`；`DPB1` → `HLA-DPB1`。
- 若 HLA-HD 输出中 **DQA1 与 DQB1 成对**出现在同一行，BioDriver 拆分为两条等位基因分别填入两个数组，并在自身 `meta` 中保留成对关系（本 JSON 不强制嵌套结构）。

### 2.3 命名与 NetMHCIIpan（后续步骤）

- NetMHCIIpan 对 **等位基因字符串格式** 有自身要求（如 `HLA-DRB1*15:01` 与 `DRB1_0101` 等），与 II 类工具对接时，请在 BioDriver 或 ImmunoGen 侧维护 **「本契约」→「工具入参」** 映射表；见 `docs/TODO.md` 中 P0 下一项（等位基因映射表）。

---

## 3. 校验行为（validate_input.py）

- 若缺少 `HLA-A` / `HLA-B` / `HLA-C` 任一键，或任一键非数组、或 I 类合计无有效等位基因 → **校验失败**。
- 若存在 `HLA-DRB1` 等可选键，则必须为**字符串数组**；若 II 类合计条数为 0，仅**提示**（不失败，便于 I-only 数据）。
- 根对象上**未知键**将打印警告并忽略，避免锁死未来扩展；生产环境建议与 BioDriver 约定白名单。

---

*文档版本：与 ImmunoGen 任务书 MHC-II / NetMHCIIpan 准备阶段一致； II 类键为新增契约。*
