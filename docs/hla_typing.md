# hla_typing.json 说明（ImmunoGen / BioDriver）

**没学过 HLA、免疫？** 先看更白话的：[`allele_naming_simple.md`](./allele_naming_simple.md)（用「手机格式统一」做比喻）。

本文只回答两件事：**文件长什么样**、**上游软件的结果怎么填进这些格子里**。  
**I 类（三节）是必填**；**II 类是选填**（没有 II 时老流程照跑，只是还做不了真实 NetMHCIIpan 对接）。

---

## 1. 文件结构（就一层 JSON）

- 最外层是对象 `{}`。
- 每个键下都是 **字符串列表** `[]`，**一条等位基因占列表中一行**（一个字符串）。

### 1.1 必须有的键（I 类，MHC-I）

| 键 | 你把它当成 |
| --- | --- |
| `HLA-A` | A 座，常见写法里会有 `*02:01` 这种带星号的号 |
| `HLA-B` | B 座 |
| `HLA-C` | C 座 |

### 1.2 可以没有的键（II 类，MHC-II，给 NetMHCIIpan 等用）

| 键 | 说明（知道名字即可，不必记免疫机制） |
| --- | --- |
| `HLA-DRB1` | 常写成 `DRB1*15:01` 或 `HLA-DRB1*15:01` 这类 |
| `HLA-DQA1` / `HLA-DQB1` / `HLA-DPA1` / `HLA-DPB1` | 同上，**各是一个列表**；谁是谁的「配对」由 BioDriver 在别处说明，本 JSON 只平铺列表 |

**可省略**某个键，或写空列表 `[]`，表示没有这条数据。但只要有这个键，**值必须是列表**，不能写 `null`。

**示例（含 II 类）**：`data/examples/hla_typing.class_ii.example.json`  
**最小示例**（仅 I 类）：各 `deliveries/R*/to_immunogen/hla_typing.json`。

---

## 2. 上游结果怎么塞进 JSON（字段对照）

### 2.1 OptiType

- 一般只出 **A / B / C 三座**（I 类），**没有** II 类。  
- 把它的 A 结果写进 `HLA-A` 的列表，B 进 `HLA-B`，C 进 `HLA-C`。  
- 原输出有时带 `HLA-` 有时不带，**保留一种风格并统一**即可（本仓库校验不强制 `HLA-` 前缀，只要 I 类能被后续工具识别）。

### 2.2 HLA-HD（举例）

- 常能出 **A、B、C 以及 DRB1、DQA1、DQB1、DPA1、DPB1** 等（以你实际报告为准）。  
- **对号入座**：`A` → `HLA-A`，`B` → `HLA-B`，`C` → `HLA-C`；`DRB1` → `HLA-DRB1`；`DQA1` → `HLA-DQA1` ……  
- 若报告里 DQA1 与 DQB1 成对出现，**拆成两个数组里的各一条**，成对关系放到 BioDriver 自己的 `meta` 里，本文件不强求嵌套。

### 2.3 和 NetMHCIIpan 的「名字对齐」

II 类名字在工具里要写成 **带 `HLA-` 的完整形式** 更常见。  
- **默认可用** `data/hla_allele_map_netmhciipan.json` + `scripts/hla_allele_to_netmhciipan.py` 做「简写补全/手工改一行」；**白话说明**在 [`allele_naming_simple.md`](./allele_naming_simple.md)。

---

## 3. 程序怎么校验（validate_input.py）

- 少 `HLA-A` / `HLA-B` / `HLA-C`，或 I 类列表里**全是空** → 报错。  
- 有 II 类键时，**必须是列表**。  
- 根上**不认识的键名**会 **警告、不失败**（方便以后加字段）。  

---

*版本：与 `scripts/hla_typing_spec.py`、P0 NetMHCIIpan 准备一致。*
