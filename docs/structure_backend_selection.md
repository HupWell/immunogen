# 结构后端选型说明（PANDORA vs ColabFold/AFM）

## 1. 结论先看（默认策略）

- 默认推荐：**组合策略**
  - 第一阶段批量：优先 `PANDORA`（模板驱动，CPU 可跑，成本低，稳定产出快）
  - 重点候选复核：再用 `ColabFold/AFM`（需要 GPU，耗时更高，但在复杂构象上更有潜力）
- 交付策略：
  - 常规流水线默认仍是 `coarse`（保证无依赖可跑）
  - 当有真实结构时，使用 `prepare_simhub_delivery.py --structure_backend pandora|afm --structure_input_pdb <pdb>`

## 2. 为什么这样选（简化版）

- `PANDORA` 更像“快而稳”的工程工具：对 peptide-MHC 模板场景友好，便于批量。
- `AFM` 更像“重武器”：对难例可能更好，但 GPU、排队和重试成本都更高。
- 本项目当前目标是先保证可复现交付，再逐步提高精度，因此采用“先 PANDORA，后 AFM 复核”的二段式更符合工程节奏。

## 3. 对比表（非生信背景可直接用）

| 维度 | PANDORA | ColabFold/AFM |
|---|---|---|
| 依赖 | 主要 CPU + 模板数据 | 通常需要 GPU（本地/云） |
| 批量能力 | 高，适合先跑全量候选 | 中，适合 Top 少量精修 |
| 单例耗时（经验） | 分钟级（取决于模板与环境） | 10 分钟到数小时（看队列/GPU） |
| 成本 | 低 | 中高 |
| 失败模式 | 模板缺失或匹配差 | 资源排队、显存/超时、随机波动 |
| 建议角色 | 默认主力 | 重点复核/疑难样本 |

> 说明：耗时是工程经验区间，不同机器和队列差异会很大。

## 4. 最小落地流程（你现在就能执行）

1. 用现有流程生成序列输入（已支持）：
   - `python scripts/prepare_mhc_chain_sequences.py --run_id R001 --strict`
2. 在结构工具侧生成真实 PDB（PANDORA 或 AFM）。
3. 回填到本项目交付：
   - `python scripts/prepare_simhub_delivery.py --run_id R001 --top_k 3 --structure_backend pandora --structure_input_pdb <你的pdb>`
   - 或 `--structure_backend afm`
4. 检查 `deliveries/<run_id>/to_simhub/<case_id>/meta.json`：
   - `structure_backend`、`structure_source`、`replaces_coarse` 是否正确。

## 5. GPU/队列建议（避免卡住主流程）

- 不要把 AFM 强行绑到默认 `run_all.py` 的必经步骤。
- 推荐独立执行结构预测任务，并把结果 PDB 回填交付脚本。
- 建议排队策略：
  - 第一轮：PANDORA 批量跑 Top 3~5
  - 第二轮：AFM 只跑 Top 1~2 或 PANDORA 低置信样本

## 6. 何时租 GPU（你的决策阈值）

- 满足任一条件再租：
  - 需要在 24 小时内完成多个病例的 AFM 结构；
  - PANDORA 无模板或结果不稳定，且该候选是关键交付对象；
  - 需要做 AFM 多次重复以评估稳定性。

## 7. 当前仓库的对应能力

- `scripts/prepare_mhc_chain_sequences.py`：按分型准备 MHC-I α + β2m 序列输入。
- `scripts/prepare_simhub_delivery.py`：支持 `coarse/pandora/afm` 交付封装与元数据审计。
- `scripts/run_all.py`：默认可不绑定重型结构步骤，保持可复现与低门槛。

