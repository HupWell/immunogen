# 如何启用 NetMHCIIpan（MHC-II 真预测）

## 1. 你需要什么

- 在 **Linux / WSL2** 中安装好 [NetMHCIIpan](https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.1) 可移植版，或把 `netMHCIIpan` 可执行文件加入 **PATH**。
- `hla_typing.json` 里要有 **能自动转换的 II 类**（至少一条 **HLA-DRB1*xx:yy** 这种；脚本会转为工具常用的 `DRB1_xxyy`）。  
  更复杂的 DQ/DP 双链名请写在 `data/hla_allele_map_netmhciipan.json` 的 `manual_overrides` 里，见 `docs/allele_naming_simple.md`。
- 候选肽在 **PEPTIDE 模式** 下对 NetMHCIIpan 常要求 **长度 ≥9**；更短的链本流程仍用 **MHC-II 代理分**。

## 2. 环境变量

| 变量 | 说明 |
|------|------|
| `NETMHCIIPAN_BIN` | 可执行文件**完整路径**（推荐），例如 `/opt/netMHCIIpan-4.1/netMHCIIpan` |
| `NETMHCIIPAN_HOME` | 工作目录；若工具必须从安装目录读数据，**必须**设置 |
| `NETMHCIIPAN_BA=0` | 不自动加 `-BA`（默认会加，以便解析 nM 列） |
| `NETMHCIIPAN_EXTRA` | 其它参数，空格分隔，例如 `length` 等（需符合你本机 `netMHCIIpan -h`） |
| `NETMHCIIPAN_TIMEOUT` | 秒，默认 600 |

## 3. 如何接到流水线

- `predict_mhc_ranking.py`：`--mhc2_backend auto`（默认，有工具则跑，否则代理） / `netmhciipan`（强制，失败就报错） / `proxy`（始终代理）。
- `run_all.py` 同样带 `--mhc2_backend`。
- 输出 CSV 中新增/使用列：`mhc2_score`、`mhc2_el_rank`、`mhc2_ba_nm`、`mhc2_class2_allele`（本工具输出中的 MHC 列名）、`mhc2_backend`（`netmhciipan` 或 `proxy`）。

## 4. Windows 说明

- 若未在 WSL/容器里装 NetMHCIIpan，**保持 `auto` 或 `proxy` 即可**，流水线照常完成。
- 在 Windows 本机**直接**跑可执行文件的情况较少；**推荐在 WSL2 里装工具并在此环境中运行** `python scripts/run_all.py`（同一路径能访问 `deliveries/`、`results/`）。

## 5. 命令探针

安装完成后可在终端执行（Linux/WSL）：

```bash
netMHCIIpan -h
```

若能看到 `NetMHCIIpan-4` 等版本信息，与 ImmunoGen 的解析器预期一致。具体参数**以你本机 `-h` 为准**；若默认 `-f/-a/-inptype` 不兼容，请发 issue 并附上 `-h` 原文。

## 6. 实现位置

- 子进程与解析：`scripts/netmhciipan_runner.py`
- 与 `rank_score` 衔接：`scripts/predict_mhc_ranking.py`
