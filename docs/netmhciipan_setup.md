# NetMHCIIpan：最简说明（conda + 分开装）

**原则：conda 只装 Python 依赖；NetMHCIIpan 是单独下载的程序，不放 conda 里也行。**

---

## 情况 A：先不装 NetMHCII（推荐新手）

- 不用做任何事。  
- 跑流程时用：`--mhc2_backend auto` 或 `proxy`（MHC-II 用**代理分**，照样出完整结果）。  
- 等需要「真 II 类」再去做情况 B。

---

## 情况 B：在 WSL2 / Linux 里接真 NetMHCII（三步）

在 **WSL2 的 Ubuntu 终端**里操作（[微软 WSL 安装说明](https://learn.microsoft.com/zh-cn/windows/wsl/install)），不要用 Windows 自带的 CMD 去跑 `netMHCIIpan`。

### 第 1 步：Python 用 conda（和平时一样）

```bash
conda create -n immunogen python=3.10 -y
conda activate immunogen
# 在仓库根目录，依赖按 README 装齐（pandas、mhcflurry、ViennaRNA 等）
```

以后跑 `python scripts/run_all.py` 时**先** `conda activate immunogen` 即可。

### 第 2 步：只装 NetMHCIIpan 本体

1. 打开：  
   [https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.1](https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.1)  
   （或更新的 4.2 页面里的 **Download**）  
2. 同意条款后下载 **可移植包** `tar.gz`，解压到固定目录，例如：  
   `~/tools/NetMHCIIpan-4.1/`
3. 在该目录下找到可执行名（常见为 `netMHCIIpan`），能跑通：  
   `~/tools/NetMHCIIpan-4.1/netMHCIIpan -h`

> 不追求「装进 conda」；**路径写对**就行。

### 第 3 步：只设 1～2 个环境变量（同一次终端里）

在**要跑 Python 的同一个 WSL 终端**里（路径改成你的）：

```bash
export NETMHCIIPAN_BIN=$HOME/tools/NetMHCIIpan-4.1/netMHCIIpan
export NETMHCIIPAN_HOME=$HOME/tools/NetMHCIIpan-4.1
```

然后进入**仓库根目录**（WSL 里若项目在 `/mnt/d/desktop/immunogen` 之类）再跑：

```bash
export PYTHONUTF8=1
python scripts/predict_mhc_ranking.py --run_id 你的run_id --mhc2_backend netmhciipan
```

**成功时**：屏幕出现 `MHC-II 层: 模式=netmhciipan`；`results/.../peptide_mhc_ranking.csv` 里 **`mhc2_backend` 列为 `netmhciipan`** 且有 `mhc2_el_rank` 等。

---

## 数据上只要多改一处（II 类）

`deliveries/<run_id>/to_immunogen/hla_typing.json` 在原有 `HLA-A/B/C` 外，**至少加**（例）：

```json
"HLA-DRB1": ["HLA-DRB1*15:01"]
```

完整多键示例见 `data/examples/hla_typing.class_ii.example.json`。没有 II 类时工具无法选链，会退回代理或报错（取决于 `--mhc2_backend`）。

---

## 想再少记一点时

- **不必记**：其它环境变量可全不设；需要时再查下表。  
- **只记两个**：`NETMHCIIPAN_BIN` + 多数情况下要 **`NETMHCIIPAN_HOME` = 安装根目录**。

| 变量 | 什么时候要设 |
|------|----------------|
| `NETMHCIIPAN_BIN` | 必设：可执行文件完整路径。 |
| `NETMHCIIPAN_HOME` | 建议设：与解压目录一致；若 `-h` 能跑、预测报找不到数据，**必设**。 |
| `NETMHCIIPAN_BA=0` | 一般**不设**；只有解析失败、且工具文档说不要 `-BA` 时再设。 |
| `NETMHCIIPAN_TIMEOUT` | 一般**不设**；肽特别多再加大。 |

---

## 和本仓库代码的对应关系

- 子进程与解析：`scripts/netmhciipan_runner.py`  
- 排名与开关：`scripts/predict_mhc_ranking.py` 的 `--mhc2_backend`  

**名字对不上、工具报错**时，用 `data/hla_allele_map_netmhciipan.json` 的 `manual_overrides` 改一行，见 `docs/allele_naming_simple.md`。
