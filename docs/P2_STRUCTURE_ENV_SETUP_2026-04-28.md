# P2 结构建模环境安装记录（2026-04-28）

## 1. 结论

本次只安装真实结构工具环境，不生成、不伪造任何 peptide-MHC 结构数据。

当前可用环境：

- `pandora_src`
  - Python: `3.11.15`
  - 已安装：`csb-pandora`
  - 已安装：`MODELLER 10.8`
  - 已安装：`MUSCLE 5.2.linux64`
  - 已配置：MODELLER 真实 license key（不在文档中记录明文）
  - 已修复：`csb-pandora 0.9` 与新版 Biopython 不兼容，已将 Biopython 降到 `1.79`
  - 已修复：PANDORA 0.9 写死 IMGT `http://` 下载地址，已改为 `https://`
  - 当前状态：正在下载 IMGT 官方 `IMGT3DFlatFiles.tgz` 结构库（约 1.6GB），用于构建 MHC-I 模板数据库
- `colabfold`
  - Python: `3.10.20`
  - 已安装：`colabfold 1.6.1`
  - 已安装：`alphafold-colabfold` 相关依赖
  - `colabfold_batch` 命令可显示帮助
  - 当前限制：`jax 0.6.2` 目前仅识别到 CPU；不匹配的 CUDA JAX 插件已清理，尚未完成 GPU 化
- 系统 BLAST
  - 已安装：`blastp 2.9.0+`
  - 已安装：`makeblastdb 2.9.0+`
  - 说明：NCBI 官方 `2.17.0+` Linux tarball 下载速度过慢，本次未完成；当前先保证 BLAST 命令入口可用

已清理环境：

- `pandora`
- `pandora_tools`

这两个环境来自前期 conda 求解失败/中断尝试，已删除，避免后续误用。

## 2. 已执行安装路径

### PANDORA 主环境

```bash
conda create -y -n pandora_src --override-channels -c defaults python=3.11 pip git
conda run -n pandora_src python -m pip install csb-pandora
conda install -y -n pandora_src --override-channels -c salilab -c conda-forge -c defaults modeller
conda install -y -n pandora_src --override-channels -c bioconda -c conda-forge -c defaults \
  https://conda.anaconda.org/bioconda/linux-64/muscle-5.1.0-h4ac6f70_1.tar.bz2
```

### ColabFold 环境

```bash
conda create -y -n colabfold --override-channels -c defaults python=3.10 pip
conda run -n colabfold python -m pip install colabfold
conda run -n colabfold python -m pip install 'colabfold[alphafold]'
```

### BLAST 命令入口

```bash
apt-get update
apt-get install -y ncbi-blast+
```

## 3. 当前验证结果

```bash
conda run -n pandora_src muscle -version
# muscle 5.2.linux64 [-]

blastp -version
# blastp: 2.9.0+

makeblastdb -version
# makeblastdb: 2.9.0+

conda run -n colabfold /root/miniconda3/envs/colabfold/bin/colabfold_batch --help
# 可正常显示帮助

conda run -n colabfold python -c "import jax; print(jax.__version__); print(jax.devices())"
# 0.6.2
# [CpuDevice(id=0)]
```

MODELLER 当前验证结果：

```text
MODELLER import OK
PANDORA import OK
```

MODELLER license 已配置并通过导入验证；文档不记录 key 明文。

## 4. 下一步

1. 等待 IMGT 官方结构库下载完成。
2. 完成 `pandora_mhci_database.pkl` 构建。
3. 数据库安装完成后，开始用 `results/<run_id>/structure_inputs/mhc_chain_sequences.fasta` 和 Top peptide-HLA 组合生成真实 PDB。

## 5. 不做的事

- 不使用粗粒度 `complex.pdb` 伪装为 PANDORA/AFM 输出。
- 不把 `meta.json` 改成 `structure_backend=pandora/afm`，除非真实 PDB 已经生成并通过质检。
- 不伪造 MODELLER license key。
