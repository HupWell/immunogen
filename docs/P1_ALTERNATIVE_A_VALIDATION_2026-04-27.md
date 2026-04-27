# P1 替代方案 A 验收归档（2026-04-27）

## 1. 执行背景

- 上游 BioDriver 暂时无法提供患者真实 HLA-II 分型。
- 按替代方案 A，使用公开文献口径（Dhall et al., 2020, TCGA-SKCM）补充 II 类分型字段，用于流程演示验证（非临床终版）。

## 2. 影响范围

- `R_public_001`
- `R002`
- `R003`

## 3. 已执行动作

1. 更新三个 run 的 `deliveries/<run_id>/to_immunogen/hla_typing.json`，新增：
   - `HLA-DRB1`
   - `HLA-DQB1`
   - `HLA-DPB1`
2. 重跑严格链路（full）并生成日志：
   - `results/<run_id>/validation_logs/strict_run_2026-04-27.log`
3. 执行严格真实性验收并生成日志：
   - `results/<run_id>/validation_logs/strict_check_2026-04-27.log`
4. 更新三个 run 的 `SELF_CHECK.md`，新增“公开队列数据、非临床终版”声明。

## 4. 执行命令

```bash
source scripts/env_netmhcpan.sh
for r in R_public_001 R002 R003; do
  python scripts/run_all.py --run_id "$r" --target full \
    --backend_mhc1_netmhcpan real_cmd --backend_mhc1_bigmhc off \
    --mhc2_backend real_tsv --require_real_mhc1_cv --require_real_mhc2 \
    --require_real_immunogenicity_prime --require_real_immunogenicity_repitope
done
```

```bash
for r in R_public_001 R002 R003; do
  python scripts/check_epitope_realization.py --run_id "$r" \
    --require_mhc1_cv_real --require_mhc2_real
done
```

## 5. 关键字段结果

### R_public_001

- `mhc1_cv_source_netmhcpan`: `real_cmd`
- `mhc2_backend`: `netmhciipan`
- `immunogenicity_source_deepimmuno`: `deepimmuno_real_cmd` = 198
- `immunogenicity_source_prime`: `prime_real_cmd` = 198
- `immunogenicity_source_repitope`: `repitope_public_dataset_knn` = 186, `repitope_public_dataset_exact` = 12

### R002

- `mhc1_cv_source_netmhcpan`: `real_cmd`
- `mhc2_backend`: `netmhciipan`
- `immunogenicity_source_deepimmuno`: `deepimmuno_real_cmd` = 18
- `immunogenicity_source_prime`: `prime_real_cmd` = 18
- `immunogenicity_source_repitope`: `repitope_public_dataset_knn` = 12, `repitope_public_dataset_exact` = 6

### R003

- `mhc1_cv_source_netmhcpan`: `real_cmd`
- `mhc2_backend`: `netmhciipan`
- `immunogenicity_source_deepimmuno`: `deepimmuno_real_cmd` = 18
- `immunogenicity_source_prime`: `prime_real_cmd` = 18
- `immunogenicity_source_repitope`: `repitope_public_dataset_knn` = 12, `repitope_public_dataset_exact` = 6

## 6. 结论与限制

- 三个 run 在替代方案 A 下已完成严格链路与验收。
- 当前结果可用于流程展示与阶段性汇报。
- 仍需在 BioDriver 提供患者真实 HLA-II 分型后，进行终版替换与全量复跑。

## 7. 补充材料下载复核记录

- 2026-04-27 重新尝试访问 Frontiers 官方页面与 XML。
- 可访问资源：
  - `https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2020.00221/xml`
  - XML 中确认 `Supplementary Table S2` 指向 `Data_Sheet_1.ZIP`，且文章说明该表包含 class I/II（DP/DQ/DR）四位 HLA 分型。
- 仍未成功获取资源：
  - `Data_Sheet_1.ZIP` 附件本体未能通过 Frontiers 常见直链路径直接下载。
  - PMC Open Access API 返回的 OA 包链接指向 `ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_package/d0/d6/PMC7113398.tar.gz`，但当前 FTP/HTTPS 路径不可取回。
  - 作者 SKCMhrp 数据页 `https://webs.iiitd.edu.in/raghava/skcmhrp/data.php` 当前显示 `Data download coming soon`。
- 因此当前 P1 替代方案 A 仍保留“公开文献口径演示级”声明，不升级为“已获得原始样本级补充表”的临床或样本级口径。
