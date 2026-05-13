# MoDES: Multi-Omics Discordance/Event State inference

多组学不一致驱动的调控事件状态推断框架。

## 概述

MoDES 是一个以 **regulatory event** 为基本分析单位的统计框架。与传统的多组学整合方法不同，MoDES 将跨组学层之间的一致与不一致作为主要信号，系统地将每个调控事件分类为可解释的状态。

MoDES 输出两层信息：
1. **biological state**: concordant, chromatin_primed, rna_only, discordant_opposite, null
2. **artifact_risk**: low, medium, high

当前版本 (v0.1.0) 为 **MoDES-RA**（RNA + ATAC）。

## Current scope of v0.1.0

- RNA + ATAC only
- binary condition only (two-group comparison)
- bulk or externally generated pseudobulk recommended
- event candidates require genome_annotation, tss_map, or external_links
- no native protein model yet
- no native spatial graph model yet

## 安装

```bash
pip install -e .
```

## 快速开始

### 方案 A：使用 external peak-gene links（推荐）

```python
import pandas as pd
from modes import MoDES, MoDEData

data = MoDEData.from_matrices(
    rna_counts="rna_counts.tsv",
    atac_counts="atac_peaks.tsv",
    metadata="sample_metadata.tsv",
    condition_col="condition",
)

links = pd.read_csv("peak_gene_links.tsv", sep="\t")
modes = MoDES(
    data=data,
    condition_col="condition",
    external_links=links,
)
result = modes.run()
```

### 方案 B：使用 GTF genome annotation

```python
modes = MoDES(
    data=data,
    condition_col="condition",
    genome_annotation="genes.gtf",
)
result = modes.run()
```

### 查看和导出结果

```python
print(result.summary())
result.to_tsv("output/")
result.to_graphml("output/event_network.graphml")
result.to_report("output/report.html")
```

### 输出文件

| 文件 | 说明 |
|------|------|
| event_table.tsv | 主输出表 |
| event_state_confidence.tsv | 状态置信度 |
| event_layer_effects.tsv | 每层效应大小 |
| event_evidence_vectors.tsv | 证据向量 |
| model_diagnostics.tsv | 模型诊断信息 |
| run_params.tsv | 运行参数 |

#### event_table.tsv 字段

| 字段 | 说明 |
|------|------|
| event_id | 事件唯一标识 |
| gene | target gene |
| peak_id | regulatory peak |
| state | biological state (concordant / chromatin_primed / rna_only / discordant_opposite / null) |
| state_confidence | 状态置信度 [0, 1] |
| artifact_risk | 技术伪影风险 (low / medium / high)，不是生物状态 |
| artifact_reason | 伪影原因（分号分隔） |
| event_pval | event-level p-value |
| event_fdr | event-level BH-corrected FDR |
| atac_coef / atac_pval / atac_fdr | ATAC 效应估计 |
| rna_coef / rna_pval / rna_fdr | RNA 效应估计 |
| rna_after_atac_coef / rna_after_atac_pval / rna_after_atac_fdr | 控制 linked ATAC peak 后的条件效应 |

## 算法流程

1. **Event Candidate Construction** — 将 peaks 链接到 target genes
2. **Effect Size Estimation** — NB GLM + 经验贝叶斯收缩估计每个模态的 condition effect
3. **Conditional Decomposition** — RNA condition effect after adjustment for the linked ATAC peak
4. **Evidence Vector Construction** — 构造 D_e = [z_ATAC, z_RNA, z_RNA|ATAC, q]
5. **State Classification** — 规则判别 + 经验贝叶斯置信度计算

## 输入格式

- **Bulk**: TSV/CSV count matrices + metadata（推荐）
- **Single-cell**: experimental in v0.1.0. 推荐先外部聚合为 pseudobulk，再用 `MoDEData.from_matrices()` 加载。`MoDEData.from_pseudobulk()` 支持简单 AnnData（RNA in X, ATAC in obsm["atac"]）
- **Spatial**: planned. Current v0.1.0 does not yet model spatial coordinates or neighborhood graphs.

## 引用

MoDES: Multi-Omics Discordance-guided decomposition of regulatory event states.

## 许可证

MIT
