# MoDES: Multi-Omics Discordance/Event State inference

多组学不一致驱动的调控事件状态推断框架。

## 概述

MoDES 是一个以 **regulatory event** 为基本分析单位的统计框架。与传统的多组学整合方法不同，MoDES 将跨组学层之间的一致与不一致作为主要信号，系统地将每个调控事件分类为可解释的状态：

- **concordant** — 完整调控链条激活（ATAC↑ RNA↑ Protein↑）
- **chromatin_primed** — 染色质已就绪，转录未启动（ATAC↑ RNA→）
- **rna_only** — RNA 变化不由局部染色质解释（trans/stability 调控）
- **protein_buffered** — 转录变化未传递到蛋白层
- **protein_memory** — 蛋白保留了过去激活状态
- **epigenetic_memory** — 表观层保留历史状态
- **spatial_niche_driven** — 空间微环境驱动
- **artifact_like** — 可能是技术伪影

当前版本 (v0.1.0) 为 **MoDES-RA**（RNA + ATAC），支持 4 种状态：concordant, chromatin_primed, rna_only, artifact_like。

## 安装

```bash
pip install -e .
```

## 快速开始

```python
from modes import MoDES, MoDEData

# 加载数据
data = MoDEData.from_matrices(
    rna_counts="rna_counts.tsv",
    atac_counts="atac_peaks.tsv",
    metadata="sample_metadata.tsv",
    condition_col="condition",
)

# 运行 MoDES 流程
modes = MoDES(data=data, condition_col="condition")
result = modes.run()

# 查看结果
print(result.summary())

# 导出结果
result.to_tsv("output/")
result.to_graphml("output/event_network.graphml")
result.to_report("output/report.html")
```

## 算法流程

1. **Event Candidate Construction** — 将 peaks 链接到 target genes
2. **Effect Size Estimation** — NB GLM + 经验贝叶斯收缩估计每个模态的 condition effect
3. **Conditional Decomposition** — 控制 ATAC 后检测 RNA 的剩余效应
4. **Evidence Vector Construction** — 构造 D_e = [z_ATAC, z_RNA, z_RNA|ATAC, q]
5. **State Classification** — 规则判别 + 经验贝叶斯概率分类

## 输入格式

- Bulk: TSV/CSV count matrices + metadata
- Single-cell: AnnData (.h5ad)，通过 pseudobulk 聚合
- Spatial: 空间坐标 + 图结构

## 引用

MoDES: Multi-Omics Discordance-guided decomposition of regulatory event states.

## 许可证

MIT
