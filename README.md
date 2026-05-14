# MoDES: Multi-Omics Discordance/Event State inference

**多组学不一致驱动的调控事件状态推断框架**

MoDES 是一个以 **regulatory event**（调控事件）为基本分析单位的统计框架。不同于传统的多组学整合方法，MoDES
将跨组学层之间的一致与不一致作为主要信号，系统地将每个调控事件分类为可解释的生物学状态，
并同时评估其技术伪影风险。

---

## 1. 基本原理

### 1.1 核心思想

传统多组学分析通常将 RNA、ATAC、Protein 等模态分别做差异分析，然后交叉验证 —
"RNA 找到了差异基因，ATAC 看看附近 peak 是否开放"。这种逻辑隐式地把一个模态当作主分析，
另一个模态当作验证层。

MoDES 的逻辑不同：

> **RNA、ATAC、protein 不是主次关系；它们是同一个调控事件的不同观测层。**

- **一致性**（concordance）说明完整调控链条激活：chromatin → transcription → protein
- **ATAC-only**（chromatin priming）说明染色质层已就绪，但转录/蛋白层尚未响应
- **RNA-only**（trans-driven）说明 RNA 变化不由局部染色质解释，可能来自 trans 调控、RNA stability 等
- **RNA-protein 不一致**说明 post-transcriptional buffering、protein memory 或技术问题
- **空间局部一致**说明 niche-driven regulation

这种视角的本质提升在于 **数据利用率**：不是把多组学数据压缩成一个 embedding，而是利用
不同组学层之间的不一致模式来推断调控事件的发生层级。

### 1.2 分析单位：regulatory event，不是 gene

一个 regulatory event 定义为：

```
e = (TF, enhancer/peak, target_gene, context)
```

例如：

> **STAT1 motif / enhancer chr1:100-200 accessibility ↑**
> → **IFIT3 RNA ↑**
> → 发生在 **disease monocyte** 中

这比 "gene-level differential expression" 或 "peak-level differential accessibility" 更接近
真实的调控生物学。

### 1.3 生物学状态分类

MoDES-RA v0.1.0 基于 RNA + ATAC 两层数据，将每个事件分类为以下生物学状态：

| 状态 | 模式 | 生物学解释 |
|---|---|---|
| `concordant` | ATAC↑ RNA↑ | 局部染色质开放驱动转录激活（完整顺式调控链） |
| `chromatin_primed` | ATAC↑ RNA→ | 染色质已就绪，转录尚未启动（epigenetic priming） |
| `rna_only` | ATAC→ RNA↑ | RNA 变化不由局部 chromatin 解释（trans 调控 / RNA stability / 未测调控层） |
| `discordant_opposite` | ATAC↑ RNA↓ 或 ATAC↓ RNA↑ | 两层方向相反，可能是复杂调控或技术问题 |
| `null` | ATAC→ RNA→ | 该事件在当前条件下没有显著变化 |

同时，每个事件附带 **artifact_risk**（技术伪影风险）：

| 风险等级 | 含义 |
|---|---|
| `low` | 数据质量良好，结果可信 |
| `medium` | 存在一定的质量问题 |
| `high` | 单模态显著 + 低质量分数，建议谨慎解读 |

这种双层设计（biological state + artifact risk）比将 "artifact_like" 作为主状态更合理：
一个事件可以同时是 `chromatin_primed` 且 `artifact_risk = high`，
而不是被简单标记为 "artifact" 而丢失生物学信息。

---

## 2. 当前版本状态

**MoDES-RA v0.1.0-alpha** — RNA + ATAC 原型。

| 能力 | 状态 |
|---|---|
| RNA + ATAC 两层分析 | ✅ 已支持 |
| 二分类条件比较 | ✅ 已支持 |
| Bulk 数据输入 | ✅ 推荐 |
| Pseudobulk 聚合 | ✅ 实验性支持 |
| Protein 层 | 🔮 计划中 (v0.2) |
| Spatial graph | 🔮 计划中 (v0.4) |
| 多分类条件 / 连续协变量 | 🔮 计划中 |
| 时间 / pseudotime 延迟 | 🔮 计划中 |

---

## 3. 安装

```bash
# 基础安装
pip install -e .

# 运行测试（可选）
pip install -r requirements-dev.txt
python -m pytest -q
```

依赖：`numpy`, `scipy`, `pandas`, `statsmodels`, `anndata`, `matplotlib`, `seaborn`, `networkx`

---

## 4. 使用方法

### 4.1 快速开始

```python
import pandas as pd
from modes import MoDES, MoDEData

# 1. 加载数据
data = MoDEData.from_matrices(
    rna_counts="rna_counts.tsv",
    atac_counts="atac_peaks.tsv",
    metadata="sample_metadata.tsv",
    condition_col="condition",
    index_col=0,
)

# 2. 加载 peak-gene 链接（推荐方式）
links = pd.read_csv("peak_gene_links.tsv", sep="\t")

# 3. 运行 MoDES
modes = MoDES(
    data=data,
    condition_col="condition",
    external_links=links,
)
result = modes.run()

# 4. 查看结果
print(result.summary())

# 5. 导出
result.to_tsv("output/")
result.to_graphml("output/network.graphml")
result.to_report("output/report.html")
```

### 4.2 分步运行

你也可以逐步执行，方便调试和交互式分析：

```python
modes = MoDES(data=data, condition_col="condition", external_links=links)

# Step 1: 构建候选事件
events = modes.build_events()

# Step 2: 估计 ATAC 和 RNA 效应
atac_effects, rna_effects = modes.estimate_effects()

# Step 3: 条件分解（RNA after ATAC）
conditional = modes.decompose()

# Step 4: 构造证据向量
evidence = modes.build_evidence()

# Step 5: 状态分类
states = modes.classify_states()

# 组装最终结果
result = modes._assemble_results()
```

### 4.3 过滤结果

```python
# 只看 concordant 事件
conc = result.filter(state="concordant")

# 排除高风险事件
clean = result.filter(exclude_high_artifact=True)

# 按 event FDR 筛选
sig = result.filter(max_event_fdr=0.1)

# 组合过滤
trusted = result.filter(
    state="concordant",
    min_confidence=0.8,
    max_event_fdr=0.05,
    exclude_high_artifact=True,
)
```

### 4.4 运行示例

```bash
python examples/minimal_bulk/run_minimal.py
```

---

## 5. 参数设置

### 5.1 MoDES 初始化参数

| 参数 | 类型 | 默认值 | 说明 |
|---|---|---|---|
| `data` | `MoDEData` | 必需 | 输入数据容器 |
| `condition_col` | `str` | 必需 | `data.obs` 中指定条件的列名（必须为二分类） |
| `covariate_cols` | `list[str]` | `[]` | 额外协变量列名 |
| `donor_col` | `str` | `None` | 捐赠者/重复标识列（作为固定效应） |
| `batch_col` | `str` | `None` | 批次标识列（作为固定效应） |
| `fdr_threshold` | `float` | `0.1` | 状态分类的 FDR 显著性阈值 |
| `genome_annotation` | `str` | `None` | GTF/GFF 文件路径，用于获取基因 TSS 坐标 |
| `external_links` | `DataFrame` | `None` | 预计算的 peak-to-gene 链接表（推荐） |
| `motif_annotation` | `DataFrame` | `None` | peak-to-TF motif 注释表 |
| `tss_map` | `dict` | `None` | 手动指定的 gene → (name, chr, tss_pos) 映射 |

### 5.2 `external_links` 格式

| 列名 | 必需 | 说明 |
|---|---|---|
| `peak_id` | ✅ | peak 标识符，格式 `chr:start-end` |
| `gene` | ✅ | 基因名称 |
| `tf_name` | 可选 | TF 注释 |
| `source` | 可选 | 链接来源（如 `scenic`, `scarlink`, `archr`） |

### 5.3 `MoDEData.from_matrices()` 参数

| 参数 | 类型 | 说明 |
|---|---|---|
| `rna_counts` | `str` 或 `DataFrame` | RNA count matrix（samples × genes） |
| `atac_counts` | `str` 或 `DataFrame` | ATAC peak count matrix（samples × peaks） |
| `metadata` | `str` 或 `DataFrame` | 样本元数据，必须包含条件列 |
| `condition_col` | `str` | 条件列名 |
| `donor_col` | `str` | 捐赠者列名（可选） |
| `batch_col` | `str` | 批次列名（可选） |
| `index_col` | `int` | TSV 文件中行索引列（默认 0） |

---

## 6. 输入文件格式

### 6.1 Bulk count matrix（TSV）

**rna_counts.tsv** — 行 = 样本，列 = 基因：

| sample | STAT1 | GZMB | IL7R |
|---|---|---|---|
| ctrl_1 | 120 | 85 | 200 |
| ctrl_2 | 98 | 92 | 180 |
| trt_1 | 350 | 90 | 220 |

**atac_counts.tsv** — 行 = 样本，列 = peaks：

| sample | chr1:100-200 | chr2:300-400 |
|---|---|---|
| ctrl_1 | 50 | 60 |
| trt_1 | 180 | 55 |

### 6.2 样本元数据

**metadata.tsv**：

| sample | condition | batch |
|---|---|---|
| ctrl_1 | control | A |
| trt_1 | treatment | B |

`condition` 列必须为二分类（如 control / treatment）。

### 6.3 Peak-gene 链接

**peak_gene_links.tsv**：

| peak_id | gene | tf_name |
|---|---|---|
| chr1:100-200 | STAT1 | IRF1 |
| chr2:300-400 | GZMB | |

> **注意**：普通 gene symbols（如 STAT1、GZMB）没有基因组坐标。如果不提供
> `external_links`、`genome_annotation` 或 `tss_map`，
> MoDES 会报错提示无法生成候选事件。

---

## 7. 算法详解

### Step 1：Event Candidate Construction（事件候选构建）

为每个 gene 寻找候选 regulatory elements：

- **promoter peaks**：TSS ± 2kb 范围内的 peaks
- **distal peaks**：TSS ± 250kb 范围内的 peaks
- **external links**：来自 SCENIC+、SCARlink、ArchR 等的预计算链接
- **motif annotation**：peak 上的 TF motif 信息

### Step 2：Effect Size Estimation（效应估计）

对每个 peak 和 gene，使用 **Negative Binomial GLM**（log link）估计 condition effect：

```
log(E[Y_u]) = α + β_cond × C_u + X_u × γ + offset(log(libsize_u))
```

其中：
- `C_u` = condition indicator（0/1）
- `X_u` = 协变量矩阵（batch、donor、age 等）
- `offset` = log library size（DESeq2-style median-of-ratios normalization）

效应估计后应用 **limma-style 经验贝叶斯方差收缩**：
- 将所有 feature 的方差 pooled 到一起，估计先验分布
- 计算后验方差，得到更稳定的 moderated t-statistic
- 对小样本场景特别有效

多重检验校正：**Benjamini-Hochberg FDR**（按 modality 分别校正）。

### Step 3：Conditional Decomposition（条件分解）

核心统计问题：**RNA condition effect 是否可以被 local chromatin accessibility 解释？**

比较两个模型：

**Model 0（RNA-only）**：
```
log(E[RNA_g]) ~ Condition + Covariates
```

**Model 1（RNA | ATAC）**：
```
log(E[RNA_g]) ~ Condition + ATAC_peak + Covariates
```

关键量：**β_cond 从 Model 0 到 Model 1 的衰减**。

| β_cond 变化 | 解释 |
|---|---|
| Model 0 显著，Model 1 不显著 | RNA 变化可被 local chromatin 解释（concordant 证据） |
| Model 0 显著，Model 1 仍显著，衰减很小 | RNA 有 ATAC 不能解释的剩余效应（rna_only 证据） |
| Model 0 不显著 | 无 RNA 层效应 |

> **注意**：条件模型中使用单个 linked peak 作为协变量，因此解释应为
> "condition effect after adjustment for the linked ATAC peak"，
> 而非 "ATAC explains RNA"。

### Step 4：Evidence Vector Construction（证据向量构造）

对每个 event `e`，构造标准化证据向量：

```
D_e = [z_ATAC, z_RNA, z_RNA|ATAC, quality_score]
```

其中：
- `z_m = coef_m / SE(coef_m)` — 标准化效应（z-score）
- `quality_score ∈ [0, 1]` — 数据质量分数（基于 sequencing depth、detection rate、batch association 等）

### Step 5：State Classification（状态分类）

两层分类：

**Stage 1 — 基于规则**：
```
if ATAC_sig and RNA_sig and same_direction → concordant
elif ATAC_sig and RNA_sig and opposite → discordant_opposite
elif ATAC_sig and not RNA_sig → chromatin_primed
elif not ATAC_sig and RNA_sig → rna_only
else → null
```

**Stage 2 — 经验贝叶斯 refinement**：
- 用 rule-based 的结果作为初始标签
- 拟合每个状态的 Gaussian 分布（在 evidence vector 空间）
- 计算每个事件的 state_confidence（后验概率近似值）

同时计算 **artifact_risk**：
```
quality_score < threshold → 检查是否为单模态显著
  yes → artifact_risk = high, artifact_reason = low_quality_score;single_modality_low_quality
  no  → artifact_risk = medium
else → artifact_risk = low
```

---

## 8. 输出文件

### 8.1 文件清单

| 文件 | 格式 | 说明 |
|---|---|---|
| `event_table.tsv` | TSV | **主输出**：每个 peak-gene 事件一行 |
| `event_state_confidence.tsv` | TSV | 每个事件的状态、置信度、artifact_risk |
| `event_layer_effects.tsv` | TSV | ATAC/RNA 层效应估计值 |
| `event_evidence_vectors.tsv` | TSV | 状态分类使用的 evidence vector |
| `model_diagnostics.tsv` | TSV | ATAC/RNA marginal GLM 模型诊断 |
| `run_params.tsv` | TSV | 运行参数和数据摘要 |
| `event_network.graphml` | GraphML | TF-peak-gene 网络（可选） |
| `report.html` | HTML | 可视化和统计摘要（可选） |

> `model_diagnostics.tsv` 当前仅包含 marginal ATAC 和 RNA GLM 的诊断信息。
> 条件分解（RNA-after-ATAC）的诊断在 `event_table.tsv` 中汇总，
> 完整的条件模型诊断将在后续版本中提供。

### 8.2 event_table.tsv 字段详解

| 字段 | 类型 | 说明 |
|---|---|---|
| `event_id` | str | 事件唯一标识符 |
| `tf_name` | str/NA | 与该 peak 关联的 TF 名称（来自 motif 注释） |
| `peak_id` | str | ATAC peak 标识符（如 `chr1:100-200`） |
| `gene` | str | target gene 名称 |
| `context` | str | 上下文标签（如 cell type + condition） |
| `state` | str | 生物学状态：`concordant` / `chromatin_primed` / `rna_only` / `discordant_opposite` / `null` |
| `state_confidence` | float | 状态置信度，范围 0–1（非严格校准的后验概率，不宜在论文中称 posterior） |
| `artifact_risk` | str | 技术风险：`low` / `medium` / `high`（不是生物学状态） |
| `artifact_reason` | str | 风险原因，分号分隔（如 `low_quality_score;single_modality_low_quality`） |
| `event_pval` | float | event-level p-value（根据 state 组合 ATAC/RNA p-value） |
| `event_fdr` | float | event-level BH-corrected FDR |
| `atac_coef` | float | condition 对 ATAC peak 的效应（log fold change 尺度） |
| `atac_se` | float | ATAC 效应的标准误 |
| `atac_pval` | float | ATAC 层 p-value |
| `atac_fdr` | float | ATAC 层 BH-corrected FDR |
| `atac_direction` | int | ATAC 效应方向（+1 / -1 / 0） |
| `rna_coef` | float | condition 对 RNA gene 的效应（log fold change 尺度） |
| `rna_se` | float | RNA 效应的标准误 |
| `rna_pval` | float | RNA 层 p-value |
| `rna_fdr` | float | RNA 层 BH-corrected FDR |
| `rna_direction` | int | RNA 效应方向（+1 / -1 / 0） |
| `rna_after_atac_coef` | float | 控制 linked ATAC peak 后的 condition 效应 |
| `rna_after_atac_se` | float | 条件效应的标准误 |
| `rna_after_atac_pval` | float | 条件效应的 p-value |
| `rna_after_atac_fdr` | float | 条件效应的 BH-corrected FDR |
| `quality_score` | float | event 质量分数（0–1，基于测序深度、检出率等） |

---

## 9. 结果解析指南

### 9.1 如何解读状态

**concordant** → 高置信度的顺式调控事件。局部 chromatin 开放伴随着转录上升。
适合用于构建 core regulatory network。

**chromatin_primed** → 染色质层已变化但转录尚未改变。可能代表：
- 发育过程中的 lineage priming
- 药物处理后早期的表观响应
- 需要二次刺激才能激活的 poised enhancer

**rna_only** → RNA 变化不由 local ATAC 解释。可能来自：
- 远端 enhancer 的 trans 调控
- RNA stability 变化
- TF 活性变化而非 chromatin 变化
- peak calling 未覆盖到关键调控区域

**discordant_opposite** → ATAC 和 RNA 方向相反。最常见的解释：
- 负反馈调控
- repressor binding
- 需要进一步实验验证

**null** → 该事件在当前条件下无显著变化。

### 9.2 如何评估可信度

1. **先看 `artifact_risk`**：优先关注 `low` 和 `medium` 的事件；`high` 的事件需要谨慎
2. **再看 `state_confidence`**：> 0.8 的事件更可靠
3. **最后看 `event_fdr`**：event-level FDR < 0.1 的事件有统计显著性保证

推荐的筛选流程：

```python
trusted = result.filter(
    min_confidence=0.7,
    max_event_fdr=0.1,
    exclude_high_artifact=True,
)
```

### 9.3 常见问题

**Q: 为什么我的 concordant 事件很少？**

A: 可能原因：
- 样本量太小（n < 5 per group），统计 power 不足
- 效应量确实很小（biological）
- 使用了严格的 FDR 阈值（尝试提高 `fdr_threshold`）
- peak-gene 链接质量不够（尝试使用高质量的 external_links）

**Q: 为什么很多事件是 null？**

A: null 是正常结果 — 大多数 peak-gene pairs 在特定条件下确实没有变化。
可以关注非 null 事件的生物学功能富集。

**Q: artifact_risk = high 的事件应该怎么处理？**

A: 不要直接删除 — 它们可能仍然包含生物学信号。建议：
- 检查原始 count 数据的质量
- 对比 `artifact_reason` 字段了解具体原因
- 在后续分析中将它们标记为 "需要验证"

---

## 10. CI / 测试

![tests](https://github.com/Lings01/MoDES/actions/workflows/tests.yml/badge.svg)

```bash
python -m pytest -q          # 运行全部测试
python -m pytest -k "states" # 运行状态分类相关测试
```

---

## 11. 引用

MoDES: Multi-Omics Discordance-guided decomposition of regulatory event states.

## 12. License

MIT
