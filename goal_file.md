Lings，下面是针对 49a8e2e 这一版的修改清单。我按优先级排好了，每条都写了：位置、问题、怎么改、验收标准。

当前 README 已经把 v0.1.0 定位成 RNA+ATAC only、binary condition、bulk 或外部 pseudobulk 推荐，并且输出表文档里写的是 state_confidence、artifact_risk、artifact_reason、event_pval、event_fdr。但源码里还有几处和这个设计不完全一致。 ￼

⸻

第一优先级：必须先改

1. 全仓库统一：confidence → state_confidence

位置：

modes/_types.py
modes/core.py
modes/report.py
modes/plotting.py
tests/
README.md

问题：

README 里的主输出字段写的是：

state_confidence

但是 core.py 里的 _event_result_columns() 仍然输出：

"state", "confidence", "quality_score"

也就是说，用户看 README 会找 state_confidence，但真实 event_table.tsv 里是 confidence。 ￼

怎么改：

1.1 修改 EventResult

在 modes/_types.py 里，把：

confidence: float

改成：

state_confidence: float

1.2 修改 core.py

把 _event_result_columns() 里的：

"state", "confidence", "quality_score",

改成：

"state", "state_confidence", "quality_score",

在 _assemble_results() 里，把类似：

confidence = state_row.iloc[0]["state_confidence"]

构造 EventResult 时改成：

state_confidence=confidence

或者直接：

state_confidence = state_row.iloc[0]["state_confidence"]

然后：

EventResult(
    ...
    state_confidence=state_confidence,
    ...
)

1.3 修改所有使用 confidence 的地方

搜索：

grep -R "\bconfidence\b" modes tests

重点改这些：

row["confidence"]
event_table["confidence"]
nlargest(..., "confidence")
df["confidence"] >= min_confidence

改成：

row["state_confidence"]
event_table["state_confidence"]
nlargest(..., "state_confidence")
df["state_confidence"] >= min_confidence

min_confidence 参数名可以保留，因为它是用户传入阈值；但内部过滤列应使用 state_confidence。

验收标准：

assert "state_confidence" in result.event_table.columns
assert "confidence" not in result.event_table.columns

README、测试、report、GraphML、plotting 都使用同一个字段名：

state_confidence

⸻

2. EB 阶段保留 artifact_reason

位置：

modes/states.py

问题：

rule-based 阶段已经生成了：

artifact_risk
artifact_reason

但是 _empirical_bayes_classify() 里只保留了 artifact_risk，没有保留 artifact_reason。所以只要 EB 生效，artifact_reason 就会丢失。states.py 里当前 biological state 和 artifact risk 的分层设计已经是对的，但 EB 分支还没有完整传递 artifact_reason。 ￼

怎么改：

2.1 空 evidence 返回时加列

把：

return pd.DataFrame(columns=[
    "event_id",
    "state",
    "state_confidence",
    "artifact_risk",
])

改成：

return pd.DataFrame(columns=[
    "event_id",
    "state",
    "state_confidence",
    "artifact_risk",
    "artifact_reason",
])

2.2 _empirical_bayes_classify() 中保留 reason

在循环里加：

current_reason = states.iloc[i].get("artifact_reason", "")

然后每个 append 分支都加入：

"artifact_reason": current_reason,

例如：

confidences_all.append({
    "event_id": row["event_id"],
    "state_confidence": float(best_prob),
    "state": best_state,
    "artifact_risk": current_risk,
    "artifact_reason": current_reason,
})

无效 evidence 分支也一样：

confidences_all.append({
    "event_id": row["event_id"],
    "state_confidence": 1.0,
    "state": current_state,
    "artifact_risk": current_risk,
    "artifact_reason": current_reason,
})

验收标准：

新增测试：

def test_empirical_bayes_preserves_artifact_reason():
    classifier = StateClassifier(use_empirical_bayes=True)
    states = classifier.classify(evidence_df_with_low_quality_event)
    assert "artifact_reason" in states.columns
    row = states.loc[states["event_id"] == "e_low_quality"].iloc[0]
    assert row["artifact_reason"] != ""

⸻

3. integration test 严格检查三种核心状态

位置：

tests/test_integration.py

问题：

现在 integration test 对 concordant 已经比较严格，但对 chromatin_primed 和 rna_only 仍然比较宽松，例如允许 primed 不是 null/rna_only，或者允许 RNA-only 被判成 discordant_opposite。这样不能证明工具真的恢复了 ground truth state。 ￼

怎么改：

把测试目标改成明确断言：

assert conc_state == "concordant"
assert primed_state == "chromatin_primed"
assert rna_state == "rna_only"

如果当前 synthetic data 不稳定，不要放宽断言，而是调整模拟数据：

增加 sample 数量
增加 effect size
降低随机噪声
固定 random seed
让 primed 的 RNA 在两组完全同分布
让 rna_only 的 ATAC 在两组完全同分布

当前测试里如果 fdr_threshold=0.5 太松，建议改回：

fdr_threshold=0.1

或在模拟数据里把效应做得更强。

验收标准：

至少有一个 integration test 同时通过：

assert conc_state == "concordant"
assert primed_state == "chromatin_primed"
assert rna_state == "rna_only"

⸻

4. 删除 artifact_like 残留

位置：

modes/report.py
tests/
README.md

问题：

当前主状态体系已经改成：

concordant
chromatin_primed
rna_only
discordant_opposite
null

artifact_like 不应该再作为主 state 出现。states.py 里 biological states 已经不含 artifact_like，这是对的；但 report 或测试里如果还有残留，会让体系混乱。 ￼

怎么改：

全仓库搜索：

grep -R "artifact_like" .

主逻辑、CSS、测试断言里都不要再出现。

report.py 里如果有：

ALLOWED_STATES = {
    "concordant",
    "chromatin_primed",
    "rna_only",
    "discordant_opposite",
    "artifact_like",
    "null",
}

改成：

ALLOWED_STATES = {
    "concordant",
    "chromatin_primed",
    "rna_only",
    "discordant_opposite",
    "null",
}

如果 CSS 里还有：

.artifact_like

删除它。

保留 artifact risk 的 CSS：

.artifact-risk-low
.artifact-risk-medium
.artifact-risk-high

验收标准：

grep -R "artifact_like" modes tests README.md

主代码和测试中应无结果。

⸻

5. min_event_fdr 改名

位置：

modes/core.py
tests/
README.md

问题：

如果函数参数叫：

min_event_fdr

但逻辑是：

df = df[df["event_fdr"] < min_event_fdr]

那它实际意思是“最大允许 FDR 阈值”，不是“最小 event FDR”。

怎么改：

把参数改成：

max_event_fdr: Optional[float] = None

或者：

event_fdr_threshold: Optional[float] = None

推荐：

max_event_fdr

示例：

def filter(
    self,
    states: Optional[List[str]] = None,
    min_confidence: Optional[float] = None,
    max_event_fdr: Optional[float] = None,
    exclude_high_artifact: bool = False,
):
    df = self.event_table.copy()
    if states is not None:
        df = df[df["state"].isin(states)]
    if min_confidence is not None and "state_confidence" in df.columns:
        df = df[df["state_confidence"] >= min_confidence]
    if max_event_fdr is not None and "event_fdr" in df.columns:
        df = df[df["event_fdr"] <= max_event_fdr]
    if exclude_high_artifact and "artifact_risk" in df.columns:
        df = df[df["artifact_risk"] != "high"]
    return MoDESResult(
        event_table=df,
        ...
    )

验收标准：

测试：

filtered = result.filter(max_event_fdr=0.1)
assert (filtered.event_table["event_fdr"] <= 0.1).all()

旧参数 min_event_fdr 不再出现。

⸻

第二优先级：建议这一轮一起改

6. events.py 加 gene / peak 坐标解析覆盖率 warning

位置：

modes/events.py
tests/test_events.py

问题：

如果没有 external_links，工具会根据 gene 坐标和 peak 坐标生成 candidate events。现在完全没有 events 时会报错，这已经很好。但还有一种情况：

1000 个 genes
只有 100 个能解析坐标
最后生成了一些 events
所以不会报错
但 900 个 genes 静默缺失

用户可能不知道大量基因或 peaks 没参与 event generation。

怎么改：

在 coordinate-based event generation 中统计：

无法解析坐标的 genes 数量
无法解析 interval 的 peaks 数量

示例：

import warnings
n_total_peaks = len(peak_names)
n_unknown_peaks = (peak_df["chr"] == "unknown").sum()
if n_unknown_peaks > 0:
    warnings.warn(
        f"{n_unknown_peaks}/{n_total_peaks} peaks could not be parsed as genomic intervals. "
        "They may be excluded from coordinate-based event generation.",
        UserWarning,
    )

gene 也类似：

n_total_genes = len(gene_names)
n_missing_genes = 0
for gene in gene_names:
    tss_info = self._tss_map.get(gene)
    if tss_info is None:
        n_missing_genes += 1
        continue
    chrom = tss_info[1]
    if chrom in {"", "unknown", None}:
        n_missing_genes += 1
        continue

最后：

if n_missing_genes > 0:
    warnings.warn(
        f"{n_missing_genes}/{n_total_genes} genes have no genomic coordinates "
        "and may be excluded from coordinate-based event generation. "
        "Provide genome_annotation, tss_map, or external_links for better coverage.",
        UserWarning,
    )

验收标准：

新增测试：

def test_unannotated_genes_warn():
    with pytest.warns(UserWarning, match="genes have no genomic coordinates"):
        builder.build(
            gene_names=["STAT1", "GZMB"],
            peak_names=["chr1:100-200"],
            external_links=None,
            genome_annotation=None,
            tss_map=None,
        )

⸻

7. Poisson fallback 不收敛时给 warning 或继续 fallback

位置：

modes/effects.py
tests/test_effects.py

问题：

现在 Poisson fallback 分支会标记：

model_used = poisson_fallback
family = poisson

这个已经比之前正确。但如果 Poisson 也没有收敛，当前代码可能仍然返回结果，只是在 summary 里显示 converged=False。这可能让用户误用不收敛的 coefficient。

怎么改，二选一：

方案 A：严格版，未收敛就继续 fallback

result3 = model3.fit(...)
if getattr(result3, "converged", False):
    result3._modes_model_used = "poisson_fallback"
    result3._modes_family = "poisson"
    result3._modes_alpha = None
    result3._modes_alpha_estimated = False
    result3._modes_dropped_covariates = False
    return result3

如果没收敛，继续走 simplified fallback。

方案 B：宽松版，保留但写 warning

result3._modes_warning = ""
if not getattr(result3, "converged", False):
    result3._modes_warning = (
        "Poisson fallback did not converge; coefficients are returned with caution."
    )

然后 model_summary 里加入：

"warning": getattr(result, "_modes_warning", "")

建议：

MVP 用方案 B 即可，因为有些真实 GLM 不收敛但仍能给近似结果。关键是不要静默。

验收标准：

model_diagnostics.tsv 中如果 Poisson fallback 未收敛，应出现：

model_used = poisson_fallback
converged = False
warning = Poisson fallback did not converge...

⸻

8. 明确 model_diagnostics.tsv 是否包含 conditional decomposition

位置：

modes/core.py
modes/decompose.py
README.md

问题：

现在 model_diagnostics.tsv 主要来自 ATAC 和 RNA 主效应模型；但 MoDES 还有 RNA_after_ATAC 条件分解。用户看到 model_diagnostics.tsv 可能会以为里面包含所有模型诊断。

两种改法：

方案 A：文档说明当前只包含 marginal GLM

README 里加一句：

model_diagnostics.tsv currently reports marginal ATAC and RNA GLM diagnostics.
Conditional RNA-after-ATAC diagnostics are summarized in event_table.tsv and will be expanded in a future version.

方案 B：把 conditional model 也加入 diagnostics

在 decompose.py 的 conditional result 中增加：

conditional_model_used
conditional_family
conditional_converged
conditional_warning

在 core.py 的 _build_model_diagnostics() 中加入：

for _, row in self.conditional_effects.iterrows():
    diagnostic_rows.append({
        "feature_id": row["event_id"],
        "modality": "RNA_after_ATAC",
        "model_used": row.get("model_used", "conditional_nb"),
        "family": row.get("family", "negative_binomial"),
        "alpha": row.get("alpha", None),
        "alpha_estimated": row.get("alpha_estimated", False),
        "converged": row.get("convergence", False),
        "dropped_covariates": row.get("dropped_covariates", False),
        "warning": row.get("warning", ""),
    })

建议：

先做方案 A，后续再补方案 B。

验收标准：

README 和真实输出一致，不夸大 model_diagnostics.tsv 的覆盖范围。

⸻

9. README 表格重新格式化成标准 Markdown

位置：

README.md

问题：

README 当前内容虽然 GitHub 页面能显示，但 raw 文件中很多表格像普通文本，不是真正 Markdown table。比如“输出文件”和“event_table.tsv 字段”现在看起来是：

文件 说明
event_table.tsv 主输出表
...
字段 说明
event_id 事件唯一标识
...

这在渲染、复制和维护时不如标准表格清晰。 ￼

怎么改：

改成标准 Markdown table。

例如：

### 输出文件
| 文件 | 说明 |
|---|---|
| `event_table.tsv` | 主输出表 |
| `event_state_confidence.tsv` | 状态置信度 |
| `event_layer_effects.tsv` | 每层效应大小 |
| `event_evidence_vectors.tsv` | 证据向量 |
| `model_diagnostics.tsv` | 模型诊断信息 |
| `run_params.tsv` | 运行参数 |

字段表：

#### `event_table.tsv` 字段
| 字段 | 说明 |
|---|---|
| `event_id` | 事件唯一标识 |
| `gene` | target gene |
| `peak_id` | regulatory peak |
| `state` | biological state |
| `state_confidence` | 状态置信度，范围 0 到 1 |
| `artifact_risk` | 技术伪影风险，取值 `low` / `medium` / `high` |
| `artifact_reason` | 伪影原因，分号分隔 |
| `event_pval` | event-level p-value |
| `event_fdr` | event-level BH-corrected FDR |
| `atac_coef` / `atac_pval` / `atac_fdr` | ATAC 效应估计 |
| `rna_coef` / `rna_pval` / `rna_fdr` | RNA 效应估计 |
| `rna_after_atac_coef` / `rna_after_atac_pval` / `rna_after_atac_fdr` | 控制 linked ATAC peak 后的条件效应 |

验收标准：

GitHub README 渲染成真正的表格，不是普通文本块。

⸻

第三优先级：可做可不做，但建议补上

10. filter() 增加 artifact risk 过滤

位置：

modes/core.py
tests/test_core.py

问题：

既然现在有 artifact_risk，用户很可能想筛掉高风险事件。

怎么改：

给 MoDESResult.filter() 增加：

exclude_high_artifact: bool = False
max_artifact_risk: Optional[str] = None

简单版：

if exclude_high_artifact and "artifact_risk" in df.columns:
    df = df[df["artifact_risk"] != "high"]

更通用版：

risk_order = {"low": 0, "medium": 1, "high": 2}
if max_artifact_risk is not None and "artifact_risk" in df.columns:
    max_rank = risk_order[max_artifact_risk]
    df = df[
        df["artifact_risk"].map(risk_order).fillna(0) <= max_rank
    ]

验收标准：

filtered = result.filter(exclude_high_artifact=True)
assert "high" not in set(filtered.event_table["artifact_risk"])

⸻

11. GraphML 里加入 artifact_risk、event_fdr

位置：

modes/core.py

问题：

如果用户把结果导出到 Cytoscape / Gephi，最好能在网络边属性里看到：

state
state_confidence
artifact_risk
event_fdr

怎么改：

在 to_graphml() 的 edge attributes 里加入：

artifact_risk=row.get("artifact_risk", "low"),
artifact_reason=row.get("artifact_reason", ""),
event_pval=float(row.get("event_pval", 1.0)),
event_fdr=float(row.get("event_fdr", 1.0)),
state_confidence=float(row.get("state_confidence", np.nan)),

验收标准：

导出的 GraphML edge 属性中包含：

artifact_risk
event_fdr
state_confidence

⸻

12. report 中按 state_confidence 排序

位置：

modes/report.py

问题：

如果你完成第 1 条，把主字段改成 state_confidence，report 里所有排序和显示也要同步。

怎么改：

top_events = results.event_table.nlargest(50, "state_confidence")

display columns：

display_cols = [
    "event_id",
    "gene",
    "peak_id",
    "state",
    "state_confidence",
    "artifact_risk",
    "event_fdr",
    "atac_coef",
    "rna_coef",
    "atac_fdr",
    "rna_fdr",
]

验收标准：

HTML report 不再引用 confidence。

⸻

建议修改顺序

你这轮按这个顺序改：

1. confidence → state_confidence，全仓库统一
2. EB 阶段保留 artifact_reason
3. integration test 严格断言 chromatin_primed 和 rna_only
4. 删除 artifact_like 残留
5. min_event_fdr → max_event_fdr
6. events.py 加 gene/peak 坐标 coverage warning
7. Poisson fallback 未收敛时写 warning
8. README 说明 model_diagnostics.tsv 当前范围
9. README 表格格式化
10. filter() 支持 artifact risk 过滤
11. GraphML 增加 artifact_risk / event_fdr / state_confidence

最少先改前 5 个：

1. state_confidence 字段统一
2. artifact_reason 在 EB 后不丢
3. integration test 严格
4. artifact_like 残留清理
5. max_event_fdr 参数命名

这 5 个改完，49a8e2e 这一版的主要不一致就基本消掉了。
