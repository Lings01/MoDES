Lings，下面是这一轮的修改 list。我按优先级排好了。你可以一条一条改。

⸻

第一优先级：必须先改

⸻

1. 把 artifact_risk 写入最终 event_table

位置：

modes/_types.py
modes/core.py
modes/report.py
tests/

问题：

现在 StateClassifier.classify() 已经生成了：

state
state_confidence
artifact_risk

但是 _assemble_results() 只把：

state
state_confidence

写入 EventResult，没有把 artifact_risk 带到最终 event_table.tsv。

这会导致 README 说有 artifact_risk，但用户实际输出里看不到。

⸻

怎么改

1.1 修改 EventResult

在 modes/_types.py 里给 EventResult 增加字段：

artifact_risk: str = "low"
artifact_reason: str = ""

例如：

@dataclass
class EventResult:
    event_id: str
    gene: str
    peak_id: str
    state: str
    state_confidence: float
    artifact_risk: str = "low"
    artifact_reason: str = ""
    atac_effect: float = np.nan
    atac_pval: float = np.nan
    atac_fdr: float = np.nan
    rna_effect: float = np.nan
    rna_pval: float = np.nan
    rna_fdr: float = np.nan
    rna_after_atac_effect: float = np.nan
    rna_after_atac_pval: float = np.nan
    rna_after_atac_fdr: float = np.nan

按你现有字段顺序合并即可，不一定照搬这个完整顺序。

⸻

1.2 修改 _event_result_columns()

在 modes/core.py 里找到事件表字段列表，加入：

"artifact_risk",
"artifact_reason",

例如：

def _event_result_columns(self) -> List[str]:
    return [
        "event_id",
        "gene",
        "peak_id",
        "state",
        "state_confidence",
        "artifact_risk",
        "artifact_reason",
        "atac_effect",
        "atac_pval",
        "atac_fdr",
        "rna_effect",
        "rna_pval",
        "rna_fdr",
        "rna_after_atac_effect",
        "rna_after_atac_pval",
        "rna_after_atac_fdr",
    ]

⸻

1.3 修改 _assemble_results()

现在你大概有类似逻辑：

state = state_row.iloc[0]["state"]
confidence = state_row.iloc[0]["state_confidence"]

改成：

state = state_row.iloc[0]["state"]
confidence = state_row.iloc[0]["state_confidence"]
artifact_risk = state_row.iloc[0].get("artifact_risk", "low")
artifact_reason = state_row.iloc[0].get("artifact_reason", "")

然后构造 EventResult 时传进去：

event_result = EventResult(
    event_id=event_id,
    gene=gene,
    peak_id=peak_id,
    state=state,
    state_confidence=confidence,
    artifact_risk=artifact_risk,
    artifact_reason=artifact_reason,
    ...
)

⸻

1.4 修改 report 显示字段

在 modes/report.py 里显示 event table 的列中加入：

"artifact_risk"

如果版面允许，也加入：

"artifact_reason"

⸻

验收标准

新增测试：

def test_event_table_contains_artifact_risk():
    result = modes.run()
    assert "artifact_risk" in result.event_table.columns
    assert "artifact_reason" in result.event_table.columns

⸻

2. 删除 artifact_like 作为主 state 的逻辑

位置：

modes/states.py
tests/test_states.py
README.md
modes/report.py

问题：

你现在的设计已经变成两层：

biological state
artifact_risk

这是对的。

但是 states.py 里仍然可能返回：

state = artifact_like

这会导致状态体系混乱。

最终建议是：

state:
  concordant
  chromatin_primed
  rna_only
  discordant_opposite
  null
artifact_risk:
  low
  medium
  high

不要再让 artifact_like 成为主状态。

⸻

怎么改

在 modes/states.py 里，把主状态分类逻辑改成：

if atac_sig and rna_sig and same_dir:
    state = "concordant"
elif atac_sig and rna_sig and not same_dir:
    state = "discordant_opposite"
elif atac_sig and not rna_sig:
    state = "chromatin_primed"
elif not atac_sig and rna_sig:
    state = "rna_only"
else:
    state = "null"

单独计算 artifact：

artifact_risk, artifact_reason = self._compute_artifact_risk(row)

例如：

def _compute_artifact_risk(self, row: pd.Series) -> tuple[str, str]:
    reasons = []
    quality_score = row.get("quality_score", 1.0)
    atac_sig = row.get("atac_fdr", 1.0) < self.fdr_threshold
    rna_sig = row.get("rna_fdr", 1.0) < self.fdr_threshold
    if quality_score < self.low_quality_threshold:
        reasons.append("low_quality_score")
    if quality_score < self.low_quality_threshold and (atac_sig ^ rna_sig):
        reasons.append("single_modality_low_quality")
    if len(reasons) == 0:
        return "low", ""
    if "single_modality_low_quality" in reasons:
        return "high", ";".join(reasons)
    return "medium", ";".join(reasons)

如果你现在已有 _compute_artifact_risk()，就不要新建，直接调整现有函数。

⸻

测试要改

原来如果测试允许：

assert state in {"artifact_like", "chromatin_primed"}

改成：

assert row["state"] == "chromatin_primed"
assert row["artifact_risk"] == "high"

或者 RNA-only 情况：

assert row["state"] == "rna_only"
assert row["artifact_risk"] == "high"

⸻

验收标准

全仓库搜索：

grep -R "artifact_like" .

理想情况下只允许在 changelog 或旧版本说明里出现，不应在主逻辑、主状态列表、测试断言、CSS class 中继续出现。

⸻

3. 修 Poisson fallback 被误标成 NB 的问题

位置：

modes/effects.py

问题：

_safe_fit_nb_glm() 里有多层 fallback：

NB default
NB alpha=1
Poisson fallback
simplified NB fallback

但是外层 model_summary 可能仍然统一标成：

family = negative_binomial
model_used = nb_fixed_alpha

如果实际用了 Poisson，却报告成 NB，这是不准确的。

⸻

怎么改

在 _safe_fit_nb_glm() 每个成功分支上给 result 挂诊断属性。

NB default 分支

result._modes_model_used = "nb_default_alpha"
result._modes_family = "negative_binomial"
result._modes_alpha = None
result._modes_alpha_estimated = False
result._modes_dropped_covariates = False

NB alpha=1 分支

result._modes_model_used = "nb_fixed_alpha"
result._modes_family = "negative_binomial"
result._modes_alpha = 1.0
result._modes_alpha_estimated = False
result._modes_dropped_covariates = False

Poisson fallback 分支

result._modes_model_used = "poisson_fallback"
result._modes_family = "poisson"
result._modes_alpha = None
result._modes_alpha_estimated = False
result._modes_dropped_covariates = False

simplified fallback 分支

result._modes_model_used = "nb_simple_fallback"
result._modes_family = "negative_binomial"
result._modes_alpha = 1.0
result._modes_alpha_estimated = False
result._modes_dropped_covariates = True

⸻

然后在 _fit_nb_glm() 里生成 summary

model_summary = {
    "family": getattr(result, "_modes_family", "unknown"),
    "model_used": getattr(result, "_modes_model_used", "unknown"),
    "alpha": getattr(result, "_modes_alpha", None),
    "alpha_estimated": getattr(result, "_modes_alpha_estimated", False),
    "dropped_covariates": getattr(result, "_modes_dropped_covariates", False),
    "converged": bool(getattr(result, "converged", True)),
}

⸻

验收标准

每个 ModalityEffect.model_summary 都能准确说明：

family
model_used
alpha
alpha_estimated
dropped_covariates
converged

尤其是 Poisson fallback 时必须显示：

model_used = poisson_fallback
family = poisson

⸻

4. 无效 coefficient 早退时不要丢 model_summary

位置：

modes/effects.py

问题：

现在如果出现：

np.isnan(coef)
np.isnan(se)
se <= 0

你会提前返回一个 failed effect，但这个分支可能没有带上 model_summary。

这样用户只知道失败了，不知道：

用的什么模型
是不是 Poisson fallback
是不是 simplified fallback
是不是 dropped covariates

⸻

怎么改

在无效 coefficient 分支里也构造 summary：

if np.isnan(coef) or np.isnan(se) or se <= 0:
    model_summary = {
        "family": getattr(result, "_modes_family", "unknown"),
        "model_used": getattr(result, "_modes_model_used", "unknown"),
        "alpha": getattr(result, "_modes_alpha", None),
        "alpha_estimated": getattr(result, "_modes_alpha_estimated", False),
        "dropped_covariates": getattr(result, "_modes_dropped_covariates", False),
        "converged": False,
        "warning": "Invalid coefficient or standard error.",
    }
    return ModalityEffect(
        feature_id=feature_id,
        coefficient=np.nan,
        standard_error=np.nan,
        pvalue=1.0,
        fdr=1.0,
        direction=0,
        convergence=False,
        model_summary=model_summary,
    )

字段名按你现有 ModalityEffect 定义调整。

⸻

验收标准

失败的 feature 也有：

model_summary["model_used"]
model_summary["warning"]
model_summary["converged"] == False

⸻

5. decompose.py 不要吞掉 NotImplementedError / ValueError

位置：

modes/decompose.py

问题：

_fit_conditional() 里有宽泛的：

except Exception:
    return self._null_conditional(...)

这可能吞掉本来应该暴露给用户的错误，例如：

multi-class condition not supported
rank deficient design
invalid design matrix

这些不应该被变成 null result。

⸻

怎么改

把 validation error 放出来：

try:
    ...
except (NotImplementedError, ValueError):
    raise
except Exception as e:
    return self._null_conditional(
        event_id=event_id,
        gene=gene,
        peak_id=peak_id,
        reason=str(e),
    )

或者更稳：

X, design_info = self._build_design_matrix(...)

放到 try 外面，避免设计矩阵错误被吞掉。

⸻

验收标准

多类别 condition 时：

with pytest.raises(NotImplementedError):
    decomposer.decompose(...)

rank deficient design 时：

with pytest.raises(ValueError, match="rank deficient"):
    decomposer.decompose(...)

不能静默返回 null conditional result。

⸻

6. 修 report summary card 的 label/value 顺序

位置：

modes/report.py

问题：

你现在可能是：

cards_data = [
    ("Total Events", str(len(self.result.event_table))),
    ...
]
for value, label in cards_data:
    ...

这会把 "Total Events" 当成 value，把数字当成 label。

⸻

怎么改

改成：

for label, value in cards_data:
    summary_cards += f"""
    <div class="card">
        <div class="value">{_esc(value)}</div>
        <div class="label">{_esc(label)}</div>
    </div>
    """

⸻

验收标准

报告里应该显示：

大数字 / 统计值
下面是 Total Events / Concordant / Primed 等标签

不能反过来。

⸻

7. integration test 改成严格断言具体状态

位置：

tests/test_integration.py

问题：

当前测试如果允许：

assert conc_state in {"concordant", "discordant_opposite"}
assert primed_state in {"chromatin_primed", "rna_only", "concordant", "discordant_opposite"}
assert rna_state in {"rna_only", "discordant_opposite"}

这个太宽松，不能证明方法真的能恢复 ground truth。

⸻

怎么改

把 synthetic data 设计得更强，然后断言写死：

assert conc_state == "concordant"
assert primed_state == "chromatin_primed"
assert rna_state == "rna_only"

如果随机性导致不稳定，优先调整模拟数据：

增加 sample size
增加 effect size
降低噪声
固定 random seed
降低 threshold 的边界敏感性

不要用过宽断言让测试通过。

⸻

验收标准

至少有一个 integration test 明确检查：

concordant → concordant
chromatin_primed → chromatin_primed
rna_only → rna_only

⸻

第二优先级：建议这一轮一起改

⸻

8. 给无法解析坐标的 genes / peaks 加 warning

位置：

modes/events.py

问题：

现在如果没有 external links，代码会尝试从 gene name 或 annotation 推坐标。完全生成不了 events 时已经会报错，这是好的。

但是还有一种情况：

1000 个 genes
只有 120 个能解析坐标
最后生成了一些 events
所以不会报错
但 880 个 genes 被静默丢掉

这会让用户以为全部基因都参与了 event generation。

⸻

怎么改

在构建 TSS map 或解析 peak interval 时统计失败数量。

例如：

n_total_genes = len(gene_names)
n_missing_genes = 0
for g in gene_names:
    if g not in tss_map or tss_map[g][1] == "":
        n_missing_genes += 1

然后：

if n_missing_genes > 0:
    warnings.warn(
        f"{n_missing_genes}/{n_total_genes} genes have no genomic coordinates "
        "and may be excluded from coordinate-based event generation. "
        "Provide genome_annotation, tss_map, or external_links for better coverage.",
        UserWarning,
    )

peak 也类似：

if n_unparsed_peaks > 0:
    warnings.warn(
        f"{n_unparsed_peaks}/{n_total_peaks} peaks could not be parsed as genomic intervals.",
        UserWarning,
    )

⸻

验收标准

新增测试：

def test_unannotated_genes_warn():
    with pytest.warns(UserWarning, match="genes have no genomic coordinates"):
        builder.build(...)

⸻

9. 清理 artifact_like 的 CSS、测试、文档残留

位置：

modes/report.py
README.md
tests/
docs if any

问题：

如果你采用：

state + artifact_risk

双层体系，就不应该再给 artifact_like 一个主状态 CSS。

⸻

怎么改

删除或停用：

.artifact_like

新增：

.artifact-risk-low {
    color: #2e7d32;
}
.artifact-risk-medium {
    color: #ef6c00;
}
.artifact-risk-high {
    color: #c62828;
    font-weight: bold;
}

报告表格中 artifact_risk 可以这样渲染：

risk = str(row.get("artifact_risk", "low"))
risk_class = f"artifact-risk-{risk}" if risk in {"low", "medium", "high"} else "artifact-risk-low"

⸻

验收标准

全仓库搜索：

grep -R "artifact_like" .

主逻辑中不再出现。

⸻

10. 增加 event-level p-value / FDR

位置：

modes/core.py
modes/states.py
modes/_types.py
tests/

问题：

现在主要有：

atac_fdr
rna_fdr
rna_after_atac_fdr
state_confidence

但最终 event state 本身没有一个统一的：

event_pval
event_fdr

一个 gene 连很多 peaks 时，RNA FDR 会被复制到很多 events，容易让 event 数量看起来膨胀。

⸻

怎么改，最小版本

在 _assemble_results() 里，根据 state 计算 event p-value：

if state == "concordant":
    event_pval = max(atac_pval, rna_pval)
elif state == "discordant_opposite":
    event_pval = max(atac_pval, rna_pval)
elif state == "chromatin_primed":
    event_pval = atac_pval
elif state == "rna_only":
    event_pval = rna_pval
else:
    event_pval = 1.0

然后对所有 events 的 event_pval 做 BH：

event_fdr = benjamini_hochberg(event_pvals)

⸻

修改字段

EventResult 增加：

event_pval: float = 1.0
event_fdr: float = 1.0

event_table 输出增加：

event_pval
event_fdr

⸻

验收标准

assert "event_pval" in result.event_table.columns
assert "event_fdr" in result.event_table.columns
assert result.event_table["event_fdr"].between(0, 1).all()

⸻

11. 增加 model_diagnostics.tsv

位置：

modes/core.py
modes/effects.py
modes/decompose.py
tests/

问题：

现在 to_tsv() 输出主要是：

event_table.tsv
event_state_confidence.tsv
event_layer_effects.tsv
event_evidence_vectors.tsv
run_params.json

但是用户还不知道每个 feature 最终用了什么模型。

⸻

怎么改

新增一个 diagnostics table：

model_diagnostics.tsv

字段建议：

feature_id
modality
model_used
family
alpha
alpha_estimated
converged
dropped_covariates
warning

⸻

数据来源

从每个 ModalityEffect.model_summary 里提取。

例如：

diagnostic_rows = []
for effect in self.atac_effects.values():
    summary = effect.model_summary or {}
    diagnostic_rows.append({
        "feature_id": effect.feature_id,
        "modality": "ATAC",
        "model_used": summary.get("model_used", "unknown"),
        "family": summary.get("family", "unknown"),
        "alpha": summary.get("alpha", None),
        "alpha_estimated": summary.get("alpha_estimated", False),
        "converged": summary.get("converged", effect.convergence),
        "dropped_covariates": summary.get("dropped_covariates", False),
        "warning": summary.get("warning", ""),
    })

RNA 同理。

⸻

修改 MoDESResult

增加属性：

model_diagnostics: pd.DataFrame

to_tsv() 里写出：

self.model_diagnostics.to_csv(outdir / "model_diagnostics.tsv", sep="\t", index=False)

⸻

验收标准

运行：

result.to_tsv("output")

后出现：

output/model_diagnostics.tsv

并包含：

model_used
family
converged
dropped_covariates

⸻

第三优先级：增强可解释性和稳定性

⸻

12. artifact_reason 做成可解释字段

位置：

modes/states.py
modes/core.py
modes/_types.py
modes/report.py

问题：

只有：

artifact_risk = high

还不够。用户会问：

为什么 high？

⸻

怎么改

让 _compute_artifact_risk() 返回：

artifact_risk, artifact_reason

artifact_reason 可以是分号分隔字符串：

low_quality_score
single_modality_low_quality
low_atac_depth
low_rna_depth
batch_associated
library_size_outlier

当前可以先实现最基础的两个：

low_quality_score
single_modality_low_quality

⸻

示例

def _compute_artifact_risk(self, row: pd.Series) -> tuple[str, str]:
    reasons = []
    quality_score = float(row.get("quality_score", 1.0))
    atac_sig = float(row.get("atac_fdr", 1.0)) < self.fdr_threshold
    rna_sig = float(row.get("rna_fdr", 1.0)) < self.fdr_threshold
    if quality_score < self.low_quality_threshold:
        reasons.append("low_quality_score")
    if quality_score < self.low_quality_threshold and (atac_sig ^ rna_sig):
        reasons.append("single_modality_low_quality")
    if not reasons:
        return "low", ""
    if "single_modality_low_quality" in reasons:
        return "high", ";".join(reasons)
    return "medium", ";".join(reasons)

⸻

验收标准

低质量单模态事件输出：

artifact_risk = high
artifact_reason = low_quality_score;single_modality_low_quality

⸻

13. 把 BIOLOGICAL_STATES 和状态分类规则完全统一

位置：

modes/states.py

问题：

如果你定义了：

BIOLOGICAL_STATES = [...]

就必须保证 _rule_based_classify() 只返回里面的状态。

⸻

怎么改

定义：

BIOLOGICAL_STATES = {
    "concordant",
    "chromatin_primed",
    "rna_only",
    "discordant_opposite",
    "null",
}

在 classify 后加防御检查：

if state not in BIOLOGICAL_STATES:
    raise ValueError(f"Invalid biological state returned: {state}")

⸻

验收标准

assert set(states["state"]).issubset(BIOLOGICAL_STATES)

⸻

14. 把 README 的输出表字段更新

位置：

README.md

问题：

你增加 artifact_risk、artifact_reason、event_pval、event_fdr 后，README 的示例输出表也要同步。

⸻

怎么改

README 里的 output table 示例改成：

event_id
gene
peak_id
state
state_confidence
artifact_risk
artifact_reason
event_pval
event_fdr
atac_effect
rna_effect
rna_after_atac_effect

并解释：

state: inferred biological event state
artifact_risk: technical-risk flag, not a biological state
artifact_reason: semicolon-separated reasons for artifact risk
event_fdr: event-level multiple-testing-adjusted significance

⸻

验收标准

README 和真实 event_table.tsv 字段一致。

⸻

15. filter() 支持 artifact risk 过滤

位置：

modes/core.py

问题：

现在 MoDESResult.filter() 支持按 confidence 和 state 过滤，但最好支持去掉高 artifact risk。

⸻

怎么改

给 filter() 增加参数：

max_artifact_risk: Optional[str] = None
exclude_high_artifact: bool = False

简单实现：

if exclude_high_artifact and "artifact_risk" in filtered.columns:
    filtered = filtered[filtered["artifact_risk"] != "high"]

或者更通用：

risk_order = {"low": 0, "medium": 1, "high": 2}
if max_artifact_risk is not None:
    max_rank = risk_order[max_artifact_risk]
    filtered = filtered[
        filtered["artifact_risk"].map(risk_order).fillna(0) <= max_rank
    ]

⸻

验收标准

filtered = result.filter(exclude_high_artifact=True)
assert "high" not in set(filtered.event_table["artifact_risk"])

⸻

16. to_report() 里增加 artifact risk summary

位置：

modes/report.py

问题：

报告现在主要总结 state distribution，最好也总结 artifact risk distribution。

⸻

怎么改

增加一个 summary section：

Artifact risk distribution:
low: N
medium: N
high: N

如果有 artifact_reason，再加 top reasons：

Top artifact reasons:
low_quality_score
single_modality_low_quality

⸻

验收标准

HTML report 能看到：

Artifact risk
low / medium / high counts

⸻

推荐修改顺序

你这一轮按这个顺序改：

1. artifact_risk 写入 EventResult 和 event_table
2. 删除 artifact_like 作为主 state
3. artifact_reason 一起写入
4. 修 Poisson fallback model_summary
5. 修 invalid coefficient 早退 model_summary
6. decompose.py 不吞 NotImplementedError / ValueError
7. 修 report summary card 顺序
8. integration test 改成严格状态断言
9. 清理 artifact_like 残留
10. README 同步输出字段

如果你时间有限，最少先改这 5 个：

1. artifact_risk 输出到 event_table
2. artifact_like 不再作为 state
3. Poisson fallback 不再误标 NB
4. decompose.py 不吞关键错误
5. integration test 严格检查 concordant / primed / rna_only

这 5 个改完，MoDES 的状态体系和输出就会一致很多。
