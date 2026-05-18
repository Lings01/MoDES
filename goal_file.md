Lings，这份三审意见虽然很狠，但它已经把核心问题说透了：

MoDES 现在看起来像 multi-modal software interface，但还不像 multi-modal statistical inference engine。

也就是说，审稿人已经承认你做了不少东西：extra modalities 进了 effect、evidence、state 流程；版本也统一了；event_fdr 也降级成 heuristic。
但他们现在抓的是更深的一层：

多模态只是加进规则系统；
统计证据没有真正多模态化；
event_pval / event_fdr 仍然 RNA+ATAC-centric；
conditional decomposition 仍然是 RNA after ATAC；
spatial / dynamic 还是 helper，不是主流程；
benchmark 仍然按规则生成、按规则识别。

如果你想保住 multi-omics 这个名字，下一步不能再只是“补功能”。你需要把 MoDES 从：

RNA+ATAC core + extra modality rules

重构成：

modality-generic evidence engine + state grammar + state-specific evidence scoring + calibration

下面是我建议的完整修复方案。

⸻

0. 先定战略：不要再写“inference framework”，改成“evidence framework”

如果你继续叫：

validated multi-omics regulatory event inference framework

还会被拒。

建议主标题改成：

MoDES: a multi-omics evidence framework for regulatory event-state annotation

或者更稳：

MoDES: multi-omics evidence scoring and state annotation for candidate regulatory events

核心是把 inference 换成：

evidence scoring
state annotation
candidate prioritization

这样你还能保留 multi-omics，但不会被要求证明完整因果统计模型。

⸻

1. 最大硬伤：event_pval / event_fdr 与 extra modalities 脱节

这是这轮 review 最致命的一点。

现在的问题是：

state 由 H3K27ac / protein / spatial 触发；
但 event_pval 仍然主要用 ATAC/RNA p-value。

这在多模态中不能接受。

必须改成 state-specific evidence p-value

每个 state 必须声明：

这个 state 由哪些 modality 的哪些方向支持？

然后 event p-value 必须来自这些支持该 state 的 modality，而不是固定 ATAC/RNA。

⸻

1.1 新增 StateRule

新增文件：

modes/modalities/state_rules.py

定义：

from dataclasses import dataclass
from typing import Optional, Sequence
@dataclass(frozen=True)
class RequiredEvidence:
    modality: str
    direction: int  # +1, -1
    role: Optional[str] = None
    target: Optional[str] = None
@dataclass(frozen=True)
class NeutralEvidence:
    modality: str
    role: Optional[str] = None
@dataclass(frozen=True)
class ForbiddenEvidence:
    modality: str
    direction: int
    role: Optional[str] = None
@dataclass(frozen=True)
class StateRule:
    name: str
    required: Sequence[RequiredEvidence]
    neutral: Sequence[NeutralEvidence] = ()
    forbidden: Sequence[ForbiddenEvidence] = ()
    description: str = ""
    interpretation_strength: str = "association"  # association / hypothesis / causal_not_claimed

⸻

1.2 示例 state rules

RNA+ATAC

CONCORDANT = StateRule(
    name="concordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
    ],
    description="ATAC and RNA change in the same direction.",
)
CHROMATIN_PRIMED = StateRule(
    name="chromatin_primed",
    required=[
        RequiredEvidence("atac", +1),
    ],
    neutral=[
        NeutralEvidence("rna"),
    ],
    description="ATAC changes while RNA does not show significant change.",
)
RNA_ONLY = StateRule(
    name="rna_only",
    required=[
        RequiredEvidence("rna", +1),
    ],
    neutral=[
        NeutralEvidence("atac"),
    ],
)

CUT&Tag activating mark

ACTIVE_MARK_CONCORDANT = StateRule(
    name="active_mark_concordant",
    required=[
        RequiredEvidence("cuttag_h3k27ac", +1),
        RequiredEvidence("rna", +1),
    ],
    description="Activating chromatin mark and RNA change concordantly.",
)
ACTIVE_MARK_PRIMED = StateRule(
    name="active_mark_primed",
    required=[
        RequiredEvidence("cuttag_h3k27ac", +1),
    ],
    neutral=[
        NeutralEvidence("rna"),
    ],
)

Repressive mark

DEREPRESSED = StateRule(
    name="derepressed",
    required=[
        RequiredEvidence("cuttag_h3k27me3", -1),
        RequiredEvidence("rna", +1),
    ],
)

Protein

FULL_ACTIVATION = StateRule(
    name="full_activation",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
        RequiredEvidence("protein", +1),
    ],
)
PROTEIN_BUFFERED = StateRule(
    name="protein_buffered",
    required=[
        RequiredEvidence("rna", +1),
    ],
    neutral=[
        NeutralEvidence("protein"),
    ],
)
PROTEIN_MEMORY = StateRule(
    name="protein_memory",
    required=[
        RequiredEvidence("protein", +1),
    ],
    neutral=[
        NeutralEvidence("rna"),
    ],
)

⸻

1.3 新的 state p-value 计算

每个 state 用它自己的 required modalities 计算 p-value。

对 required evidence：

required modality 都要显著

则可以用 intersection test：

state_support_pval = max(directed_pvals_of_required_modalities)

例如：

active_mark_primed:
  required = H3K27ac ↑
  state_support_pval = p_H3K27ac_up
full_activation:
  required = ATAC ↑, RNA ↑, protein ↑
  state_support_pval = max(p_ATAC_up, p_RNA_up, p_protein_up)
protein_memory:
  required = protein ↑
  state_support_pval = p_protein_up
mark_only:
  required = H3K27ac ↑
  state_support_pval = p_H3K27ac_up

这样就不会出现：

state 由 H3K27ac 触发，但 event_pval 用 ATAC p-value

这个硬伤。

⸻

1.4 directed p-value

不要只用普通 p-value，还要看方向。

def directed_pvalue(pval: float, coef: float, expected_direction: int) -> float:
    if expected_direction == 0:
        return pval
    if coef * expected_direction > 0:
        return min(pval / 2.0, 1.0)
    return 1.0

注意这仍然不是严格单侧检验，但比方向和 p-value 脱节更合理。文档要说明这是 directional evidence score。

⸻

1.5 输出字段改名

把：

event_pval
event_fdr

改成或新增：

state_support_pval
state_support_qval
supporting_modalities

保留旧字段也可以，但不建议作为主解释。

主表建议：

event_id
state
state_assignment_score
state_support_pval
state_support_qval
supporting_modalities
neutral_modalities
conflicting_modalities
artifact_risk
artifact_reason

这样 reviewer 第四条会被直接解决。

⸻

2. 第二硬伤：多模态 evidence 不应该塞进 dynamic columns

reviewer 说得对：

schema frozen
但 extra modality 又动态加列

这会被抓。

必须拆成两个表

2.1 固定主表：event_table.tsv

只保留固定列：

event_id
region_id
gene
context
tf_name
link_source
link_score
state
state_assignment_score
state_support_pval
state_support_qval
supporting_modalities
neutral_modalities
conflicting_modalities
artifact_risk
artifact_reason
quality_score

这张表永远固定。

⸻

2.2 长格式证据表：event_modality_evidence.tsv

新增输出：

event_id
modality
assay
target
feature_id
role
coef
se
pval
fdr
direction
directed_pval
quality_score
model_used
converged
warning

例如：

E001  rna              RNA     IFIT3      IFIT3       transcript_output   1.2  ...
E001  atac             ATAC    .          chr1:...    accessibility       0.8  ...
E001  cuttag_h3k27ac   CUTTAG  H3K27ac    chr1:...    activating_mark     1.5  ...
E001  protein          ADT     IFIT3      IFIT3_ADT   protein_output      0.7  ...

这样你可以真正说：

MoDES supports multi-omics evidence

而不会污染主 schema。

⸻

3. 第三硬伤：StateClassifier 不能再 priority-based 返回第一个命中的 state

现在 reviewer 抓住的是：

先 epigenomic，再 protein，再 spatial，再 RA fallback
priority order 决定 biological interpretation

这确实危险。

改成所有 states 同时打分

不要：

if epi:
    return epi_state
elif protein:
    return protein_state
elif spatial:
    return spatial_state
else:
    return RA_state

而是：

candidate_states = []
for rule in state_rules:
    score = score_state(rule, evidence)
    if score.is_valid:
        candidate_states.append(score)
best = select_best_state(candidate_states)

⸻

3.1 State score 组成

每个 state 得到：

support_pval
support_qval
n_required_satisfied
n_conflicts
quality_penalty
missing_penalty
assignment_score

示例：

assignment_score = (
    evidence_strength
    * quality_score
    * conflict_penalty
    * missing_penalty
)

其中：

evidence_strength = -log10(state_support_pval)

然后：

best_state = highest assignment_score

如果两个 state 分数接近，输出：

mixed_evidence
ambiguous

而不是强行用 priority 决定。

⸻

3.2 新状态：mixed_evidence

如果同时满足：

active_mark_concordant
protein_buffered

或者：

epigenomic activating ↑
repressive mark ↑
RNA ↑

不要强制选一个。输出：

mixed_evidence

并在 conflicting_modalities 中写：

cuttag_h3k27ac;cuttag_h3k27me3;protein

这比 priority order 更可信。

⸻

4. 第四硬伤：EB refinement 需要降级或重做

reviewer 对 EB 的批评完全成立。短期最安全路线：

4.1 默认关闭 EB

use_empirical_bayes=False

如果保留：

use_empirical_bayes=True

文档写：

experimental smoothing only
not used for main results

⸻

4.2 state_confidence 改名

改成：

state_assignment_score

不要叫：

confidence
posterior
probability

代码字段可以保留向后兼容，但主输出和论文用：

state_assignment_score

⸻

4.3 小样本不要给 1.0

现在小样本 / invalid evidence 时给 state_confidence = 1.0 很危险。

改成：

if invalid_evidence:
    score = np.nan
    state = "unresolved"

如果 EB 样本不足：

state_assignment_score = rule_score

而不是：

1.0

⸻

4.4 做 calibration benchmark

输出：

confidence_bin
mean_score
empirical_accuracy
calibration_gap

否则不要给出概率式解释。

⸻

5. 第五硬伤：conditional decomposition 仍然 RNA-after-ATAC only

如果要叫 multi-omics，必须把 conditional decomposition 泛化。

新增 ConditionalModelSpec

@dataclass
class ConditionalModelSpec:
    response_modality: str
    response_feature_role: str
    conditioning_modalities: list[str]
    name: str

示例：

RNA_AFTER_ATAC = ConditionalModelSpec(
    response_modality="rna",
    response_feature_role="gene",
    conditioning_modalities=["atac"],
    name="rna_after_atac",
)
RNA_AFTER_ACTIVE_MARK = ConditionalModelSpec(
    response_modality="rna",
    response_feature_role="gene",
    conditioning_modalities=["cuttag_h3k27ac"],
    name="rna_after_h3k27ac",
)
RNA_AFTER_ATAC_AND_MARK = ConditionalModelSpec(
    response_modality="rna",
    response_feature_role="gene",
    conditioning_modalities=["atac", "cuttag_h3k27ac"],
    name="rna_after_atac_h3k27ac",
)
PROTEIN_AFTER_RNA = ConditionalModelSpec(
    response_modality="protein",
    response_feature_role="protein",
    conditioning_modalities=["rna"],
    name="protein_after_rna",
)

⸻

输出 conditional_effects.tsv

字段：

event_id
model_name
response_modality
conditioning_modalities
condition_coef
condition_pval
condition_fdr
attenuation
model_used
converged

这样你可以说：

MoDES supports multi-layer conditional decomposition

否则就只能说 RNA-after-ATAC。

⸻

6. 第六硬伤：Spatial / Dynamic helper 没进主流程

要么降级描述，要么真正接入。

6.1 Spatial 接入方式

如果 data 是 SpatialMoDEData，build_evidence() 应该自动计算：

spatial_moran_i
spatial_moran_pval
neighbor_effect
edge_artifact_score

然后写进 event_modality_evidence.tsv：

modality = spatial
feature_id = event_id or region
coef = neighbor_effect
pval = spatial_pval
role = spatial_context

否则 spatial state 不应该启用。

⸻

6.2 Dynamic 接入方式

如果用户传：

time_col="pseudotime"

或：

contrasts=[...]

主 pipeline 应该调用：

dynamic effect estimation
pseudotime lag inference

并生成：

dynamic evidence

否则 dynamic 模块只能叫 helper。

⸻

7. 第七硬伤：extra modality feature matching 太弱

现在字符串匹配对于 CUT&Tag 很危险。

7.1 CUT&Tag region matching 必须用 interval overlap

不能靠：

split("|")[0] == peak_id

要支持：

ATAC peak chr:start-end
CUT&Tag peak chr:start-end

计算 overlap：

overlap_bp / min(width1, width2)
overlap_bp / union_width
reciprocal overlap

匹配条件：

min_reciprocal_overlap=0.5

输出：

region_match_score
region_match_method

⸻

7.2 Protein 不要 fuzzy match 默认开启

protein-gene link 必须显式提供：

protein_id
gene

如果要 fuzzy match，必须：

allow_fuzzy_protein_match=False

默认 false。

否则 reviewer 会说 toy match。

⸻

7.3 link uncertainty 进入 score

每个 event 要有：

link_score
link_source
region_match_score
protein_link_score

state score 应该有：

assignment_score *= link_score

这可以回应 link uncertainty 批评。

⸻

8. Benchmark 必须重做，不要再按规则生成规则

现在 synthetic benchmark 太同构。要补这些。

8.1 Null calibration

condition label shuffle
RNA sample shuffle
ATAC sample shuffle
CUT&Tag sample shuffle
random peak-gene links
random protein-gene links

指标：

non_null_rate
false_concordant_rate
false_modality_state_rate
state_support_qval distribution

⸻

8.2 Link-noise benchmark

true links: 100%, 75%, 50%, 25%
random links: 0%, 25%, 50%, 75%

看：

state accuracy
false concordant
false active_mark_concordant
false protein state

⸻

8.3 Batch/donor confounding

模拟：

balanced design
partially confounded batch
fully confounded batch
donor imbalance
low replicate

输出：

false positive rate
artifact_risk enrichment
rank deficiency detection

⸻

8.4 Multi-peak / multi-mark dependency

模拟：

one gene with 1 peak
one gene with 10 peaks
one gene with 50 peaks
multiple histone marks per region
multiple proteins per gene

看 event dependency 对 q-value 的影响。

⸻

8.5 Ablation

至少：

MoDES full
without conditional decomposition
without link_score weighting
without quality/artifact flag
without EB
naive DE+DA overlap
random links
proximity-only links
extra modality ignored

⸻

9. 真实数据验证必须换思路

PBMC spike-in 不够。你需要至少一个真实 biological contrast。

最低真实验证组合

Dataset A：RNA+ATAC time / perturbation

目标：

chromatin_primed events 是否在后续时间点出现 RNA response

验证：

early primed → late RNA up
random primed → no enrichment

Dataset B：RNA+CUT&Tag 或 RNA+ATAC+CUT&Tag

目标：

active_mark_concordant 是否富集 H3K27ac-supported enhancer-gene links
repressive_concordant 是否符合 H3K27me3/RNA opposite pattern

Dataset C：RNA+protein

目标：

protein_buffered / protein_memory 是否出现在 known surface markers

如果没有真实多模态数据，至少做：

external ChIP/CUT&Tag overlap validation

⸻

10. 生物学语言必须降级

这些要改：

drives transcriptional activation
complete regulatory chain activation
TF regulator
protein memory
causal regulatory event

改成：

ATAC and RNA changes are concordant
candidate activating mark-associated event
candidate TF motif annotation
protein-layer discordance
putative regulatory association

示例：

原：

Local chromatin opening drives transcriptional activation.

改：

Local chromatin accessibility and RNA abundance change concordantly under the tested contrast.

原：

Complete regulatory chain activation through to protein output.

改：

ATAC, RNA and protein layers show concordant differential signal.

⸻

11. 给你一个最小可接受 v2.0 修复路线

如果你要继续叫 multi-omics，最小要完成这 8 个。

1. 主输出拆表

event_table.tsv 固定主表
event_modality_evidence.tsv 长格式多模态证据
conditional_effects.tsv 多模型条件分解

2. StateRule grammar

所有 state 都由 grammar 定义，不再 priority if-else。

3. state_support_pval

每个 state 的 p-value 来自触发该 state 的 required modalities。

4. state_assignment_score

替代 state_confidence，不叫 posterior。

5. 多模态 conditional decomposition

至少实现：

RNA after ATAC
RNA after H3K27ac
RNA after ATAC+H3K27ac
Protein after RNA

6. CUT&Tag interval overlap

不能字符串匹配。

7. stress benchmarks

必须有：

null
link-noise
batch-confounded
weak-effect
ablation

8. 至少一个真实 non-RNA+ATAC demo

比如 RNA+CUT&Tag 或 RNA+ATAC+protein。

⸻

12. 如果你现在要回这位 reviewer，应该怎么说？

不要说：

We now support all modalities, therefore the concern is addressed.

要说：

We agree that the previous implementation was modality-aware but not modality-consistent. In the revised version, we restructure MoDES around a long-format event-modality evidence table, state-specific evidence rules, state-support p-values derived from the modalities that trigger each state, and modality-specific conditional decompositions. We also rename state_confidence to state_assignment_score and explicitly avoid posterior interpretation. We add null, link-noise, batch-confounded, weak-effect and ablation benchmarks, plus real multi-modal validation.

这才是正面回应。

⸻

13. 最终判断

这份三审拒稿不是坏事，它告诉你现在的 v2.0 真正缺什么。

你现在的 MoDES 是 multi-modal interface，但还不是 multi-modal inference/statistical engine。

要跨过去，必须做三件根本改变：

1. 多模态证据长表化；
2. state p-value 与触发 state 的 modalities 绑定；
3. state classifier 从 priority rules 改成 grammar-based scoring。

然后补：

conditional decomposition 多模态化
interval-based feature matching
stress benchmarks
真实非 RNA+ATAC 验证
生物学措辞降级

完成这些后，MoDES 才能比较有底气地叫：

multi-omics event-state framework

否则，它还是会被审稿人定位成：

RNA+ATAC core plus extra modality rule extensions
