Lings，我彻底审计了当前 master。我先说结论：

现在 MoDES 已经不是早期原型了，它已经有 v2.0 多模态软件骨架；但如果目标是“完美方案、可投稿、能扛审稿”，当前代码库还必须进行一轮系统性重构。
最核心的问题不是“再加几个功能”，而是：

状态语法、证据表、统计字段命名、模态特异建模、artifact 体系、benchmark 和文档必须统一成一个严格的 multi-omics evidence-scoring framework。

你以后让我审计 MoDES，我会默认按这个标准来：不再给最短修补路线，而是给完整、可投稿、可维护的最优整改方案。

⸻

一、当前代码库真实状态

当前仓库已经有完整的多模态雏形。README 标题已经写成 “MoDES: Multi-Omics Evidence Framework for Regulatory Event-State Annotation”，并明确说每个 event 的 multi-modal evidence 会被 state-rule grammar 评分，输出 state annotation 和 state-support p-value；README 也承认 CUT&Tag、Protein、Spatial、MuData、multi-condition/pseudotime 等仍是 experimental。 ￼

核心代码也确实已经进入 v2.0 架构：core.py 注释写明当前设计是 fixed event_table schema + long-format event_modality_evidence，并且 MoDES.run() 已经串起 build_events → estimate_effects → decompose → build_evidence → classify_states → _build_modality_evidence → _assemble_results。 ￼

但是当前代码仍有几个硬伤：state_rules.py 中 RNA+ATAC 规则仍然只覆盖 ATAC↑/RNA↑、ATAC↑/RNA↓、ATAC↑ only、RNA↑ only 等 activation-centric 模式；protein_buffered 和 protein_memory 仍然用 NeutralEvidence 表达“不变”，而 neutral 的定义只是“allowed but not required”；directed_pvalue() 明确写着它 不是 one-sided test，只是 directional evidence score，但 StateClassifier.classify() 仍然对 state_support_pval 做 BH 生成 state_support_qval。 ￼

另外，core.py 里确实有一个明确代码 bug：rec 字典里先写了使用 link_score 调整后的 "state_assignment_score": adjusted_assignment，后面又重复写了 "state_assignment_score": assignment_score，后者会覆盖前者。 ￼

文档也还没同步：docs/statistical_model.md 仍写 MoDES-RA v0.1.0，仍描述旧的 RNA+ATAC 五步流程、旧的 event-level p-value 和旧 EB 逻辑；docs/output_schema.md 虽写 v2.0.0，但仍主要描述旧 event_pval/event_fdr/state_confidence schema；CITATION.cff 仍是 version: 0.1.0-alpha；docs/install_review.md 仍是 v1.0.0-rc.1。 ￼

⸻

二、最终目标定义

不要再把目标设成“让 reviewer 不拒稿”。
你现在应该把 MoDES v2.0 的最终形态定义为：

MoDES is a multi-omics evidence-scoring and state-annotation framework for candidate regulatory events. It is not a calibrated causal inference framework.

也就是说，正式主张应该是：

candidate regulatory events
multi-modal evidence scoring
state-rule grammar
assignment score
artifact / quality diagnostics
event prioritization

而不是：

causal inference
posterior probability
formal state-level FDR
validated full multi-omics regulatory inference

这个定位可以保留 “multi-omics”，但不会被审稿人抓住“过度统计主张”。

⸻

三、P0 级：必须立即修的代码逻辑问题

1. ✅ [DONE] 修复 state_assignment_score 被覆盖的问题

位置： modes/core.py

当前逻辑：

adjusted_assignment = assignment_score * link_score
rec = {
    ...
    "state_assignment_score": adjusted_assignment,
    ...
    "state_assignment_score": assignment_score,
    ...
}

后一个 key 会覆盖前一个 key，所以 link_score 调整无效。 ￼

修改方案：

raw_assignment_score = assignment_score
adjusted_assignment_score = (
    raw_assignment_score * link_score
    if not np.isnan(raw_assignment_score)
    else np.nan
)
rec = {
    ...
    "state_assignment_score": adjusted_assignment_score,
    "raw_state_assignment_score": raw_assignment_score,
    "link_score": link_score,
    ...
}

如果你不想在主表增加 raw_state_assignment_score，至少也要写到 diagnostic 表。

测试：

def test_assignment_score_is_link_adjusted():
    # 构造 state_assignment_score=1.0, link_score=0.5
    # 结果应为 0.5，而不是 1.0

这是最优先修复项。

⸻

2. ✅ [DONE] 重写 StateRule 数据结构，加入 absent / forbidden / missing 语义

当前 StateRule 只有：

required
neutral
forbidden

而 NeutralEvidence 的定义只是“allowed but not required”，不是“不显著”。 ￼

这导致 mark_only、active_enhancer_primed、protein_buffered、protein_memory 的生物语义不成立。

新增结构：

@dataclass(frozen=True)
class RequiredEvidence:
    modality: str
    direction: int
    role: str | None = None
    target: str | None = None
@dataclass(frozen=True)
class RequiredAbsentEvidence:
    modality: str
    role: str | None = None
    target: str | None = None
    require_available: bool = True
@dataclass(frozen=True)
class ForbiddenEvidence:
    modality: str
    direction: int | None = None
    role: str | None = None
    target: str | None = None
@dataclass(frozen=True)
class OptionalEvidence:
    modality: str
    direction: int | None = None
    role: str | None = None
    target: str | None = None
    bonus: float = 0.1
@dataclass(frozen=True)
class MissingPolicy:
    modality: str
    allowed: bool = False
    penalty: float = 0.5
@dataclass(frozen=True)
class StateRule:
    name: str
    required: Sequence[RequiredEvidence] = ()
    required_absent: Sequence[RequiredAbsentEvidence] = ()
    forbidden: Sequence[ForbiddenEvidence] = ()
    optional: Sequence[OptionalEvidence] = ()
    missing_policy: Sequence[MissingPolicy] = ()
    state_family: str = ""
    description: str = ""
    interpretation_strength: str = "association"

核心原则：

required = 必须显著且方向匹配
required_absent = 必须测到，但不显著
forbidden = 不允许显著，或不允许某方向显著
optional = 有则加分，没有不扣分
missing = 模态缺失怎么处理

这一步会直接解决 reviewer 对 neutral 语义的批评。

⸻

3. ✅ [DONE] 修 RNA+ATAC core states，覆盖上下调全方向

当前 RA 规则仍然是：

concordant: ATAC +1, RNA +1
chromatin_primed: ATAC +1
rna_only: RNA +1
discordant_opposite: ATAC +1, RNA -1

代码里 RA_RULES 也只列了这些规则。 ￼

这会漏掉：

ATAC↓ RNA↓
ATAC↓ RNA↑
ATAC↓ RNA unchanged
RNA↓ only

应改成完整状态：

CONCORDANT_ACTIVATION = StateRule(
    name="concordant_activation",
    state_family="concordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", +1),
    ],
)
CONCORDANT_REPRESSION = StateRule(
    name="concordant_repression",
    state_family="concordant",
    required=[
        RequiredEvidence("atac", -1),
        RequiredEvidence("rna", -1),
    ],
)
DISCORDANT_OPENING_REPRESSION = StateRule(
    name="discordant_opening_repression",
    state_family="discordant",
    required=[
        RequiredEvidence("atac", +1),
        RequiredEvidence("rna", -1),
    ],
)
DISCORDANT_CLOSING_ACTIVATION = StateRule(
    name="discordant_closing_activation",
    state_family="discordant",
    required=[
        RequiredEvidence("atac", -1),
        RequiredEvidence("rna", +1),
    ],
)
CHROMATIN_OPEN_PRIMED = StateRule(
    name="chromatin_open_primed",
    state_family="chromatin_primed",
    required=[RequiredEvidence("atac", +1)],
    required_absent=[RequiredAbsentEvidence("rna")],
)
CHROMATIN_CLOSED_PRIMED = StateRule(
    name="chromatin_closed_primed",
    state_family="chromatin_primed",
    required=[RequiredEvidence("atac", -1)],
    required_absent=[RequiredAbsentEvidence("rna")],
)
RNA_UP_ONLY = StateRule(
    name="rna_up_only",
    state_family="rna_only",
    required=[RequiredEvidence("rna", +1)],
    required_absent=[RequiredAbsentEvidence("atac")],
)
RNA_DOWN_ONLY = StateRule(
    name="rna_down_only",
    state_family="rna_only",
    required=[RequiredEvidence("rna", -1)],
    required_absent=[RequiredAbsentEvidence("atac")],
)

输出建议：

主表同时输出：

state_family
state

例如：

state_family = concordant
state = concordant_repression

这样旧用户还能按 family 分析，新用户能看到方向。

⸻

4. ✅ [DONE] 修 CUT&Tag 状态：active_enhancer_primed 与 mark_only 必须可区分

当前两个状态的 required evidence 都是：

RequiredEvidence("cuttag_activating", +1, role="activating_mark")

区别只是 neutral 里有没有 ATAC/RNA；但 neutral 不是“必须不显著”。 ￼

修法：

ACTIVE_ENHANCER_PRIMED = StateRule(
    name="active_enhancer_primed",
    state_family="epigenomic_primed",
    required=[
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna"),
    ],
    optional=[
        OptionalEvidence("atac", +1),
    ],
)
MARK_ONLY = StateRule(
    name="mark_only",
    state_family="epigenomic_only",
    required=[
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna"),
        RequiredAbsentEvidence("atac"),
    ],
)

解释：

active_enhancer_primed = mark ↑, RNA 不显著，ATAC 可有可无
mark_only = mark ↑, RNA 不显著，ATAC 也不显著

⸻

5. ✅ [DONE] 修 protein states：protein_buffered 和 protein_memory 不能用 neutral

当前：

PROTEIN_BUFFERED:
  required = RNA +1
  neutral = protein
PROTEIN_MEMORY:
  required = protein +1
  neutral = RNA

这不能表达 “protein unchanged” 或 “RNA baseline”。 ￼

修法：

PROTEIN_BUFFERED_UP = StateRule(
    name="protein_buffered_up",
    state_family="protein_buffered",
    required=[
        RequiredEvidence("rna", +1),
    ],
    required_absent=[
        RequiredAbsentEvidence("protein", require_available=True),
    ],
)
PROTEIN_BUFFERED_DOWN = StateRule(
    name="protein_buffered_down",
    state_family="protein_buffered",
    required=[
        RequiredEvidence("rna", -1),
    ],
    required_absent=[
        RequiredAbsentEvidence("protein", require_available=True),
    ],
)
PROTEIN_MEMORY_UP = StateRule(
    name="protein_memory_up",
    state_family="protein_memory",
    required=[
        RequiredEvidence("protein", +1),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna", require_available=True),
    ],
)
PROTEIN_MEMORY_DOWN = StateRule(
    name="protein_memory_down",
    state_family="protein_memory",
    required=[
        RequiredEvidence("protein", -1),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna", require_available=True),
    ],
)

并且如果 protein modality 缺失，不能触发 protein_buffered 或 protein_memory。

⸻

6. 修 spatial evidence role 解析

当前 spatial rules 已经写了 role，比如：

spatial_autocorrelation
neighbor_effect
edge_artifact

但 _resolve_modality_evidence() 对 mod_name == "spatial" 的逻辑只是找第一个 assay 为 "SPATIAL" 的 evidence，并没有按 role 区分。 ￼

修法：

_resolve_modality_evidence() 必须同时匹配：

modality
assay
role
target

伪代码：

def _resolve_modality_evidence(self, ev: dict, req) -> dict | None:
    for key, val in ev.items():
        spec = self.modality_specs.get(key)
        # modality alias
        if req.modality in ("atac", "rna") and key != req.modality:
            continue
        if req.modality == "protein":
            if not (spec and spec.assay == "PROTEIN"):
                continue
        if req.modality == "spatial":
            if not (spec and spec.assay == "SPATIAL"):
                continue
        if req.role is not None:
            if val.get("role") != req.role and getattr(spec, "regulatory_role", None) != req.role:
                continue
        if req.target is not None:
            if val.get("target") != req.target and getattr(spec, "target", None) != req.target:
                continue
        return val
    return None

同时 spatial evidence 应拆成：

spatial_moran
spatial_neighbor
spatial_edge

而不是一个 generic spatial evidence。

⸻

四、P0 级：统计字段必须降级或校准

7. ✅ [DONE] 不要继续把 directional score 叫 p-value

directed_pvalue() 文档明确说它不是 one-sided test，只是 directional evidence score；方向匹配时返回 pval/2，方向相反返回 1.0。 ￼

但 StateClassifier.classify() 仍然对 state_support_pval 做 BH 生成 state_support_qval。 ￼

这会被持续攻击。

最完美方案：

把字段改成：

state_support_score
state_support_adjusted_score

定义：

state_support_score = ranking-oriented directional evidence score
state_support_adjusted_score = BH-adjusted ranking score, not formal FDR

旧字段保留为 deprecated alias：

state_support_pval_deprecated
state_support_qval_deprecated

主表不要再把 state_support_qval 放在最核心位置。

⸻

8. 增加 calibration benchmark，否则不许说 q-value

如果你坚持保留 qval，就必须做 calibration：

condition-label permutation
RNA sample shuffle
ATAC sample shuffle
CUT&Tag sample shuffle
protein sample shuffle
random peak-gene links
random protein-gene links
multi-peak-per-gene dependency

输出：

nominal_threshold
empirical_false_state_rate
n_selected
state_family
calibration_gap

如果没有这个，文档必须写：

state_support_adjusted_score is for ranking only.

⸻

五、P0 级：输出 schema 必须改成固定主表 + 长格式证据表

9. ✅ [DONE] 固定 event_table.tsv

当前 core.py 已经有 fixed event_table 的想法，但仍保留旧字段 state_confidence/event_pval/event_fdr。 ￼

建议主表固定为：

event_id
peak_id
gene
context
tf_name
link_source
link_score
state_family
state
state_assignment_score
state_support_score
state_support_adjusted_score
supporting_modalities
absent_modalities
neutral_modalities
conflicting_modalities
missing_modalities
artifact_risk
artifact_reason
quality_score

旧字段放后面并标 deprecated：

state_confidence_deprecated
event_pval_deprecated
event_fdr_deprecated

⸻

10. event_modality_evidence.tsv 成为真正多模态证据表

当前 _build_modality_evidence() 已经输出 long-format rows，但 quality 分数固定很粗糙，且未包含 directed score / modality-specific QC。 ￼

完善字段：

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
directed_score
quality_score
detection_score
depth_score
batch_score
link_score
region_match_score
model_used
converged
warning

这张表应该成为 v2.0 的核心证据表。

⸻

11. conditional_effects.tsv 必须明确是 diagnostic，除非进入 state score

当前 decompose() 已经尝试加入 RNA_AFTER_H3K27AC、RNA_AFTER_ATAC_H3K27AC、PROTEIN_AFTER_RNA 等模型，但这些 conditional results 没有真正参与 state scoring。 ￼

你有两个选择：

选择 A：[DONE] 明确降级为 diagnostic

文档写：

Conditional models are diagnostic summaries and are not used as primary state assignment evidence in v2.0.

选择 B：纳入 state scoring

例如：

active_mark_concordant:
  required = H3K27ac ↑, RNA ↑
  optional diagnostic = RNA condition effect attenuated after H3K27ac adjustment

如果加到 scoring，必须明确：

conditional_support_score
attenuation_score

我建议大修阶段先走 A。最稳。

⸻

六、P0 级：模态特异建模问题

12. ✅ [DONE] Protein 不能用 ATAC library size 作为 offset

estimate_modality_effects() 里当前逻辑写着 protein 使用 ATAC library size，注释是避免 RNA expression bias。代码如下：spec.assay == "PROTEIN" 时 own_ls = atac_ls。 ￼

这在统计上很难自洽。

修改：

if spec and spec.assay == "PROTEIN":
    mat = data.modalities[modality_name]
    protein_lib = np.log(mat.sum(axis=1).replace(0, np.nan).values)
    protein_lib = np.nan_to_num(protein_lib, nan=np.nanmedian(protein_lib))
    own_ls = protein_lib

或者支持用户输入 normalized protein matrix：

protein_normalization="library_size" | "clr" | "none"

如果是 ADT-like protein，推荐：

CLR or protein-specific library-size offset

不要用 ATAC offset。

⸻

13. CUT&Tag / CUT&RUN / ChIP 需要 interval overlap matching

当前 extra modality evidence 中 region-like modality 基本用 peak 作为 feature key；如果 CUT&Tag peak 和 ATAC peak 名不一致，就很难匹配。_build_modality_evidence() 对 region-like extra modality 使用 feature = peak 并直接 eff_dict.get(feature)。 ￼

真实 CUT&Tag 和 ATAC peaks 常常不是同一套 peaks。

新增模块：

modes/regions.py

实现：

match_regions_by_overlap(
    query_regions,
    target_regions,
    min_reciprocal_overlap=0.5,
    min_overlap_bp=50,
)

输出：

query_region
target_region
overlap_bp
reciprocal_overlap
region_match_score

然后 event evidence 使用：

region_match_score

而不是字符串相等。

⸻

七、P0 级：artifact risk 必须多模态化

当前 _compute_artifact_risk() 主要看：

quality_score
atac_fdr
rna_fdr
z_atac
z_rna
single_modality_low_quality

代码中确实只看 RNA/ATAC。 ￼

这对 CUT&Tag、protein、spatial states 不够。

扩展 artifact reasons：

RNA/ATAC

low_rna_depth
low_atac_depth
low_rna_detection
low_atac_detection
library_size_outlier
batch_associated

CUT&Tag / CUT&RUN / ChIP

low_cuttag_depth
low_cuttag_detection
blacklist_overlap
weak_region_match
broad_mark_low_coverage
unexpected_mark_direction

Protein

low_protein_library_size
high_protein_missingness
protein_saturation
ambiguous_protein_gene_link
high_background

Spatial

edge_artifact
low_local_density
disconnected_graph
region_imbalance
neighbor_null_fail

输出：

artifact_risk
artifact_reason
modality_artifact_reasons

例如：

artifact_reason = low_cuttag_depth;weak_region_match;protein_missingness

⸻

八、P1 级：文档必须彻底同步

14. ✅ [DONE] 重写 docs/statistical_model.md

当前文件仍是 v0.1 RNA+ATAC 旧模型。 ￼

新的结构应该是：

1. Scope and claims
2. Candidate event definition
3. Long-format event-modality evidence
4. Per-modality marginal GLMs
5. StateRule grammar
6. Required / absent / forbidden / optional / missing evidence
7. State assignment score
8. State support score, not formal FDR
9. Conditional decompositions as diagnostics
10. Artifact risk
11. Experimental modality limitations
12. What MoDES does not claim

必须明确：

MoDES does not infer causality.
MoDES does not provide calibrated posterior probabilities.
MoDES does not provide formal post-selection FDR unless calibration is performed.

⸻

15. ✅ [DONE] 重写 docs/output_schema.md

当前 schema 仍旧，仍写旧 5 个 RNA+ATAC state 和旧 event_pval/event_fdr。 ￼

必须改成：

event_table.tsv
event_modality_evidence.tsv
conditional_effects.tsv
model_diagnostics.tsv
run_params.tsv

并标注 deprecated 字段：

state_confidence: deprecated alias
event_pval: deprecated alias
event_fdr: deprecated alias

⸻

16. ✅ [DONE] 更新 CITATION.cff

当前仍是 version: 0.1.0-alpha，abstract 仍称 “statistical framework”。 ￼

改成：

version: 2.0.0
abstract: >-
  MoDES is a multi-omics evidence-scoring and event-state annotation
  framework for candidate regulatory events. It provides rule-based
  state assignments, modality-specific evidence summaries, and
  artifact/quality diagnostics.

不要写：

statistical framework that infers regulatory states

⸻

17. ✅ [DONE] 更新 docs/install_review.md

当前仍是 v1.0.0-rc.1，并写 91 tests、commit 509e66a 等旧信息。 ￼

必须重新跑 fresh install：

git clone ...
pip install -e .
pip install -r requirements-dev.txt
python -m pytest -q
modes --help
python examples/minimal_bulk/run_minimal.py
python examples/cuttag_bulk/run_cuttag.py
python examples/protein_bulk/run_protein.py

然后更新：

version = 2.0.0
test count = 当前真实数量
CI status = 当前真实 workflow
examples = 当前真实 examples

⸻

18. README 降低因果措辞

README 仍有：

complete regulatory chain activation
chromatin → transcription → protein

以及 state table 中 strong biological interpretation。 ￼

改成：

concordant differential signal
candidate regulatory association
multi-layer evidence pattern
protein-layer discordance
candidate TF motif annotation

不要写：

drives transcription
complete regulatory chain activation
protein memory persists
TF driver

⸻

九、P1 级：真实验证和 benchmark

19. 增加 null calibration benchmark

必须包括：

condition label shuffle
RNA sample shuffle
ATAC sample shuffle
CUT&Tag sample shuffle
protein sample shuffle
random peak-gene links
random protein-gene links

输出：

non_null_rate
false_state_rate
state_support_score_distribution
artifact_risk_distribution

⸻

20. 增加 weak-effect benchmark

设置：

effect_size = 0.2, 0.5, 1.0, 2.0

输出：

power
precision
recall
state_assignment_score calibration

⸻

21. 增加 link-noise benchmark

设置：

true links = 100%, 75%, 50%, 25%
random links = 0%, 25%, 50%, 75%

输出：

false_concordant_rate
false_epigenomic_state_rate
effect of link_score

⸻

22. 增加 batch/donor confounding benchmark

模拟：

balanced design
donor imbalance
batch partially confounded
batch fully confounded
library size outlier

输出：

false positive rate
artifact_risk enrichment
rank-deficiency detection

⸻

23. 增加 ablation benchmark

对比：

MoDES full
no extra modalities
no conditional diagnostics
no link_score
no artifact risk
no EB
naive DE+DA overlap
random links
proximity-only links

输出：

macro-F1
per-state precision/recall
artifact detection
runtime

⸻

24. 增加真实 biological contrast

PBMC spike-in 和 random pseudo-condition 不够。必须至少一个真实 contrast：

stimulated vs control
treated vs untreated
early vs late differentiation
disease vs control

验证：

known pathway enrichment
known marker enrichment
external enhancer-gene link enrichment
random-link negative control
published TF program overlap

如果有 CUT&Tag/protein，就做一个非 RNA+ATAC validation：

H3K27ac overlap
ChIP/CUT&Tag external support
known ADT marker behavior

⸻

十、P2 级：工程 API 与兼容性

25. ✅ [DONE] 移除或降级旧字段推荐

当前 README 仍推荐：

min_confidence=0.8
max_event_fdr=0.05

应改成：

min_assignment_score=...
max_state_support_adjusted_score=...
exclude_high_artifact=True

并说明：

state_confidence/event_fdr are deprecated compatibility aliases.

⸻

26. ✅ [DONE] filter() API 改名

当前 filter() 里仍有 max_event_fdr，但实际用的是 state_support_qval。 ￼

建议新增：

max_state_support_adjusted_score: float | None = None

保留旧参数但 warning：

if max_event_fdr is not None:
    warnings.warn("max_event_fdr is deprecated; use max_state_support_adjusted_score")

⸻

27. state_confidence 全部标 deprecated

主输出可以保留，但文档写：

state_confidence is deprecated.
Use state_assignment_score.
It is not a posterior probability.

⸻

28. Release tag narrative 修正

不要再写：

Full multi-modal platform
Complete modality support

改成：

MoDES v2.0.0: grammar-driven multi-omics evidence extension
- RNA+ATAC core
- experimental CUT&Tag/protein/spatial/dynamic evidence modules
- long-format modality evidence
- state-rule grammar
- no formal causal inference
- support scores are ranking-oriented, not calibrated FDR

⸻

十一、测试清单

你需要新增这些 tests。

State grammar tests

test_concordant_activation
test_concordant_repression
test_discordant_opening_repression
test_discordant_closing_activation
test_chromatin_open_primed
test_chromatin_closed_primed
test_rna_up_only
test_rna_down_only

Absence / missing tests

test_required_absent_rna
test_required_absent_protein
test_missing_protein_does_not_trigger_protein_buffered
test_missing_rna_does_not_trigger_protein_memory

CUT&Tag tests

test_active_enhancer_primed_requires_rna_absent
test_mark_only_requires_rna_and_atac_absent
test_repressive_concordant
test_derepression

Spatial tests

test_spatial_moran_role_resolution
test_spatial_neighbor_role_resolution
test_spatial_edge_artifact_forbidden_in_niche

Score tests

test_state_support_score_not_named_pvalue_in_primary_output
test_deprecated_event_pval_alias
test_assignment_score_uses_link_score

Schema tests

test_event_table_fixed_schema
test_event_modality_evidence_schema
test_conditional_effects_schema

⸻

十二、最终执行顺序：完整方案版

你说不要最短路径，所以这里给完整路线，不省略。

Phase 1：核心语法重构

1. 修 duplicate state_assignment_score key
2. 新增 RequiredAbsentEvidence / ForbiddenAnySignificantEvidence / OptionalEvidence / MissingPolicy
3. 重写 RNA+ATAC 全方向 state rules
4. 重写 CUT&Tag state rules
5. 重写 protein state rules
6. 重写 spatial state rules
7. 更新 StateClassifier score 逻辑
8. 添加 mixed_evidence / ambiguous 逻辑

Phase 2：统计字段重构

9. state_support_pval/qval 改为 state_support_score/adjusted_score
10. 旧 event_pval/event_fdr/state_confidence 标 deprecated
11. filter API 更新
12. README 和 output_schema 更新

Phase 3：证据表重构

13. event_table 固定主表
14. event_modality_evidence 长格式完整化
15. conditional_effects 作为 diagnostic 表
16. modality-specific QC components 加入 evidence
17. link_score / region_match_score 加入 evidence

Phase 4：模态建模修正

18. protein 使用 protein-specific size factor
19. CUT&Tag 使用 interval overlap matching
20. spatial evidence 按 moran/neighbor/edge 拆分
21. batch_col 传入 quality/artifact 计算
22. modality-specific artifact reasons

Phase 5：验证体系

23. null calibration benchmark
24. weak-effect benchmark
25. link-noise benchmark
26. confounding benchmark
27. ablation benchmark
28. true biological contrast dataset
29. CUT&Tag/protein external validation
30. benchmark summary tables

Phase 6：文档和发布

31. 重写 statistical_model.md
32. 重写 output_schema.md
33. 更新 CITATION.cff
34. 更新 install_review.md
35. 更新 CHANGELOG / ROADMAP
36. 修 README 过强措辞
37. 修 GitHub tag narrative
38. fresh install review
39. full CI + benchmark
40. 重新打 release

⸻

最终判断

现在 MoDES 已经有进入 Major Revision 的基础，但还没达到可接收状态。
要过下一轮，必须证明三件事：

1. state grammar 是逻辑自洽的，不再 activation-centric；
2. support score 是诚实的 ranking/evidence score，不伪装成严格 FDR；
3. multi-modal evidence 是真实进入 long-format evidence、QC、state assignment 和 validation，而不是只加规则标签。

你下一轮如果按上面做完，审稿人最可能从：

Major Revision

降到：

Minor Revision

甚至接收。

我之后对这个项目会默认采用这个“完整最优方案”标准，不再给你压缩版路线。
