Lings，我检查了当前 master，这轮审稿人的判断基本是对的：现在不能接收，也不能说已经小修；但已经不该 Reject，确实是 targeted Major Revision。

你已经把很多大方向修对了：state_rules.py 现在有 RequiredAbsentEvidence / ForbiddenEvidence / OptionalEvidence / MissingPolicy，RNA+ATAC 状态也扩展到了 concordant_activation / concordant_repression / discordant_opening_repression / discordant_closing_activation / chromatin_open_primed / chromatin_closed_primed / rna_up_only / rna_down_only / null；StateClassifier 已经输出 state_assignment_score / state_support_score / state_support_adjusted_score / supporting_modalities / absent_modalities / conflicting_modalities / missing_modalities；core.py 也有固定 event_table 加长格式 event_modality_evidence 的结构。(raw.githubusercontent.com￼) (raw.githubusercontent.com￼) (raw.githubusercontent.com￼)

但是审稿人列的几个硬 bug 仍然成立，尤其是：mark_only 和 active_enhancer_primed 的规则仍会打平；required_absent(require_available=True) 仍没有真正区分“模态缺失”和“测到了但不显著”；filter() 里 state_support_adjusted_score 的过滤方向是反的；README 还在推荐不存在的 min_confidence 和旧的 max_event_fdr 用法；tag 仍写着 “Full multi-modal MoDES platform / Complete modality support”，这和 README/ROADMAP 的 experimental 表述冲突。(github.com￼) (raw.githubusercontent.com￼) (github.com￼)

下面是完整修改 list。

⸻

P0：会直接影响结果或用户运行的硬修复

1. ✅ [DONE] 修 mark_only 被 active_enhancer_primed 抢走的问题

当前问题

现在 ACTIVE_ENHANCER_PRIMED 和 MARK_ONLY 都只要求：

RequiredEvidence("cuttag_activating", +1, role="activating_mark")

区别是：

active_enhancer_primed:
  required_absent = RNA
  optional = ATAC +
mark_only:
  required_absent = RNA + ATAC

但 _score_rule() 里的 specificity bonus 只看 n_required，不看 n_absent_ok，所以在 CUT&Tag ↑, RNA absent, ATAC absent 这种真正 mark_only 场景下，两者分数很可能一样；因为 EPIGENOMIC_RULES 里 ACTIVE_ENHANCER_PRIMED 排在 MARK_ONLY 前面，排序稳定时前者会赢。(raw.githubusercontent.com￼) (raw.githubusercontent.com￼)

修改

在 _score_rule() 中把 required_absent 的满足数量纳入 specificity：

specificity = (
    1.0
    + 0.30 * max(n_required - 1, 0)
    + 0.20 * n_absent_ok
    + 0.10 * len(rule.forbidden)
)

这样：

mark_only:
  required = 1
  required_absent satisfied = 2
active_enhancer_primed:
  required = 1
  required_absent satisfied = 1

在 ATAC absent 时，mark_only 会自然得分更高。

同时建议把 ACTIVE_ENHANCER_PRIMED 改成更明确的逻辑：

ACTIVE_ENHANCER_PRIMED = StateRule(
    name="active_enhancer_primed",
    state_family="epigenomic_primed",
    required=[
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna", require_available=True),
    ],
    optional=[
        OptionalEvidence("atac", +1, bonus=0.25),
    ],
    description="Activating mark increases while RNA is measured but not significant; ATAC support is optional.",
)

MARK_ONLY 保持：

MARK_ONLY = StateRule(
    name="mark_only",
    state_family="epigenomic_only",
    required=[
        RequiredEvidence("cuttag_activating", +1, role="activating_mark"),
    ],
    required_absent=[
        RequiredAbsentEvidence("rna", require_available=True),
        RequiredAbsentEvidence("atac", require_available=True),
    ],
    description="Activating mark changes while RNA and ATAC are measured but not significant.",
)

新测试

def test_mark_only_wins_when_mark_up_rna_absent_atac_absent():
    evidence = pd.DataFrame({
        "event_id": ["e1"],
        "atac_fdr": [1.0],
        "atac_direction": [0],
        "atac_pval": [1.0],
        "atac_coef": [0.0],
        "rna_fdr": [1.0],
        "rna_direction": [0],
        "rna_pval": [1.0],
        "rna_coef": [0.0],
        "h3k27ac_z": [8.0],
        "h3k27ac_fdr": [1e-8],
        "h3k27ac_pval": [1e-8],
        "h3k27ac_direction": [1],
        "h3k27ac_coef": [2.0],
        "quality_score": [1.0],
    })
    clf = StateClassifier(
        modality_specs={
            "h3k27ac": ModalitySpec(
                name="h3k27ac",
                assay="CUTTAG",
                feature_type="region",
                target="H3K27ac",
                regulatory_role="activating_mark",
                expected_rna_direction=1,
            )
        }
    )
    out = clf.classify(evidence)
    assert out.loc[0, "state"] == "mark_only"

⸻

2. ✅ [DONE] 修 required_absent(require_available=True) 的 missing 逻辑

当前问题

required_absent 的语义应该是：

该模态必须被测量 / 匹配到该 event，并且不显著。

但现在如果 extra modality 没匹配到 feature，EvidenceBuilder 会填：

fdr = 1.0
direction = 0
coef = nan

随后 _modality_evidence_map() 会把它当成存在的 evidence。这样“没有 protein evidence / 没匹配到 protein feature”会被误认为“protein measured but not significant”。这会错误触发 protein_buffered、protein_memory、mark_only 等 required_absent 状态。(raw.githubusercontent.com￼)

修改

所有 modality evidence 必须带：

available
matched
measured
missing_reason

在 EvidenceBuilder 中，extra modality 没匹配到 feature 时，不要只填 fdr=1.0。应同时填：

record[f"{mod_name}_available"] = True
record[f"{mod_name}_matched"] = False
record[f"{mod_name}_measured"] = False
record[f"{mod_name}_missing_reason"] = "no_feature_match"
record[f"{mod_name}_fdr"] = np.nan
record[f"{mod_name}_pval"] = np.nan
record[f"{mod_name}_direction"] = 0
record[f"{mod_name}_coef"] = np.nan

如果 modality 根本不存在：

available = False
matched = False
measured = False
missing_reason = "modality_absent"

如果 feature 匹配到了但不显著：

available = True
matched = True
measured = True
fdr >= threshold
direction = 0 or sign(coef)

然后 _modality_evidence_map() 返回：

ev[mod_name] = {
    "fdr": ...,
    "direction": ...,
    "pval": ...,
    "coef": ...,
    "available": bool(row.get(f"{mod_name}_available", False)),
    "matched": bool(row.get(f"{mod_name}_matched", False)),
    "measured": bool(row.get(f"{mod_name}_measured", False)),
    "missing_reason": row.get(f"{mod_name}_missing_reason", ""),
}

在 required_absent 检查中：

mod_ev = self._resolve_modality_evidence(ev, ra)
if mod_ev is None or not mod_ev.get("available", False) or not mod_ev.get("matched", False):
    if ra.require_available:
        n_missing += 1
        missing_mods.append(f"{ra.modality}(unavailable_or_unmatched)")
        n_absent_fail += 1
    else:
        n_absent_ok += 1
        absent.append(ra.modality)
elif not mod_ev.get("measured", True):
    if ra.require_available:
        n_missing += 1
        missing_mods.append(f"{ra.modality}(not_measured)")
        n_absent_fail += 1
    else:
        n_absent_ok += 1
else:
    sig = mod_ev["fdr"] < self.fdr_threshold
    if not sig:
        n_absent_ok += 1
        absent.append(ra.modality)
    else:
        n_absent_fail += 1
        conflicting.append(f"{ra.modality}(should_be_absent)")

新测试

def test_protein_buffered_requires_matched_protein_evidence():
    # RNA up, protein modality exists but no protein feature matched
    # should NOT trigger protein_buffered_up
def test_protein_buffered_triggers_when_protein_measured_and_not_significant():
    # RNA up, protein matched, measured, fdr=1.0
    # should trigger protein_buffered_up

⸻

3. ✅ [DONE] 修 state_support_adjusted_score 过滤方向反了的问题

当前问题

state_support_adjusted_score 是：

越大 = 证据越强

因为它是：

-log10(adjusted pseudo-p)

但 MoDESResult.filter() 现在有：

max_state_support_adjusted_score
df = df[df["state_support_adjusted_score"] <= max_state_support_adjusted_score]

这会保留弱证据事件，过滤掉强证据事件。旧参数 max_event_fdr 也被映射到 state_support_adjusted_score，同样用 <=。(raw.githubusercontent.com￼)

修改

新增正确参数：

min_state_support_adjusted_score: float | None = None

逻辑：

if min_state_support_adjusted_score is not None:
    df = df[df["state_support_adjusted_score"] >= min_state_support_adjusted_score]

保留旧参数但不要直接错误映射。对 max_event_fdr 给 warning，并转换成 score 阈值：

if max_event_fdr is not None:
    warnings.warn(
        "max_event_fdr is deprecated. "
        "Use min_state_support_adjusted_score or max_state_support_pseudo_q.",
        DeprecationWarning,
    )
    # 如果仍然允许旧用法，把 q 阈值转换成 -log10(q)
    score_thr = -np.log10(max(max_event_fdr, 1e-15))
    df = df[df["state_support_adjusted_score"] >= score_thr]

建议再加一个更明确的参数：

max_state_support_pseudo_q: float | None = None

但如果不想增加太多 API，先做 min_state_support_adjusted_score 即可。

README 修改

把：

trusted = result.filter(
    state="concordant",
    min_confidence=0.8,
    max_event_fdr=0.05,
    exclude_high_artifact=True,
)

改成：

trusted = result.filter(
    state_family="concordant",
    min_assignment_score=0.8,
    min_state_support_adjusted_score=2.0,
    exclude_high_artifact=True,
)

这里 2.0 相当于 10^-2 量级的 ranking score，但文档必须说它不是 formal FDR。

新测试

def test_filter_keeps_high_support_scores():
    df = pd.DataFrame({
        "state_support_adjusted_score": [0.5, 2.0, 5.0],
        ...
    })
    result = MoDESResult(event_table=df)
    out = result.filter(min_state_support_adjusted_score=2.0)
    assert set(out["state_support_adjusted_score"]) == {2.0, 5.0}

⸻

4. ✅ [DONE] 修 README/API 不一致

当前问题

README 仍推荐不存在的参数：

min_confidence=0.8

但 MoDESResult.filter() 没有这个参数。README 还推荐 state="concordant"，但当前具体状态名已经是：

concordant_activation
concordant_repression

concordant 更像 state_family。(github.com￼) (raw.githubusercontent.com￼)

修改

README 的过滤示例改为：

# View concordant family events
conc = result.filter(state_family="concordant")
# Or specific states
conc_act = result.filter(state="concordant_activation")
conc_rep = result.filter(state="concordant_repression")
# Exclude high artifact-risk events
clean = result.filter(exclude_high_artifact=True)
# Use assignment/support scores for ranking-oriented filtering
trusted = result.filter(
    state_family="concordant",
    min_assignment_score=0.8,
    min_state_support_adjusted_score=2.0,
    exclude_high_artifact=True,
)

如果保留旧字段，README 应写：

Deprecated compatibility fields:
- state_confidence_deprecated
- event_pval_deprecated
- event_fdr_deprecated
Do not use them for new analyses.

⸻

5. ✅ [DONE] 修 README 的状态表

当前问题

README 仍写 RNA+ATAC core 五状态：

concordant
chromatin_primed
rna_only
discordant_opposite
null

但代码现在输出具体状态：

concordant_activation
concordant_repression
discordant_opening_repression
discordant_closing_activation
chromatin_open_primed
chromatin_closed_primed
rna_up_only
rna_down_only
null

README 必须同步。(github.com￼) (raw.githubusercontent.com￼)

修改

README 状态表改成：

state_family	state	Pattern	Interpretation
concordant	concordant_activation	ATAC↑ RNA↑	Accessibility and RNA increase concordantly
concordant	concordant_repression	ATAC↓ RNA↓	Accessibility and RNA decrease concordantly
discordant	discordant_opening_repression	ATAC↑ RNA↓	Accessibility increase with RNA decrease
discordant	discordant_closing_activation	ATAC↓ RNA↑	Accessibility decrease with RNA increase
chromatin_primed	chromatin_open_primed	ATAC↑ RNA not significant	Open chromatin candidate without RNA response
chromatin_primed	chromatin_closed_primed	ATAC↓ RNA not significant	Closing chromatin candidate without RNA response
rna_only	rna_up_only	RNA↑ ATAC not significant	RNA up without linked ATAC change
rna_only	rna_down_only	RNA↓ ATAC not significant	RNA down without linked ATAC change
null	null	no significant required evidence	No assigned event-state signal

注意措辞要用：

concordant differential signal
candidate event

不要用：

drives transcription

⸻

6. ✅ [DONE] 修 GitHub tag / release narrative 过度主张

当前问题

GitHub v2.0.0 tag 仍写：

Full multi-modal MoDES platform
Complete modality support
MoDES-Spatial ✅
MoDES-Dynamic ✅

但 ROADMAP 仍把 multi-condition、pseudotime、real multi-modal validation、full multi-layer event model 放到 hardening 或 v2.1+。(github.com￼) (raw.githubusercontent.com￼)

修改

如果能改 release notes，就改成：

v2.0.0: Grammar-driven multi-modal evidence-scoring prototype
- RNA+ATAC core supported
- CUT&Tag/protein/spatial evidence modules are experimental
- StateRule grammar added
- Long-format event_modality_evidence.tsv added
- State support scores are ranking-oriented, not formal FDR
- Spatial and dynamic modules are helper APIs under validation

不要写：

Full multi-modal platform
Complete modality support

如果 tag message 已经不能改，建议新建：

v2.0.1

并在 release notes 中明确修正 v2.0.0 的 overclaim。

⸻

7. ✅ [DONE] 修 multi-condition / dynamic 叙事

当前问题

effects.py 仍然在非二分类 condition 时抛出：

MoDES v0.1 supports only binary condition.

这说明 multi-condition 仍不是主 pipeline 支持。审稿人说得对：不能说 MoDES-Dynamic complete。(raw.githubusercontent.com￼)

修改

把错误信息改成当前版本并诚实描述：

raise NotImplementedError(
    "MoDES v2.0 main effect-estimation pipeline supports binary contrasts only. "
    "For multi-condition designs, run explicit pairwise contrasts or use the "
    "experimental dynamic helper APIs. Full multi-condition modeling is planned."
)

README/ROADMAP 中写：

Multi-condition / pseudotime: experimental helper APIs; not yet validated as end-to-end main pipeline.

⸻

8. ✅ [DONE] 明确 artifact risk 对 experimental modalities 是 heuristic

当前问题

_compute_artifact_risk() 已经遍历 extra modalities，检查 extra modality z-score 和 region match score，这是进步；但还没有 protein antibody background、CUT&Tag FRiP/blacklist、spatial density/edge 等完整 QC。(raw.githubusercontent.com￼)

修改

README / statistical_model 加：

artifact_risk is a heuristic quality flag. For experimental modalities,
it currently uses available depth, detection, z-score and matching diagnostics;
it is not a validated assay-specific QC classifier.

代码里可以增加更具体字段，但最少先让文档不夸大。

如果继续完善代码，建议：

protein_missingness_score
protein_background_score
cuttag_region_match_score
cuttag_blacklist_flag
spatial_edge_score
spatial_density_score

先输出到 event_modality_evidence.tsv，不用都进入主表。

⸻

P1：必须补的 targeted tests

审稿人说 integration tests 仍像 smoke test，这点成立。你需要加专门测试，不要只测 pipeline 能跑。

1. ✅ [DONE] RNA+ATAC 方向测试
2. ✅ [DONE] CUT&Tag 状态测试
3. ✅ [DONE] Protein absent/missing 测试
4. ✅ [DONE] Filter 测试
5. [NOTED] Spatial role 测试 — requires spatial evidence infrastructure not yet wired into main EvidenceBuilder

⸻

P2：文档同步清单

README

必须改：

✅ [DONE] Filtering example
✅ [DONE] RNA+ATAC state table
✅ [DONE] Deprecated fields explanation
✅ [DONE] Multi-condition/pseudotime status
✅ [DONE] Experimental modality caveat
✅ [DONE] Causal language

⸻

docs/statistical_model.md

已经比以前好，但要再补：

1. required_absent / missing 的完整语义
2. state_support_adjusted_score 的过滤方向
3. experimental modalities artifact risk 是 heuristic
4. conditional_effects 是 diagnostic
5. binary contrast limitation

⸻

docs/output_schema.md

必须确认：

event_table.tsv 包含 state_family 和具体 state
deprecated fields 名字和代码一致：
  state_confidence_deprecated
  event_pval_deprecated
  event_fdr_deprecated

⸻

docs/install_review.md

当前写：

119 tests pass
1 pre-existing test failure

这两个不能同时作为“通过”叙事。(raw.githubusercontent.com￼)

改成二选一：

All tests passed: 119/119

或者：

118/119 passed; known failing test is ...

不要模糊。

⸻

P3：release narrative 修复

Tag / Release

因为 GitHub tags 页面仍显示 v2.0.0 过度宣称，你需要新 tag：

v2.0.1

release notes：

v2.0.1: Targeted Major Revision Fixes
- Fix state-rule specificity for mark_only vs active_enhancer_primed
- Distinguish missing/unmatched modalities from measured non-significant modalities
- Fix state_support_adjusted_score filtering direction
- Update README and output schema to v2.0 grammar states
- Downgrade support scores to ranking-oriented evidence scores
- Clarify experimental status of spatial/dynamic/protein/CUT&Tag modules

这样可以压过旧 tag 的 overclaim。

⸻

最终判断

审稿人这轮意见基本正确。你现在离 Minor Revision 很近，但有 3 个必须立即修的硬点：

1. mark_only vs active_enhancer_primed tie / priority bug
2. required_absent missing/unlinked modality 被当成 not significant
3. state_support_adjusted_score 过滤方向反了

再加 2 个用户层面的硬点：

4. README 推荐不存在的 min_confidence 参数
5. README 用旧 state 名，和 v2.0 具体 state 不一致

如果你修掉这五个，再补 targeted tests，下一轮就可以从 Major Revision 降到 Minor Revision。当前版本还不能接收，但已经不是推倒重来的大修了。
