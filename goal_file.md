Lings，这次从 Reject 变成 Major Revision 是一个关键信号：审稿人已经承认 MoDES 不是“空壳多模态”了。现在主要矛盾不是“你有没有写代码”，而是：

你已经把多模态接进了软件架构，但还没有把多模态接进统计解释、证据评分、验证体系和文档叙事。

我刚看了当前仓库 README，它确实已经把 MoDES 写成 v2.0 的 multi-modal regulatory event-state inference，并把 CUT&Tag、Protein、Spatial、MuData 等标为 experimental；同时 README 仍有一些较强的生物学解释，比如 “complete regulatory chain activation” 和 “chromatin opening drives transcriptional activation” 这类语言，容易被审稿人抓住。 ￼

下面是我建议你按 “能过大修” 的标准来做的完整方案。

⸻

总体策略

这轮不要再继续盲目加功能。你已经有足够多功能了。现在要做的是：

1. 修 state grammar 的逻辑漏洞；
2. 修 state_support_pval/qval 的统计表述；
3. 把文档彻底更新到 v2.0；
4. 把真实验证和 benchmark 做到能说服人；
5. 把因果措辞全部降级成 association / annotation / prioritization。

你现在最危险的不是代码不够，而是 claim 仍然比证据大。

建议把论文定位改成：

MoDES is a multi-omics evidence-scoring and event-state annotation framework for candidate regulatory links.

不要写：

MoDES is a validated multi-omics regulatory inference framework.

这两个句子差别非常大。前者能过大修，后者很容易继续被拒。

⸻

一、必须重写 state grammar

审稿人指出的第一个大问题很准确：你现在的 RNA+ATAC core states 偏 activation-centric。

如果你的 concordant 只覆盖：

ATAC ↑, RNA ↑

那它漏掉了：

ATAC ↓, RNA ↓

而后者在很多真实数据中也很常见，比如 enhancer closing + gene downregulation。

1.1 修 RNA+ATAC core state

建议不要只用 concordant 一个名字。改成方向明确的 state：

concordant_activation:
  ATAC ↑, RNA ↑
concordant_repression:
  ATAC ↓, RNA ↓
discordant_opening_repression:
  ATAC ↑, RNA ↓
discordant_closing_activation:
  ATAC ↓, RNA ↑
chromatin_open_primed:
  ATAC ↑, RNA not significant
chromatin_closed_primed:
  ATAC ↓, RNA not significant
rna_up_only:
  RNA ↑, ATAC not significant
rna_down_only:
  RNA ↓, ATAC not significant
null:
  neither significant

如果你想保留旧名字，也可以在输出里加一个更高层字段：

state_family = concordant / discordant / chromatin_primed / rna_only / null
state = concordant_activation / concordant_repression / ...

这样兼容旧版本，也解决审稿人的方向性批评。

⸻

1.2 给 StateRule 加 direction relation

现在每个 rule 如果都是固定 +1 或 -1，会让状态表爆炸。更优雅的是支持：

same_direction
opposite_direction
any_direction

例如：

StateRule(
    name="concordant",
    required_relation=[
        DirectionRelation("atac", "rna", relation="same")
    ],
)

但从工程速度考虑，我建议先用显式 state：

concordant_activation
concordant_repression
...

更容易写测试，也更容易解释。

⸻

二、修 neutral / absent / forbidden / missing 的语义

审稿人指出 protein_buffered 和 mark_only 的逻辑太松，这是非常重要的问题。

现在如果一个 rule 写：

required: RNA ↑
neutral: protein

这不等于：

RNA ↑ and protein unchanged

它可能只是：

protein missing
protein noisy
protein weak
protein significant but not used

所以必须把 evidence 状态拆开。

2.1 新增四类 evidence condition

建议 StateRule 支持这几种条件：

required_significant:
  某模态必须显著，且方向匹配
required_absent:
  某模态必须不显著
forbidden_significant:
  某模态不能显著，或者不能以某方向显著
optional_support:
  有则加分，没有不扣分
missing_allowed:
  该模态缺失时仍可分类，但降低 assignment score

例如：

protein_buffered:
  required_significant: RNA ↑
  required_absent: protein

而不是：

required: RNA ↑
neutral: protein

mark_only 应该是：

required_significant: active_mark ↑
required_absent: RNA
required_absent: ATAC

active_enhancer_primed 应该是：

required_significant: active_mark ↑
optional_support: ATAC ↑
required_absent: RNA

这样 active_enhancer_primed 和 mark_only 才能区分。

⸻

2.2 缺失模态不能等同于“不显著”

这是多组学里非常关键的点。

protein not measured

不能被解释成：

protein unchanged

所以状态规则里必须区分：

not significant because measured and p >= threshold
missing because modality absent

建议 evidence 表里增加：

modality_available
measured
significant
missing_reason

如果 protein 不存在，那么不应触发 protein_buffered 或 protein_memory。

⸻

三、state_support_pval / qval 要降级或重新定义

审稿人现在已经从“event_fdr 完全错”退一步，承认 state_support_pval 是进步。但他们继续抓 post-selection 和 directed p-value。这个必须正面处理。

3.1 不要叫 p-value / q-value，除非你证明它们有校准意义

如果 directed_pvalue() 本身文档写了不是严格 one-sided test，那你不能再把它放进 BH 后叫 q-value。

建议改名：

state_support_score
state_support_rank
state_support_bh_score

或者更折中：

state_support_pseudo_p
state_support_pseudo_q

但最好不要出现 pval/qval 这种强统计词。

我最推荐：

state_support_score
state_support_adjusted_score

文档写：

State support scores are ranking-oriented evidence scores. They are not formal post-selection p-values and should not be interpreted as calibrated FDR.

⸻

3.2 如果你坚持保留 qval，就必须做 calibration

那就需要 benchmark 证明：

under null:
  selected state false discovery rate is controlled or at least empirically bounded

最低限度要做：

condition label permutation
random peak-gene links
random protein-gene links
random CUT&Tag-region links
sample label shuffle per modality

然后报告：

fraction of selected states at q < 0.05 / 0.1

如果不做这个，就别叫 qval。

⸻

四、主输出 schema 必须拆成固定主表 + 长格式证据表

审稿人说得很对：你不能一边说 event_table schema frozen，一边动态追加 {mod}_coef/{mod}_pval。

4.1 固定 event_table.tsv

主表只保留固定列：

event_id
region_id / peak_id
gene
context
tf_name
link_source
link_score
state_family
state
state_assignment_score
state_support_score
supporting_modalities
neutral_modalities
conflicting_modalities
artifact_risk
artifact_reason
quality_score

不要再在主表里动态追加：

cuttag_h3k27ac_coef
protein_pval
spatial_neighbor_effect
...

⸻

4.2 所有模态证据进入 event_modality_evidence.tsv

这个表应该是 v2.0 的核心。

字段：

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
model_used
converged
warning

这样 reviewer 看到就会明白：
MoDES 真正的 multi-omics 部分是长格式 evidence，而不是一堆动态列拼到主表里。

⸻

4.3 conditional 结果单独放 conditional_effects.tsv

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
warning

这样你可以诚实地说：

Conditional models are diagnostic layers. They are reported separately and optionally used in state scoring.

不要把它们暗示成所有状态的核心 inferential engine，除非你真的把它们纳入 scoring。

⸻

五、conditional decomposition 必须进入 state score，或明确降级成 diagnostic

现在审稿人说：

multi-modal conditional models 被输出，但没有真正影响 state assignment。

这个批评很合理。

你有两个选择。

选择 A：降级

文档写：

Conditional decomposition is an optional diagnostic module.
State assignment is based on marginal cross-modality evidence.

这很诚实，也容易过审。

选择 B：纳入 scoring

例如：

active_mark_concordant:
  required: H3K27ac ↑, RNA ↑
  optional diagnostic: RNA condition effect attenuated after H3K27ac adjustment

assignment score 可以加一项：

if attenuation_after_mark > threshold:
    score *= 1.2

但要小心，不要说这是因果证据。只能说：

conditional attenuation-supported

我建议大修阶段用选择 A。先过稿，不要把模型搞太复杂。

⸻

六、修 protein 和 CUT&Tag 的建模细节

6.1 Protein 不能用 ATAC library size 做 offset

审稿人抓得很准。Protein/ADT 应该有自己的 size factor。

建议：

protein_size_factor = protein total counts per sample

或者：

CLR / centered log ratio for ADT-like data

至少你不能默认用 ATAC library size。

实现上：

if modality.assay == "PROTEIN":
    offset = log(protein_counts.sum(axis=1))
elif modality.feature_type == "region":
    offset = log(modality_counts.sum(axis=1))
elif modality.assay == "RNA":
    offset = log(rna_counts.sum(axis=1))

同时文档写：

Protein modeling is experimental and uses protein-specific library-size normalization.

⸻

6.2 CUT&Tag / CUT&RUN / ChIP-seq 不能完全复用 ATAC 假设

CUT&Tag/CUT&RUN/ChIP 的 region signal 和 ATAC 不同：

histone mark 有 broad/narrow
background 不同
peak width 不同
target specificity 不同

短期可以仍用 NB GLM，但必须写：

generic count-based marginal model

不要写：

assay-optimized model

同时至少要支持：

peak width covariate / offset
blacklist flag
region_match_score
target-specific role

⸻

七、spatial 证据必须拆成不同 modality/role

审稿人指出 spatial evidence 可能同一列重复满足多个 required evidence。这个要修。

不要只有：

modality = spatial

而应该有：

spatial_moran
spatial_neighbor
spatial_edge

或者：

modality = spatial
role = spatial_autocorrelation
modality = spatial
role = neighbor_effect
modality = spatial
role = edge_artifact

_resolve_modality_evidence() 必须按：

modality + role

匹配，而不是只按 modality。

spatial_niche_driven 应该要求：

required_significant: spatial_neighbor
optional_support: spatial_moran
forbidden_significant: spatial_edge_artifact

而不是两个 generic spatial evidence。

⸻

八、benchmark 必须从“规则同构”升级成 stress + validation

你现在 benchmark 已经比以前多，但 reviewer 仍说“按规则生成，再按规则识别”。要破这个问题，必须加入和规则不完全同构的 benchmark。

8.1 必须增加的 benchmark

Null benchmark

condition label shuffle
modality sample shuffle
random links

目标：

non-null rate 低
support_score 不应过高
artifact flag 不应乱触发

Link-noise benchmark

true links: 100 / 75 / 50 / 25%
random links: 0 / 25 / 50 / 75%

目标：

state accuracy 随 link noise 合理下降
random links 不应产生强 biological state enrichment

Confounded benchmark

batch partially confounded
donor imbalance
low replicate
library size outlier

目标：

false positive rate
artifact_risk enrichment
model warnings

Weak-effect benchmark

effect size small / medium / strong

目标：

power curve
state recovery curve

Ablation benchmark

full MoDES
no conditional
no extra modality
no artifact
no EB
naive overlap
random links

目标：

证明每个组件有价值

⸻

8.2 必须补真实数据验证

你现在不能再用 PBMC pseudo-condition + spike-in 当 real validation。

最小需要一个真实 contrast：

stimulation vs control
treated vs untreated
time-course early vs late
differentiation stage A vs B

如果暂时找不到多模态数据，至少要有：

RNA+ATAC 真实 contrast
external validation:
  known pathways
  known marker genes
  known enhancer-gene links
  ChIP/CUT&Tag overlap
  published TF programs

CUT&Tag/protein 的真实验证可以是 supplement，但 RNA+ATAC core 的真实验证必须有。

⸻

九、文档必须彻底更新到 v2.0

这不是小事。现在 Major Revision 明确要求文档同步。

你要重写这些文件：

docs/statistical_model.md
docs/output_schema.md
docs/install_review.md
CITATION.cff
README.md
ROADMAP.md
CHANGELOG.md

9.1 statistical_model.md

必须覆盖：

1. v2.0 architecture
2. event_modality_evidence long table
3. StateRule grammar
4. required / absent / forbidden / missing semantics
5. state_assignment_score
6. state_support_score
7. why support score is not formal FDR
8. conditional_effects as diagnostic
9. experimental modalities
10. limitations

不要继续写 v0.1 RNA+ATAC 五步流程作为主文档。

⸻

9.2 output_schema.md

必须写：

event_table.tsv fixed schema
event_modality_evidence.tsv long schema
conditional_effects.tsv schema
deprecated compatibility fields:
  state_confidence
  event_pval
  event_fdr

并写清楚：

state_confidence is deprecated alias of state_assignment_score.
event_pval/event_fdr are deprecated aliases or legacy ranking fields.

⸻

9.3 CITATION.cff

改成当前版本：

version: 2.0.0

或者即将重投版本。

⸻

9.4 install_review.md

重新跑 fresh install，更新：

version
test count
CI status
examples
benchmarks

⸻

十、生物学措辞必须降级

这类句子要改：

Local chromatin opening drives transcriptional activation
Complete regulatory chain activation through to protein output
Protein memory persists
TF regulator

改成：

Local chromatin accessibility and RNA abundance change concordantly
ATAC, RNA, and protein layers show concordant differential signal
Protein-layer signal persists relative to RNA signal
candidate TF motif annotation

这一步很重要。Reviewer 已经明确说这不是文字洁癖，是机制主张过强。

⸻

十一、给你一份大修任务清单

按优先级执行。

P0：必须修，不修下一轮仍会被拒

1. 重写 StateRule：
   - 支持 up/down 双方向
   - 支持 required_significant / required_absent / forbidden / missing
   - 修 active_enhancer_primed vs mark_only
   - 修 protein_buffered / protein_memory 逻辑
2. 把 event_table 固定化，extra modalities 移到 event_modality_evidence.tsv
3. 把 state_support_pval/qval 改名或降级：
   - state_support_score
   - state_support_adjusted_score
   或明确 qval only for ranking
4. 不再把 directed_pvalue 当正式 p-value 使用，或完成 calibration
5. 更新 statistical_model.md / output_schema.md 到 v2.0
6. 修 protein normalization，不再用 ATAC library size
7. spatial evidence 拆成 moran / neighbor / edge roles
8. 降级 README 所有因果措辞

⸻

P1：大修验证必须补

9. Null calibration benchmark
10. Link-noise benchmark
11. Batch/donor-confounded benchmark
12. Weak-effect benchmark
13. Ablation benchmark
14. 真实 biological contrast 数据
15. external validation：
    - known pathway
    - known markers
    - published enhancer-gene links
    - ChIP/CUT&Tag overlap or TF program

⸻

P2：文档和 release 统一

16. CITATION.cff version 更新
17. install_review 更新
18. ROADMAP 改掉 overclaim
19. GitHub tag message 改掉 “Full multi-modal platform”
20. README dependency 加 scikit-learn
21. README clearly says experimental modalities are not validated inference layers

⸻

十二、如果你要写 rebuttal，可以这么写

你可以在 response letter 里这样承认并回应：

We agree that the previous v2.0 revision made the software modality-aware but did not sufficiently clarify the statistical interpretation of state support values. In the revised version, we separate fixed event-level summaries from long-format modality evidence, replace priority-based rules with explicit StateRule conditions including required, absent, forbidden and missing evidence, and rename state confidence fields to assignment scores. We also remove claims of formal post-selection FDR and provide calibration and null benchmarks. Conditional decompositions are now presented as diagnostic models rather than causal evidence unless explicitly validated.

这段比“我们已经支持多模态了”有力得多。

⸻

十三、最终判断

这次 Major Revision 是好消息。它说明：

审稿人已经承认 MoDES 有软件架构价值；
现在卡的是统计解释、规则逻辑、文档同步和真实验证。

你要做的不是再堆更多模态，而是：

把已有多模态从 interface 变成 evidence framework；
把 state rules 从优先级 if/else 变成严谨 grammar；
把 p/q/confidence 全部降级或校准；
把 benchmark 从 toy recovery 变成 stress validation；
把真实数据从 spike-in 变成 biological contrast。

如果这些完成，下一轮有机会从 Major Revision 走向接收；如果只是继续加功能，下一轮仍然会被拒。
