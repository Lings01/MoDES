Lings，如果你要一口气把 MoDES 推到 1.0，我建议不要把 1.0 理解成“多加几个功能”，而要理解成：

API 稳定、结果可信、真实数据验证充分、bulk / single-cell pseudobulk / spatial 至少有清楚边界、文档可复现、CI/测试/benchmark 完整。

当前 f9fdb89 已经是一个不错的 v0.1.0-alpha 候选：commit 说明里写了 strict integration test、CI、README、requirements split、StateClassifier 单元测试、import smoke test 等；README 当前也把范围写成 RNA+ATAC、binary condition、bulk / pseudobulk prototype，protein 和 spatial graph 仍是 planned。也就是说，现在的核心不是继续修 alpha，而是进入 1.0 产品化 + 方法学验证 阶段。 ￼

下面我按非常细的任务清单给你列。

⸻

0. 先定义 MoDES 1.0 到底是什么

0.1 1.0 的最低定义

MoDES 1.0 应该至少做到：

1. RNA+ATAC bulk 可稳定使用
2. RNA+ATAC single-cell pseudobulk 可稳定使用
3. spatial 数据有正式支持方式：
   - 要么 native spatial graph
   - 要么明确支持 region-level / spot-group pseudobulk
4. 输入输出格式稳定
5. state classification 稳定
6. artifact risk 体系稳定
7. 至少 2 个真实数据 demo
8. 至少 1 个系统 benchmark
9. CI 全绿
10. 文档完整
11. API 不随便破坏

按照 Semantic Versioning，1.0.0 通常表示 public API 已经稳定，之后的补丁和小版本更新不能随意破坏已有用户代码。 ￼

所以你现在要先确定：

MoDES 1.0 是 RNA+ATAC 稳定版？
还是 RNA+ATAC+Protein+Spatial 都稳定？

我建议现实一点：

MoDES 1.0 = RNA+ATAC bulk + single-cell pseudobulk + spatial region-pseudobulk 稳定版
MoDES 1.1 / 1.2 再加 protein 和 native spatial graph

但如果你坚持 1.0 必须包含 protein 和 native spatial，我也在下面列了完整任务。

⸻

1. Release / 版本管理任务

1.1 立刻打 alpha tag

git tag v0.1.0-alpha
git push origin v0.1.0-alpha

1.2 新建 CHANGELOG.md

内容结构：

# Changelog
## v0.1.0-alpha
### Added
- RNA+ATAC regulatory event state inference
- Biological states: concordant, chromatin_primed, rna_only, discordant_opposite, null
- artifact_risk / artifact_reason
- event_pval / event_fdr
- model_diagnostics.tsv
- GraphML export
- HTML report
- minimal bulk example
- CI workflow
### Limitations
- binary condition only
- RNA+ATAC only
- bulk / pseudobulk recommended
- no native protein layer
- no native spatial graph

1.3 新建 ROADMAP.md

建议版本路线：

v0.1.0-alpha: RNA+ATAC bulk / pseudobulk MVP
v0.2.0-alpha: real-data smoke tests + benchmark
v0.3.0-alpha: single-cell pseudobulk workflow hardening
v0.4.0-alpha: spatial region-pseudobulk support
v0.5.0-beta: documentation + API freeze
v0.8.0-beta: real benchmark + method report
v1.0.0: stable RNA+ATAC release

1.4 明确 public API

1.0 前要冻结这些 API：

from modes import MoDES, MoDEData
MoDEData.from_matrices(...)
MoDEData.from_pseudobulk(...)
MoDES(...)
MoDES.run()
MoDES.build_events()
MoDES.estimate_effects()
MoDES.decompose()
MoDES.build_evidence()
MoDES.classify_states()
MoDESResult.to_tsv()
MoDESResult.to_graphml()
MoDESResult.to_report()
MoDESResult.filter()

1.5 写 API stability policy

文档中写：

From v1.0.0 onward, MoDES will follow semantic versioning.
Breaking API changes require a major version bump.

⸻

2. 当前代码最后一轮硬化

你现在的 alpha 基本能用，但推到 1.0 前还要补这些。

2.1 修 requirements.txt 文件格式

我这里重新打开 f9fdb89 的 raw requirements.txt，它仍显示为一行空格分隔依赖，而不是一行一个依赖。 ￼

必须改成：

numpy>=1.21
scipy>=1.7
pandas>=1.3
statsmodels>=0.13
anndata>=0.8
matplotlib>=3.5
seaborn>=0.11
networkx>=2.6

requirements-dev.txt：

pytest>=7.0
pytest-cov>=4.0
ruff>=0.5
black>=24.0
mypy>=1.8
build>=1.0
twine>=5.0

2.2 修 README raw 格式

我打开 f9fdb89 的 raw README，仍然只有 8 行，很多 Markdown 内容挤在长行里。 ￼

需要真正保存为多行 Markdown，尤其是这些部分：

Overview
Current status
Installation
Quick start
Output files
event_table.tsv fields
Input formats
Algorithm
Limitations
Roadmap

2.3 检查 CI 是否真的跑当前最新 commit

你说 CI 已通过，我按你的结果接受。仓库 Actions 页面我这里仍能看到旧的 9981707 workflow run 记录，页面加载不完整，不能作为我这边的确认依据。 ￼

你本地确认后，把 README 顶部加 badge：

![tests](https://github.com/Lings01/MoDES/actions/workflows/tests.yml/badge.svg)

2.4 加 coverage

CI 里加：

python -m pytest --cov=modes --cov-report=term-missing

目标：

v0.1-alpha: >=70%
v0.5-beta: >=80%
v1.0: >=85%

2.5 加 lint

新增 .ruff.toml 或 pyproject.toml：

[tool.ruff]
line-length = 100
target-version = "py310"
[tool.ruff.lint]
select = ["E", "F", "I", "B", "UP"]

CI 加：

ruff check modes tests

2.6 加格式化检查

black --check modes tests

2.7 加 typing 检查

mypy modes

先允许宽松：

[tool.mypy]
ignore_missing_imports = true
disallow_untyped_defs = false

2.8 加 package build 检查

python -m build
twine check dist/*

2.9 加 fresh install smoke test

CI 里加：

python -m venv /tmp/modes-test
source /tmp/modes-test/bin/activate
pip install dist/*.whl
python -c "from modes import MoDES, MoDEData; print('ok')"

⸻

3. 数据输入层要做到 1.0 稳定

3.1 MoDEData.from_matrices() 稳定化

必须测试这些输入：

TSV path
CSV path
pandas DataFrame
sample index 顺序不同
RNA 有 sample，ATAC 缺 sample
metadata 缺 sample
重复 sample ID
重复 gene name
重复 peak ID
空矩阵
全 0 sample
非整数 count
NaN count
负数 count

每个情况都要有明确报错。

3.2 增加输入验证函数

新增：

MoDEData.validate()

检查：

RNA 和 ATAC sample 完全对齐
metadata index 完全对齐
condition_col 存在
condition 正好 2 类
RNA counts 非负
ATAC counts 非负
无 NaN
无 infinite
每个 sample RNA library size > 0
每个 sample ATAC library size > 0
peak IDs 唯一
gene IDs 唯一

3.3 支持 sparse matrix

单细胞/空间数据非常稀疏，必须支持：

scipy.sparse.csr_matrix
scipy.sparse.csc_matrix

当前如果全部转 dense，大数据会爆内存。

需要：

MoDEData.rna
MoDEData.atac

允许：

pd.DataFrame
scipy sparse + var_names + obs_names

或者先内部标准化为 AnnData-like。

3.4 支持 MuData / h5mu

如果你要认真做 single-cell multiome，最终要支持：

MuData
AnnData RNA modality
AnnData ATAC modality

新增：

MoDEData.from_mudata(
    mdata,
    rna_mod="rna",
    atac_mod="atac",
    condition_col="condition",
    donor_col=None,
    batch_col=None,
)

3.5 支持 10x Multiome 文件夹

新增：

MoDEData.from_10x_multiome(
    filtered_feature_bc_matrix_h5,
    fragments_tsv=None,
    metadata=None,
    ...
)

第一版可以只支持已经处理好的 RNA/ATAC matrix，不碰 fragments。

⸻

4. Event candidate construction 要稳定

4.1 external links schema 固化

要求列：

peak_id
gene

可选列：

tf_name
source
score
distance
link_type

新增 validator：

validate_external_links(links)

检查：

peak_id 列存在
gene 列存在
peak_id 在 ATAC columns 中
gene 在 RNA columns 中
重复 peak-gene link 如何处理
score 是否在 0-1

4.2 external links 过滤策略

参数：

min_link_score=None
allowed_sources=None
drop_missing_features=True

4.3 peak ID 解析更稳

支持：

chr1:100-200
chr1-100-200
chr1_100_200
chr1:100..200

建议只官方支持一种，其他给 warning。

4.4 GTF / GFF 解析

现在如果用户不用 external links，就靠 gene 坐标。1.0 要让 GTF 解析可靠：

gene_name
gene_id
transcript_id
chr
start
end
strand
TSS

参数：

gene_id_col="gene_name"
feature_type="gene"

4.5 TSS map 输入格式文档化

允许：

tss_map = {
    "STAT1": ("STAT1", "chr2", 191000000)
}

但要文档写清楚。

4.6 candidate event 去重

同一个：

peak_id + gene + tf_name

重复时处理：

保留最高 score
合并 source
记录 n_sources

4.7 event_id 稳定

不要用随机 ID。建议：

event_id = "{peak_id}|{gene}|{tf_name or NA}"

或者 hash：

sha1(f"{peak_id}|{gene}|{tf_name}").hexdigest()[:12]

1.0 后 event_id 格式不要再变。

⸻

5. 统计模型要到 1.0 可信

5.1 binary condition 明确

1.0 如果仍然只支持 binary condition，也可以，但必须写清楚。

参数：

reference_condition="control"
target_condition="treatment"

不要默认按 sorted categories，因为：

case/control
treated/untreated
disease/healthy

排序可能不符合用户预期。

新增：

MoDES(..., contrast=("treatment", "control"))

输出里记录：

coef = treatment vs control

5.2 direction 定义固定

你现在已改成 coefficient sign，这个对。1.0 文档写：

direction is sign(coefficient), independent of FDR.
significance is determined separately by FDR.

5.3 NB GLM backend 说明准确

当前是：

statsmodels NB GLM with fixed/default alpha

不要夸大为：

fully estimated gene-wise dispersion

文档写：

MoDES uses a practical NB GLM fallback hierarchy:
1. negative binomial GLM
2. fixed-alpha NB
3. Poisson fallback
4. simplified fallback if enabled

5.4 fallback 策略参数化

新增：

allow_poisson_fallback=True
allow_simplified_fallback=False

默认建议：

poisson fallback: True
simplified fallback: False

因为 simplified fallback 会丢 covariates，风险更大。

5.5 rank deficient design 处理

当前已经会报错。1.0 还要输出更好的错误：

Design matrix is rank deficient.
Possible causes:
- condition fully confounded with donor
- condition fully confounded with batch
- too many covariates for sample size

5.6 donor/batch 建模说明

当前 donor/batch 是 fixed effects。文档写清楚：

donor_col and batch_col are currently modeled as fixed-effect covariates.
Mixed-effect models are planned.

5.7 event-level FDR 定义固定

文档写清楚：

concordant: max(atac_pval, rna_pval)
discordant_opposite: max(atac_pval, rna_pval)
chromatin_primed: atac_pval
rna_only: rna_pval
null: 1.0

5.8 state_confidence 定义固定

不能叫 posterior。文档写：

state_confidence is an empirical confidence score derived from rule-based initialization and optional empirical Bayes refinement.
It should not be interpreted as a calibrated posterior probability unless calibration is performed.

5.9 artifact_risk 定义增强

当前：

low / medium / high

1.0 要支持原因：

low_quality_score
single_modality_low_quality
low_atac_depth
low_rna_depth
batch_associated
library_size_outlier
weak_link_score

可以先实现前 2-3 个，其他 planned。

5.10 quality_score 分解

不要只有一个 opaque 分数。输出：

quality_score
depth_score
detection_score
batch_score
link_score

⸻

6. Single-cell pseudobulk 推到 1.0

6.1 明确 1.0 支持方式

写成：

MoDES 1.0 supports single-cell multiome through donor-aware pseudobulk analysis.
Native cell-level mixed models are not part of 1.0.

6.2 from_pseudobulk() 稳定

输入：

MoDEData.from_pseudobulk(
    adata,
    groupby=["donor", "condition", "cell_type"],
    condition_col="condition",
    donor_col="donor",
    batch_col="batch",
    min_cells_per_group=20,
    rna_layer=None,
    atac_layer=None,
)

测试：

groupby 列缺失
condition_col 缺失
donor_col 缺失
某些 group cell 数不足
ATAC layer 不存在
RNA/ATAC 细胞顺序不一致
empty group

6.3 pseudobulk metadata 输出

每个 pseudobulk sample 应有：

group_id
donor
condition
cell_type
batch
n_cells
rna_total_counts
atac_total_counts

6.4 min_cells_per_group 策略

默认：

min_cells_per_group=20

低于阈值：

drop + warning

输出：

dropped_groups.tsv

6.5 多 cell type 批量分析

新增 helper：

run_by_context(
    data,
    context_col="cell_type",
    ...
)

输出：

context/event_table.tsv
all_contexts_event_table.tsv

6.6 context 字段稳定

event_table 中：

context = cell_type or user-defined context

6.7 单细胞 example

新增：

examples/singlecell_pseudobulk/
  README.md
  run_pseudobulk.py
  make_pseudobulk_from_anndata.py

6.8 单细胞文档

写：

Do not treat cells as independent biological replicates.
Use donor × condition × cell_type pseudobulk.

⸻

7. Spatial 支持到 1.0 的两个选择

这里你必须做决定。

选择 A：1.0 只支持 spatial region-pseudobulk

这是现实路线。

定义：

Spatial support in 1.0 means:
spatial spots/cells can be aggregated by sample × region/niche × condition,
then analyzed as bulk-like RNA+ATAC matrices.
Native spatial graph modeling is planned for 1.1+.

需要做：

MoDEData.from_spatial_pseudobulk(...)

输入：

RNA matrix
ATAC matrix
metadata with region/niche
optional coordinates

但 v1.0 不使用 coordinates 建模，只保存 metadata。

输出 context：

context = spatial_region

选择 B：1.0 支持 native spatial graph

如果你坚持 1.0 必须 native spatial，需要做很多。

7B.1 新增 SpatialMoDEData

class SpatialMoDEData(MoDEData):
    coords: pd.DataFrame
    spatial_graph: scipy.sparse.csr_matrix
    region_labels: Optional[pd.Series]

7B.2 坐标输入

支持：

spot_id, x, y

7B.3 自动建图

参数：

graph_method="knn" or "radius"
n_neighbors=6
radius=None

7B.4 空间 evidence

新增：

spatial_autocorrelation
region_effect
neighbor_effect
edge_artifact_score

7B.5 空间状态

新增：

spatial_region_specific
spatial_niche_driven
spatial_artifact_edge

7B.6 空间输出

spatial_event_table.tsv
spatial_event_maps.h5ad
spatial_graph.graphml

7B.7 空间测试

known spatial cluster event recovery
randomized coordinates negative control
edge artifact detection

我的建议：1.0 用选择 A，native spatial 放 1.1 或 2.0。

⸻

8. Protein 层任务

如果 1.0 要包含 protein，就做这些。否则放到 1.1。

8.1 数据结构扩展

MoDEData 增加：

protein: Optional[pd.DataFrame]
protein_names

8.2 输入

MoDEData.from_matrices(
    rna_counts=...,
    atac_counts=...,
    protein_counts=...,
)

8.3 protein feature linking

protein 和 gene 的 link：

protein_id
gene
protein_name

例：

CD3D_ADT -> CD3D
CD4_ADT -> CD4

8.4 protein effect estimation

新增：

protein_effects = estimator.estimate_protein_effects(...)

模型：

ADT/protein count: NB GLM
bulk proteomics: Gaussian linear model

8.5 conditional protein decomposition

Protein ~ Condition + RNA + ATAC + Covariates

8.6 新 evidence vector

D_e = [z_ATAC, z_RNA, z_RNA|ATAC, z_Protein, z_Protein|RNA, quality]

8.7 新状态

full_activation: ATAC↑ RNA↑ Protein↑
protein_buffered: RNA↑ Protein→
protein_memory: RNA→ Protein↑
protein_opposite: RNA↑ Protein↓

8.8 输出字段

protein_coef
protein_pval
protein_fdr
protein_direction
protein_after_rna_coef
protein_after_rna_pval
protein_after_rna_fdr

8.9 文档

写清楚：

Protein layer is optional.
If protein is absent, MoDES-RA states are used.
If protein is present, MoDES-RAP states are used.

我的建议：不要把 protein 放进 1.0；放到 1.1。

⸻

9. Benchmark 必须做

1.0 不能只有单元测试。必须有 benchmark。

9.1 synthetic benchmark

目录：

benchmarks/simulated_event_states/

模拟状态：

concordant
chromatin_primed
rna_only
discordant_opposite
null
high_artifact_risk

输出：

truth.tsv
predicted_event_table.tsv
metrics.tsv
confusion_matrix.png

指标：

accuracy
macro_F1
per-state precision
per-state recall
artifact_risk precision
artifact_risk recall
event_fdr calibration
runtime
peak_memory

9.2 semi-real benchmark

用真实 count matrix 的 library size / sparsity，然后 spike-in known effects。

步骤：

1. 读真实 RNA+ATAC matrix
2. 抽取 peak-gene links
3. 随机指定 truth states
4. 对 counts 做 controlled perturbation
5. 跑 MoDES
6. 评估 state recovery

9.3 negative control benchmark

shuffle condition labels
shuffle peak-gene links
shuffle ATAC matrix sample labels

预期：

多数 event → null
event_fdr 受控
artifact_risk 不应乱升

9.4 baseline 对照

至少做：

baseline_overlap:
  ATAC significant + RNA significant = concordant
  ATAC only = primed
  RNA only = rna_only
baseline_correlation:
  peak-gene correlation + DE/DA
MoDES:
  current method

9.5 benchmark 脚本

python benchmarks/simulated_event_states/run_benchmark.py
python benchmarks/simulated_event_states/summarize.py

9.6 benchmark 文档

benchmarks/README.md

⸻

10. 真实数据验证

10.1 minimal example 保持

现在已有：

examples/minimal_bulk/

1.0 要加入 expected outputs：

examples/minimal_bulk/expected/event_table.tsv
examples/minimal_bulk/expected/report.html

10.2 real PBMC smoke test

目的：

真实 10x Multiome 数据能跑通

不要过度解释生物学。

输出：

notebooks/01_pbmc_multiome_smoke_test.ipynb

10.3 real biological demo

需要真正有 condition / time / differentiation。

选择一个：

SHARE-seq mouse skin
10x multiome stimulation dataset
hematopoiesis multiome
tumor treatment multiome

目标：

chromatin_primed events 是否富集 lineage TF
concordant events 是否富集成熟 marker
rna_only events 是否富集 trans/stress pathways

10.4 reproducibility demo

如果有 donor：

leave-one-donor-out
run MoDES
看 top event overlap

10.5 marker sanity check

输出：

known_marker_event_enrichment.tsv
motif_enrichment.tsv

⸻

11. 文档体系

11.1 README 保持简洁

README 只放：

what is MoDES
installation
quickstart
input/output
current limitations
citation

11.2 docs/

新增：

docs/
  installation.md
  input_formats.md
  output_schema.md
  singlecell_pseudobulk.md
  spatial_region_pseudobulk.md
  statistical_model.md
  artifact_risk.md
  benchmark.md
  faq.md

11.3 output schema 文档

每个输出文件都列：

字段名
类型
范围
含义
是否可空

11.4 FAQ

必须回答：

Q: Can I run MoDES on cells directly?
A: Not recommended. Use donor-aware pseudobulk.
Q: Can MoDES run on RNA-only spatial transcriptomics?
A: Not in v1.0. It requires RNA+ATAC or region-level paired inputs.
Q: Is state_confidence a posterior probability?
A: No, it is an empirical confidence score.
Q: Does MoDES infer peak-gene links?
A: It can generate proximity candidates, but external links are recommended.
Q: Is MoDES a multiome integration method?
A: No, it is an event-state decomposition method.

⸻

12. CLI 命令行接口

1.0 应该有 CLI，否则只有 Python API 不够像工具。

12.1 新增 entry point

setup.py：

entry_points={
    "console_scripts": [
        "modes=modes.cli:main",
    ],
}

12.2 CLI 子命令

modes run
modes validate-input
modes build-events
modes benchmark
modes report

12.3 modes run

modes run \
  --rna rna_counts.tsv \
  --atac atac_counts.tsv \
  --metadata metadata.tsv \
  --condition condition \
  --external-links peak_gene_links.tsv \
  --out output/

可选：

--donor donor
--batch batch
--covariates age,sex
--fdr-threshold 0.1
--max-event-fdr 0.1
--exclude-high-artifact

12.4 modes validate-input

输出：

input_validation_report.txt

12.5 CLI 测试

subprocess.run(["modes", "--help"])
subprocess.run(["modes", "run", ...])

⸻

13. Packaging

13.1 迁移到 pyproject.toml

setup.py 可以保留，但 1.0 推荐 pyproject.toml。

内容：

[project]
name = "modes-bio"
version = "1.0.0"
description = "Multi-omics discordance-guided regulatory event state inference"
requires-python = ">=3.10"
dependencies = [...]

13.2 包名检查

modes 这个名字可能太泛，PyPI 可能冲突。建议包名：

modes-bio

但 import 仍然：

import modes

13.3 添加 license

README 写 MIT，但仓库里要有：

LICENSE

13.4 添加 citation

CITATION.cff

13.5 添加 manifest

MANIFEST.in

包含 examples、README、LICENSE。

13.6 构建测试

python -m build
twine check dist/*

⸻

14. 性能与可扩展性

14.1 runtime profiling

数据规模：

small: 10 samples, 100 events
medium: 50 samples, 10k events
large: 200 samples, 100k events

记录：

runtime
peak memory
n_events/sec

14.2 并行化

可加：

n_jobs=1

用于：

effect estimation
conditional decomposition

用 joblib 或 multiprocessing。

14.3 cache

中间结果：

events.tsv
atac_effects.pkl
rna_effects.pkl
conditional.tsv
evidence.tsv

允许 resume：

modes.run(resume=True)

14.4 大矩阵处理

避免 .values 全部 dense 化。

⸻

15. 错误信息与用户体验

15.1 所有常见失败要有友好错误

例如：

No candidate events generated
Condition is not binary
Design matrix rank deficient
Zero library size
External links missing required columns
Peak IDs not found in ATAC matrix
Gene IDs not found in RNA matrix

15.2 所有 warning 要可追踪

输出：

warnings.log

15.3 report 里显示 warnings

HTML report 增加：

Input warnings
Model warnings
Dropped groups
Fallback models

⸻

16. 测试矩阵

16.1 单元测试

覆盖：

data loading
input validation
event construction
external links
TSS parsing
effect estimation
design matrix
rank deficiency
conditional decomposition
evidence vector
state classification
artifact risk
event FDR
GraphML export
HTML report
CLI

16.2 集成测试

minimal bulk
single-cell pseudobulk
external links
GTF links
with covariates
with donor
with batch

16.3 回归测试

保存 expected outputs。

tests/expected/minimal_bulk/event_table.tsv

不要要求浮点完全一致，用 tolerance。

16.4 CI matrix

Python 3.10
Python 3.11
Python 3.12
Ubuntu
macOS optional

16.5 dependency lower bound test

确保最低依赖版本可用。

16.6 dependency latest test

确保最新依赖版本可用。

⸻

17. 方法学报告 / manuscript

17.1 技术报告标题

MoDES: Discordance-aware regulatory event state inference from RNA+ATAC multiome data

17.2 核心图

Figure 1: Method overview
Figure 2: State definitions
Figure 3: Simulation benchmark
Figure 4: Real data smoke test
Figure 5: Biological demo
Figure 6: Runtime / scalability

17.3 Methods 写清楚

event candidate construction
marginal effect estimation
conditional decomposition
evidence vector
state classification
artifact risk
event-level FDR

17.4 Limitations 写清楚

binary condition only
single linked peak conditional model
fixed-effect donor/batch
RNA+ATAC only in 1.0
spatial graph not native if选择 region-pseudobulk
state_confidence not calibrated posterior

⸻

18. 1.0 Release Gate

你可以设一个 checklist。全部通过才打 1.0。

18.1 功能 gate

[ ] bulk RNA+ATAC works
[ ] single-cell pseudobulk works
[ ] spatial region-pseudobulk works, or explicitly out of scope
[ ] external_links works
[ ] genome_annotation works
[ ] output schema stable
[ ] CLI works

18.2 统计 gate

[ ] state definitions documented
[ ] event_fdr documented
[ ] state_confidence documented
[ ] artifact_risk documented
[ ] fallback model diagnostics exported
[ ] rank deficiency handled

18.3 测试 gate

[ ] pytest passes
[ ] coverage >=85%
[ ] CI passes Python 3.10/3.11/3.12
[ ] minimal example test passes
[ ] benchmark script passes

18.4 数据 gate

[ ] synthetic benchmark done
[ ] semi-real spike-in benchmark done
[ ] real PBMC smoke test done
[ ] biological demo done

18.5 文档 gate

[ ] README complete
[ ] docs/input_formats.md
[ ] docs/output_schema.md
[ ] docs/singlecell_pseudobulk.md
[ ] docs/statistical_model.md
[ ] docs/faq.md
[ ] CHANGELOG.md
[ ] LICENSE
[ ] CITATION.cff

18.6 packaging gate

[ ] pyproject.toml
[ ] build succeeds
[ ] twine check succeeds
[ ] fresh install succeeds

⸻

19. 我建议的 1.0 范围

我强烈建议 MoDES 1.0 不要贪。

1.0 包含

RNA+ATAC bulk
RNA+ATAC single-cell pseudobulk
spatial region-pseudobulk as bulk-like analysis
CLI
benchmark
real data demo
stable output schema

1.0 不包含

native cell-level mixed model
native spatial graph model
protein layer
multi-class condition
continuous pseudotime lag
calibrated Bayesian posterior
best-in-class peak-gene linking

这些放到：

v1.1 protein
v1.2 native spatial graph
v1.3 multi-condition / pseudotime
v2.0 full multi-layer event model

⸻

20. 最短执行路线

如果你想最快推进：

Week 1:
  fix requirements raw format
  fix README raw format
  tag v0.1.0-alpha
  add CHANGELOG/LICENSE/CITATION
  add CLI skeleton
Week 2:
  input validation hardening
  single-cell pseudobulk example
  docs/output_schema.md
  docs/singlecell_pseudobulk.md
Week 3:
  synthetic + semi-real benchmark
  confusion matrix / metrics
Week 4:
  real PBMC smoke test
  one biological demo
Week 5:
  package build
  coverage/lint/type checks
  v0.5-beta
Week 6:
  polish docs
  run all benchmarks
  freeze API
  tag v1.0.0

⸻

一句话总结

你现在不要再把精力放在“这个 alpha 还有没有小 bug”。从现在开始，MoDES 到 1.0 的核心工作是：

稳定 API
硬化输入输出
做 single-cell pseudobulk 正式支持
给 spatial 一个明确 1.0 范围
补 benchmark 和真实数据 demo
加 CLI / docs / packaging
冻结输出 schema

最推荐的 1.0 定义是：

MoDES 1.0 = 稳定的 RNA+ATAC regulatory event-state inference 工具，支持 bulk、single-cell pseudobulk、spatial region-pseudobulk，原生 protein 和 native spatial graph 放到 1.x 后续版本。
