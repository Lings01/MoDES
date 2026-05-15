Lings，我看了当前 master 的仓库完整性。结论是：

仓库已经从 MVP 原型进入“接近 1.0 工程化”的阶段，但现在不能打 1.0，因为当前最新 CI 是失败的。 GitHub Actions 最新 tests #14 对应 commit 3857c3c，状态是 Failure，并且有 lint、lower-bound、pytest 3.10/3.11/3.12 相关失败或取消记录。这个必须作为最高优先级修复。 ￼

好的一面是，你现在仓库结构已经明显完整了：根目录已有 .github/workflows、benchmarks、docs、examples、modes、notebooks、tests、CHANGELOG.md、CITATION.cff、LICENSE、MANIFEST.in、ROADMAP.md 等；docs/ 里也已经有 benchmark.md、faq.md、input_formats.md、installation.md、output_schema.md、singlecell_pseudobulk.md、statistical_model.md；benchmarks/ 里已有 simulated、semi-real、negative control、baseline comparison 四类目录；examples/ 里已有 minimal_bulk 和 singlecell_pseudobulk。 ￼

下面是我建议你接下来按优先级执行的详细修改列表。

⸻

0. 当前总体判断

当前 MoDES 已经有这些优势：

1. 核心代码模块完整：
   data / events / effects / decompose / states / core / report / plotting / cli
2. 测试目录完整：
   test_cli / test_core / test_data / test_decompose / test_effects /
   test_events / test_import / test_integration / test_states / test_utils
3. 文档目录已经建立：
   installation / input_formats / output_schema / statistical_model /
   singlecell_pseudobulk / benchmark / FAQ
4. benchmark 目录已经开始成型：
   simulated_event_states / semi_real_spikein / negative_control /
   baseline_comparison
5. examples 已经覆盖：
   minimal_bulk
   singlecell_pseudobulk
6. release 基础文件已有：
   LICENSE / CITATION.cff / CHANGELOG.md / ROADMAP.md / MANIFEST.in

但是现在最关键的问题是：

当前 master 最新 CI 失败；
benchmark/CI/lint 还没有稳定；
1.0 release gate 还没过；
spatial 仍然不是 native spatial graph；
protein layer 还没实现；
multi-condition / pseudotime 还没实现；
README 和部分文档仍有“愿景大于当前能力”的地方需要压实。

README 当前也仍把 MoDES-RA 定位为 v0.1.0-alpha — RNA + ATAC prototype，并明确 protein layer、spatial graph、multi-class condition、pseudotime delay 都是 planned，这个定位是准确的。 ￼

⸻

P0：立刻修，阻塞 1.0

1. 修复当前 GitHub Actions 失败

位置：

.github/workflows/tests.yml
tests/
benchmarks/
requirements*.txt

现状：

最新 Actions run tests #14 是 Failure；summary 里显示：

Status Failure
lint exit code 1
lower-bound exit code 1
pytest 3.10 exit code 4
pytest 3.11 canceled / failed
pytest 3.12 canceled

我看不到完整日志，因为 GitHub 要求登录查看 logs，但失败状态和失败 job 已经明确。 ￼

要做：

1. 先在本地跑：

python -m pytest -q
python -m pytest -q --cov=modes --cov-report=term-missing
ruff check modes/ tests/

2. 逐个修失败：

lint 失败 → 先修 ruff
lower-bound 失败 → 检查最低依赖版本是否真的支持当前代码
pytest exit code 4 → 通常是 pytest 参数、collection、配置或依赖问题

3. GitHub Actions 里不要同时加太多 gate。建议拆成：

jobs:
  unit-tests:
    python 3.10 / 3.11 / 3.12
  lint:
    ruff only
  lower-bound:
    暂时允许失败，等依赖边界稳定后再设为 required

4. 在 CI 稳定前，不要再追加新功能。

验收标准：

GitHub Actions tests #latest 全部绿色
Python 3.10 / 3.11 / 3.12 通过
lint 通过
lower-bound 要么通过，要么暂时标记 continue-on-error

⸻

2. 修 requirements.txt / requirements-dev.txt 的真实换行格式

位置：

requirements.txt
requirements-dev.txt

现状：

Web raw 视图显示 requirements.txt 和 requirements-dev.txt 都是 Total lines: 1。我不能排除是 GitHub/raw 显示折叠问题，但从审计角度看，这很危险；如果真实文件是一行空格分隔，pip install -r requirements.txt 很容易出问题。 ￼

要做：

requirements.txt 应该是：

numpy>=1.21
scipy>=1.7
pandas>=1.3
statsmodels>=0.13
anndata>=0.8
matplotlib>=3.5
seaborn>=0.11
networkx>=2.6

requirements-dev.txt 应该是：

pytest>=7.0
pytest-cov>=4.0
ruff>=0.5
black>=24.0
mypy>=1.8
build>=1.0
twine>=5.0

如果你暂时不用 black/mypy/build/twine，可以先只保留：

pytest>=7.0
pytest-cov>=4.0
ruff>=0.5

验收标准：

pip install -r requirements.txt
pip install -r requirements-dev.txt

在 fresh venv 里能通过。

⸻

3. 暂时降低 CI 的 lower-bound 风险

位置：

.github/workflows/tests.yml

问题：

当前 workflow 里有 lower-bound job；Actions summary 显示 lower-bound job 失败。 ￼

对于科学计算包，numpy/scipy/pandas/statsmodels/anndata 的低版本组合很容易互相不兼容。你现在还在快速开发阶段，不应该让 lower-bound 阻塞主开发。

要做：

先改成：

lower-bound:
  continue-on-error: true

或者暂时删除 lower-bound job。

等 1.0 前再做：

最低支持版本测试
最新依赖测试
Python 版本矩阵

验收标准：

主 pytest + lint 必须通过；
lower-bound 可以作为 warning，不阻塞 merge。

⸻

4. 修 lint 问题，建立风格基线

位置：

modes/
tests/
benchmarks/
examples/

现状：

当前最新 Actions summary 显示 lint job exit code 1。 ￼

要做：

1. 本地跑：

ruff check modes/ tests/ benchmarks/ examples/

2. 先自动修：

ruff check modes/ tests/ benchmarks/ examples/ --fix

3. 再手动修不能自动修的部分。
4. 新建 pyproject.toml 或 .ruff.toml：

[tool.ruff]
line-length = 100
target-version = "py310"
[tool.ruff.lint]
select = ["E", "F", "I", "B", "UP"]
ignore = [
  "E501",  # 如果暂时不想处理长行
]

验收标准：

ruff check modes/ tests/ benchmarks/ examples/

返回 0。

⸻

P1：1.0 前必须补齐的工程完整性

5. 迁移到 pyproject.toml

位置：

pyproject.toml
setup.py

现状：

仓库有 setup.py，但 pyproject.toml raw 打开失败或未正常存在；根目录列表没有显示 pyproject.toml。 ￼

要做：

新增 pyproject.toml：

[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"
[project]
name = "modes-bio"
version = "0.1.0-alpha"
description = "Multi-omics discordance-guided regulatory event state inference"
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
authors = [
  {name = "Lings"}
]
dependencies = [
  "numpy>=1.21",
  "scipy>=1.7",
  "pandas>=1.3",
  "statsmodels>=0.13",
  "anndata>=0.8",
  "matplotlib>=3.5",
  "seaborn>=0.11",
  "networkx>=2.6",
]
[project.scripts]
modes = "modes.cli:main"

setup.py 可以保留一段时间，但最终推荐由 pyproject.toml 作为主入口。

验收标准：

python -m build
pip install dist/*.whl
modes --help

都通过。

⸻

6. CLI 与 CHANGELOG 的矛盾要修

位置：

modes/cli.py
CHANGELOG.md
README.md
ROADMAP.md

问题：

仓库现在有 modes/cli.py，而且 cli.py 已经实现了 modes run 和 modes validate-input 的基本逻辑；但是 CHANGELOG.md 的 limitation 里仍写着 “No CLI interface”。这已经不一致。 ￼

要做：

1. CHANGELOG.md 改成：

CLI interface is experimental:
- modes run
- modes validate-input

2. ROADMAP.md 里 v0.5.0-beta - CLI interface 改成：

CLI hardening and API freeze

3. README 增加：

modes run \
  --rna rna_counts.tsv \
  --atac atac_counts.tsv \
  --metadata metadata.tsv \
  --condition condition \
  --external-links peak_gene_links.tsv \
  --out output \
  --report \
  --network

4. 给 CLI 加测试：

tests/test_cli.py

覆盖：

modes --help
modes run on minimal example
modes validate-input
missing required args
invalid external links

验收标准：

modes --help
modes validate-input --help
modes run --help

全部可用，README/CHANGELOG/ROADMAP 不冲突。

⸻

7. 为 validate-input 输出机器可读报告

位置：

modes/cli.py
modes/data.py

现状：

validate-input 现在主要 print 文本，并可写一个简单 report。 ￼

要做：

新增：

modes validate-input \
  --rna rna.tsv \
  --atac atac.tsv \
  --metadata metadata.tsv \
  --condition condition \
  --external-links links.tsv \
  --out validation.json

输出 JSON：

{
  "ok": false,
  "errors": [],
  "warnings": [],
  "n_samples": 20,
  "n_genes": 1000,
  "n_peaks": 5000,
  "n_links": 10000,
  "n_links_matched": 9500
}

验收标准：

CI 里跑一次：

modes validate-input ... --out validation.json
python -c "import json; json.load(open('validation.json'))"

⸻

8. 增加 fresh wheel install test

位置：

.github/workflows/tests.yml

要做：

新增 job：

build:
  runs-on: ubuntu-latest
  steps:
    - uses: actions/checkout@v4
    - uses: actions/setup-python@v5
      with:
        python-version: "3.11"
    - run: python -m pip install --upgrade pip build twine
    - run: python -m build
    - run: twine check dist/*
    - run: pip install dist/*.whl
    - run: python -c "from modes import MoDES, MoDEData; print('ok')"
    - run: modes --help

验收标准：

打包成功
wheel 可安装
console_scripts 可用

⸻

9. 明确 Python 支持范围

位置：

README.md
pyproject.toml
.github/workflows/tests.yml

现状：

Actions 里尝试 3.10/3.11/3.12，最新失败 run 中 3.10/3.11/3.12 有失败/取消记录。 ￼

要做：

如果 3.12 暂时麻烦，就先官方支持：

Python 3.10 / 3.11

README 写清楚：

Tested on Python 3.10 and 3.11.
Python 3.12 support is planned / experimental.

等 3.12 CI 绿了再加入 official support。

⸻

P2：代码质量与 API 稳定性

10. 统一所有输出 schema

位置：

modes/core.py
docs/output_schema.md
README.md
tests/test_core.py

要做：

冻结 event_table.tsv 字段：

event_id
tf_name
peak_id
gene
context
state
state_confidence
artifact_risk
artifact_reason
event_pval
event_fdr
quality_score
atac_coef
atac_se
atac_pval
atac_fdr
atac_direction
rna_coef
rna_se
rna_pval
rna_fdr
rna_direction
rna_after_atac_coef
rna_after_atac_se
rna_after_atac_pval
rna_after_atac_fdr

加测试：

def test_event_table_schema_exact():
    expected = [...]
    assert list(result.event_table.columns) == expected

验收标准：

输出列顺序稳定，文档和代码完全一致。

⸻

11. 给 MoDESResult.filter() 增加更多组合过滤

位置：

modes/core.py
tests/test_core.py

建议参数：

result.filter(
    state=None,
    states=None,
    min_confidence=None,
    max_event_fdr=None,
    max_atac_fdr=None,
    max_rna_fdr=None,
    max_rna_after_atac_fdr=None,
    exclude_high_artifact=False,
    max_artifact_risk=None,
    min_quality_score=None,
    context=None,
    genes=None,
    peaks=None,
)

要做：

1. state 支持单个字符串；
2. states 支持列表；
3. max_artifact_risk 支持：

low / medium / high

4. genes / peaks 支持 list。

验收标准：

每个过滤条件都有单元测试。

⸻

12. 增加 MoDESResult.save() / load()

位置：

modes/core.py

要做：

新增：

result.save("output/modes_result.pkl")
result = MoDESResult.load("output/modes_result.pkl")

或者更稳：

result.to_dir("output/")
result = MoDESResult.from_dir("output/")

验收标准：

保存再读取后：

pd.testing.assert_frame_equal(result.event_table, loaded.event_table)

⸻

13. 增加运行参数完整记录

位置：

modes/core.py
run_params.tsv/json

要记录：

MoDES version
git commit, if available
Python version
numpy/pandas/statsmodels/anndata versions
condition_col
contrast
fdr_threshold
donor_col
batch_col
covariate_cols
n_samples
n_genes
n_peaks
n_events
n_external_links
n_dropped_links
runtime_seconds

验收标准：

run_params.json 可追溯整次分析。

⸻

P3：数据输入层增强

14. MoDEData.validate() 强化

位置：

modes/data.py
tests/test_data.py

现状：

已有输入验证雏形，但需要产品级强化。

要检查：

RNA/ATAC/metadata index 完全一致
index 是否重复
gene columns 是否重复
peak columns 是否重复
是否有 NaN
是否有 inf
是否有负数
是否全 0 sample
condition 是否存在
condition 是否二分类
donor_col/batch_col/covariate_cols 是否存在
是否有过少 replicate

建议返回结构：

@dataclass
class ValidationReport:
    errors: list[str]
    warnings: list[str]
    summary: dict

而不是只返回字符串 list。

验收标准：

modes validate-input 和 Python API 共用同一个 validator。

⸻

15. 增加重复 feature name 处理策略

位置：

modes/data.py

要做：

遇到重复 gene/peak：

默认报错
可选 aggregate="sum" 合并
可选 make_unique=True 自动重命名

建议默认：

raise ValueError("Duplicate gene names found")

验收标准：

重复 gene/peak 都有测试。

⸻

16. 增加 sparse matrix 支持

位置：

modes/data.py
effects.py
decompose.py

问题：

单细胞和空间矩阵很稀疏。当前 MoDEData 主要用 pd.DataFrame，后续大数据会内存爆。

要做：

第一步不要全改架构，先支持：

MoDEData.from_anndata(...)

内部保存：

AnnData / scipy sparse

或者在 pseudobulk 后再转 DataFrame。

更现实的 1.0 做法：

明确 v1.0 的 MoDES 统计模型使用 pseudobulk dense matrix；
cell-level sparse matrix 只用于聚合，不用于直接 GLM。

如果要支持稀疏聚合：

scipy.sparse.csr_matrix groupby sum

验收标准：

from_pseudobulk() 输入 sparse AnnData 能聚合出 pseudobulk。

⸻

17. MuData 支持

位置：

modes/data.py
docs/singlecell_pseudobulk.md

要做：

新增：

MoDEData.from_mudata(
    mdata,
    rna_mod="rna",
    atac_mod="atac",
    groupby=None,
    condition_col="condition",
    donor_col=None,
    batch_col=None,
)

两种模式：

groupby=None → 直接读取 paired cells，不推荐做 DE
groupby=[...] → 聚合 pseudobulk，推荐

验收标准：

用一个 toy MuData 测试。

⸻

P4：event construction 与 peak-gene links

18. external links schema validator

位置：

modes/events.py
tests/test_events.py
docs/input_formats.md

要做：

新增：

validate_external_links(links, gene_names, peak_names)

检查：

required columns: peak_id, gene
optional: tf_name, source, score, distance
peak_id 是否在 ATAC matrix 中
gene 是否在 RNA matrix 中
score 是否数值
重复 peak-gene 怎么处理

输出：

n_links
n_links_matched
n_links_dropped_peak_missing
n_links_dropped_gene_missing

验收标准：

不匹配链接不会静默进入结果。

⸻

19. event_id 稳定化

位置：

modes/events.py

建议：

现在 event ID 类似：

gene_peak_source

更稳定的做法：

event_id = sha1(f"{peak_id}|{gene}|{tf_name or 'NA'}|{source}").hexdigest()[:12]

同时输出：

event_key = peak_id|gene|tf_name|source

验收标准：

同样输入每次生成同样 event_id。

⸻

20. peak-gene link source 归一化

要做：

如果来自：

SCENIC+
SCARlink
ArchR
Signac
Cicero
user
promoter
distal_250kb

输出统一：

link_source
link_score
link_distance

验收标准：

event_table 里能追踪 link 由谁提供。

⸻

P5：统计模型增强

21. 明确 contrast，不要隐式按 sorted condition

位置：

modes/effects.py
modes/decompose.py
modes/core.py
docs/statistical_model.md

问题：

二分类条件如果按 sorted categories，可能方向不符合用户预期。

要做：

新增参数：

contrast=("treatment", "control")

或者：

reference_condition="control"
target_condition="treatment"

输出记录：

coef = target vs reference

验收标准：

测试：

control/treatment
treatment/control
case/healthy

方向符合用户指定 contrast。

⸻

22. fallback 策略参数化

位置：

modes/effects.py

要做：

新增：

allow_poisson_fallback=True
allow_simplified_fallback=False

默认不允许 simplified fallback 丢 covariates，除非用户明确打开。

验收标准：

如果完整模型失败：

allow_simplified_fallback=False → 报错或返回 failed effect
allow_simplified_fallback=True → 使用简化模型并写 warning

⸻

23. cis_ATAC_score 聚合模型

位置：

modes/decompose.py

当前问题：

现在 conditional decomposition 是：

RNA_g ~ Condition + linked_peak_ATAC + covariates

这对单个 peak-gene event 可以，但一个 gene 通常有多个 cis peaks。

要做：

增加 gene-level cis score：

cis_ATAC_score(g) = weighted sum of all linked peaks for gene g

权重：

external link score
distance decay
uniform

模型：

RNA_g ~ Condition + cis_ATAC_score(g) + covariates

输出：

rna_after_peak_coef
rna_after_cis_atac_coef
cis_atac_score_method

验收标准：

用户可选：

conditional_mode="single_peak"
conditional_mode="cis_score"

⸻

24. artifact risk 增强

位置：

modes/states.py
modes/data.py

当前 artifact risk 主要来自 quality score 和单模态信号。要增强为：

low_quality_score
single_modality_low_quality
low_atac_depth
low_rna_depth
library_size_outlier
batch_associated
weak_link_score
low_group_replicates

输出：

artifact_risk
artifact_reason
depth_score
batch_score
link_score

验收标准：

每个 artifact reason 都有构造测试。

⸻

P6：single-cell pseudobulk 做到稳定

25. from_pseudobulk() 变成正式支持

位置：

modes/data.py
examples/singlecell_pseudobulk/
docs/singlecell_pseudobulk.md

要做：

稳定支持：

groupby=["donor", "condition", "cell_type"]
min_cells_per_group=20
rna_layer=None
atac_layer=None

输出 obs 增加：

group_id
donor
condition
cell_type
n_cells
rna_total_counts
atac_total_counts

验收标准：

toy AnnData / sparse AnnData / group drop / missing layer 都有测试。

⸻

26. 多 cell type 批量运行

要做：

新增：

run_by_context(
    data,
    context_col="cell_type",
    min_samples_per_context=4,
)

输出：

all_contexts_event_table.tsv
context_summary.tsv

验收标准：

多个 cell type 可一键跑，不需要用户写 for-loop。

⸻

P7：spatial 能力边界与增强

27. v1.0 先支持 spatial region-pseudobulk

不要现在直接做 native spatial graph。

新增：

MoDEData.from_spatial_pseudobulk(
    rna_counts,
    atac_counts,
    metadata,
    region_col="region",
    sample_col="sample",
    condition_col="condition",
)

定义：

spatial support in v1.0 = region/sample-level pseudobulk, not native graph.

输出：

context = spatial_region

验收标准：

空间 region 示例跑通。

⸻

28. native spatial graph 放到 v1.2+

未来设计：

SpatialMoDEData(
    rna,
    atac,
    obs,
    coords,
    spatial_graph,
)

新增状态：

spatial_region_specific
spatial_niche_driven
spatial_artifact_edge

但不要现在混进 v1.0。

⸻

P8：benchmarks 完整化

29. 修当前 benchmark 与 CI 失败

位置：

benchmarks/
.github/workflows/tests.yml

现状：

最新 commit 说明说 synthetic benchmark 100% accuracy，但 Actions 失败。 ￼

要做：

1. 先把 benchmark 从 CI 主测试里拆出来：

unit tests: 每次跑
benchmarks: 手动或 nightly

2. benchmark 脚本加 --quick 模式：

python benchmarks/simulated_event_states/run_benchmark.py --quick

3. CI 只跑 quick benchmark。

验收标准：

unit tests 绿色
quick benchmark 绿色
full benchmark 可本地/手动运行

⸻

30. benchmark 产物标准化

每个 benchmark 输出：

truth.tsv
predicted_event_table.tsv
metrics.tsv
confusion_matrix.tsv
confusion_matrix.png
runtime.tsv

每个 benchmark README 写：

what it tests
expected behavior
how to run
how to interpret

⸻

31. negative control 强化

要做：

负控包括：

shuffle condition labels
shuffle peak-gene links
shuffle ATAC sample labels
random external links

预期：

majority null
low false concordant rate
event_fdr roughly controlled

⸻

32. baseline comparison 真实可解释

要做：

baseline 至少有：

naive_overlap:
  ATAC sig + RNA sig = concordant
  ATAC only = primed
  RNA only = rna_only
correlation_baseline:
  peak-gene correlation + DA/DE
MoDES:
  full method

输出：

MoDES vs naive overlap macro-F1
per-state precision/recall
artifact risk specificity

⸻

P9：真实数据 demo

33. PBMC multiome smoke test

位置：

notebooks/
examples/
docs/benchmark.md

要做：

添加：

notebooks/01_pbmc_multiome_smoke_test.ipynb

目的：

真实 10x multiome 数据可以跑完整 pipeline

不要声称发现疾病机制。

输出：

event_table head
state distribution
artifact risk distribution
runtime
memory

⸻

34. biological demo

选择一个有真实状态变化的数据：

differentiation
stimulation
treatment
time course

目标：

chromatin_primed events 富集 lineage TF
concordant events 富集成熟 marker
rna_only events 富集 stress/trans response

输出：

demo_event_table.tsv
motif_enrichment.tsv
marker_enrichment.tsv

⸻

P10：文档与发布

35. 文档表格格式再整理

README 当前在 GitHub 页面能读，但仍有一些表格在 web 解析中显示为纯文本，如 “State Pattern Biological Interpretation” 和 “Capability Status”。建议全部改成标准 Markdown table。 ￼

要做：

全部检查：

README.md
docs/*.md
benchmarks/README.md

把：

Metric Description
Accuracy Fraction...

改成：

| Metric | Description |
|---|---|
| Accuracy | Fraction of events correctly classified |

⸻

36. CITATION.cff 改成真实作者信息

位置：

CITATION.cff

现状：

作者现在是：

MoDES contributors

比较粗糙。 ￼

要做：

改成：

authors:
  - family-names: Ling
    given-names: Rongsong

如果你不想写全名，可以保持，但 1.0 建议正式化。

⸻

37. README badge

加：

![tests](https://github.com/Lings01/MoDES/actions/workflows/tests.yml/badge.svg)

但只有在 CI 绿色后再加，否则 badge 会红。

⸻

38. Release checklist

新增：

docs/release_checklist.md

内容：

[ ] CI green
[ ] unit tests pass
[ ] lint pass
[ ] quick benchmarks pass
[ ] example runs
[ ] README updated
[ ] CHANGELOG updated
[ ] version bumped
[ ] tag pushed

⸻

P11：1.0 后功能路线

39. protein layer，v1.1

新增：

MoDES-RAP = RNA + ATAC + Protein

新增状态：

full_activation
protein_buffered
protein_memory
protein_opposite

新增输入：

protein_counts
protein_gene_links

新增输出：

protein_coef
protein_fdr
protein_after_rna_coef

⸻

40. native spatial graph，v1.2

新增：

SpatialMoDEData
coords
spatial_graph
region labels
neighbor effect
spatial autocorrelation

新增状态：

spatial_region_specific
spatial_niche_driven
spatial_artifact_edge

⸻

41. multi-condition / pseudotime，v1.3

新增：

contrast matrix
multi-class condition
time-course
pseudotime lag
ATAC→RNA delay

⸻

最推荐的执行顺序

你现在不要再开新功能。先按下面顺序做：

1. 修当前 CI failure
2. 修 lint / pytest / lower-bound
3. 确认 requirements/dev-requirements 真实多行可安装
4. 把 CLI / CHANGELOG / ROADMAP 的矛盾修掉
5. 加 pyproject.toml 和 build job
6. 冻结 output schema，并加 schema exact test
7. 强化 input validation
8. 强化 from_pseudobulk
9. benchmark quick/full 分离
