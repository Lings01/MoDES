Lings，我按英文关键词做了文献和代码库检索。结论是：这个方向有新意，但新意不能写成“多组学整合”本身；真正的新意应该写成“跨组学不一致模式的调控事件状态分解”。

也就是说，不要提出一个“新的 integration/embedding 方法”。这个方向已经很拥挤。你要提出的是：

一个把 RNA、ATAC、蛋白、甲基化、代谢物、空间位置等多层信息统一到“生物事件”层面，并系统判定该事件处于 chromatin priming、transcriptional activation、protein buffering、epigenetic memory、spatial niche effect 或 technical artifact 的方法。

我建议暂名：

MoDES：Multi-Omics Discordance/Event State inference

中文名：

多组学不一致驱动的调控事件状态推断框架

下面是完整技术报告。

⸻

1. 执行摘要

目前多组学领域已经有很多方法能做 joint embedding、batch correction、cell clustering、feature selection、cross-modality imputation、spatial registration 等任务。2025 年 Nature Methods 的 benchmark 总结了 40 个单细胞多模态整合方法，覆盖 vertical、diagonal、mosaic、cross integration 四类输入结构，以及 dimension reduction、clustering、batch correction、classification、feature selection、imputation、spatial registration 七类任务。这个结果说明：再做一个普通的“多组学整合算法”很难有明显新意。  ￼

但是，我没有检索到一个成熟代码库同时做到下面四件事：第一，把多组学结果的“不一致”作为主要信号；第二，把不一致分解成可解释的调控状态；第三，同时适配 bulk、单细胞和空间转录组/空间多组学；第四，以 regulatory event 而不是 gene、peak、cluster、latent factor 作为核心分析单位。现有最接近的工具分别覆盖部分问题，例如 SCENIC+ 和 SCARlink 做 RNA+ATAC 调控链接，MOFA/MEFISTO 做多组学潜变量，DRAGON 做 paired omics 图模型，MISTy/SIMO/MultiGATE 做空间或多组学整合，但它们不是以“多组学不一致状态分类”为核心输出。 ￼

所以这个方向有新意，但要避免写成：

我们提出一种新的多组学整合方法。

应该写成：

我们提出一种 discordance-aware regulatory event model，
不是把多组学压缩成一个 embedding，
而是利用不同组学层之间的一致与不一致，
推断调控事件发生在哪一层、是否有时间延迟、是否有空间 niche 依赖、是否可能是技术伪影。

⸻

2. 检索范围与关键词

我主要按这些英文关键词检索：

single-cell multiome RNA ATAC protein discordance
multi-omics discordance regulatory inference
RNA ATAC protein regulatory event inference
single-cell multiome enhancer gene protein integration
bulk multi-omics transcriptome proteome discordance
spatial multi-omics integration regulatory inference
multi-omics evidence graph
multi-omics differential correlation
RNA protein discordance CITE-seq
SCENIC+ SCARlink GLUE MOFA MEFISTO MISTy SIMO MultiGATE code

检索重点不是列全所有多组学软件，而是判断你这个想法在现有方法中的空位。

⸻

3. 现有方法地图

3.1 普通多组学整合：已经很拥挤

现有主流目标通常是：

joint embedding
dimension reduction
batch correction
clustering
classification
imputation
feature selection
spatial registration

2025 年的单细胞多模态 benchmark 明确把当前方法划成四类 integration 结构，并评估 40 个方法在 64 个真实数据集和 22 个模拟数据集上的表现。这说明现有方法生态已经非常复杂，而且大多数方法的核心输出仍是 embedding、cell label、imputed modality 或 selected feature。 ￼

典型例子包括 Seurat WNN。WNN 的目标是学习每个细胞中不同模态的相对信息量，并基于多模态数据定义 cellular state；这是非常实用的单细胞 multimodal integration 方法，但它的主要输出仍是 cell-level neighbor graph / embedding / clustering，而不是解释一个具体调控事件为何在 ATAC、RNA、protein 层不一致。 ￼

MOFA/MOFA+ 是经典多组学因子分析框架，能从多个 omics matrix 中学习 latent factors，并处理 missing values、shared/omics-specific variation；MEFISTO 在 MOFA 基础上加入时间或空间连续结构，用于建模随时间/空间平滑变化的 factors。它们很适合发现全局 variation source，但不直接把一个 enhancer→gene→protein 事件分类为 primed、buffered、memory 或 artifact。 ￼

DIABLO/mixOmics 是 bulk 多组学里常用的 supervised integration 思路，目标是最大化多个数据集之间的 common/correlated information，并选择能区分 phenotype group 的特征。它非常适合 biomarker panel discovery，但核心思想仍偏“共变特征选择”，不是“调控链条层级分解”。 ￼

⸻

3.2 RNA+ATAC 调控推断：SCENIC+、SCARlink 很强，但仍有空位

SCENIC+ 是很接近你目标的现有方法。它结合 single-cell chromatin accessibility、gene expression 和 motif discovery，推断 enhancer-driven gene regulatory networks，也就是 eGRN。SCENIC+ 的三步流程包括 candidate enhancer identification、TF motif enrichment、TF-enhancer-target gene linking，输出 eRegulons。 ￼

SCARlink 也非常接近。它用 single-cell ATAC+RNA multiome 数据，从 gene 周围 ±250 kb 的 500 bp tiles 的 chromatin accessibility 预测 gene expression，使用 regularized Poisson regression，并用 learned coefficients / Shapley values 推断 enhancer-gene links。SCARlink 的优点是绕开简单 pairwise peak-gene correlation，也不强依赖 peak calling。 ￼

但是，这两个方法的主要问题是：它们更偏 RNA+ATAC 的 regulatory link / eGRN inference。它们不以“多组学不一致状态”为主输出，也不天然覆盖 RNA→protein buffering、protein memory、bulk proteogenomics、空间 niche effect、technical discordance adjudication 这些问题。你的方法不能和 SCENIC+/SCARlink 正面抢“谁更会做 peak-gene link”，而应该把它们作为 candidate event generator 或 baseline。

⸻

3.3 多组学图模型：DRAGON 有启发，但不是这个问题

DRAGON 是一个很有参考价值的 paired omics network inference 方法。它用 Gaussian Graphical Model 估计多组学 partial correlation network，并且显式考虑不同 omics layer 的变量数、噪音和边密度差异；作者用 TCGA breast cancer 的 transcriptome–methylome 数据展示了 promoter methylation 与 gene expression 的 regulatory relationship。 ￼

但 DRAGON 的核心是 conditional dependence network。它能告诉你两个 feature 是否在控制其他变量后仍相关，却不会直接回答：

这个事件是 chromatin priming？
是 RNA-only regulation？
是 protein buffering？
是空间 niche 驱动？
还是技术伪影？

所以 DRAGON 可作为 MoDES 的网络 baseline，而不是替代品。

⸻

3.4 空间多组学：很多方法在做空间映射或 embedding，但机制分解仍不足

MISTy 是空间分析中很重要的 explainable framework。它可以处理 highly multiplexed spatial data，通过不同 spatial views 建模 intracellular、juxtacrine、paracrine 或 broader tissue view 的 marker relationships。这个思想对 MoDES 的空间扩展非常有价值，因为我们也需要把 event effect 分成 cell-intrinsic、local niche 和 tissue-scale 三层。 ￼

SIMO 解决的是 spatial transcriptomics 与多个 non-spatial single-cell omics 数据的空间映射和整合，支持 RNA、ATAC、DNA methylation 等数据，并用于构建 spatial patterns 和 spatial GRN。它的价值是 spatial mapping，而不是专门解释跨组学不一致。 ￼

MultiGATE 是 2025 年 Nature Communications 的 spatial multi-omics integration/regulatory inference 方法，使用 two-level graph attention auto-encoder 整合 spatial multi-omics，并推断 cross-modality regulatory relationships。它代表了空间多组学深度学习方向，但也进一步说明：如果你不想卷 AI，就应该避开“更好的 embedding/autoencoder”，转向可解释统计建模和事件状态推断。 ￼

空间多组学综述也明确指出，跨层解释必须考虑不同 molecular layer 的调控关系：chromatin accessibility 允许 TF binding，mRNA abundance 影响但不决定 protein level，protein 还受 translation、stability、technical artifacts 等影响；偏离预期相关性的模式既可能是生物调控，也可能是抗体交叉反应等技术问题。这个观点正好支持 MoDES 的问题设定。 ￼

⸻

3.5 “不一致”相关工具：已有，但还不够机制化

OmicsTIDE 是一个直接提到 concordance/discordance 的工具。它是 web-based application，用于研究 omics datasets 之间 regulatory trends 的一致和不一致；GitHub 说明它用于分析 multi-omics data sets 中 regulatory trends 的 concordance 和 discordance。 ￼

discordant R 包用于检测 molecular feature pairs 的 differential correlation，基于 mixture models 判断不同条件下 pairwise association 是否改变。这个包关注的是 differential correlation，而不是 enhancer→RNA→protein 的调控状态分类。 ￼

所以，不一致分析不是完全没人做；但现有工具大多停留在 trend comparison、pairwise correlation、visualization 或 differential correlation。你的新意应该是：把不一致模式升格为 regulatory event state，并且给出统计检验、状态概率、空间/细胞类型/条件分解。

⸻

4. 新意判断

4.1 不新颖的版本

下面这些说法不够新：

我们提出一个多组学整合方法。
我们把 RNA 和 ATAC 结合起来找调控关系。
我们把 RNA 和 protein 做相关性分析。
我们找 concordant / discordant genes。
我们做一个多组学网络。
我们做一个 bulk/sc/spatial 都能用的统一可视化工具。

这些方向都已有很多方法或相近工具。

⸻

4.2 有新意的版本

真正有新意的是这个：

我们提出一个 regulatory-event-centric framework，
把多组学数据从 feature-level / cell-level / embedding-level
转换到 event-level evidence space。
每个 event 不再只有“显著/不显著”，
而是被判定为一种可解释状态：
chromatin-primed
cis-driven concordant activation
RNA-only / trans-driven
protein-buffered
protein-memory
epigenetic-memory
spatial-niche-driven
artifact-like

这和现有方法的区别是：

维度	现有主流方法	MoDES 的目标
分析单位	cell、gene、peak、factor、embedding	regulatory event
核心目标	integration / clustering / imputation / GRN	事件状态分解
对不一致的态度	噪音、低相关、可视化结果	生物信号或技术伪影
输出	latent space、cluster、network edge	state probability + layer-specific effect
bulk/sc/spatial	通常分别处理	同一 evidence abstraction 下处理
解释	后验解释	模型内解释

这就是可以写成方法学论文的核心。

⸻

5. 方法定义

5.1 核心对象：regulatory event，而不是 gene

一个 event 定义为：

e = (r, c, g, p, h, q)

其中：

符号	含义
r	regulator，例如 TF、miRNA、kinase、pathway regulator
c	cis-regulatory element，例如 enhancer、promoter、methylated region、ATAC peak
g	target gene
p	protein 或 ADT marker
h	pathway / metabolite / phosphosite，可选
q	context，例如 cell type、condition、spatial niche、time point

例如：

STAT1 motif / enhancer accessibility ↑
→ IFIT3 RNA ↑
→ IFIT3 protein 或 IFN pathway marker ↑
→ 只发生在 disease monocytes 的 inflammatory niche

这才是一个多组学事件。

⸻

5.2 输入数据

MoDES 应该支持三类数据结构。

A. Bulk multi-omics

单位是 sample：

sample × RNA
sample × ATAC / ChIP / methylation
sample × proteomics / phosphoproteomics
sample × metabolomics
sample metadata: condition, batch, donor, tissue, purity, cell composition

适用场景：

TCGA/CPTAC-like cohort
bulk RNA + ATAC
bulk RNA + proteomics
bulk RNA + methylation + protein
time-course perturbation bulk multi-omics

⸻

B. Single-cell / single-nucleus multiome

单位可以是 cell，也可以是 pseudobulk：

cell × RNA
cell × ATAC
cell × ADT protein, optional
cell × TCR/BCR, optional
cell metadata: donor, batch, condition, cell type, pseudotime

对于 condition effect，我建议默认走 pseudobulk：

donor × condition × cell type

原因是 differential analysis 在单细胞里如果直接把 cell 当独立 replicate，很容易出现 donor-level pseudoreplication。2024 年 Nature Communications 的 scATAC differential accessibility benchmark 也显示，在单细胞多组学中，把细胞聚合成 pseudobulk 的 DA 方法在 RNA/ATAC concordance 评估中整体优于非 pseudobulk 方法。 ￼

⸻

C. Spatial transcriptomics / spatial multi-omics

单位可以是 spot、bin、segmented cell 或 pixel：

spot/cell × RNA
spot/cell × protein, optional
spot/cell × ATAC / histone mark / methylation, optional
spatial coordinates
neighborhood graph
image-derived region labels, optional
cell composition, optional

适用场景：

Visium RNA + protein
spatial ATAC-RNA
spatial RNA + non-spatial ATAC reference
Xenium/CosMx/MERFISH + protein imaging
spatial transcriptomics + inferred cell composition

⸻

6. 统一统计框架

6.1 观测矩阵

设分析单位为 u。

在 bulk 中：

u = \text{sample}

在单细胞中：

u = \text{cell}

或：

u = \text{donor} \times \text{condition} \times \text{cell type}

在空间中：

u = \text{spot / cell / spatial bin}

每个模态为：

m \in \{RNA, ATAC, Protein, Methylation, Metabolite, Phospho\}

观测为：

Y^{(m)}_{u,f}

其中 f 是 gene、peak、protein、CpG region、metabolite 等 feature。

⸻

6.2 每个事件的 evidence vector

对每个 event e，构造多组学证据向量：

D_e =
[
\delta^{ATAC}_e,
\delta^{RNA}_e,
\delta^{RNA|ATAC}_e,
\delta^{Protein}_e,
\delta^{Protein|RNA,ATAC}_e,
\delta^{Methylation}_e,
\delta^{Spatial}_e,
\tau_e
]

其中：

量	含义
\delta^{ATAC}_e	condition / pseudotime / spatial 对 enhancer accessibility 的 effect
\delta^{RNA}_e	对 target gene RNA 的 effect
\delta^{RNA|ATAC}_e	控制 cis-ATAC 后 RNA 仍剩余的 effect
\delta^{Protein}_e	对 protein abundance 的 effect
\delta^{Protein|RNA,ATAC}_e	控制 RNA/ATAC 后 protein 层仍剩余的 effect
\delta^{Methylation}_e	methylation 层 effect
\delta^{Spatial}_e	spatial niche 或空间自相关 effect
\tau_e	ATAC→RNA 或 RNA→protein 的时间/伪时间延迟

这个 evidence vector 是 MoDES 的核心。
所有 bulk、single-cell、spatial 数据最终都映射到这个 event-level evidence space。

⸻

7. 分层模型

7.1 ATAC 层

对 enhancer / promoter / chromatin region：

Y^{ATAC}_{u,c}
\sim
\text{NB}(\mu^{ATAC}_{u,c}, \theta_c)

或对 binary accessibility：

Y^{ATAC}_{u,c}
\sim
\text{Binomial}(n_u, q_{u,c})

模型：

g_A(\mu^{ATAC}_{u,c})
=
\alpha^A_e
+
\delta^A_e C_u
+
X_u \beta^A
+
b^A_{donor}
+
b^A_{batch}

其中 C_u 是 condition、pseudotime、spatial region 或 perturbation。

⸻

7.2 RNA 层：先看总效应

Y^{RNA}_{u,g}
\sim
\text{NB}(\mu^{RNA}_{u,g}, \theta_g)

g_R(\mu^{RNA}_{u,g})
=
\alpha^R_e
+
\delta^R_e C_u
+
X_u \beta^R
+
b^R_{donor}
+
b^R_{batch}

这里 \delta^R_e 是 RNA 层总 effect。

⸻

7.3 RNA 层：再控制 ATAC

g_R(\mu^{RNA}_{u,g})
=
\alpha^{R|A}_e
+
\eta_e A_{u,c}
+
\rho^{R|A}_e C_u
+
X_u \beta^{R|A}
+
b_{donor}
+
b_{batch}

解释：

* 如果 \delta^R_e 显著，但 \rho^{R|A}_e 不显著，说明 RNA 变化可以被 local chromatin accessibility 解释。
* 如果 \rho^{R|A}_e 仍显著，说明 RNA 层有 ATAC 不能解释的剩余效应，可能是 trans-regulation、RNA stability、batch、stress response 或未测量调控层。

⸻

7.4 Protein 层

Y^{Protein}_{u,p}
\sim
\text{NB}(\mu^{Protein}_{u,p}, \phi_p)

或 bulk proteomics 中：

Y^{Protein}_{u,p}
\sim
N(\mu^{Protein}_{u,p}, \sigma_p^2)

条件模型：

g_P(\mu^{Protein}_{u,p})
=
\alpha^P_e
+
\kappa_R R_{u,g}
+
\kappa_A A_{u,c}
+
\rho^{P|R,A}_e C_u
+
X_u \beta^P
+
b_{donor}
+
b_{batch}

解释：

* RNA 变、protein 不变：可能是 protein buffering。
* RNA 不变、protein 变：可能是 protein memory、protein stability、surface marker retention、technical antibody issue。
* 控制 RNA 和 ATAC 后 protein 仍有 condition effect：说明 protein 层有独立调控。

RNA-protein 关系本来就不是简单线性。2025 年关于 gene expression 与 protein abundance 关系的综述指出，mRNA-protein 相关性会受统计方法、实验设计、translation rate、protein degradation、post-transcriptional modification 等因素影响，并且单细胞数据正在成为研究这个问题的新来源。 ￼

⸻

7.5 Methylation / epigenetic memory 层

对 methylation：

Y^{Meth}_{u,c}
\sim
\text{BetaBinomial}(n_{u,c}, q_{u,c}, \phi_c)

模型：

\text{logit}(q_{u,c})
=
\alpha^M_e
+
\delta^M_e C_u
+
X_u\beta^M
+
b_{donor}
+
b_{batch}

典型解释：

methylation / ATAC 改变
但 RNA 和 protein 暂时没变
→ epigenetic memory / priming

这在 bulk methylome + RNA 或 single-cell multiome 中都能用。

⸻

8. 事件状态分类

MoDES 的最终输出不是“相关/不相关”，而是：

P(Z_e = z | D_e)

其中：

z \in
\{
null,
concordant,
chromatin\ primed,
RNA\ only,
protein\ buffered,
protein\ memory,
epigenetic\ memory,
spatial\ niche,
artifact
\}

⸻

8.1 状态 1：cis-driven concordant activation

模式：

ATAC ↑
RNA ↑
Protein ↑

解释：

局部染色质开放、转录上升、蛋白执行上升。
这是完整调控链条激活。

数学判据：

sign(\delta^A_e)
=
sign(\delta^R_e)
=
sign(\delta^P_e)

且：

|\rho^{R|A}_e| < |\delta^R_e|

说明 RNA effect 可部分被 ATAC 解释。

⸻

8.2 状态 2：chromatin primed

模式：

ATAC ↑
RNA unchanged
Protein unchanged

解释：

染色质层已经准备好，但转录/蛋白层尚未响应。

如果有 time-course 或 pseudotime，则进一步要求：

ATAC(t) \rightarrow RNA(t + \tau)

且：

\tau > 0

SHARE-seq 提出的 chromatin potential 思想已经证明，chromatin accessibility 可以作为 lineage priming 和未来 cell fate 的信号；这可以作为 MoDES 中 primed state 的理论基础之一。 ￼

⸻

8.3 状态 3：RNA-only / trans-driven event

模式：

ATAC unchanged
RNA ↑
Protein weak / unchanged

解释：

RNA 变化不能由局部 chromatin accessibility 解释，
可能来自 trans TF activity、RNA stability、stress response、
enhancer 未覆盖、peak calling 问题，或其他调控层。

数学判据：

\delta^R_e \neq 0

\delta^A_e \approx 0

且：

\rho^{R|A}_e \approx \delta^R_e

⸻

8.4 状态 4：protein buffered

模式：

ATAC ↑ 或 RNA ↑
Protein unchanged / ↓

解释：

转录变化没有传递到蛋白层。
可能是 translation control、protein degradation、protein half-life、post-transcriptional buffering。

数学判据：

\delta^R_e \neq 0

但：

\delta^P_e \approx 0

或：

sign(\delta^P_e) \neq sign(\delta^R_e)

⸻

8.5 状态 5：protein memory

模式：

RNA unchanged
Protein ↑

解释：

蛋白保留了过去激活状态；
RNA 已经回落，但 protein 或 surface marker 仍存在。

这对 CITE-seq、TEA-seq、spatial proteomics 特别重要。CITE-seq 最初就是为了弥补 scRNA-seq 不能直接测量 cell-surface protein levels 的问题，而 TEA-seq 进一步把 transcriptome、epitope 和 chromatin accessibility 同时测量到同一批单细胞中。 ￼

⸻

8.6 状态 6：epigenetic memory

模式：

ATAC / methylation changed
RNA unchanged
Protein unchanged

解释：

当前转录状态不明显，但表观层保留历史状态。

适用：

trained immunity
chronic inflammation
tumor resistance
developmental memory
treatment-exposed residual cells

⸻

8.7 状态 7：spatial-niche-driven event

模式：

事件在某个空间邻域显著
但在全局 cell type / sample 层不显著

解释：

这个 event 不是单纯 cell-intrinsic，
而是由 spatial niche、邻近细胞、组织区域、氧/营养梯度或局部炎症信号驱动。

空间模型中加入图正则：

\lambda
\sum_{(u,v)\in E}
w_{uv}
(\phi_u - \phi_v)^2

或 CAR prior：

\phi_u | \phi_{-u}
\sim
N
\left(
\frac{\sum_v w_{uv}\phi_v}{\sum_v w_{uv}},
\frac{\sigma^2}{\sum_v w_{uv}}
\right)

⸻

8.8 状态 8：artifact-like discordance

模式：

只有一个模态异常
且该异常与 QC、depth、ambient、batch、antibody background、spatial edge effect 强相关

解释：

可能不是生物学，而是技术伪影。

这很重要，因为多组学不一致不一定都是“高级生物学”。空间多组学综述也强调，mRNA-protein discordance 既可能来自 post-transcriptional regulation，也可能来自 technical artifacts，例如 antibody cross-reactivity，因此需要 uncertainty quantification 和机制知识来判断。 ￼

⸻

9. MoDES 如何同时适配 bulk、单细胞、空间

9.1 Bulk 版本

输入单位

u = sample

推荐模型

limma / voom
negative binomial GLM
linear mixed model
Gaussian graphical model
Bayesian hierarchical model

特点

Bulk 版本最适合：

RNA + proteomics
RNA + methylation
RNA + ATAC
RNA + phosphoproteomics
RNA + metabolomics
time-course perturbation
large clinical cohort

必须处理的问题

Bulk tissue 最大问题是 cell composition confounding。因此模型中要加入：

CellComposition_u

例如：

Y^{RNA}_{u,g}
\sim
Condition_u
+
CellComposition_u
+
Batch_u
+
Purity_u
+
Sex_u
+
Age_u

如果有空间或单细胞 reference，可以先估计 bulk cell composition，然后在 MoDES 中作为 covariate。

Bulk 的核心输出

event_id
condition_effect_ATAC
condition_effect_RNA
condition_effect_protein
RNA_residual_after_ATAC
protein_residual_after_RNA
state_probability
state_label

Bulk 的优势

Bulk multi-omics 的 replicate 通常更明确，适合临床 cohort 和 time-course。缺点是无法直接区分 cell type，所以最好和 deconvolution 或 single-cell reference 配合。

⸻

9.2 单细胞版本

输入单位

两层：

cell-level for state heterogeneity / pseudotime / trajectory
pseudobulk-level for condition effect

也就是：

u_1 = cell

u_2 = donor \times condition \times celltype

推荐策略

单细胞中不要直接对每个 cell 做 condition test，然后把 cell 当独立样本。推荐：

1. 先做 cell type / state annotation
2. 每个 donor × condition × celltype 聚合 pseudobulk
3. 在 pseudobulk 上估计 ATAC/RNA/protein effect
4. 在 cell-level 上估计连续轨迹、pseudotime、event activity distribution

单细胞的核心优势

单细胞可以回答 bulk 不能回答的问题：

事件发生在哪个 cell type？
是否只发生在 rare state？
ATAC 是否先于 RNA？
RNA-protein discordance 是否只存在于某个 activation state？
同一 clone 是否有不同 multi-omics state？

单细胞输出

event × celltype × condition state table
event activity per cell
event pseudotime lag
cell-level discordance score
donor-replicated FDR

⸻

9.3 空间版本

输入单位

u = spot / cell / spatial\ bin

空间图

构造：

G = (V, E)

其中节点是 spot/cell，边来自：

physical adjacency
kNN in spatial coordinates
histology-based neighborhood
cell-cell contact graph
region graph

空间效应分解

把 event effect 分成三部分：

Effect_e
=
Intrinsic_e
+
LocalNiche_e
+
TissueRegion_e
+
Residual_e

对应：

cell-intrinsic regulatory event
near-neighbor niche effect
large-scale anatomical region effect
unexplained / artifact

空间模型

Y_{u,e}
=
\alpha_e
+
\beta^{intrinsic}_e X_{u,e}
+
\beta^{niche}_e N_{u,e}
+
\phi_{region(u),e}
+
\epsilon_{u,e}

其中：

N_{u,e}
=
\sum_{v \in Neighbor(u)}
w_{uv} X_{v,e}

空间版本的核心输出

spatially localized primed event
spatially localized protein-memory event
niche-driven enhancer activation
region-specific RNA-only event
artifact-like edge / low-quality region

⸻

10. 代码库检索结果与可复用模块

下面是和 MoDES 关系最密切的代码库。

代码库	已有功能	和 MoDES 的关系
SCENIC+	Python 包，用 scRNA-seq + scATAC-seq 构建 enhancer-driven GRN	可作为 event candidate generator / baseline，不建议正面替代它  ￼
SCARlink	用 tile-based ATAC 和 regularized Poisson regression 预测 gene expression，推断 regulatory regions	可借鉴 RNA-from-ATAC conditional model  ￼
MOFA2 / MEFISTO	多组学 latent factor；MEFISTO 支持 temporal/spatial covariates	可作为 latent factor baseline，而不是事件状态模型  ￼
GLUE	Graph-linked unified embedding，用于 single-cell multi-omics integration 和 regulatory inference	可作为 unpaired/mosaic integration baseline  ￼
mixOmics / DIABLO	R 包，supervised multi-omics integration 和 biomarker selection	bulk baseline，特别是 phenotype discrimination  ￼
MISTy	多视角空间建模，分析 intra/intercellular relationships	空间扩展可借鉴 view decomposition  ￼
OmicsTIDE2.0	研究 multi-omics regulatory trend 的 concordance / discordance	最接近“不一致”概念的可视化工具，但不是事件状态推断  ￼
discordant	R 包，用 mixture models 检测 omics feature pairs 的 differential correlation	可作为 pairwise discordance baseline  ￼
scMultiBench	2025 benchmark 的代码，覆盖多种 multimodal integration tasks	可用来设计 benchmark 和 baseline runner  ￼
scMultiSim	可模拟 gene expression、chromatin accessibility、RNA velocity、spatial locations，并包含 GRN/CCI/technical noise	可用于制造 ground-truth event-state benchmark  ￼

⸻

11. 与现有方法的差异矩阵

方法	RNA+ATAC	RNA+protein	spatial	bulk	event-state	discordance as signal	technical artifact state
WNN/Seurat	可用	可用	间接	否	否	否	否
MOFA/MOFA+	可用	可用	MEFISTO 可用	可用	否	间接	否
DIABLO	可用	可用	否	强	否	间接	否
SCENIC+	强	否	间接	否	部分 eGRN	部分	否
SCARlink	强	否	否	否	部分 peak-gene	部分	否
GLUE	强	部分	部分	否	否	间接	否
DRAGON	可用	可用	否	强	否	pairwise/network	否
MISTy	间接	可用	强	否	否	间接	否
OmicsTIDE	可用	可用	否	强	否	是，但偏趋势	否
MoDES	强	强	强	强	核心输出	核心思想	显式建模

这张表就是 proposal 的核心逻辑。

⸻

12. 最小可行版本设计

我建议第一版不要吃太多数据类型。最现实的 MVP 是：

MoDES-RA:
RNA + ATAC regulatory event state inference

先支持：

bulk RNA + ATAC
single-cell multiome RNA + ATAC
spatial RNA + ATAC / spatial RNA + inferred ATAC reference

第二版再加入 protein：

MoDES-RAP:
RNA + ATAC + Protein

⸻

12.1 MVP 输入

RNA count matrix
ATAC peak matrix
sample/cell/spot metadata
condition labels
donor labels
batch labels
cell type labels, optional
spatial coordinates, optional
genome annotation
motif annotation
peak-to-gene candidate links, optional

支持格式：

AnnData
MuData
Seurat object exported h5ad
pseudobulk count tables
plain TSV matrices

⸻

12.2 MVP 输出

event_table.tsv
event_state_probability.tsv
event_layer_effects.tsv
event_network.graphml
event_spatial_map.h5ad
report.html

核心输出表：

event_id	TF	peak	gene	context	ATAC_effect	RNA_effect	RNA_after_ATAC	state	posterior
E001	STAT1	chr1:…	IFIT3	disease monocyte	+2.1	+1.7	+0.2	concordant	0.91
E002	RUNX3	chr7:…	GZMB	CD8 T	+1.8	+0.1	+0.0	chromatin_primed	0.84
E003	CEBPB	chr8:…	IL1B	macrophage niche	+0.2	+2.3	+2.1	RNA_only	0.79
E004	FOXP3	chr10:…	IL2RA	Treg	+0.4	-0.1	+0.0	protein_memory	0.76

⸻

13. 算法流程

Step 1：生成 candidate events

对每个 gene g，找候选 regulatory elements：

promoter peaks
gene body peaks
±250 kb distal peaks
co-accessible peaks
motif-supported peaks
known enhancer annotations
methylation regions near promoter/enhancer
protein encoded by gene, if available

可以复用 SCENIC+、SCARlink、ArchR peak-to-gene links 或 Signac links 作为候选，不需要第一版自己从零构造所有 enhancer-gene links。SCARlink 已经证明，用 gene 周围 chromatin accessibility 预测 gene expression 是可行路线；SCENIC+ 已经证明 motif + accessibility + expression 可用于 enhancer-driven GRN inference。 ￼

⸻

Step 2：每个模态估计 effect size

对每个 event e，估计：

ATAC effect
RNA effect
protein effect
methylation effect
spatial effect

每个 effect 都需要：

estimate
standard error
p value
FDR
direction
quality score

⸻

Step 3：条件分解

核心是比较以下模型。

RNA-only model

RNA_g \sim Condition + Covariates

RNA conditioned on ATAC

RNA_g \sim Condition + ATAC_{cis(g)} + Covariates

Protein conditioned on RNA and ATAC

Protein_p \sim Condition + RNA_g + ATAC_{cis(g)} + Covariates

通过 condition coefficient 的衰减判断 effect 发生在哪一层。

⸻

Step 4：构造 evidence vector

D_e =
[
z_A,
z_R,
z_{R|A},
z_P,
z_{P|R,A},
z_M,
z_S,
q_e
]

其中：

z_m = \frac{\hat{\delta}_m}{SE(\hat{\delta}_m)}

q_e 是质量分数，例如：

read depth
dropout
ambient score
ADT background
ATAC fragment depth
spatial edge score
batch association
donor reproducibility

⸻

Step 5：状态判别

第一版可以用 rule-based + empirical Bayes。

例如：

if ATAC_sig and RNA_sig and same_direction:
    concordant
if ATAC_sig and not RNA_sig:
    chromatin_primed
if RNA_sig and not ATAC_sig:
    RNA_only
if RNA_sig and not Protein_sig:
    protein_buffered
if Protein_sig and not RNA_sig:
    protein_memory
if spatial_sig and not global_sig:
    spatial_niche
if modality_sig but QC_associated:
    artifact_like

第二版再改成 probabilistic finite mixture：

D_e | Z_e = k
\sim
N(\mu_k, \Sigma_k)

P(Z_e = k | D_e)
=
\frac{
\pi_k N(D_e;\mu_k,\Sigma_k)
}{
\sum_l \pi_l N(D_e;\mu_l,\Sigma_l)
}

输出：

most_likely_state
posterior_probability
local_FDR = 1 - posterior_probability

⸻

14. 空间扩展

空间版本的关键不是简单把 spot 加进去，而是判断 event 的驱动来源。

14.1 三层空间分解

对每个 event：

intrinsic view:
  当前 spot/cell 自己的 RNA/ATAC/protein
local niche view:
  邻近 spot/cell 的 ligand、cell type composition、event activity
tissue-scale view:
  anatomical region、distance to boundary、tumor core/invasive edge、vascular niche

这和 MISTy 的 multi-view 思想兼容，但 MoDES 的目标不是预测 marker，而是给每个 regulatory event 判定状态。MISTy 已经展示了 spatial views 可以用于理解 marker interactions 和 intra-/intercellular relationships。 ￼

⸻

14.2 空间 event 模型

Y_{u,e}
=
\alpha_e
+
\beta^{I}_e I_{u,e}
+
\beta^{N}_e N_{u,e}
+
\beta^{R}_e Region_u
+
\phi_{u,e}
+
\epsilon_{u,e}

其中：

* I_{u,e}：cell-intrinsic event evidence；
* N_{u,e}：neighbor event evidence；
* Region_u：组织区域；
* \phi_{u,e}：空间平滑 random effect。

如果：

\beta^N_e \gg \beta^I_e

则判为：

spatial-niche-driven event

如果：

\beta^I_e \gg \beta^N_e

则判为：

cell-intrinsic event

⸻

15. Bulk 扩展

Bulk 不是单细胞的低分辨率替代品，而是 MoDES 很重要的应用场景。

15.1 为什么 bulk 仍然重要

Bulk multi-omics 通常有：

更多 biological replicates
更多临床 metadata
更成熟的 proteomics / phosphoproteomics / metabolomics
更容易做 time-course
更容易连接 phenotype

很多 single-cell 数据没有 proteome 或 phosphoproteome，而 bulk 有。

⸻

15.2 Bulk 中的 event state

Bulk 中一样可以判定：

methylation ↑, RNA ↓
→ methylation-suppressed transcription
RNA ↑, protein unchanged
→ protein buffering
RNA unchanged, protein ↑
→ protein-level regulation / stability
ATAC ↑, RNA unchanged
→ epigenetic priming
RNA/protein effect disappears after cell composition covariate
→ composition-driven rather than molecular regulation

⸻

15.3 Bulk 必须加入 composition correction

bulk tissue 模型：

Y_{s,m,f}
=
\alpha
+
\beta Condition_s
+
\gamma CellComposition_s
+
\eta Purity_s
+
\lambda Batch_s
+
\epsilon_s

如果加入 cell composition 后 \beta 大幅下降，说明原始 effect 很可能是 composition-driven。

⸻

16. 验证策略

16.1 模拟数据

可以用 scMultiSim 构造 ground truth。scMultiSim 能生成 gene expression、chromatin accessibility、RNA velocity、spatial locations，并考虑 GRN、cell-cell interaction 和 technical noise；这非常适合构造 MoDES 的 benchmark，因为你可以人为设定哪些 event 是 primed、concordant、buffered 或 spatial-niche-driven。 ￼

模拟任务：

Task 1: recover concordant ATAC→RNA events
Task 2: recover chromatin-primed events
Task 3: recover RNA-only events
Task 4: recover protein-buffered events
Task 5: recover spatial-niche events
Task 6: distinguish artifact-like single-modality noise

⸻

16.2 真实单细胞数据

适合用：

10x Multiome RNA+ATAC
SHARE-seq RNA+ATAC
TEA-seq RNA+ATAC+ADT
CITE-seq RNA+ADT
perturbation multiome
stimulation time-course multiome

TEA-seq 尤其适合验证 RNA+ATAC+protein 三层状态，因为它同时测 transcriptomics、epitopes 和 chromatin accessibility。 ￼

⸻

16.3 真实空间数据

适合用：

spatial RNA+protein
spatial ATAC-RNA
spatial transcriptomics + single-cell ATAC/RNA reference
spatial transcriptomics + histology region labels

SIMO、MultiGATE、MISTy 可以作为 spatial baseline 或 comparison framework。SIMO 明确面向 spatial transcriptomics 与多种 non-spatial single-cell omics 的整合；MultiGATE 面向 spatial multi-omics integration 与 regulatory inference；MISTy 面向空间多视角 marker relationship 建模。 ￼

⸻

16.4 评价指标

不要只用 embedding 指标。MoDES 应该用 event-level 指标：

event state recovery accuracy
known perturbation target enrichment
motif enrichment
eQTL / GWAS enrichment for enhancer events
replication across donors
cross-cohort reproducibility
protein/RNA validation consistency
spatial autocorrelation of spatial-niche events
artifact detection specificity

SCARlink 使用 promoter capture Hi-C、fine-mapped eQTL 和 GWAS enrichment 验证 enhancer-gene links，这可以作为 MoDES 中 regulatory event validation 的参考。 ￼

⸻

17. 方法论文可以怎么写

17.1 标题

英文：

Discordance-aware regulatory event inference from bulk, single-cell and spatial multi-omics data

或：

MoDES: Multi-omics discordance-guided decomposition of regulatory event states across bulk, single-cell and spatial profiles

中文：

一种基于多组学不一致模式的调控事件状态分解方法

⸻

17.2 摘要逻辑

可以这样写：

现有多组学方法主要关注联合表示、聚类、补全和调控网络推断。
然而，不同组学层之间的不一致通常被视为噪音或后验解释对象，
而很少被系统建模为调控状态。
我们提出 MoDES，
一个以 regulatory event 为基本单位的统计框架。
MoDES 将 RNA、chromatin accessibility、protein、methylation 和 spatial context
映射到统一的 event-level evidence space，
并通过条件回归、层级模型和空间图正则，
将事件分类为 concordant activation、chromatin priming、RNA-only regulation、
protein buffering、protein memory、epigenetic memory、spatial niche effect
或 artifact-like discordance。
MoDES 可用于 bulk、single-cell 和 spatial multi-omics 数据，
并输出每个事件的层级效应、状态概率和可解释证据链。

⸻

18. 和“测一个组学，另一个组学验证”的区别

传统逻辑：

RNA 找差异基因
ATAC 看附近 peak 是否开放
protein 看 marker 是否变化

MoDES 逻辑：

RNA、ATAC、protein 不是主次关系；
它们是同一个事件的不同观测层。
一致性说明完整链条激活；
ATAC-only 说明 priming；
RNA-only 说明 trans/post-transcriptional/未测层调控；
RNA-protein 不一致说明 buffering/memory/technical issue；
空间局部一致说明 niche-driven regulation。

这就是数据利用率提升的本质。

⸻

19. 推荐开发路线

Phase 1：RNA+ATAC，bulk + single-cell pseudobulk

目标：

实现 event candidate construction
实现 ATAC/RNA effect estimation
实现 RNA-after-ATAC conditional decomposition
实现 4 个状态：
  concordant
  chromatin_primed
  RNA_only
  artifact_like

输入：

bulk RNA+ATAC
10x multiome pseudobulk

输出：

event_state.tsv
event_report.html

这是最小论文雏形。

⸻

Phase 2：加入 ADT/proteomics

新增状态：

protein_buffered
protein_memory
full_activation

输入：

CITE-seq
TEA-seq
bulk RNA+proteomics
bulk RNA+ATAC+proteomics

这一步会让方法从 RNA+ATAC regulatory inference 扩展到真正的 multi-layer regulation。

⸻

Phase 3：加入 spatial graph

新增状态：

spatial_niche_driven
region_specific_priming
spatial_artifact_like

输入：

Visium RNA+protein
spatial ATAC-RNA
spatial transcriptomics + single-cell reference

输出：

spatial event maps
niche-driven event table

⸻

Phase 4：加入 time / pseudotime lag

新增：

ATAC→RNA delay
RNA→protein delay
event ordering
early-driver score

这会让 MoDES 从静态状态分类变成动态调控分析。

⸻

20. 技术风险

风险 1：不一致可能是技术噪音

解决：

加入 QC covariates
加入 modality-specific reliability score
加入 donor-level replication
加入 artifact-like state
不要把所有 discordance 都解释成生物学

⸻

风险 2：peak-gene link 本身不可靠

解决：

第一版不要声称自己发明最好的 peak-gene link
复用 SCENIC+ / SCARlink / ArchR / Cicero 作为候选
MoDES 关注 event state decomposition

⸻

风险 3：bulk 和 single-cell 的统计单位不同

解决：

统一到 event-level evidence vector
bulk: sample-level effect
single-cell: pseudobulk effect + cell-level distribution
spatial: spot/cell effect + graph smoothness

⸻

风险 4：protein 数据覆盖有限

解决：

protein 层可选
RNA+ATAC 是核心版本
protein 用于扩展到 CITE-seq / TEA-seq / proteogenomics

⸻

21. 最终判断

这个方向值得做，但题目必须收紧。

最好的定位不是：

多组学整合方法

而是：

多组学调控事件状态分解方法

最好的第一版不是一口气支持所有组学，而是：

RNA + ATAC 的 condition effect decomposition

然后扩展到：

RNA + ATAC + protein
bulk + single-cell + spatial

我对新意的判断是：

版本	新意	风险
普通 multi-omics integration	低	已经很卷
RNA+ATAC peak-gene link	中低	SCENIC+/SCARlink 很强
多组学 concordance/discordance visualization	中	OmicsTIDE 已有相近概念
多组学 differential correlation	中	discordant/DRAGON 等已有相近方法
event-level discordance state inference	高	需要清楚 benchmark
bulk/sc/spatial 统一 event evidence framework	高	工程和模型设计要稳
加入 artifact-like state 与 spatial niche decomposition	更高	需要好数据验证

最实际的项目核心可以写成一句话：

MoDES 把 bulk、单细胞和空间多组学数据统一为 regulatory-event-level evidence，并把跨组学一致/不一致模式解释为可检验的调控状态，而不是把另一个组学仅作为验证层。

这条路线有方法学空间，也符合多组学趋势，而且能避开“又一个 embedding / AI integration tool”的拥挤赛道。
