---
title: "T4 Checkpoint 1: How each method calculates its result"
subtitle: "Group T4, Tumor surface target prioritization (CPSC 4430). Code repository: https://github.com/JacksonMahr782/Machine-Learning"
---

This sheet gives the equation behind every feature, label, model and metric in our plan. Part A covers the LUAD gene-dependency model. Part B covers the statistics in the TROP2 paper we reproduce.

**Notation.** $g$ is a gene, $s$ is a TCGA sample, $\ell$ is a DepMap cell line, $n$ is the number of genes in a split, $y_g$ is the true label and $\hat y_g$ is the model's prediction.

# Part A. LUAD gene-dependency model (one row per gene)

## A1. Features (from public data, never from DepMap)

Expression from Xena is $x_{gs} = \log_2(\text{TPM}_{gs} + 1)$. $T$ is the set of LUAD tumour samples and $N$ is the set of normal lung samples.

**Mean tumour expression** (`expr_mean_tumor`):
$$\mu_g^{T} = \frac{1}{|T|}\sum_{s \in T} x_{gs}$$

**Tumour-vs-normal enrichment** (`expr_enrichment`):
$$e_g = \mu_g^{T} - \mu_g^{N}, \qquad \mu_g^{N} = \frac{1}{|N|}\sum_{s \in N} x_{gs}$$

**Variability** (`expr_sd_tumor`):
$$\sigma_g = \sqrt{\frac{1}{|T|-1}\sum_{s \in T}\left(x_{gs} - \mu_g^{T}\right)^2}$$

The amplification features come from cBioPortal GISTIC calls, $c_{gs} \in \{-2,-1,0,+1,+2\}$, over the $S$ tumour samples. $\mathbf{1}[\cdot]$ is 1 when the condition is true and 0 otherwise.

**Amplification frequency** (`gistic_amp_freq`):
$$a_g = \frac{1}{|S|}\sum_{s \in S}\mathbf{1}[c_{gs} \ge 1]$$

**High-level amplification frequency** (`gistic_highamp_freq`):
$$h_g = \frac{1}{|S|}\sum_{s \in S}\mathbf{1}[c_{gs} = 2]$$

**Surface flag** (`is_surface`, from the Human Protein Atlas). It is used for eligibility and stratification only, never as a model input:
$$u_g = \mathbf{1}[\text{HPA main or additional location includes "Plasma membrane"}]$$

The feature vector for gene $g$ is $\mathbf{x}_g = (\mu_g^{T}, e_g, \sigma_g, a_g, h_g)$.

## A2. Label (from DepMap only)

$L_g$ is the set of LUAD cell lines (Model.csv `OncotreeCode = LUAD`) with a non-missing Chronos effect $c_{g\ell}$ for gene $g$:
$$y_g = -\operatorname{median}_{\ell \in L_g}\, c_{g\ell} \quad \text{only if } |L_g| \ge 3; \text{ otherwise the gene is unlabeled.}$$

A higher $y_g$ means LUAD cells depend on the gene more. Here $0$ means no effect and about $1$ is a typical essential gene. To prevent leakage, no other DepMap-derived value enters $\mathbf{x}_g$.

## A3. Split and scaling

Labeled genes are split 70/15/15 into train, validation and test, stratified on $u_g$ with seed 42. That gives every split the same share of surface genes. Each feature is standardized with the training mean and SD only:
$$z_{gj} = \frac{x_{gj} - \bar x_j^{\text{train}}}{s_j^{\text{train}}}$$

## A4. Models

**Constant baseline.** It predicts the training mean for every gene:
$$\hat y_g = \bar y^{\text{train}} = \frac{1}{n_{\text{train}}}\sum_{g \in \text{train}} y_g$$

**Ridge regression.** It is a linear model with an L2 penalty. $\lambda$ is chosen on the validation set:
$$\hat y_g = \beta_0 + \mathbf{z}_g^{\top}\boldsymbol\beta, \qquad \hat{\boldsymbol\beta} = \arg\min_{\boldsymbol\beta}\sum_{g \in \text{train}}\left(y_g - \beta_0 - \mathbf{z}_g^{\top}\boldsymbol\beta\right)^2 + \lambda\lVert\boldsymbol\beta\rVert_2^2$$

The penalized coefficients have a closed form, with $y$ and $Z$ centred:
$$\hat{\boldsymbol\beta} = (Z^{\top}Z + \lambda I)^{-1}Z^{\top}\mathbf{y}$$

**Random forest.** It averages $B$ regression trees. Each tree $T_b$ is grown on a bootstrap sample of training genes. At every split it picks the feature and threshold that most reduce the squared error inside the two child nodes:
$$\hat y_g = \frac{1}{B}\sum_{b=1}^{B} T_b(\mathbf{x}_g)$$

**Gradient boosting.** It adds small trees one at a time. Each new tree $f_m$ is fitted to the current residuals, which are the negative gradient of the squared loss. $\nu$ is the learning rate:
$$F_0 = \bar y^{\text{train}}, \qquad r_g^{(m)} = y_g - F_{m-1}(\mathbf{x}_g), \qquad F_m(\mathbf{x}) = F_{m-1}(\mathbf{x}) + \nu f_m(\mathbf{x})$$

**Training loss (all models).** It is the mean squared error:
$$\text{MSE} = \frac{1}{n}\sum_{g}\left(y_g - \hat y_g\right)^2$$

## A5. Evaluation metrics

Each metric is reported twice on the held-out test genes: once over all genes, and once over surface genes only ($u_g = 1$).

**Mean absolute error** (lower is better):
$$\text{MAE} = \frac{1}{n}\sum_{g}\left|y_g - \hat y_g\right|$$

**Spearman rank correlation** (higher is better; it checks whether the model puts genes in the right order). With $d_g$ the difference between gene $g$'s rank by $y$ and its rank by $\hat y$, and no ties:
$$\rho = 1 - \frac{6\sum_g d_g^2}{n(n^2-1)}$$

With ties, $\rho$ is the Pearson correlation of the ranks.

## A6. Improvement and ablation

**Improvement.** Expression × amplification interaction terms are added, so the model can learn that an amplified gene that is also highly expressed matters more:
$$\mathbf{x}_g^{+} = \left(\mathbf{x}_g,\; \mu_g^{T} a_g,\; e_g a_g\right)$$

**Amplification ablation.** All amplification-derived columns ($a_g$, $h_g$ and the interaction terms) are removed and the model is refitted on the same split. The effect of amplification is the change in test performance:
$$\Delta\text{MAE} = \text{MAE}_{\text{no amp}} - \text{MAE}_{\text{full}}, \qquad \Delta\rho = \rho_{\text{full}} - \rho_{\text{no amp}}$$

Positive values mean amplification helps.

## A7. Ranking surface candidates

Test genes with $u_g = 1$ (at least 20) are ranked with the course public-data guide's locked prioritization score. Predicted dependency $\hat y_g$ is one input and counts as supporting evidence only. Nonessential surface antigens stay eligible. The exact score formula is taken from the T4 project overview and is fixed before test genes are scored.

# Part B. TROP2 paper reproduction (E-MTAB-16585, Figure 1)

## B1. What we reproduce

This part reproduces two things. The first is TACSTD2 (TROP2) and LGR5 expression in each cell state. The second is the proportion of each cell state in each donor. The donor, not the cell, is the unit of analysis.

**Per-cell normalized expression** (typical single-cell pipeline; confirm the exact steps in Figure_01.ipynb). Here $k_{ig}$ is the raw count of gene $g$ in cell $i$:
$$\tilde x_{ig} = \log\left(1 + 10^4 \cdot \frac{k_{ig}}{\sum_{g'} k_{ig'}}\right)$$

**Donor cell-state proportion.** For donor $d$ and state $k$:
$$p_{d,k} = \frac{n_{d,k}}{\sum_{k'} n_{d,k'}}$$

Here $n_{d,k}$ is the number of donor $d$'s cells annotated as state $k$.

## B2. Statistical tests the paper uses

**Wilcoxon rank-sum test** (comparing two groups). Rank all values together. $R_1$ is the sum of group 1's ranks:
$$U = R_1 - \frac{n_1(n_1+1)}{2}$$

**Benjamini–Hochberg correction** (many genes tested at once). Sort the $m$ p-values $p_{(1)} \le \dots \le p_{(m)}$. Then reject every hypothesis up to the largest $k$ with:
$$p_{(k)} \le \frac{k}{m}\,q$$

**Log-rank test** (survival curves in the 1,336-patient cohorts). At each event time $j$, $O_{1j}$ is the observed number of deaths in group 1. The expected number is $E_{1j} = d_j\, n_{1j}/n_j$, and $V_j$ is its variance:
$$\chi^2_{\text{LR}} = \frac{\left(\sum_j (O_{1j} - E_{1j})\right)^2}{\sum_j V_j}, \qquad 1 \text{ degree of freedom}$$

**Chi-squared test** (counts in categories):
$$\chi^2 = \sum_{\text{cells}}\frac{(O - E)^2}{E}$$

**Welch's t-test** (two means):
$$t = \frac{\bar x_1 - \bar x_2}{\sqrt{s_1^2/n_1 + s_2^2/n_2}}$$

**Limiting-dilution assay** (stem-cell frequency $f$). Under a single-hit Poisson model, the chance a well seeded with $c$ cells shows no growth is:
$$P(\text{no growth}) = e^{-f c}$$

$f$ is estimated by maximum likelihood. Groups are compared with a likelihood-ratio test.
