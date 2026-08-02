# Getting Started with genpca

## Why generalized PCA?

Standard PCA ([`prcomp()`](https://rdrr.io/r/stats/prcomp.html)) treats
every observation and every variable equally. That is fine when the data
are homogeneous, but in many applications the assumption is wrong: some
observations are noisier than others, some variables are correlated by
design, or there is known spatial or temporal structure you want the
decomposition to respect. Generalized PCA encodes that prior knowledge
through a row metric `M` and a column metric `A`. When both are
identity, you recover ordinary PCA.

## Quick start

``` r

set.seed(1)
X <- matrix(rnorm(80 * 30), 80, 30)
fit <- genpca(X, ncomp = 3, preproc = multivarious::center())
fit$sdev
#> [1] 14.29106 13.45933 12.99494
```

![Singular values from the 3-component
fit.](genpca_files/figure-html/quick-plot-1.png)

Singular values from the 3-component fit.

With identity metrics,
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
reproduces standard PCA. The interesting cases are when `M` or `A` is
non-trivial.

## The decomposition

### Two metrics, one inner product

Let $`X`$ be the $`n \times p`$ data matrix. Generalized PCA asks you to
supply two symmetric positive semi-definite matrices:

- $`M`$ ($`n \times n`$), the **row metric**, which puts a geometry on
  the space of observations;
- $`A`$ ($`p \times p`$), the **column metric**, which puts a geometry
  on the space of variables.

Together they define an inner product on $`n \times p`$ matrices,

``` math
\langle Y, Z \rangle_{M,A} \;=\; \operatorname{tr}\!\left(Y^{\top} M Z A\right),
\qquad
\|Y\|_{M,A}^{2} \;=\; \operatorname{tr}\!\left(Y^{\top} M Y A\right).
```

When $`M = I_n`$ and $`A = I_p`$ this is the ordinary Frobenius inner
product, and everything below collapses to the familiar SVD and PCA.

### The problem it solves

The **generalized least-squares matrix decomposition** (GMD) of Allen,
Grosenick & Taylor (2014) is the best rank-$`K`$ approximation of
$`X`$*in that norm*:

``` math
\min_{U, D, V} \;\bigl\|X - U D V^{\top}\bigr\|_{M,A}^{2}
\qquad \text{subject to} \qquad
U^{\top} M U = I_K, \quad V^{\top} A V = I_K,
```

with $`D = \operatorname{diag}(d_1 \ge d_2 \ge \cdots \ge d_K \ge 0)`$.

Only two things changed relative to the SVD, and they are the same
change twice: the discrepancy is measured in the $`(M,A)`$ norm rather
than the Frobenius norm, and the factors are orthonormal *in the
metrics* rather than Euclidean-orthonormal.

### The solution: whiten, decompose, unwhiten

The problem has a closed form. Take symmetric square roots
$`M = M^{1/2}M^{1/2}`$ and $`A = A^{1/2}A^{1/2}`$, whiten the data,

``` math
\widetilde{X} \;=\; M^{1/2} X A^{1/2},
```

and compute its *ordinary* SVD,
$`\widetilde{X} = \widetilde{U} D \widetilde{V}^{\top}`$. Then

``` math
U = M^{-1/2}\widetilde{U}, \qquad V = A^{-1/2}\widetilde{V}, \qquad D \text{ unchanged}
```

solves the problem above (with pseudo-inverses when a metric is
singular). That is the whole idea: **GMD is an SVD performed in a
different coordinate system.** Everything else is bookkeeping and
computation.

Eliminating the whitening gives the form the solvers actually use — an
eigenproblem for a single operator,

``` math
X^{\top} M X A \, v_k \;=\; d_k^{2}\, v_k, \qquad v_k^{\top} A v_k = 1,
```

and the $`k`$-th generalized principal component (the score vector) is

``` math
z_k \;=\; X A v_k \;=\; d_k\, u_k .
```

In the package these are `components(fit)` $`= A\,`$`ov` for the
loadings and `scores(fit)` $`=`$`ou %*% diag(sdev)` for the scores, with
`ou` and `ov` holding the metric-orthonormal factors $`U`$ and $`V`$.

### What the metrics mean

There are two readings of $`M`$ and $`A`$, and both are worth carrying.

**As geometry.** The metrics define what “distance” means, and they act
crosswise. Observations are rows of $`X`$, living in $`\mathbb{R}^p`$,
so the *column* metric sets the distance between two observations,
$`\|x_i - x_j\|_A^2 = (x_i - x_j)^{\top} A\,(x_i - x_j)`$. Variables are
columns, living in $`\mathbb{R}^n`$, so the *row* metric sets the
distance between two variables — and, when $`M`$ is diagonal, simply
weights the observations. Choose $`A`$ to say which variables should be
treated as similar; choose $`M`$ to say which observations count more.

**As statistics.** Suppose the noise is separable, meaning the errors
follow a matrix-normal law
$`E \sim \mathcal{MN}(0, \Sigma_{\text{row}}, \Sigma_{\text{col}})`$,
equivalently $`\operatorname{vec}(E)`$ has covariance
$`\Sigma_{\text{col}} \otimes \Sigma_{\text{row}}`$. Then the
log-likelihood of a rank-$`K`$ mean is, up to constants, exactly
$`-\tfrac12\|X - UDV^{\top}\|^2_{M,A}`$ with

``` math
M = \Sigma_{\text{row}}^{-1}, \qquad A = \Sigma_{\text{col}}^{-1}.
```

So the metrics are not an arbitrary regularizer: setting them to the
inverse noise covariances makes GMD the maximum-likelihood low-rank fit.
This is why the package ships
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md),
which learns $`(M, A)`$ by penalized maximum likelihood, and
[`mnpca_mrl()`](https://bbuchsbaum.github.io/genpca/reference/mnpca_mrl.md),
which does the same with sparse precision matrices.

## Relation to the French school

Readers coming from *analyse des données* will recognize all of this.
The French tradition — Benzécri’s correspondence analysis, formalized by
Escoufier as the **duality diagram** — describes an analysis by a
*triplet* $`(X, Q, D)`$: a data table, a metric on the variable space,
and a set of row masses. The analysis is then the eigendecomposition of
the operator $`X^{\top} D X Q`$. That is, term for term, the
eigenproblem two sections above. Correspondence analysis, multiple
correspondence analysis, principal coordinates, discriminant analysis,
canonical correlation, and PLS are all recovered by choosing the triplet
appropriately; the `ade4` package is built directly on this formalism,
and Abdi’s papers on the generalized SVD — cited in this package’s
`DESCRIPTION` — are the same construction in English.

The one genuine trap is notation: **$`Q`$ means opposite things in the
two literatures.**

| Role | `genpca` | Allen et al. (2014) | Duality diagram |
|:---|:---|:---|:---|
| Row metric / observation weights | `M` | $`Q`$ | $`D`$ |
| Column metric / variable geometry | `A` | $`R`$ | $`Q`$ |

What Allen, Grosenick & Taylor contribute on top of the classical theory
is less the geometry than the framing and the machinery: the
decomposition is posed as a *best rank-$`K`$ approximation* with an
optimality theorem rather than only an eigenanalysis; the metrics are
allowed to be large, dense, and structured (spatial kernels, temporal
covariances) instead of the diagonal weights typical of the classical
treatments; and the formulation extends cleanly to regularized and
sparse variants, which is what
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
implements.

## Where this pays off

**Genomics and microbiome data.** Genes in a pathway or taxa on a
phylogeny are not exchangeable. A column metric built from a
gene-network Laplacian or a phylogenetic distance encodes that, so
components respect known biology; double principal coordinate analysis
is precisely this choice of triplet. On the row side, a diagonal $`M`$
down-weights samples with low sequencing depth, and a kinship matrix as
$`M`$ absorbs relatedness.

**Imaging.** In fMRI the natural table is time $`\times`$ voxels.
Neighbouring voxels are correlated, so an $`A`$ built from the voxel
adjacency graph or a spatial smoothing operator yields spatially
coherent components instead of salt-and-pepper maps — the motivating
example in Allen et al. On the row side, $`M`$ can carry AR(1) temporal
whitening or simply down-weight motion-corrupted frames.

**Time series and functional data.** When the variables are ordered —
time points, wavelengths, positions along a transect — you usually want
loadings that vary smoothly rather than jumping between neighbours. GPCA
delivers that through the metric, which is Silverman’s “smoothing by
choice of norm” construction: components become smooth functions of the
index. This is the foundation
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
builds on, adding sparsity to smoothness.

One direction trap is worth stating explicitly, because the two
interfaces in this package take *opposite* inputs. Write
$`\Omega = D_2^{\top}D_2`$ for a second-difference roughness operator,
so that $`v^{\top}\Omega v`$ is large for wiggly $`v`$.

- In
  [`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md),
  $`A`$ is a **metric**, and a metric amplifies its own dominant
  directions: `components(fit)` $`= A\,`$`ov`. To get smooth loadings,
  pass a *smoothing* operator — $`A = (I + \alpha\Omega)^{-1}`$, or a
  kernel or adjacency matrix, as in the banded example above. Passing
  $`\Omega`$ itself makes loadings **rougher**, not smoother.
- In
  [`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md),
  the same information enters as a **constraint**,
  $`v^{\top}(I + \alpha\Omega)\,v \le 1`$, which charges wiggly $`v`$
  against a fixed budget. There you supply the roughness operator
  directly (via `alpha_v` and `spat_cds`), and larger `alpha_v` means
  smoother.

Metric form and constraint form are inverse to one another; mixing them
up silently produces the opposite of what you wanted, so it is worth
checking the smoothness of your loadings against a small sweep of
$`\alpha`$.

**Space $`\times`$ time.** The case where both metrics earn their keep
at once: $`M`$ from the temporal covariance, $`A`$ from the spatial
covariance, on a single space-by-time table. Under a separable noise
model this is exactly the maximum-likelihood fit from the statistics
reading above, which makes it the principled choice rather than a
convenient one. The same structure appears in EEG and MEG (sensor
covariance estimated from empty-room noise), in climate fields, and
generally in any table whose two margins each carry their own
dependence.

## How the package computes it

The literal recipe — form $`M^{1/2}`$ and $`A^{1/2}`$, whiten, take an
SVD — costs $`O(n^3 + p^3)`$ and densifies everything in sight. That is
fine for small dense problems and hopeless for a spatial metric over
200,000 voxels. So
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
offers several backends that compute the same decomposition by different
routes:

| `method` | Approach | Best for |
|:---|:---|:---|
| `"eigen"` | Works on the smaller side: a $`p \times p`$ Gram eigenproblem when $`n \ge p`$, an $`n \times n`$ one when $`n < p`$. Metric square roots are cached. | Small to medium dense problems (the default) |
| `"spectra"` | Matrix-free. Supplies only the *action* $`w \mapsto X^{\top}MXA\,w`$ to a C++ Krylov eigensolver; neither the whitened matrix nor the Gram matrix is ever formed. | Large problems with sparse structured metrics |
| `"randomized"` | Randomized range finder with power iterations and a polish step, all in the metric geometry. | Wide ($`p \gg n`$), low-rank problems |
| `"deflation"` | Generalized power iteration one component at a time, with the residual applied implicitly so sparse inputs stay sparse. | Tight memory budgets |

Underneath, metrics are handled as operator closures rather than
materialized matrices, with element-wise fast paths for diagonal metrics
and Cholesky factors for general SPD ones.
[`genpca_cov()`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md)
covers the case where you only have the cross-product $`C = X^{\top}MX`$
rather than $`X`$ itself. The backends are mathematically equivalent and
the test suite asserts they agree; they differ only in cost and in the
accuracy trade-offs of iterative solvers.

### References

Allen, G. I., Grosenick, L., & Taylor, J. (2014). A generalized
least-square matrix decomposition. *Journal of the American Statistical
Association*, 109(505), 145–159.

Escoufier, Y. (1987). The duality diagram: a means for better practical
applications. In *Development in Numerical Ecology*. Springer.

Holmes, S. (2008). Multivariate data analysis: the French way. In
*Probability and Statistics: Essays in Honor of David A. Freedman*. IMS.

Dray, S., & Dufour, A.-B. (2007). The ade4 package: implementing the
duality diagram for ecologists. *Journal of Statistical Software*,
22(4).

Silverman, B. W. (1996). Smoothed functional principal components
analysis by choice of norm. *The Annals of Statistics*, 24(1), 1–24.

## Row weights (diagonal metric)

A diagonal `M` says “weight observation `i` by `m_ii`.” Useful when rows
represent groups of different sizes, or when some rows are more reliable
than others.

``` r

row_wt <- runif(nrow(X), 0.5, 1.5)
M <- Diagonal(x = row_wt)
fit_row <- genpca(X, M = M, ncomp = 3, preproc = multivarious::center())
head(multivarious::scores(fit_row), 4)
#>             PC1         PC2        PC3
#> Obs1  0.7120376 -0.66586763  0.3811552
#> Obs2  1.8471851 -1.09234435 -0.5769734
#> Obs3  0.1803497  3.08148414  1.3099005
#> Obs4 -1.8055237 -0.08106447 -2.1575624
```

![Identity (left) vs row-weighted (right) component scores. Symbol size
in the right panel scales with row
weight.](genpca_files/figure-html/row-compare-1.png)

Identity (left) vs row-weighted (right) component scores. Symbol size in
the right panel scales with row weight.

## Column structure (smoothness)

A non-identity column metric encodes prior structure on the variables.
Here we add a small AR-style coupling between adjacent columns:

``` r

p <- ncol(X)
A <- Diagonal(p)
for (i in seq_len(p - 1)) {
  A[i, i + 1] <- 0.15
  A[i + 1, i] <- 0.15
}
A <- A + 0.05 * Diagonal(p)
fit_col <- genpca(X, A = A, ncomp = 3, preproc = multivarious::center())
fit_col$sdev
#> [1] 14.79937 14.00442 13.28632
```

![Banded column metric A (left) and the resulting PC1 loadings, identity
vs A (right). The banded A pulls adjacent loadings together, suppressing
high-frequency wiggle.](genpca_files/figure-html/col-metric-plot-1.png)

Banded column metric A (left) and the resulting PC1 loadings, identity
vs A (right). The banded A pulls adjacent loadings together, suppressing
high-frequency wiggle.

## Choosing the number of components

`fit$sdev` holds the generalized singular values `d_k` of the
metric-whitened data, not standard deviations: with identity metrics and
centering, `fit$sdev` equals `prcomp(X)$sdev * sqrt(nrow(X) - 1)`. A
scree plot of the squared values is the quickest first look:

``` r

barplot(fit$sdev^2,
        names.arg = paste0("PC", seq_along(fit$sdev)),
        xlab = "Component", ylab = "Squared generalized singular value",
        col = "grey60", border = NA)
```

![Squared generalized singular value by
component.](genpca_files/figure-html/scree-1.png)

Squared generalized singular value by component.

For a more principled choice, split the rows into train and holdout, fit
on the training set, and compare holdout reconstruction error across
values of `ncomp`. With informative metrics, a small `ncomp` often beats
a larger unweighted PCA.

## Real data: USArrests

``` r

data("USArrests")
X_real  <- as.matrix(USArrests[, c("Murder", "Assault", "Rape")])
pop_wt  <- USArrests$UrbanPop / mean(USArrests$UrbanPop)
M_real  <- Diagonal(x = pop_wt)
col_sd  <- apply(X_real, 2, sd)
A_real  <- Diagonal(x = 1 / col_sd^2)

fit_real <- genpca(X_real, M = M_real, A = A_real, ncomp = 2,
                   preproc = multivarious::center())
scores2d <- multivarious::scores(fit_real)
```

![GPCA on USArrests with row weights from urbanisation and column
weights inversely proportional to
variance.](genpca_files/figure-html/usarrests-plot-1.png)

GPCA on USArrests with row weights from urbanisation and column weights
inversely proportional to variance.

PC1 separates states by overall violent-crime rate, with
high-urbanisation states pulled toward the high end through the row
weights. PC2 picks up the residual contrast between assault-heavy and
rape-heavy states once the dominant rate axis is removed.

## Object verbs

[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
returns a `bi_projector` from the `multivarious` ecosystem. The verbs
you reach for most:

| Verb | Returns | What it gives you |
|:---|:---|:---|
| `multivarious::scores(fit)` | `n x k` | Observation coordinates in component space |
| `multivarious::components(fit)` | `p x k` | Loadings (variable directions) |
| `multivarious::project(fit, X_new)` | `n_new x k` | Out-of-sample scores |
| `reconstruct(fit, comp = 1:k)` | `n x p` | Rank-`k` reconstruction of `X` |

These verbs are shared across
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md),
[`genpls()`](https://bbuchsbaum.github.io/genpca/reference/genpls.md),
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md), and
[`rpls()`](https://bbuchsbaum.github.io/genpca/reference/rpls.md), so
once you learn them on `genpca` they carry over.

## What are `ou`, `ov`, `u`, and `v`?

Under the hood,
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
stores the raw factors of the decomposition alongside the verbs above.
`ou` and `ov` are the metric-orthonormal factors of the paper (Allen,
Grosenick & Taylor 2014): `ou` is orthonormal in the row metric
(`t(ou) %*% M %*% ou` is the identity) and `ov` is orthonormal in the
column metric (`t(ov) %*% A %*% ov` is the identity). `u` and `v` are
the same factors after applying the metric (`u = M %*% ou`,
`v = A %*% ov`), which is the form used for
[`components()`](https://bbuchsbaum.github.io/multivarious/reference/components.html).
The table below ties these to the verbs:

| Quantity | Definition | Relation to verbs |
|:---|:---|:---|
| `ou` | Row-metric-orthonormal left factor | `t(ou) %*% M %*% ou = I` |
| `ov` | Column-metric-orthonormal right factor | `t(ov) %*% A %*% ov = I` |
| `v` (`components(fit)`) | `A %*% ov` | Loadings in the original variable space |
| `s` (`scores(fit)`) | `X %*% A %*% ov = ou %*% diag(sdev)` | Equal to `project(fit, X)` on the training data |

In other words, `scores(fit)` and `project(fit, X)` compute the same
quantity two ways – one directly from `X` and the column loadings, the
other from the stored orthonormal factors and the singular values – and
they agree exactly on the training data.

## Where next

- [`vignette("gpca-metrics")`](https://bbuchsbaum.github.io/genpca/articles/gpca-metrics.md)
  covers metric recipes (AR(1), spatial kernels, graph Laplacians), the
  smoother-versus-precision orientation rule, and the
  [`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
  learner.
- [`vignette("structured-noise")`](https://bbuchsbaum.github.io/genpca/articles/structured-noise.md)
  works out which noise a metric can and cannot remove, with worked
  simulations – start here if your data has several kinds of structure
  at once.
- [`vignette("gpca-scale")`](https://bbuchsbaum.github.io/genpca/articles/gpca-scale.md)
  covers backends, sparse workflows, and covariance-only GPCA.
- [`vignette("gplssvd-reference")`](https://bbuchsbaum.github.io/genpca/articles/gplssvd-reference.md)
  covers the generalized PLS / PLS-SVD interface and includes a
  contributor-oriented reference implementation.
