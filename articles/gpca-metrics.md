# GPCA Metrics: Building M and A

This vignette collects practical recipes for row and column metrics,
plus notes on SPD remedies and the experimental
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
learner.

## Why metrics

Metrics encode weighting and correlation. The row metric `M` changes how
observations are compared; the column metric `A` changes how variables
are compared. Setting either to something other than the identity is how
you tell the decomposition what you already know about the data: that
the samples are a time series, that the variables sit on a spatial grid
or fall into groups, that some measurements are noisier than others.
Ordinary PCA has no way to accept that information and treats every row
and column alike. GPCA writes it into the objective being optimised.

Structure alone is not enough, though: you also have to get the
*direction* right, and supplying a matrix versus its inverse produces
opposite results. The next section is about that, and it is worth
reading before the recipes.

Both metrics must also be symmetric and positive semi-definite. That
rarely gets in the way, but it does constrain what you can pass; see
[SPD requirements and remedies](#spd-requirements-and-remedies) below
for what the requirement means and what
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
does when a metric falls short of it.

## Heteroscedastic diagonals

The simplest non-trivial metric is a diagonal that down-weights noisy
rows or columns:

``` r

set.seed(42)
n <- 60; p <- 20
X <- matrix(rnorm(n * p), n, p)

col_noise_sd <- runif(p, 0.5, 2)
A <- Diagonal(x = 1 / col_noise_sd^2)
row_noise_sd <- runif(n, 0.7, 1.3)
M <- Diagonal(x = 1 / row_noise_sd^2)

fit <- genpca(X, M = M, A = A, ncomp = 3,
              preproc = multivarious::center())
fit$sdev
#> [1] 14.69818 11.24551 10.96887
```

![Inverse-variance weights on columns (top) and rows (bottom). Noisier
dimensions get smaller
weights.](gpca-metrics_files/figure-html/hetero-plot-1.png)

Inverse-variance weights on columns (top) and rows (bottom). Noisier
dimensions get smaller weights.

### Inverse column variance *is* scaled PCA

Inverse-variance weighting on the columns is not a new idea in disguise:
it is exactly the standardisation that `prcomp(scale. = TRUE)` performs.
GPCA whitens with the square root $`A^{1/2}`$, so setting
$`A = \operatorname{diag}(1/s_j^2)`$ makes
$`X A^{1/2} = X \operatorname{diag}(1/s_j)`$ — the column-scaled matrix
that correlation-matrix PCA decomposes.

``` r

set.seed(42)
Xv <- matrix(rnorm(60 * 20), 60, 20) %*% diag(runif(20, 0.5, 3))
sds <- apply(Xv, 2, sd)

g  <- genpca(Xv, A = Diagonal(x = 1 / sds^2), ncomp = 5,
             preproc = multivarious::center())
pr <- prcomp(Xv, scale. = TRUE)

# scores agree component by component
sapply(1:3, function(k) cor(multivarious::scores(g)[, k], pr$x[, k]))
#> [1] -1 -1 -1
```

The scores are identical (up to sign). The reported `sdev` values differ
by a single constant, because
[`prcomp()`](https://rdrr.io/r/stats/prcomp.html) divides its singular
values by $`\sqrt{n-1}`$ and
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
reports them unnormalised:

``` r

rbind(genpca = g$sdev[1:5],
      prcomp = pr$sdev[1:5],
      ratio  = g$sdev[1:5] / pr$sdev[1:5])
#>             [,1]      [,2]      [,3]     [,4]     [,5]
#> genpca 11.527935 10.943470 10.310907 9.651052 9.172200
#> prcomp  1.500809  1.424718  1.342366 1.256460 1.194119
#> ratio   7.681146  7.681146  7.681146 7.681146 7.681146
sqrt(nrow(Xv) - 1)
#> [1] 7.681146
```

The constant ratio is the whole difference. This is worth internalising
as the baseline: **a diagonal `A` generalises column scaling**, and
everything else in this vignette — kernels, Laplacians, AR(1) structure
— generalises it further by letting the metric go off-diagonal.

### Weighting both margins at once

The example above sets `M` and `A` simultaneously, which is a reasonable
thing to want when both observations and variables are heteroscedastic.
Two consequences are easy to miss.

**The two weightings interact.** GPCA works with $`M^{1/2} X A^{1/2}`$,
so the row weights change each column’s effective variance and the
column weights change each row’s. Estimating row and column standard
deviations from the raw data and applying both at once therefore
standardises *neither* margin:

``` r

Xc <- scale(Xv, center = TRUE, scale = FALSE)
W  <- diag(1 / apply(Xc, 1, sd)) %*% Xc %*% diag(1 / apply(Xc, 2, sd))

range(apply(W, 1, sd))   # row SDs, would be constant if standardised
#> [1] 0.4075883 0.6182529
range(apply(W, 2, sd))   # column SDs
#> [1] 0.4740343 0.5317369
```

It moves in the right direction — the row SDs are markedly more uniform
than in `Xc` — but one pass does not land on unit variance either way,
because each rescaling perturbs the other’s normalisation. Getting both
margins standardised requires alternating between them until they
settle, which is precisely what
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
does when it iterates between `M` and `A`.

**Only the product of the two scales is identified.** Replacing
$`(M, A)`$ with $`(cM, A/c)`$ leaves $`M^{1/2} X A^{1/2}`$ untouched, so
the fit cannot distinguish them:

``` r

M0 <- Diagonal(x = 1 / apply(Xv, 1, sd)^2)
A0 <- Diagonal(x = 1 / sds^2)

f1 <- genpca(Xv, M = M0,     A = A0,     ncomp = 4, preproc = multivarious::center())
f2 <- genpca(Xv, M = 7 * M0, A = A0 / 7, ncomp = 4, preproc = multivarious::center())

max(abs(f1$sdev - f2$sdev))
#> [1] 3.552714e-15
```

The singular values and the component subspace are identical; only the
scores pick up a constant factor, since their normalisation is tied to
the scale of `M`. The practical upshot is that there is no point tuning
the overall magnitude of `M` against that of `A` — it is the *relative*
weighting within each metric that changes the answer. This indeterminacy
is also why
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
takes a `scale_fix` argument: when both metrics are learned, the split
has to be pinned down by convention rather than by the data.

## Which way does a metric point?

Before the recipes, the single most important thing to get right — and
the easiest to get backwards. **A metric amplifies its own dominant
eigendirections.**

Seeing why means being precise about what comes back from a fit. GPCA
factorises $`X \approx U D V^{\top}`$, where $`V`$ is orthonormal *in
the column metric* rather than in the ordinary sense:
$`V^{\top} A V = I`$. The loadings returned by `components(fit)` are not
$`V`$ but $`AV`$:

``` r

set.seed(1)
Xd <- matrix(rnorm(400), 40, 10)
Ad <- crossprod(matrix(rnorm(100), 10, 10)) / 10 + diag(10)
fd <- genpca(Xd, A = Ad, ncomp = 3, preproc = multivarious::center())

max(abs(multivarious::components(fd) - as.matrix(Ad %*% fd$ov)))
#> [1] 0
```

That multiplication by $`A`$ is the whole story. It stretches every
direction in proportion to the eigenvalue $`A`$ assigns it, so the
patterns $`A`$ scores highly are the patterns that dominate the loadings
you read off. (The bare factor $`V`$ is kept in the `ov` slot if you
ever need it, but
[`components()`](https://bbuchsbaum.github.io/multivarious/reference/components.html)
is what you should normally interpret.)

For a spatial or temporal structure there are two natural matrices, and
they point in *opposite* directions:

| You supply | $`v^{\top}Av`$ measures | Large eigenvalues on | Components come out |
|:---|:---|:---|:---|
| A **smoother**: kernel $`K`$, adjacency $`I + \alpha W`$, $`(I+\alpha L)^{-1}`$, heat kernel $`e^{-tL}`$ | agreement between neighbours | smooth patterns | **smoother** |
| A **precision**: Laplacian $`L`$, $`I + \alpha L`$, $`K^{-1}`$, inverse AR(1) | disagreement across edges (Dirichlet energy $`\sum_{i\sim j}(v_i - v_j)^2`$) | rough patterns | **rougher** |

Both are legitimate, because they encode different beliefs about *where
the noise lives*. The bridge is $`A = \Sigma_{\text{col}}^{-1}`$, and
the step that is easy to skip is the inversion: the metric and the noise
covariance share eigenvectors but carry **reciprocal** weights, so a
direction the metric scores highly is a direction the noise model calls
quiet.

That reciprocal is worth seeing rather than taking on faith. On a cycle
graph the Laplacian and the adjacency share Fourier eigenvectors
exactly, so “roughness” is unambiguously frequency:

``` r

p <- 32
Wc <- matrix(0, p, p)
for (i in 1:p) { Wc[i, i %% p + 1] <- 1; Wc[i %% p + 1, i] <- 1 }
Lc <- diag(rowSums(Wc)) - Wc

ec <- eigen(Lc, symmetric = TRUE)
o  <- order(ec$values)
Vc <- ec$vectors[, o]          # smoothest first
lc <- ec$values[o]             # Dirichlet energy = roughness

Ac  <- diag(p) + 0.45 * Wc     # a smoother (PSD)
Sig <- solve(Ac)               # the noise covariance it implies

modes <- c(1, 16, 32)          # smoothest, middling, roughest
data.frame(
  roughness   = round(lc[modes], 3),
  metric_wt   = round(sapply(modes, function(k) t(Vc[, k]) %*% Ac  %*% Vc[, k]), 3),
  noise_var   = round(sapply(modes, function(k) t(Vc[, k]) %*% Sig %*% Vc[, k]), 3)
)
#>   roughness metric_wt noise_var
#> 1         0       1.9     0.526
#> 2         2       1.0     1.000
#> 3         4       0.1    10.000
```

The metric weight *falls* with roughness while the implied noise
variance *rises*, each the reciprocal of the other. Note this is a
statement about $`\Sigma_{\text{col}} = A^{-1}`$, not about the
adjacency matrix you supplied: $`A`$ itself has its large eigenvalues on
the **smooth** directions. Inverting is what moves the variance to the
rough end.

- **Smoother as metric** $`\Rightarrow`$$`A`$ is large on smooth
  directions $`\Rightarrow`$$`\Sigma_{\text{col}} = A^{-1}`$ is large on
  the rough ones $`\Rightarrow`$ “the noise is high-frequency speckle;
  the signal is smooth.” This is denoising, and it is what you want when
  you are after spatially coherent maps.
- **Precision as metric** $`\Rightarrow`$$`A`$ is large on rough
  directions $`\Rightarrow`$$`\Sigma_{\text{col}}`$ is large on the
  smooth ones $`\Rightarrow`$ “a smooth field — drift, a scanner
  gradient, a global trend — is the nuisance.” Whitening it away lets
  fine-scale structure surface. This is ordinary generalized least
  squares against correlated noise.

One caveat on how literally to read this.
$`A = \Sigma_{\text{col}}^{-1}`$ is an interpretive frame, not a
constraint
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
enforces — the algorithm only ever whitens with $`A^{1/2}`$. The
covariance reading is how you should *choose* a metric, and it is what
makes “smoother $`\Rightarrow`$ denoising” more than a slogan.

So the question is never “adjacency or Laplacian?” in the abstract. It
is: **is the smooth thing my signal, or my nuisance?**

![The same graph, two metrics. Left: a smoother concentrates PC1 on the
smooth blob. Right: the Laplacian whitens the smooth field away, so PC1
locks onto the fine-scale checkerboard that was buried underneath
it.](gpca-metrics_files/figure-html/orientation-demo-1.png)

The same graph, two metrics. Left: a smoother concentrates PC1 on the
smooth blob. Right: the Laplacian whitens the smooth field away, so PC1
locks onto the fine-scale checkerboard that was buried underneath it.

Two practical notes on the graph matrices themselves. A raw adjacency
$`W`$ is **indefinite** (its eigenvalues sum to zero), so it is not a
valid metric on its own — shift it, as in $`I + \alpha W`$ with
$`\alpha`$ small enough to keep it PSD. A raw Laplacian $`L`$ is PSD but
**singular**: $`L\mathbf{1} = 0`$, so the spatially constant pattern has
zero length under it.
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
accepts that (the whitening uses a pseudo-inverse), but adding a small
ridge, $`L + \varepsilon I`$, is usually what you want.

Finally, the strength matters and is not cosmetic. In the example above
the first component only switched from the smooth blob to the
checkerboard once $`\alpha`$ passed roughly 50; at $`\alpha = 2`$ the
blob still dominated. Sweep it and look at your loadings.

## Recipes

Each recipe below is labelled with the direction it produces. All three
are written here as **precision** matrices, i.e. the noise-whitening
stance; drop the
[`solve()`](https://rdrr.io/pkg/Matrix/man/solve-methods.html) to get
the smoothing version instead.

### AR(1) row metric (whitens temporal autocorrelation)

Standard GLS treatment of serially correlated observations, as in fMRI
prewhitening: it removes temporal autocorrelation rather than imposing
temporal smoothness.

``` r

rho     <- 0.7
n_t     <- 60
idx     <- 0:(n_t - 1)
Sigma_r <- outer(idx, idx, function(i, j) rho^abs(i - j))
M_ar1   <- solve(Sigma_r + 1e-3 * diag(n_t))
```

### Spatial RBF kernel (whitens smooth spatial noise)

``` r

coords <- as.matrix(expand.grid(x = 1:8, y = 1:8))
d2     <- as.matrix(dist(coords))^2
ell    <- 2
K      <- exp(-d2 / (2 * ell^2))
A_rbf  <- solve(K + 1e-3 * diag(nrow(K)))   # precision: emphasises fine scale
# A_smooth <- K + 1e-3 * diag(nrow(K))      # kernel itself: smooth loadings
```

### Graph Laplacian (emphasises contrast across edges)

``` r

W <- bandSparse(30, k = c(-1, 0, 1),
                diagonals = list(rep(0.2, 29),
                                 rep(1, 30),
                                 rep(0.2, 29)))
D     <- Diagonal(x = rowSums(W))
A_lap <- (D - W) + 1e-2 * Diagonal(nrow(W))
# A_smooth <- solve(A_lap)                  # smoother: spatially coherent loadings
```

![Three structured metrics, all shown in the precision (noise-whitening)
orientation. Off-diagonal banding is what couples nearby rows or
variables -- it tells GPCA 'treat these dimensions as related, not
independent'.](gpca-metrics_files/figure-html/recipe-plots-1.png)

Three structured metrics, all shown in the precision (noise-whitening)
orientation. Off-diagonal banding is what couples nearby rows or
variables – it tells GPCA ‘treat these dimensions as related, not
independent’.

### A note on `sfpca()`

[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
takes the *opposite* input for the same intent. There the structure
enters as a constraint, $`v^{\top}(I + \alpha\Omega)v \le 1`$, which
charges rough $`v`$ against a fixed budget — so you pass the roughness
operator (Laplacian, second differences) directly via `alpha_v`, and
larger `alpha_v` means **smoother**. Metric form and constraint form are
inverse to one another: the same Laplacian smooths in
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md) and
roughens in
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md).

## Learning metrics with `gpca_mle()`

When you suspect structured noise but lack good priors,
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md)
learns SPD metrics by alternating GPCA with matrix-normal maximum
likelihood:

``` r

set.seed(1)
n_m <- 40; p_m <- 10
X_mle <- matrix(rnorm(n_m * p_m), n_m, p_m)
fit_mle <- gpca_mle(X_mle, ncomp = 2, max_iter = 6,
                    lambda = 1e-3, scale_fix = "trace",
                    method = "eigen", verbose = FALSE)
range(diag(as.matrix(fit_mle$M)))
#> [1] 156827.4 249196.5
range(diag(as.matrix(fit_mle$A)))
#> [1] 2.583950 3.330251
```

![Learned row metric M (left) and column metric A (right) on i.i.d.
Gaussian data. With identity-noise input the learner stays close to
identity, with mild row- and column-specific
damping.](gpca-metrics_files/figure-html/mle-plot-1.png)

Learned row metric M (left) and column metric A (right) on i.i.d.
Gaussian data. With identity-noise input the learner stays close to
identity, with mild row- and column-specific damping.

Keep `lambda` non-zero, start with a small `ncomp`, and watch warnings –
they signal a metric was repaired during the iteration.

## SPD requirements and remedies

### Why the requirement exists

GPCA measures variance in the inner products
$`\langle u, v\rangle_M = u^\top M
v`$ and $`\langle x, y\rangle_A = x^\top A y`$, and the solvers whiten
the data with the square roots $`M^{1/2}`$ and $`A^{1/2}`$. Both steps
need the metrics to be symmetric positive semi-definite (PSD). If a
metric has a negative eigenvalue, vectors in that direction have
negative squared length, the square root is not real, and “maximise
variance” no longer picks out anything meaningful.

Note that *semi*-definite is enough. **Singular metrics are perfectly
legal.** A graph Laplacian is rank-deficient by construction — it has a
zero eigenvalue on the constant vector — and
[`genpca()`](https://bbuchsbaum.github.io/genpca/reference/genpca.md)
takes it without complaint. A zero eigenvalue simply means that
direction is given no weight. Only *negative* eigenvalues are a problem.

### What happens when a metric falls short

Metrics usually fail the test for dull reasons rather than modelling
mistakes: a covariance estimated from fewer samples than variables, a
kernel matrix with round-off in the last few digits, a matrix that is
asymmetric by $`10^{-16}`$ because of the order in which it was
assembled. The check itself is tolerant — eigenvalues down to
`-tol * max(abs(diag()))` with `tol = 1e-6` count as non-negative — so
ordinary floating-point noise never triggers a remedy, and a metric that
passes is used exactly as supplied.

When a metric does fail, `constraints_remedy` decides what happens next:

| Value | What it does | What it costs |
|----|----|----|
| `"ridge"` (default) | Adds a diagonal shift (from the Gershgorin bound, with a [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html) fallback for small dense matrices) just large enough to make the matrix positive definite. Preserves sparsity. | The shift pulls the metric toward a multiple of the identity, diluting the structure you supplied. A *large* shift means the input was badly indefinite — diagnose it rather than absorb it. |
| `"clip"` | Eigendecomposes and sets the negative eigenvalues to zero, leaving the rest of the spectrum exactly as it was. | Densifies the matrix, so it refuses sparse input larger than 2000×2000. Use `"ridge"` at that size. |
| `"identity"` | Replaces the offending metric with the identity. | Discards your structure and quietly hands back ordinary PCA, with the run still reporting success. Pass `verbose = TRUE` to be told when this fires. |
| `"error"` | Refuses the input. | Nothing — this is the right setting when the metric comes from a pipeline that ought to be producing a valid one. [`genpca_cov()`](https://bbuchsbaum.github.io/genpca/reference/genpca_cov.md) defaults to it for that reason. |

Two behaviours are worth knowing about because they are silent:

- `"ridge"` and `"clip"` both symmetrise the input by mirroring the
  **upper** triangle onto the lower. If you pass a genuinely asymmetric
  matrix, the lower triangle is discarded without a word. Only `"error"`
  reports it.
- `"identity"` and a large `"ridge"` shift both produce a fit that looks
  healthy while encoding much less structure than you intended.

So if a metric matters to the result, verify it rather than trusting the
remedy to have done something sensible:

``` r

# Is the metric usable as-is? (small problems)
range(eigen(as.matrix(A), symmetric = TRUE, only.values = TRUE)$values)

# Did a remedy fire, and how hard?
fit <- genpca(X, A = A, M = M, ncomp = 3, verbose = TRUE)
```

Two habits that avoid the problem to begin with: scale metrics before
use (dividing by the mean diagonal, say) so that `M` and `A` are
comparably conditioned, and prefer building a metric that is PSD by
construction — a kernel, a Laplacian, an inverse-variance diagonal —
over one estimated and then patched. With
[`gpca_mle()`](https://bbuchsbaum.github.io/genpca/reference/gpca_mle.md),
warnings during the iteration mean a metric was repaired along the way
and the learned result should be inspected.

## Where next

See
[`vignette("structured-noise")`](https://bbuchsbaum.github.io/genpca/articles/structured-noise.md)
for how to choose the transfer function when several kinds of structure
are present at once, and
[`vignette("gpca-scale")`](https://bbuchsbaum.github.io/genpca/articles/gpca-scale.md)
for backend choices, sparse workflows, and covariance-only GPCA.
