library(genpca)
library(Matrix)
X <- diag(c(1, 1e6)); W <- diag(c(1, 1e-10))
for (side in c('rows', 'columns')) for (seed in c(1234L, 1:12)) {
 M <- if (side == 'rows') W else diag(2)
 A <- if (side == 'columns') W else diag(2)
 f <- genpca(X, M=M, A=A, ncomp=2, method='randomized',
             jitter_metric=0, seed_randomized=seed)
 cat(side, seed, 'd:', format(f$sdev, digits=8),
     'reconstruction error:', norm(as.matrix(multivarious::reconstruct(f))-X,'F')/norm(X,'F'), '\n')
}
set.seed(2); X <- matrix(rnorm(180), 30, 6)
for (s in c(1,100)) {
 r <- gpca_mle(s*X, ncomp=2, max_iter=8)
 cat('MLE scale', s, 'loglik', r$loglik, 'rescale delta', r$loglik_rescale_delta, 'refit delta', r$loglik_refit_delta, '\n')
}
sessionInfo()
