library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- diag(n)  # Use identity for simplicity
R <- diag(p)  # Use identity for simplicity
X_centered <- scale(X, center=TRUE, scale=FALSE)

# With identity matrices, this should be standard SVD
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

# Standard SVD for comparison
svd_res <- svd(X_centered, nu=k, nv=k)

cat("With identity constraint matrices:\n")
cat("gmd_fast_cpp d:", res_cpp$d, "\n")
cat("svd d:", svd_res$d[1], "\n")
cat("Match?", abs(res_cpp$d - svd_res$d[1]) < 1e-10, "\n\n")

# Check eigenvector alignment
u_corr <- abs(cor(res_cpp$u, svd_res$u))
v_corr <- abs(cor(res_cpp$v, svd_res$v))

cat("U correlation:", u_corr, "\n")
cat("V correlation:", v_corr, "\n")
