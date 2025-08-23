library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- diag(n)
R <- diag(p)
X_centered <- scale(X, center=TRUE, scale=FALSE)

res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)
res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())

cat("With identity constraints (should be standard SVD):\n\n")

cat("What res_cpp returns:\n")
cat("u dimensions:", dim(res_cpp$u), "\n")
cat("v dimensions:", dim(res_cpp$v), "\n")
cat("u column norms:", apply(res_cpp$u, 2, function(x) sqrt(sum(x^2))), "\n")
cat("v column norms:", apply(res_cpp$v, 2, function(x) sqrt(sum(x^2))), "\n\n")

cat("What genpca stores internally:\n")
cat("ou (should be orthonormal U):\n")
ou <- res_r$ou
cat("  dimensions:", dim(ou), "\n")
cat("  column norms:", apply(ou, 2, function(x) sqrt(sum(x^2))), "\n\n")

cat("ov (should be orthonormal V):\n")
ov <- res_r$ov
cat("  dimensions:", dim(ov), "\n")
cat("  column norms:", apply(ov, 2, function(x) sqrt(sum(x^2))), "\n\n")

cat("What multivarious returns:\n")
cat("scores dimensions:", dim(multivarious::scores(res_r)), "\n")
cat("components dimensions:", dim(multivarious::components(res_r)), "\n\n")

cat("Direct comparison:\n")
cat("res_cpp$u == ou?", all.equal(res_cpp$u, ou), "\n")
cat("res_cpp$v == ov?", all.equal(res_cpp$v, ov), "\n")
