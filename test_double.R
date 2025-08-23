library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Get raw result from C++
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

# What genpca does with it
# From gpca.R line 383-385:
# M_ou <- M %*% svdfit$u
# scores_mat <- sweep(M_ou, 2, svdfit$d, `*`)

M_ou <- Q %*% res_cpp$u
scores_computed <- sweep(M_ou, 2, res_cpp$d, `*`)

# Get actual scores from genpca
res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())
scores_actual <- multivarious::scores(res_r)

cat("Score norms:\n")
cat("Computed from res_cpp:", sqrt(sum(scores_computed^2)), "\n")
cat("Actual from genpca:", sqrt(sum(scores_actual^2)), "\n")
cat("Ratio:", sqrt(sum(scores_actual^2)) / sqrt(sum(scores_computed^2)), "\n\n")

# The issue: res_cpp$u is now orthonormal (after my fix)
# But the C++ code actually returns scores = M*U*D, not orthonormal U!
# My normalization is breaking it!

cat("res_cpp$u column norm (should be 1 if orthonormal):", sqrt(sum(res_cpp$u^2)), "\n")
