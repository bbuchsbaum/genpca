library(genpca)
library(Matrix)
library(testthat)

# Reproduce the test case
set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1 # Dense SPD
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1 # Dense SPD
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Get results from both methods
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)
res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())

cat("Testing k=1 case:\n")
cat("res_cpp$d:", res_cpp$d, "\n")
cat("sdev(res_r):", multivarious::sdev(res_r), "\n")
cat("Match?", all.equal(res_cpp$d, multivarious::sdev(res_r)), "\n\n")

# Get the scores
U1 <- res_cpp$u
U2 <- multivarious::scores(res_r)

cat("res_cpp$u dimensions:", dim(U1), "\n")
cat("scores(res_r) dimensions:", dim(U2), "\n")

cat("\nBefore normalization:\n")
cat("U1 column norm:", sqrt(sum(U1^2)), "\n")
cat("U2 column norm:", sqrt(sum(U2^2)), "\n")

# Now compare if they point to the same direction (up to sign)
cat("\nDirect comparison (should be same up to scale):\n")
scale_factor <- U2[1] / U1[1]
cat("Scale factor from first element:", scale_factor, "\n")
cat("All elements have same scale?", all(abs(U2/U1 - scale_factor) < 1e-10, na.rm=TRUE), "\n")

# Compare using the test's logic
U1_norm <- U1 / sqrt(sum(U1^2))
U2_norm <- U2 / sqrt(sum(U2^2))

corr <- abs(sum(U1_norm * U2_norm))
cat("\nCorrelation between normalized vectors:", corr, "\n")
cat("Close to 1?", abs(corr - 1) < 1e-5, "\n")
