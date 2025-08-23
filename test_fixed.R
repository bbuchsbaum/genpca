library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Get results from both methods
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)
res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())

cat("Testing the fix:\n")
cat("res_cpp$d:", res_cpp$d, "\n")
cat("sdev(res_r):", multivarious::sdev(res_r), "\n")
cat("Match?", all.equal(res_cpp$d, multivarious::sdev(res_r), check.attributes=FALSE), "\n\n")

# Get the scores - now they should match!
scores_cpp <- res_cpp$u * res_cpp$d  # res_cpp$u is scores/d, multiply by d to get scores
scores_r <- multivarious::scores(res_r)

cat("Score comparison:\n")
cat("scores_cpp norm:", sqrt(sum(scores_cpp^2)), "\n")
cat("scores_r norm:", sqrt(sum(scores_r^2)), "\n")

# Normalize and check correlation
normalize <- function(x) x / sqrt(sum(x^2))
corr <- abs(sum(normalize(scores_cpp) * normalize(scores_r)))
cat("Correlation:", corr, "\n")
cat("Close to 1?", abs(corr - 1) < 1e-5, "\n")
