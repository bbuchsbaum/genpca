library(genpca)
library(testthat)

# Test the normalization fix
set.seed(123)
n <- 20; p <- 10; k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- diag(n)
R <- diag(p)
X_centered <- scale(X, center=TRUE, scale=FALSE)

res <- genpca:::gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)

cat("After normalization fix:\n")
cat("res$u column norms:", apply(res$u, 2, function(x) sqrt(sum(x^2))), "\n")
cat("res$v column norms:", apply(res$v, 2, function(x) sqrt(sum(x^2))), "\n")
cat("res$d:", res$d, "\n")

# Test with genpca
res_gpca <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())
cat("\ngenpca scores column norms:", apply(multivarious::scores(res_gpca), 2, function(x) sqrt(sum(x^2))), "\n")
cat("genpca sdev:", multivarious::sdev(res_gpca), "\n")

# Check if they match
cat("\nDo sdev match?", all.equal(res$d, multivarious::sdev(res_gpca)), "\n")
