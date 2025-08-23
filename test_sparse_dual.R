library(genpca)
library(Matrix)

set.seed(4)
n <- 15
p <- 25
k <- 4
X <- matrix(rnorm(n*p), n, p)
Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
R <- Diagonal(p, x = runif(p, 0.5, 1.5))

X_centered <- scale(X, center=TRUE, scale=FALSE)

res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)

cat("Test p > n, sparse constraints:\n")
cat("n =", n, ", p =", p, ", k =", k, "\n")
cat("res_cpp$d:", res_cpp$d, "\n")
cat("sdev(res_r):", multivarious::sdev(res_r), "\n")
cat("d match?", all.equal(res_cpp$d, multivarious::sdev(res_r), check.attributes=FALSE), "\n\n")

# Check components
V1 <- res_cpp$v
V2 <- multivarious::components(res_r)

cat("Component dimensions:\n")
cat("res_cpp$v:", dim(V1), "\n")
cat("components(res_r):", dim(V2), "\n\n")

# Normalize and compute correlation matrix
V1_norm <- apply(V1, 2, function(x) x / sqrt(sum(x^2)))
V2_norm <- apply(V2, 2, function(x) x / sqrt(sum(x^2)))

corr_mat <- abs(t(V1_norm) %*% V2_norm)

cat("Correlation matrix between normalized components:\n")
print(round(corr_mat, 4))

cat("\nSum of squared correlations:", sum(corr_mat^2), "\n")
cat("Should be close to k =", k, "\n")
cat("Difference from k:", abs(sum(corr_mat^2) - k), "\n")

# Check diagonal elements
cat("\nDiagonal of correlation matrix (should be close to 1):\n")
cat(diag(corr_mat), "\n")

# Check if it's a permutation issue
cat("\nMax correlation per row:", apply(corr_mat, 1, max), "\n")
cat("Max correlation per col:", apply(corr_mat, 2, max), "\n")
