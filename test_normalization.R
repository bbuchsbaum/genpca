library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- diag(n)
R <- diag(p)
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Call the function
res <- genpca:::gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)

cat("Testing normalization:\n")
cat("res$d:", res$d, "\n")
cat("length(res$d):", length(res$d), "\n")
cat("Column norms of res$u:\n")
for(i in 1:ncol(res$u)) {
  cat("  Column", i, ":", sqrt(sum(res$u[,i]^2)), "\n")
}
