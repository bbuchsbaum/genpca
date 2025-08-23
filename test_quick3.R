library(genpca)
set.seed(123)
n <- 20; p <- 10; k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- diag(n)
R <- diag(p)

# Center the data
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Call the function
res <- genpca:::gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)

cat("u'u:\n")
print(round(crossprod(res$u), 3))

cat("\nv'v:\n")
print(round(crossprod(res$v), 3))

cat("\nCheck if u'u = diag(d²):\n")
cat("diag(u'u):", diag(crossprod(res$u)), "\n")
cat("d²:", res$d^2, "\n")

cat("\nSo the scores are scaled by d (not normalized).\n")
