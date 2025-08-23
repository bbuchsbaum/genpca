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

cat("\nShould u be scaled by d? Check u'u/d²:\n")
scaled_u <- res$u %*% diag(1/res$d)
print(round(crossprod(scaled_u), 3))
