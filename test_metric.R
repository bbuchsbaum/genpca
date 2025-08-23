library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1
X_centered <- scale(X, center=TRUE, scale=FALSE)

res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

# Check orthonormality in the Q-metric (for U) and R-metric (for V)
u1 <- res_cpp$u[,1, drop=FALSE]
v1 <- res_cpp$v[,1, drop=FALSE]

# U should be orthonormal in Q-metric: u'*Q*u = 1
uQu <- t(u1) %*% Q %*% u1
cat("u'*Q*u (should be 1 if Q-orthonormal):", uQu[1,1], "\n")

# V should be orthonormal in R-metric: v'*R*v = 1  
vRv <- t(v1) %*% R %*% v1
cat("v'*R*v (should be 1 if R-orthonormal):", vRv[1,1], "\n")

# Standard norms
cat("\nStandard Euclidean norms:\n")
cat("||u||:", sqrt(sum(u1^2)), "\n")
cat("||v||:", sqrt(sum(v1^2)), "\n")
