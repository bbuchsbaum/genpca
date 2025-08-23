library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1
X_centered <- scale(X, center=TRUE, scale=FALSE)

# What problem does each method solve?
# gmd_fast_cpp solves: (X'*M*X*A)v = d^2*v and returns u = X*A*v/d
# eigen method solves the same in theory but through different approach

# Let's verify the eigenvalue problem
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

# Manually verify the eigenvalue problem
XtMX <- t(X_centered) %*% Q %*% X_centered
XtMXA <- XtMX %*% R

# Check if v satisfies (X'MXA)v = d^2*v
v1 <- res_cpp$v[,1]
d1 <- res_cpp$d[1]

lhs <- XtMXA %*% v1
rhs <- d1^2 * v1

cat("Checking eigenvalue equation for first component:\n")
cat("||XtMXA*v - d^2*v||:", sqrt(sum((lhs - rhs)^2)), "\n")
cat("d^2:", d1^2, "\n\n")

# Check the relationship u = X*A*v/d
u1_computed <- (X_centered %*% R %*% v1) / d1
u1_returned <- res_cpp$u[,1]

cat("Checking u = X*A*v/d:\n")
cat("||u_computed - u_returned||:", sqrt(sum((u1_computed - u1_returned)^2)), "\n")
