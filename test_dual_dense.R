library(genpca)
library(Matrix)

# Test dual path with DENSE constraints
set.seed(2)
n <- 15
p <- 20
k <- 5
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1 # Dense SPD
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1 # Dense SPD

X_centered <- scale(X, center=TRUE, scale=FALSE)

res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

V1 <- res_cpp$v
V2 <- multivarious::components(res_r)

# Normalize and compute correlation matrix
V1_norm <- apply(V1, 2, function(x) x / sqrt(sum(x^2)))
V2_norm <- apply(V2, 2, function(x) x / sqrt(sum(x^2)))

corr_mat <- abs(t(V1_norm) %*% V2_norm)

cat("Test p > n, DENSE constraints:\n")
cat("Sum of squared correlations:", sum(corr_mat^2), "\n")
cat("Should be close to k =", k, "\n")
cat("Difference from k:", abs(sum(corr_mat^2) - k), "\n")
cat("Pass with tol=1e-5?", abs(sum(corr_mat^2) - k) < 1e-5 * k, "\n\n")

# Test primal path with SPARSE constraints
set.seed(1)
n <- 20
p <- 10
k <- 3
X <- matrix(rnorm(n*p), n, p)
Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
R <- Diagonal(p, x = runif(p, 0.5, 1.5))

X_centered <- scale(X, center=TRUE, scale=FALSE)

res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)

V1 <- res_cpp$v
V2 <- multivarious::components(res_r)

# Normalize and compute correlation matrix
V1_norm <- apply(V1, 2, function(x) x / sqrt(sum(x^2)))
V2_norm <- apply(V2, 2, function(x) x / sqrt(sum(x^2)))

corr_mat <- abs(t(V1_norm) %*% V2_norm)

cat("Test p <= n, SPARSE constraints:\n")
cat("Sum of squared correlations:", sum(corr_mat^2), "\n")
cat("Should be close to k =", k, "\n")
cat("Difference from k:", abs(sum(corr_mat^2) - k), "\n")
cat("Pass with tol=1e-5?", abs(sum(corr_mat^2) - k) < 1e-5 * k, "\n")
