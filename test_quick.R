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

cat("Results structure:\n")
str(res)

cat("\nDimensions:\n")
cat("u:", dim(res$u), "\n")
cat("v:", dim(res$v), "\n")
cat("d:", length(res$d), "\n")

cat("\nNames:\n")
cat("colnames(u):", colnames(res$u), "\n")
cat("colnames(v):", colnames(res$v), "\n")
cat("names(d):", names(res$d), "\n")

# Check that scores and components are orthonormal
cat("\nOrthonormality checks:\n")
cat("u'u diagonal?:", all(abs(diag(crossprod(res$u)) - 1) < 1e-10), "\n")
cat("v'v diagonal?:", all(abs(diag(crossprod(res$v)) - 1) < 1e-10), "\n")
