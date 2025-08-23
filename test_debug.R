library(genpca)
library(Matrix)

set.seed(1)
n <- 20; p <- 10; k <- 1
X <- matrix(rnorm(n*p), n, p)
Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1
R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1
X_centered <- scale(X, center=TRUE, scale=FALSE)

# Call gmd_fast_cpp directly
res_cpp <- genpca:::gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)

# Call genpca with use_cpp=FALSE to use the R implementation
res_r_nocpp <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=FALSE, method="eigen", preproc=multivarious::pass())

# Call genpca with use_cpp=TRUE
res_r_cpp <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, method="spectra", preproc=multivarious::pass())

cat("Comparing eigenvalues:\n")
cat("gmd_fast_cpp d:", res_cpp$d, "\n")
cat("genpca (use_cpp=FALSE) sdev:", multivarious::sdev(res_r_nocpp), "\n")  
cat("genpca (use_cpp=TRUE) sdev:", multivarious::sdev(res_r_cpp), "\n\n")

# Check if the scores are consistent
scores_cpp_direct <- res_cpp$u * res_cpp$d  # Convert back to scores
scores_r_nocpp <- multivarious::scores(res_r_nocpp)
scores_r_cpp <- multivarious::scores(res_r_cpp)

cat("Score norms:\n")
cat("gmd_fast_cpp (scaled back):", sqrt(sum(scores_cpp_direct^2)), "\n")
cat("genpca (use_cpp=FALSE):", sqrt(sum(scores_r_nocpp^2)), "\n")
cat("genpca (use_cpp=TRUE):", sqrt(sum(scores_r_cpp^2)), "\n\n")

# Check correlation between different methods
normalize <- function(x) x / sqrt(sum(x^2))
cat("Correlations between normalized scores:\n")
cat("cpp_direct vs nocpp:", abs(sum(normalize(scores_cpp_direct) * normalize(scores_r_nocpp))), "\n")
cat("cpp_direct vs cpp:", abs(sum(normalize(scores_cpp_direct) * normalize(scores_r_cpp))), "\n")
cat("nocpp vs cpp:", abs(sum(normalize(scores_r_nocpp) * normalize(scores_r_cpp))), "\n")
