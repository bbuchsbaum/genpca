set.seed(1)

n <- 6
p <- 4

# (1) Generate data X (n x p)
X <- matrix(rnorm(n * p), nrow = n)

# (2) Generate PD row constraint M (n x n)
tmpM <- crossprod(matrix(rnorm(n * n), n, n))
M    <- tmpM + diag(n)*0.5  # ensure positive definiteness

# (3) Generate PD column constraint A (p x p)
tmpA <- crossprod(matrix(rnorm(p * p), p, p))
A    <- tmpA + diag(p)*0.5

# Load your GMD code or library, then:
res_defl <- genpca(
  X,
  A = A,
  M = M,
  ncomp     = 2,
  deflation = FALSE,   # iterative deflation
  threshold = 1e-9,
  geneig    = FALSE,  # specifically NOT direct generalized eigen
  use_cpp   = FALSE
)

# Load your GMD code or library, then:
res_g <- genpca(
  X,
  A = Matrix(A),
  M = Matrix(M),
  ncomp     = 2,
  deflation = FALSE,   # iterative deflation
  threshold = 1e-9,
  geneig    = TRUE,  # specifically NOT direct generalized eigen
  use_cpp   = FALSE
)

cat("=== Deflation approach returned ===\n")
cat("res_defl$d is:\n"); print(res_defl$d)          # Might be NULL
cat("res_defl$v (loadings):\n"); print(res_defl$v)
cat("res_defl$u (scores):\n");   print(res_defl$u)

# ------------------------------------------------------------------
#  A) Manually compute deflation 'singular values' from (u_i, v_i)
# ------------------------------------------------------------------
# In the deflation code, the singular value for component i is typically:
#     d[i] = u_i^T * M * X * (A * v_i)
#  or equivalently:
#     d[i] = (u_i^T M) X (A v_i)
#
# We'll compute that for each extracted component i.

my_d <- numeric(2)
for (i in seq_len(2)) {
  u_i <- res_defl$u[, i, drop=FALSE]  # (n x 1)
  v_i <- res_defl$v[, i, drop=FALSE]  # (p x 1)
  # The GMD "singular value" is the scalar:
  #   d[i] = (u_i^T M) X (A v_i)
  num_i <- t(u_i) %*% (M %*% (X %*% (A %*% v_i)))  # 1x1
  my_d[i] <- as.numeric(num_i)
}

cat("\nManually computed deflation singular values:\n")
print(my_d)

# ------------------------------------------------------------------
#  B) Solve generalized eigenproblem:
#     (X' M X) v = lambda A v
# ------------------------------------------------------------------
Mx <- crossprod(X, M) %*% X  # p x p
A_inv_Mx <- solve(A, Mx)     # p x p

res_eig <- eigen(A_inv_Mx, symmetric=FALSE)
lambda_vals <- Re(res_eig$values)   # largest are the top GMD
v_eigs      <- Re(res_eig$vectors)

# We'll look at the top 2
lambda_top2 <- lambda_vals[1:2]
v_top2      <- v_eigs[,1:2, drop=FALSE]
d_top2      <- sqrt(lambda_top2)

cat("\n=== Manual generalized eigen ===\n")
cat("Top 2 lambda:\n"); print(lambda_top2)
cat("sqrt(lambda):\n"); print(d_top2)
cat("Eigenvectors (top 2):\n"); print(v_top2)

# ------------------------------------------------------------------
#  C) Compare
# ------------------------------------------------------------------
df_comp <- data.frame(
  d_deflation_manual = my_d,
  d_eig_sqrt         = d_top2
)
cat("\nCompare manually computed deflation singular values vs sqrt(lambda):\n")
print(df_comp)

# Also check correlation for the top loadings:
v_defl_1 <- res_defl$v[,1]; v_defl_2 <- res_defl$v[,2]
v_eig_1  <- v_top2[,1];      v_eig_2  <- v_top2[,2]

cat("\nCorrelation of 1st loadings:", cor(v_defl_1, v_eig_1), "\n")
cat("Correlation of 2nd loadings:", cor(v_defl_2, v_eig_2), "\n")