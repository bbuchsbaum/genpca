# Implement the correct generalized PLS formulation
# Based on expert feedback

library(Matrix)
library(RSpectra)

# Test case with non-identity metrics
set.seed(123)
n <- 10; p <- 5; q <- 4
X <- matrix(rnorm(n * p), n, p)
Y <- matrix(rnorm(n * q), n, q)

# Create some non-identity metrics
Mx <- diag(runif(n, 0.5, 2))
My <- diag(runif(n, 0.5, 2))
Ax <- diag(runif(p, 0.5, 2))
Ay <- diag(runif(q, 0.5, 2))

# Step 1: Compute metric square roots
Mx_half <- diag(sqrt(diag(Mx)))
My_half <- diag(sqrt(diag(My)))

# Step 2: Transform data
X_tilde <- Mx_half %*% X  # n x p
Y_tilde <- My_half %*% Y  # n x q

# Step 3: Build K_X operator
Ax_reg <- Ax + diag(p) * 1e-8
inv_Ax <- solve(Ax_reg)
K_X <- X_tilde %*% inv_Ax %*% t(X_tilde)  # n x n

# Step 4: Compute K_X^{1/2} via eigendecomposition
eig_KX <- eigen(K_X, symmetric = TRUE)
pos_idx <- eig_KX$values > 1e-10
KX_half <- eig_KX$vectors[, pos_idx] %*% 
           diag(sqrt(eig_KX$values[pos_idx])) %*% 
           t(eig_KX$vectors[, pos_idx])
KX_half_inv <- eig_KX$vectors[, pos_idx] %*% 
               diag(1/sqrt(eig_KX$values[pos_idx])) %*% 
               t(eig_KX$vectors[, pos_idx])

# Step 5: Build K_Y operator (for testing, compute explicitly)
Ay_reg <- Ay + diag(q) * 1e-8
inv_Ay <- solve(Ay_reg)
K_Y <- Y_tilde %*% inv_Ay %*% t(Y_tilde)  # n x n

# Step 6: The key insight - we need the correct coupling
# For generalized PLS, the coupling is through the cross-covariance
# The operator should be symmetric in the K_X^{1/2} space

# Build the symmetric operator
# A_sym = K_X^{1/2} * K_Y * K_X^{1/2}
A_sym <- KX_half %*% K_Y %*% KX_half

# Get eigenvalues/vectors
eig_sym <- eigen(A_sym, symmetric = TRUE)
k <- 2
lambda <- eig_sym$values[1:k]
u_sym <- eig_sym$vectors[, 1:k]

cat("Eigenvalues of symmetric operator:", lambda, "\n")
cat("Square roots:", sqrt(lambda), "\n\n")

# Check orthogonality
cat("Orthogonality of symmetric eigenvectors:\n")
print(round(t(u_sym) %*% u_sym, 6))

# Transform back to dual space
U_dual <- KX_half_inv %*% u_sym
S_dual <- K_Y %*% U_dual

# Check the key relationship
cat("\nChecking t(U_dual) * S_dual:\n")
print(round(t(U_dual) %*% S_dual, 6))

# The diagonal should be the eigenvalues
cat("\nDiagonal values:", diag(t(U_dual) %*% S_dual), "\n")
cat("Expected (eigenvalues):", lambda, "\n")

# Transform to primal scores
# Key insight: we need to scale by sqrt(lambda) to get the right normalization
scores_x <- solve(Mx_half, U_dual %*% diag(sqrt(lambda)))
scores_y <- solve(My_half, S_dual %*% diag(1/sqrt(lambda)))

# Check cross-block reconstruction
cross_cov <- t(scores_x) %*% My %*% scores_y
cat("\nCross-block covariance:\n")
print(round(cross_cov, 6))
cat("\nDiagonal:", diag(cross_cov), "\n")
cat("Expected (singular values):", sqrt(lambda), "\n")

# Compute loadings
load_x <- inv_Ax %*% t(X) %*% My %*% scores_y %*% diag(1/sqrt(lambda))
load_y <- inv_Ay %*% t(Y) %*% Mx %*% scores_x %*% diag(1/sqrt(lambda))

# Check Omega-orthogonality
cat("\nChecking Ax-orthogonality of loadings:\n")
print(round(t(load_x) %*% Ax %*% load_x, 6))

cat("\nChecking Ay-orthogonality of loadings:\n")
print(round(t(load_y) %*% Ay %*% load_y, 6))