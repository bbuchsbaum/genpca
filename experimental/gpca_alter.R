# # library(Matrix)
# # library(PRIMME)  # Assuming 'PRIMME' is the package containing 'eigs_sym'
# # 
# # gmdLA2 <- function(X, Q, R, k = min(n, p), n, p) {
# #   # Ensure Q and R are symmetric and sparse
# #   Q <- as(Q, "symmetricMatrix")
# #   R <- as(R, "symmetricMatrix")
# #   
# #   # Compute M = X' * Q * X
# #   M <- crossprod(X, Q %*% X)
# #   
# #   # Define functions for matrix-vector multiplication
# #   M_mult <- function(x) as.vector(M %*% x)
# #   R_mult <- function(x) as.vector(R %*% x)
# #   
# #   # Use eigs_sym from PRIMME to solve the generalized eigenvalue problem
# #   eig_result <- eigs_sym(A = M_mult, B = R_mult, NEig = k, which = "LA", n = p)
# #   
# #   # Extract eigenvalues and eigenvectors
# #   lambda <- eig_result$values
# #   V <- eig_result$vectors
# #   
# #   # Compute singular values (GMD values)
# #   D <- diag(sqrt(abs(lambda)), nrow = k, ncol = k)
# #   
# #   # Compute U
# #   D_inv <- diag(1 / sqrt(abs(lambda)), nrow = k, ncol = k)
# #   U <- X %*% (R %*% V %*% D_inv)
# #   
# #   # Normalize U and V with respect to Q and R
# #   for (i in 1:k) {
# #     # Normalize U with respect to Q
# #     norm_U <- sqrt(as.numeric(crossprod(U[, i], Q %*% U[, i])))
# #     if (norm_U != 0) {
# #       U[, i] <- U[, i] / norm_U
# #     }
# #     
# #     # Normalize V with respect to R
# #     norm_V <- sqrt(as.numeric(crossprod(V[, i], R %*% V[, i])))
# #     if (norm_V != 0) {
# #       V[, i] <- V[, i] / norm_V
# #     }
# #   }
# #   
# #   # Return results
# #   list(u = U, v = V, d = diag(D), k = k)
# # }
# 
# # gmdLA_hybrid <- function(X, Q, R, k = min(n, p), n, p) {
# #   # Determine sparsity and size
# #   is_sparse_Q <- inherits(Q, "sparseMatrix")
# #   is_sparse_R <- inherits(R, "sparseMatrix")
# #   size_threshold <- 1e6  # Adjust based on system memory
# #   
# #   # Estimate the number of non-zero elements
# #   nnz_Q <- length(Q@x)
# #   nnz_R <- length(R@x)
# #   
# #   # Decide whether to use function-based approach
# #   use_function_interface <- (nnz_Q + nnz_R) > size_threshold
# #   
# #   if (use_function_interface) {
# #     # Function-based approach
# #     Q <- as(Q, "symmetricMatrix")
# #     R <- as(R, "symmetricMatrix")
# #     
# #     # Compute M = X' * Q * X
# #     M <- crossprod(X, Q %*% X)
# #     
# #     # Define functions for matrix-vector multiplication
# #     M_mult <- function(x) as.vector(M %*% x)
# #     R_mult <- function(x) as.vector(R %*% x)
# #     
# #     # Solve generalized eigenvalue problem using PRIMME
# #     eig_result <- eigs_sym(A = M_mult, B = R_mult, NEig = k, which = "LA", n = p)
# #     
# #     # Extract eigenvalues and eigenvectors
# #     lambda <- eig_result$values
# #     V <- eig_result$vectors
# #     
# #     # Compute singular values
# #     D <- diag(sqrt(abs(lambda)), nrow = k, ncol = k)
# #     
# #     # Compute U
# #     D_inv <- diag(1 / sqrt(abs(lambda)), nrow = k, ncol = k)
# #     U <- X %*% (R %*% V %*% D_inv)
# #     
# #     # Normalize U and V
# #     for (i in 1:k) {
# #       # Normalize U with respect to Q
# #       norm_U <- sqrt(as.numeric(t(U[, i]) %*% Q %*% U[, i]))
# #       if (norm_U != 0) {
# #         U[, i] <- U[, i] / norm_U
# #       }
# #       
# #       # Normalize V with respect to R
# #       norm_V <- sqrt(as.numeric(t(V[, i]) %*% R %*% V[, i]))
# #       if (norm_V != 0) {
# #         V[, i] <- V[, i] / norm_V
# #       }
# #     }
# #     
# #     # Return GMD factors
# #     list(u = U, v = V, d = diag(D), k = k)
# #     
# #   } else {
# #     # Direct matrix passing approach
# #     Q <- as.matrix(Q)
# #     R <- as.matrix(R)
# #     
# #     # Compute M = X' * Q * X
# #     M <- crossprod(X, Q %*% X)
# #     
# #     # Solve the generalized eigenvalue problem M v = lambda R v
# #     eig_result <- eigen(solve(R) %*% M, symmetric = FALSE)
# #     
# #     # Select the top k eigenvalues and eigenvectors
# #     idx <- order(Re(eig_result$values), decreasing = TRUE)[1:k]
# #     lambda <- Re(eig_result$values[idx])
# #     V <- eig_result$vectors[, idx]
# #     
# #     # Compute singular values
# #     D <- diag(sqrt(abs(lambda)), nrow = k, ncol = k)
# #     
# #     # Compute U
# #     D_inv <- diag(1 / sqrt(abs(lambda)), nrow = k, ncol = k)
# #     U <- X %*% (R %*% V %*% D_inv)
# #     
# #     # Normalize U and V with respect to Q and R
# #     for (i in 1:k) {
# #       # Normalize U with respect to Q
# #       norm_U <- sqrt(as.numeric(t(U[, i]) %*% Q %*% U[, i]))
# #       if (norm_U != 0) {
# #         U[, i] <- U[, i] / norm_U
# #       }
# #       
# #       # Normalize V with respect to R
# #       norm_V <- sqrt(as.numeric(t(V[, i]) %*% R %*% V[, i]))
# #       if (norm_V != 0) {
# #         V[, i] <- V[, i] / norm_V
# #       }
# #     }
# #     
# #     # Return GMD factors
# #     list(u = U, v = V, d = diag(D), k = k)
# #   }
# # }
