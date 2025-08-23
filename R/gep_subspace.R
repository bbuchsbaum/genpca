

# Factor a matrix once with regularization
#' @keywords internal
#' @noRd
factor_mat <- function(M, reg = 1e-3, max_tries = 5) {
  d <- nrow(M)
  for (i in seq_len(max_tries)) {
    M_reg <- M + Diagonal(d, reg)
    ch <- try(Cholesky(M_reg, LDL = FALSE), silent = TRUE)
    if (!inherits(ch, "try-error")) {
      return(list(ch=ch, final_reg=reg))
    }
    reg <- reg*10
  }
  stop("Unable to factor matrix even after multiple attempts.")
}


#' @keywords internal
#' @noRd
orthonormalize <- function(X) {
  # thin QR; avoids allocating the full Q when d >> q
  QR <- qr(X)
  qr.Q(QR, complete = FALSE)
}


#' Fast sub-space solver for a small block of generalized eigen-pairs
#'
#' Uses pre-conditioned sub-space iteration on the operator
#' \eqn{S_2^{-1} S_1} (or its inverse) to obtain the `q`
#' largest or smallest generalized eigen-values/vectors of
#' \eqn{S_1 v = \lambda S_2 v}.
#'
#' @param S1,S2 Symmetric positive-(semi)definite `dgCMatrix`
#'   (or dense) matrices of the same dimension \eqn{d\times d}.
#' @param q Number of eigen-pairs required (`q << d`).
#' @param which `"largest"` or `"smallest"`.
#' @param max_iter,tol Stopping rule - iteration stops when
#'   `max(abs(lambda_new - lambda_old)/abs(lambda_old)) < tol`.
#' @param V0 Optional `d x q` initial block (will be orthonormalised).
#' @param seed Optional integer seed for reproducible random initialisation.
#' @param reg_S,reg_T Ridge terms added to `S1`/`S2` and the small
#'   `q x q` Gram matrix to guarantee invertibility.
#' @param verbose Logical - print convergence info.
#'
#' @return A list with components
#'   \describe{
#'     \item{values}{length-`q` numeric vector of Ritz eigen-values.}
#'     \item{vectors}{`d x q` matrix, columns are orthonormal eigen-vectors
#'                    in the *original* S-inner-product.}
#'   }
#' @keywords internal
solve_gep_subspace <- function(S1, S2, q = 2, which = c("largest", "smallest"),
                               max_iter = 100, tol = 1e-6, V0 = NULL, seed = NULL,
                               reg_S = 1e-3, reg_T = 1e-6, verbose = FALSE) {
  which <- match.arg(which)
  d <- nrow(S1)
  
  # Choose operator (avoid explicit inverse)
  # If largest: we iterate with S2^-1 S1 (i.e. solve S2 * V_hat = S1 * V)
  # If smallest: we iterate with S1^-1 S2 (i.e. solve S1 * V_hat = S2 * V)
  
  if (which == "largest") {
    # Factor S2 once
    s_fact <- factor_mat(S2, reg=reg_S)
    # S2^{-1} S1 V via triangular solves on the Cholesky of S2
    solve_step <- function(V) solve(s_fact$ch, S1 %*% V)
    final_reg_S <- s_fact$final_reg # Store the final reg used
  } else {
    # smallest
    s_fact <- factor_mat(S1, reg=reg_S)
    # S1^{-1} S2 V via triangular solves on the Cholesky of S1
    solve_step <- function(V) solve(s_fact$ch, S2 %*% V)
    final_reg_S <- s_fact$final_reg # Store the final reg used
  }
  
  if (is.null(V0)) {
    if (!is.null(seed)) set.seed(seed)
    V <- matrix(rnorm(d*q), d, q)
  } else {
    V <- V0
  }
  
  V <- orthonormalize(V)
  
  lambda_old <- NULL
  
  for (iter in seq_len(max_iter)) {
    Y <- solve_step(V)
    
    # Orthonormalize Y
    Y <- orthonormalize(Y)
    
    # Form S and T efficiently (avoid d x d intermediate products)
    S_mat <- t(Y) %*% (S1 %*% Y)  # q x q
    T_mat <- t(Y) %*% (S2 %*% Y)  # q x q (always S2 here)
    
    # Symmetrize numerically
    S_mat <- Matrix::forceSymmetric(S_mat)
    T_mat <- Matrix::forceSymmetric(T_mat)
    
    # Robust Cholesky-whitening of T
    reg <- reg_T
    R <- NULL
    for (tries in 0:5) {
      T_reg <- as.matrix(T_mat) + diag(reg, ncol(T_mat))
      ok <- TRUE
      R <- try(chol(T_reg), silent = TRUE)
      if (inherits(R, "try-error")) { ok <- FALSE }
      if (ok) break
      reg <- reg * 10
    }
    if (!is.numeric(R)) stop("Unable to chol() the T matrix even after regularization.")
    
    # C = R^{-T} S R^{-1} (symmetric)
    C <- backsolve(R, t(backsolve(R, t(as.matrix(S_mat)), transpose = TRUE)))
    C <- (C + t(C))/2  # Symmetrize
    
    # Eigen decomposition of whitened matrix
    E <- eigen(C, symmetric = TRUE)
    
    # Select eigenvalues/vectors based on which
    ord <- if (which == "largest") {
      order(E$values, decreasing = TRUE)
    } else {
      order(E$values, decreasing = FALSE)
    }
    ord <- ord[seq_len(q)]
    
    lambda <- E$values[ord]
    U <- E$vectors[, ord, drop = FALSE]
    
    # Transform back: W = R^{-1} U
    W <- backsolve(R, U)
    
    V_new <- Y %*% W
    
    if (!is.null(lambda_old)) {
      rel_change <- max(abs(lambda - lambda_old) / pmax(abs(lambda), 1e-12))
      if (rel_change < tol) {
        V <- V_new
        if (verbose) message("Converged in ", iter, " iterations.")
        break
      }
    }
    
    V <- V_new
    lambda_old <- lambda
  }
  
  if (iter == max_iter) {
    warning("Reached max_iter without convergence (delta_lambda > tol)")
  }
  
  # Optional but recommended: enforce S2-orthonormality again
  G <- crossprod(V, S2 %*% V)
  G <- (G + t(G))/2  # Symmetrize
  Rg <- chol(as.matrix(G) + diag(reg_T, q))
  V <- V %*% solve(Rg)
  
  # Recompute lambda as Rayleigh quotients in the now S2-orthonormal basis
  lambda <- diag(crossprod(V, S1 %*% V))
  
  list(values = as.numeric(lambda), vectors = V, final_reg_S = final_reg_S)
}

##############################################
# Example usage:
# d <- 50
# q <- 3
# 
# set.seed(1)
# A <- matrix(rnorm(d*d), d, d)
# S2 <- crossprod(A) + diag(d)*0.1
# B <- matrix(rnorm(d*d), d, d)
# S1 <- crossprod(B) + diag(d)*0.1
# 
# res_largest <- solve_gep_subspace(S1, S2, q = q, which = "largest", reg_S = 1e-4, reg_T = 1e-6)
# cat("Largest 3 eigenvalues:\n", res_largest$values, "\n")
# 
# res_smallest <- solve_gep_subspace(S1, S2, q = q, which = "smallest", reg_S = 1e-3, reg_T = 1e-6)
# cat("Smallest 3 eigenvalues:\n", res_smallest$values, "\n")
