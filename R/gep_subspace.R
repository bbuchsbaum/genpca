

# Factor a matrix once with regularization
#' @keywords internal
#' @noRd
factor_mat <- function(M, reg = 1e-3, max_tries = 5) {
  d <- nrow(M)
  M_reg <- M
  for (i in seq_len(max_tries)) {
    M_reg <- M_reg + Diagonal(d, reg)
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
#'   `max(abs(λ_new - λ_old)/abs(λ_old)) < tol`.
#' @param V0 Optional `d×q` initial block (will be orthonormalised).
#' @param seed Optional integer seed for reproducible random initialisation.
#' @param reg_S,reg_T Ridge terms added to `S1`/`S2` and the small
#'   `q×q` Gram matrix to guarantee invertibility.
#' @param verbose Logical – print convergence info.
#'
#' @return A list with components
#'   \describe{
#'     \item{values}{length-`q` numeric vector of Ritz eigen-values.}
#'     \item{vectors}{`d×q` matrix, columns are orthonormal eigen-vectors
#'                    in the *original* S-inner-product.}
#'   }
#' @keywords internal
#' @examples
#' d <- 100; q <- 4
#' S2 <- crossprod(matrix(rnorm(d*d), d, d)) + Diagonal(d)*0.1
#' S1 <- crossprod(matrix(rnorm(d*d), d, d)) + Diagonal(d)*0.1
#' solve_gep_subspace(S1, S2, q)
solve_gep_subspace <- function(S1, S2, q = 2, which = c("largest", "smallest"),
                               max_iter = 100, tol = 1e-6, V0 = NULL, seed = NULL,
                               reg_S = 1e-3, reg_T = 1e-6, verbose = FALSE) {
  which <- match.arg(which)
  d <- nrow(S1)
  
  # Depending on which eigenvalues we want:
  # If largest: we iterate with S2^-1 S1 (i.e. solve S2 * V_hat = S1 * V)
  # If smallest: we iterate with S1^-1 S2 (i.e. solve S1 * V_hat = S2 * V)
  
  if (which == "largest") {
    # Factor S2 once
    s_fact <- factor_mat(S2, reg=reg_S)
    # Pre-compute inverse factor L^-1 where S2 ~ LL'
    Linv <- solve(s_fact$ch, Diagonal(d))
    solve_step <- function(V) Linv %*% (S1 %*% V)
    final_reg_S <- s_fact$final_reg # Store the final reg used
  } else {
    # smallest
    s_fact <- factor_mat(S1, reg=reg_S)
    # Pre-compute inverse factor L^-1 where S1 ~ LL'
    Linv <- solve(s_fact$ch, Diagonal(d))
    solve_step <- function(V) Linv %*% (S2 %*% V)
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
    V_hat <- solve_step(V)
    
    # Orthonormalize V_hat
    V_hat <- orthonormalize(V_hat)
    
    # Form S and T efficiently (avoid d x d intermediate products)
    S1_Vhat <- S1 %*% V_hat         # d x q
    S2_Vhat <- S2 %*% V_hat         # d x q
    S_mat <- t(V_hat) %*% S1_Vhat # q x q
    T_mat <- t(V_hat) %*% S2_Vhat # q x q
    
    # Regularize T if needed
    T_mat_reg <- T_mat + Diagonal(q, reg_T)
    
    # Ensure invertibility
    tries <- 0
    success <- FALSE
    while (tries < 5 && !success) {
      cond_num <- try({
        T_inv <- solve(as.matrix(T_mat_reg))
        TRUE
      }, silent=TRUE)
      
      if (inherits(cond_num, "try-error")) {
        T_mat_reg <- T_mat_reg + Diagonal(q, reg_T * 10^(tries+1))
        tries <- tries + 1
      } else {
        success <- TRUE
      }
    }
    if (!success) {
      stop("Unable to invert T_mat even after increasing reg. Possibly q too large or data ill-conditioned.")
    }
    
    T_inv <- solve(T_mat_reg)
    M <- T_inv %*% S_mat
    
    # The eigenvalues of M are the generalized eigenvalues
    # largest/smallest:
    # If largest: just eigen(M)
    # If smallest: eigen(M) but we want the smallest generalized eigenvalues.
    # Actually, since eigen(...) returns descending order:
    # For largest: top q from eigen(M) are largest
    # For smallest: top q from eigen(-M) are smallest eigenvals
    if (which == "smallest") {
      eig_res <- eigen(-M, symmetric = TRUE)
      lambda <- -eig_res$values[1:q]
      W <- eig_res$vectors[,1:q]
    } else {
      eig_res <- eigen(M, symmetric = TRUE)
      lambda <- eig_res$values[1:q]
      W <- eig_res$vectors[,1:q]
    }
    
    V_new <- V_hat %*% W
    V_new <- orthonormalize(V_new)
    
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
    warning("Reached max_iter without convergence (Δλ > tol)")
  }
  
  list(values = lambda, vectors = V, final_reg_S = final_reg_S)
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