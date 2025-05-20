library(Matrix)

orthonormalize_S2 <- function(B, S2) {
  M <- t(B) %*% S2 %*% B
  # Add slightly more stable regularization if needed
  if (rcond(M) < 1e-14) {
    M <- M + Diagonal(nrow(M), 1e-10)
  }
  ch <- chol(M)
  B_new <- B %*% chol2inv(ch)
  B_new
}

gep_objective <- function(B, S1, S2) {
  BtS2B <- t(B) %*% S2 %*% B
  BtS1B <- t(B) %*% S1 %*% B
  M_inv <- solve(BtS2B)
  sum(diag(M_inv %*% BtS1B))
}

gep_gradient <- function(B, S1, S2) {
  BtS2B <- t(B) %*% S2 %*% B
  BtS1B <- t(B) %*% S1 %*% B
  M_inv <- solve(BtS2B)
  M_inv_S1 <- M_inv %*% BtS1B
  # gradient of trace((BtS2B)^{-1} BtS1B)
  2 * (S1 %*% B %*% M_inv - S2 %*% B %*% M_inv_S1 %*% M_inv)
}

solve_gep_gradient <- function(S1, S2, q = 2, which = c("largest","smallest"),
                               max_iter = 1000, tol = 1e-6, V0 = NULL,
                               step_size = 1e-3, beta = 0.5, c_armijo = 1e-4,
                               max_line_search = 20) {
  which <- match.arg(which)
  d <- nrow(S1)
  
  if (is.null(V0)) {
    B <- matrix(rnorm(d*q), d, q)
  } else {
    B <- V0
  }
  B <- orthonormalize_S2(B, S2)
  
  J_old <- gep_objective(B, S1, S2)
  lambda_old <- NULL
  
  for (iter in seq_len(max_iter)) {
    grad <- gep_gradient(B, S1, S2)
    direction <- if (which == "smallest") -grad else grad
    
    # Armijo backtracking line search
    alpha <- step_size
    # Armijo condition: J(B + alpha * dir) >= J(B) + c_armijo * alpha * <grad, dir>
    # for largest: want increase in J
    # for smallest: we minimized J, want decrease
    # unify by considering the objective difference:
    
    # directional derivative approximation:
    dir_dot_grad <- sum(direction * grad)
    # For largest, we want J_new > J_old + c_armijo * alpha * dir_dot_grad
    # Note that dir_dot_grad in largest case should be positive if correct direction
    # For smallest, since direction = -grad, dir_dot_grad < 0, J_new < J_old + ...
    # Actually for smallest: want J_new <= J_old + c_armijo * alpha * dir_dot_grad (dir_dot_grad<0)
    
    # We'll implement conditions accordingly:
    improvement_condition <- function(J_new, J_old, alpha, dir_dot_grad) {
      if (which == "largest") {
        return(J_new >= J_old + c_armijo * alpha * dir_dot_grad)
      } else {
        # smallest
        return(J_new <= J_old + c_armijo * alpha * dir_dot_grad)
      }
    }
    
    success <- FALSE
    for (ls_iter in seq_len(max_line_search)) {
      B_new <- B + alpha * direction
      B_new <- orthonormalize_S2(B_new, S2)
      J_new <- gep_objective(B_new, S1, S2)
      
      if (improvement_condition(J_new, J_old, alpha, dir_dot_grad)) {
        success <- TRUE
        break
      } else {
        alpha <- alpha * beta
      }
    }
    
    if (!success) {
      # no improvement
      break
    }
    
    B <- B_new
    J_old <- J_new
    
    BtS2B_new <- t(B) %*% S2 %*% B
    BtS1B_new <- t(B) %*% S1 %*% B
    eig_q <- eigen(solve(BtS2B_new) %*% BtS1B_new, symmetric=TRUE)
    lambdas <- eig_q$values[1:q]
    
    if (!is.null(lambda_old)) {
      rel_change <- max(abs(lambdas - lambda_old) / pmax(abs(lambdas), 1e-12))
      if (rel_change < tol) {
        break
      }
    }
    lambda_old <- lambdas
  }
  
  list(values = lambda_old, vectors = B)
}

##############################################
# Example usage:
d <- 50
q <- 3
set.seed(1)
A <- matrix(rnorm(d*d), d, d)
S2 <- crossprod(A) + diag(d)*0.1
B_mat <- matrix(rnorm(d*d), d, d)
S1 <- crossprod(B_mat) + diag(d)*0.1

res_gradient_largest <- solve_gep_gradient(S1, S2, q=q, which="largest", max_iter=2000, step_size=1e-4)
cat("Approx largest eigenvalues:\n", res_gradient_largest$values, "\n")


res_gradient_smallest <- solve_gep_gradient(S1, S2, q=q, which="smallest", max_iter=200, step_size=1e-4)
cat("Approx smallest eigenvalues:\n", res_gradient_smallest$values, "\n")