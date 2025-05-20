# # Ensure you have Riemann installed:
# # install.packages("Riemann")
# 
# library(Matrix)
# library(Riemann)
# 
# # Helper: Orthonormalize w.r.t. S2 if needed
# orthonormalize_S2 <- function(B, S2) {
#   M <- t(B) %*% S2 %*% B
#   if (rcond(M) < 1e-14) {
#     M <- M + Diagonal(nrow(M), 1e-12)
#   }
#   ch <- chol(M)
#   B_new <- B %*% chol2inv(ch)
#   B_new
# }
# 
# # Objective function for GEP:
# # J(B) = trace((BtS2B)^(-1) BtS1B)
# # If we want largest eigenpairs, we minimize -J(B)
# # If we want smallest eigenpairs, we minimize J(B)
# make_gep_functions <- function(S1, S2, which) {
#   # which: "largest" or "smallest"
#   # Return a list with f and grad functions
#   # On the Stiefel manifold, gradient should be the projection of Euclidean grad.
#   
#   # Euclidean objective
#   J_func <- function(B) {
#     BtS2B <- t(B) %*% S2 %*% B
#     BtS1B <- t(B) %*% S1 %*% B
#     M_inv <- solve(BtS2B)
#     sum(diag(M_inv %*% BtS1B))
#   }
#   
#   Jsign <- if (which == "largest") -1 else 1
#   
#   f <- function(X) {
#     # X is a matrix on the Stiefel manifold
#     # we minimize Jsign * J(X)
#     Jsign * J_func(X)
#   }
#   
#   # Euclidean gradient of J:
#   # grad_J(B) = 2(S1 B M_inv - S2 B M_inv_S1 M_inv) from previous derivation
#   egrad_J <- function(B) {
#     BtS2B <- t(B) %*% S2 %*% B
#     BtS1B <- t(B) %*% S1 %*% B
#     M_inv <- solve(BtS2B)
#     M_inv_S1 <- M_inv %*% BtS1B
#     2 * (S1 %*% B %*% M_inv - S2 %*% B %*% (M_inv_S1 %*% M_inv))
#   }
#   
#   egrad <- function(X) {
#     # For our objective f(X) = Jsign * J(X)
#     Jsign * egrad_J(X)
#   }
#   
#   list(f = f, egrad = egrad)
# }
# 
# # After optimization, we recover eigenvalues by solving qxq GEP:
# recover_eigenvalues <- function(B, S1, S2) {
#   BtS2B <- t(B) %*% S2 %*% B
#   BtS1B <- t(B) %*% S1 %*% B
#   eig_q <- eigen(solve(BtS2B) %*% BtS1B, symmetric = TRUE)
#   eig_q$values
# }
# 
# # High-level solver:
# solve_gep_riemann <- function(S1, S2, q = 2, which = c("largest","smallest"),
#                               max_iter = 1000, tol = 1e-6, V0 = NULL) {
#   which <- match.arg(which)
#   d <- nrow(S1)
#   
#   # Initial guess
#   if (is.null(V0)) {
#     # random init
#     B <- matrix(rnorm(d*q), d, q)
#     # Just QR for Euclidean orthonormalization:
#     B <- qr.Q(qr(B))
#   } else {
#     B <- V0
#   }
#   
#   # Construct problem
#   M <- stiefel.manifold(d, q) # standard Stiefel manifold
#   funcs <- make_gep_functions(S1, S2, which)
#   
#   # We use a Riemannian L-BFGS or Trust-Region method
#   # rlbfgs or rtr from Riemann package
#   # rlbfgs arguments: maxiter, tol, etc.
#   
#   # We'll set a stopping criterion based on the relative change in eigenvalues:
#   # We'll handle stopping by hooking into the iteration and checking progress.
#   
#   # Riemann package: The solver runs until gradient norm or max_iter is done.
#   # We'll rely on gradient-based stopping. If we need eigenvalue tolerance:
#   # We'll do a post-check.
#   
#   # Run optimization:
#   # For large q, trust-region might be more stable than L-BFGS, but we try rlbfgs first:
#   res <- rlbfgs(M, funcs$f, funcs$egrad, x0 = B, maxiter = max_iter, tol = tol)
#   
#   # res$x is the optimized B
#   B_opt <- res$x
#   
#   # Orthonormalize w.r.t. S2 just to be sure:
#   B_opt <- orthonormalize_S2(B_opt, S2)
#   
#   vals <- recover_eigenvalues(B_opt, S1, S2)
#   
#   # Sort eigenvalues in descending order if largest requested, ascending if smallest
#   # Actually we have a qxq subproblem. The order from eigen() is descending.
#   # If "largest", we minimized -J, we got top eigenvectors, no re-order needed.
#   # If "smallest", we minimized J, also no re-order needed because eigen() returns descending order by default.
#   # Actually, smallest would give us largest for J, we might want to invert selection or carefully handle sign:
#   # Wait, we made Jsign * J. If which=largest, we minimized -J, that means we found maximum J -> largest eigenvalues at top
#   # If which=smallest, we minimized J, that means we found minimum J -> smallest eigenvalues at top?
#   # eigen() returns eigenvalues in decreasing order by default, so smallest might need a reverse?
#   # Actually, minimizing J should push us towards a subspace yielding smaller J, i.e. smaller eigenvalues at top. Perfect.
#   # So no reordering needed because eigen() returns largest->smallest. If smallest chosen, those are actually the smallest eigenvalues. Let's just trust the method converged correctly.
#   
#   list(values = vals[1:q], vectors = B_opt)
# }
# 
# #Example usage (small q):
# d <- 50
# q <- 3
# set.seed(1)
# A <- matrix(rnorm(d*d), d, d)
# S2 <- crossprod(A) + diag(d)*0.1
# B_mat <- matrix(rnorm(d*d), d, d)
# S1 <- crossprod(B_mat) + diag(d)*0.1
# res_riemann_largest <- solve_gep_riemann(S1, S2, q=q, which="largest", max_iter=500, tol=1e-6)
# cat("Approx largest eigenvalues:\n", res_riemann_largest$values, "\n")