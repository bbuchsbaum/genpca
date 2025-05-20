# Import necessary functions from packages
#' @importFrom Matrix Matrix Diagonal crossprod tcrossprod t solve bandSparse sparseMatrix Cholesky rowSums diag<-
#' @importFrom RSpectra svds eigs
#' @importFrom FNN get.knn
#' @importFrom stats quantile var median mad lm predict sd
NULL

# Function to construct the second differences matrix
#' @rdname sfpca
#' @keywords internal
second_diff_matrix <- function(n) {
  # n: Length of the time series
  if (n < 3) stop("n must be at least 3 to compute second differences.")
  
  # Number of second differences is (n - 2)
  m <- n - 2
  
  # Create the diagonals for the second difference matrix
  # Diagonals for D2: (1, -2, 1) pattern for n-2 rows, n cols
  # k=0: 1s at (1,1), (2,2), ..., (m,m)
  # k=1: -2s at (1,2), (2,3), ..., (m,m+1)
  # k=2: 1s at (1,3), (2,4), ..., (m,n)
  diagonals <- list(
    rep(1, m),   # for k=0 offset
    rep(-2, m),  # for k=1 offset
    rep(1, m)    # for k=2 offset
  )

  # Build the (n - 2) x n second difference matrix using bandSparse
  # Correct offsets k = c(0, 1, 2) relative to main diagonal
  D2 <- Matrix::bandSparse(m, n, k = 0:2, diagonals = diagonals, symmetric = FALSE)
  
  return(D2)
}




#' Sparse and Functional Principal Components Analysis (SFPCA) with Spatial Coordinates
#'
#' Performs Sparse and Functional PCA on a data matrix, allowing for both sparsity and smoothness
#' in the estimated principal components. Includes heuristic estimation of penalty parameters.
#' The spatial smoothness penalty is constructed based on provided spatial coordinates.
#'
#' @param X A numeric data matrix of dimensions n (observations/time points) by p (variables/space).
#' @param K The number of principal components to estimate.
#' @param spat_cds A matrix of spatial coordinates for each column of X (variables). Each row
#'   corresponds to a spatial dimension (e.g., x, y, z), and each column corresponds to a variable.
#' @param lambda_u Sparsity penalty parameter for u. If NULL, estimated heuristically.
#' @param lambda_v Sparsity penalty parameter for v. If NULL, estimated heuristically.
#' @param alpha_u Smoothness penalty parameter for u. If NULL, estimated heuristically.
#' @param alpha_v Smoothness penalty parameter for v. If NULL, estimated heuristically.
#' @param Omega_u A positive semi-definite matrix for smoothness penalty on u. If NULL, defaults to
#'   second differences penalty (sparse matrix).
#' @param penalty_u The penalty function for u. Can be "l1" (lasso), "scad", or a custom function.
#' @param penalty_v The penalty function for v. Can be "l1" (lasso), "scad", or a custom function.
#' @param uthresh Quantile for selecting `lambda_u` when estimated.
#' @param vthresh Quantile for selecting `lambda_v` when estimated.
#' @param knn Number of nearest neighbours for constructing `Omega_v`.
#' @param max_iter Maximum number of iterations for the alternating optimization.
#' @param tol Tolerance for convergence.
#' @param verbose Logical; if TRUE, prints progress messages.
#' @return A list containing the estimated singular values `d`, left singular vectors `u`, right
#'   singular vectors `v`, and the penalty parameters used.
#' @examples
#' library(Matrix)
#' set.seed(123)
#' # Simulate a small example due to resource constraints
#' n <- 100  # Number of time points
#' p <- 50   # Number of spatial locations
#' X <- matrix(rnorm(n * p), n, p)
#' # Simulate spatial coordinates
#' spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)  # 3D coordinates
#' result <- sfpca(X, K = 2, spat_cds = spat_cds)
#' @export
sfpca <- function(X, K, spat_cds,
                  lambda_u = NULL, lambda_v = NULL,
                  alpha_u = NULL, alpha_v = NULL,
                  Omega_u = NULL,
                  penalty_u = "l1", penalty_v = "l1",
                  uthresh=.9, vthresh=.9,
                  knn=min(6, ncol(X) - 1),  # Number of nearest neighbors for spatial penalty
                  max_iter = 100, tol = 1e-6, verbose = FALSE) {
 
  n <- nrow(X)
  p <- ncol(X)
  d_list <- numeric(K)
  u_list <- matrix(0, n, K)
  v_list <- matrix(0, p, K)
  lambda_u_list <- numeric(K)
  lambda_v_list <- numeric(K)
  alpha_u_list <- numeric(K)
  alpha_v_list <- numeric(K)
  
  # Convert X to sparse matrix if not already
  X <- Matrix(X)
  
  # Default Omega_u if not provided (second differences)
  if (is.null(Omega_u)) {
    D_n <- second_diff_matrix(n)
    Omega_u <- Matrix::crossprod(Matrix(D_n, sparse = TRUE))
    #D_n <- diff(diag(n), differences = 2)
    #Omega_u <- crossprod(Matrix(D_n, sparse = TRUE))
  }
  
  # Ensure Omega_u is sparse
  #Omega_u <- as(Omega_u, "dgCMatrix")
  
  # Construct Omega_v based on spatial coordinates
  # We will use a weighted graph Laplacian where weights are based on distances
  if (is.null(spat_cds)) {
    stop("spat_cds must be provided for constructing the spatial penalty matrix Omega_v.")
  }
  
  if (verbose) cat("Constructing spatial penalty matrix Omega_v based on spat_cds...\n")
  
  Omega_v <- construct_spatial_penalty(spat_cds, k=knn)
  
  # Deflation method
  X_residual <- X
  for (k in 1:K) {
    if (verbose) cat("Component", k, "\n")
    # Estimate initial u and v using SVD
    svd_res <- svds(X_residual, k = 1, nu = 1, nv = 1)
    u_init <- as.numeric(svd_res$u)
    v_init <- as.numeric(svd_res$v)
    # Heuristic estimation of penalty parameters
    # Sparsity penalties based on quantiles of initial singular vectors
    if (is.null(lambda_u)) {
      lambda_u_k <- quantile(abs(u_init), uthresh)
    } else {
      lambda_u_k <- lambda_u
    }
    if (is.null(lambda_v)) {
      lambda_v_k <- quantile(abs(v_init), vthresh)
    } else {
      lambda_v_k <- lambda_v
    }
    
    # Smoothness penalties based on variance of second differences or spatial distances
    if (is.null(alpha_u)) {
      u_diff2 <- diff(u_init, differences = 2)
      alpha_u_k <- 1 / (var(u_diff2) + 1e-6)
    } else {
      alpha_u_k <- alpha_u
    }
    if (is.null(alpha_v)) {
      # For spatial smoothness, we can set alpha_v based on the average distance between points
      avg_dist <- mean(Omega_v@x)
      alpha_v_k <- 1 / (avg_dist + 1e-6)
    } else {
      alpha_v_k <- alpha_v
    }
    
    # Store penalty parameters
    lambda_u_list[k] <- lambda_u_k
    lambda_v_list[k] <- lambda_v_k
    alpha_u_list[k] <- alpha_u_k
    alpha_v_list[k] <- alpha_v_k
    
    if (verbose) {
      cat("Heuristic penalties for component", k, ":\n")
      cat("lambda_u =", lambda_u_k, "lambda_v =", lambda_v_k, "\n")
      cat("alpha_u =", alpha_u_k, "alpha_v =", alpha_v_k, "\n")
    }
    
    # Run rank-1 SFPCA with estimated penalties
    result <- sfpca_rank1(X_residual,
                          lambda_u = lambda_u_k, lambda_v = lambda_v_k,
                          alpha_u = alpha_u_k, alpha_v = alpha_v_k,
                          Omega_u = Omega_u, Omega_v = Omega_v,
                          penalty_u = penalty_u, penalty_v = penalty_v,
                          max_iter = max_iter, tol = tol, verbose = verbose,
                                               u_init = u_init, v_init = v_init)
    
    
    d_list[k] <- result$d
    u_list[, k] <- as.vector(result$u)
    v_list[, k] <- as.vector(result$v)
    # Deflate X
    X_residual <- X_residual - d_list[k] * result$u %o% result$v
  }
  
  return(list(d = d_list, u = u_list, v = v_list,
              lambda_u = lambda_u_list, lambda_v = lambda_v_list,
              alpha_u = alpha_u_list, alpha_v = alpha_v_list))
}

# Function to construct the spatial penalty matrix Omega_v based on spat_cds
#' @rdname sfpca
#' @keywords internal
construct_spatial_penalty <- function(spat_cds, method = "distance", k = 6L) {
  p <- ncol(spat_cds)
  if (p <= k) {
     warning(paste0("Number of spatial locations (p=", p, ") is less than or equal to knn (k=", k, "). ",
                    "Reducing k to p-1 = ", p-1, "."))
     k <- p - 1
  }
  if (k < 1) {
      stop("knn (k) must be at least 1.")
  }

  if (method == "distance") {
    # Compute pairwise distances between spatial coordinates is too slow for large p
    # Build adjacency matrix based on k-nearest neighbors
    knn_res <- FNN::get.knn(t(spat_cds), k = k)
    indices <- knn_res$nn.index
    distances_knn <- knn_res$nn.dist

    # Initialize sparse matrix for weights
    # Use 1 / (dist + epsilon) for weights
    weights <- 1 / (as.vector(t(distances_knn)) + 1e-6)
    
    # Create sparse matrix W with weights
    W <- Matrix::sparseMatrix(i = rep(1:p, each = k),
                              j = as.vector(t(indices)),
                              x = weights,
                              dims = c(p, p), index1 = TRUE) # Ensure index1=TRUE
    
    # Make W symmetric: W = (W + W^T) / 2
    # This ensures that if j is a neighbor of i, i is also considered a neighbor of j
    # with potentially averaged weight.
    W <- (W + Matrix::t(W)) / 2
    # Ensure diagonal is zero after symmetrization
    Matrix::diag(W) <- 0
    
  } else if (method == "neighbor") {
    # If a neighbor graph is provided, construct W based on it
    # This part can be customized based on available neighbor information
    stop("Neighbor graph method not implemented in this example.")
  } else {
    stop("Invalid method for constructing spatial penalty.")
  }
  
  # Construct graph Laplacian L = D - W
  # D is the diagonal matrix of degrees (row sums of W)
  D <- Matrix::Diagonal(x = Matrix::rowSums(W))
  L <- D - W
  
  # Omega_v is the graph Laplacian
  Omega_v <- L
  
  return(Omega_v)
}

# Helper function for rank-1 SFPCA (same as before)
#' @noRd
sfpca_rank1 <- function(X,
                        lambda_u, lambda_v,
                        alpha_u, alpha_v,
                        Omega_u, Omega_v,
                        penalty_u, penalty_v,
                        max_iter, tol, verbose,
                        u_init = NULL, v_init = NULL) {
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required but not installed.")
  }
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("Package 'RSpectra' is required but not installed.")
  }
  n <- nrow(X)
  p <- ncol(X)
  # Initialize u and v
  if (is.null(u_init) || is.null(v_init)) {
    svd_res <- RSpectra::svds(X, k = 1, nu = 1, nv = 1)
    u <- as.numeric(svd_res$u)
    v <- as.numeric(svd_res$v)
    d <- svd_res$d[1]
  } else {
    u <- u_init
    v <- v_init
  }
  
  # Precompute S matrices and Lipschitz constants
  S_u <- Matrix::Diagonal(n) + alpha_u * Omega_u
  S_v <- Matrix::Diagonal(p) + alpha_v * Omega_v
  
  # Compute the largest eigenvalues using eigs_sym
  L_u <- eigs_sym(S_u, 1, which = "LM")$values  
  L_v <- eigs_sym(S_v, 1, which = "LM")$values
  
  # Ensure L_u and L_v are positive scalars
  L_u <- as.numeric(L_u)
  L_v <- as.numeric(L_v)
  
  if (is.na(L_u) || L_u <= 0) {
    stop("Invalid Lipschitz constant L_u. Check the construction of S_u.")
  }
  if (is.na(L_v) || L_v <= 0) {
    stop("Invalid Lipschitz constant L_v. Check the construction of S_v.")
  }
  
  iter <- 0
  obj_old <- -Inf
  repeat {
    iter <- iter + 1
    # Update u with fixed v
    u <- sfpca_proximal_operator(X %*% v, S_u, lambda_u, L_u, penalty_u)
    u_norm <- sqrt(as.numeric(Matrix::crossprod(u, S_u %*% u)))
    if (u_norm > 0) u <- u / u_norm else u <- rep(0, n)
    
    # Update v with fixed u
    v <- sfpca_proximal_operator(Matrix::crossprod(X, u), S_v, lambda_v, L_v, penalty_v)
    v_norm <- sqrt(as.numeric(Matrix::crossprod(v, S_v %*% v)))
    if (v_norm > 0) v <- v / v_norm else v <- rep(0, p)
    
    # Compute objective value
    obj <- as.numeric(
      Matrix::crossprod(u, X %*% v) -
        lambda_u * penalty_function(u, penalty_u, lambda_u) -
        lambda_v * penalty_function(v, penalty_v, lambda_v)
    )
    if (verbose) cat("Iteration", iter, "Objective:", obj, "\n")
    
    # Check convergence
    if (abs(obj - obj_old) < tol || iter >= max_iter) break
    obj_old <- obj
  }
  # Final scaling
  d <- as.numeric(Matrix::crossprod(u, X %*% v))
  return(list(u = u, v = v, d = d))
}

#' @noRd
sfpca_proximal_operator <- function(z, S, lambda, L, penalty) {
  # Gradient step
  y <- z / L
  
  # Check for valid y
  if (any(is.na(y) | is.infinite(y))) {
    stop("Invalid values in y during proximal operator computation.")
  }
  
  # Adjust for S matrix
  x <- Matrix::solve(S + Matrix::Diagonal(n = length(y), x = 1e-6), y)  # Regularization for numerical stability
  
  # Proximal operator depending on penalty
  if (penalty == "l1") {
    x <- soft_threshold(x, lambda / L)
  } else if (penalty == "scad") {
    x <- scad_proximal(x, lambda / L)
  } else {
    stop("Unsupported penalty function.")
  }
  return(x)
}

# Soft-thresholding operator (same as before)
#' @noRd
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(0, abs(x) - lambda)
}

# SCAD proximal operator (same as before)
#' @noRd
scad_proximal <- function(x, lambda, a = 3.7) {
  s <- sign(x)
  abs_x <- abs(x)
  x1 <- s * pmax(0, abs_x - lambda)
  idx <- abs_x > lambda & abs_x <= a * lambda
  x1[idx] <- s[idx] * ((a - 1) * abs_x[idx] - a * lambda) / (a - 2)
  x1[abs_x > a * lambda] <- x[abs_x > a * lambda]
  return(x1)
}

# Penalty function value (same as before)
#' @noRd
penalty_function <- function(x, penalty, lambda = 1) {
  if (penalty == "l1") {
    sum(abs(x))
  } else if (penalty == "scad") {
    a <- 3.7
    abs_x <- abs(x)
    penalty_values <- ifelse(
      abs_x <= lambda,
      lambda * abs_x,
      ifelse(
        abs_x <= a * lambda,
        (-abs_x^2 + 2 * a * lambda * abs_x - lambda^2) /
          (2 * (a - 1)) + lambda^2,
        ((a + 1) * lambda^2) / 2
      )
    )
    sum(penalty_values)
  } else {
    stop("Unsupported penalty function.")
  }
}

# Efficient eigenvalue computations for symmetric sparse matrices (same as before)
#' @noRd
eigs_sym <- function(A, k = 1, which = "LM") {
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("Package 'RSpectra' is required but not installed.")
  }
  RSpectra::eigs(A, k = k, which = which, opts = list(tol = 1e-3), symmetric = TRUE)
}

