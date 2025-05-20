#' @keywords internal
#' @import assertthat
prep_constraints <- function(X, A, M, tol = 1e-8, remedy = c("error", "ridge", "clip", "identity"), verbose = FALSE) {
  n <- nrow(X)
  p <- ncol(X)
  
  remedy <- match.arg(remedy)
  
  # --- Process A (Column constraints) ---
  if (is.null(A)) {
    A <- Matrix::sparseMatrix(i=1:p, j=1:p, x=1.0, dims=c(p,p))
  } else if (is.vector(A)) {
    assert_that(length(A) == p, msg="Length of vector A must equal ncol(X)")
    assert_that(all(A >= -tol), msg="Diagonal elements of A (from vector) must be non-negative.")
    A <- Matrix::Diagonal(n=p, x=A)
  } else {
    # Ensure it's a Matrix object first for isSymmetric etc.
    if (!is(A, "Matrix")) { A <- Matrix::Matrix(A, sparse=TRUE) }
    assert_that(nrow(A) == p, msg=paste("nrow(A) != ncol(X) -- ", nrow(A), " != ", p))
    assert_that(ncol(A) == p, msg=paste("ncol(A) != ncol(X) -- ", ncol(A), " != ", p))
    if (!Matrix::isSymmetric(A)) {
        warning("Matrix A is not symmetric, taking the upper triangle.")
        A <- Matrix::forceSymmetric(A, uplo="U")
    }
    # Ensure sparse format (forceSymmetric might return dsyMatrix)
    if (!is(A, "sparseMatrix")) { A <- as(A, "sparseMatrix") }
    # Check PSD (only eigenvalues for efficiency)
    min_eig_A <- tryCatch(RSpectra::eigs_sym(A, k=1, which="SA")$values,
                          error=function(e) {NA}) # Handle potential errors
    if (is.na(min_eig_A) || min_eig_A < -tol) {
        if (verbose) message("Matrix A is not PSD (min eig: ", signif(min_eig_A, 4), "). Applying remedy: ", remedy)
        if (remedy == "error") {
            stop("Matrix A must be positive semi-definite (smallest eigenvalue >= -tol). Error was: ", min_eig_A)
        } else if (remedy == "ridge") {
            ridge_val <- tol - min_eig_A # Amount to add to diagonal
            A <- A + Matrix::Diagonal(p, x = ridge_val)
            A <- (A + Matrix::t(A)) / 2 # Re-symmetrize
        } else if (remedy == "clip") {
            A_dense <- as.matrix((A + Matrix::t(A))/2) # Ensure symmetry for eigen
            E <- eigen(A_dense, symmetric = TRUE)
            vals_clipped <- pmax(E$values, tol) # Clip negative eigenvalues
            A <- Matrix::Matrix(E$vectors %*% Matrix::diag(vals_clipped) %*% Matrix::t(E$vectors), sparse=TRUE)
            A <- (A + Matrix::t(A)) / 2 # Re-symmetrize after reconstruction
        } else if (remedy == "identity") {
            A <- Matrix::Diagonal(p)
        }
    }
  }
  
  # --- Process M (Row constraints) ---
  if (is.null(M)) {
    M <- Matrix::sparseMatrix(i=1:n, j=1:n, x=1.0, dims=c(n,n))
  } else if (is.vector(M)) {
    assert_that(length(M) == n, msg="Length of vector M must equal nrow(X)")
    assert_that(all(M >= -tol), msg="Diagonal elements of M (from vector) must be non-negative.")
    M <- Matrix::Diagonal(n=n, x=M)
  } else {
    if (!is(M, "Matrix")) { M <- Matrix::Matrix(M, sparse=TRUE) }
    assert_that(nrow(M) == n, msg=paste("nrow(M) != nrow(X) -- ", nrow(M), " != ", n))
    assert_that(ncol(M) == n, msg=paste("ncol(M) != nrow(X) -- ", ncol(M), " != ", n))
    if (!Matrix::isSymmetric(M)) {
        warning("Matrix M is not symmetric, taking the upper triangle.")
        M <- Matrix::forceSymmetric(M, uplo="U")
    }
    if (!is(M, "sparseMatrix")) { M <- as(M, "sparseMatrix") }
    # Check PSD
    min_eig_M <- tryCatch(RSpectra::eigs_sym(M, k=1, which="SA")$values,
                          error=function(e) {NA})
    if (is.na(min_eig_M) || min_eig_M < -tol) {
        if (verbose) message("Matrix M is not PSD (min eig: ", signif(min_eig_M, 4), "). Applying remedy: ", remedy)
        if (remedy == "error") {
            stop("Matrix M must be positive semi-definite (smallest eigenvalue >= -tol). Error was: ", min_eig_M)
        } else if (remedy == "ridge") {
            ridge_val <- tol - min_eig_M
            M <- M + Matrix::Diagonal(n, x = ridge_val)
            M <- (M + Matrix::t(M)) / 2
        } else if (remedy == "clip") {
            M_dense <- as.matrix((M + Matrix::t(M))/2)
            E <- eigen(M_dense, symmetric = TRUE)
            vals_clipped <- pmax(E$values, tol)
            M <- Matrix::Matrix(E$vectors %*% Matrix::diag(vals_clipped) %*% Matrix::t(E$vectors), sparse=TRUE)
            M <- (M + Matrix::t(M)) / 2
        } else if (remedy == "identity") {
            M <- Matrix::Diagonal(n)
        }
    }
  }
  
  # Ensure consistent sparse format (dgCMatrix is common)
  list(A = as(A, "dgCMatrix"), M = as(M, "dgCMatrix"))
}

#' Generalised Principal Components Analysis (GPCA)
#'
#' Implements the Generalised Least‑Squares Matrix Decomposition of
#' Allen, Grosenick & Taylor (2014) for data observed in a **row**
#' inner‑product space M and a **column** inner‑product space A.
#' Setting M = I_n, A = I_p recovers ordinary PCA.
#'
#' @section Method:
#' We compute the rank‑ncomp factors UDVT that minimise
#' \deqn{ \|X - UDV^\top\|_{M,A}^2
#'       = \mathrm{tr}\!\bigl(M\, (X-UDV^\top)\,A\,(X-UDV^\top)^\top\bigr) }
#' subject to UT M U = I, VT AV = I. (Allen et al., 2014).
#' Three methods are available via the `method` argument:
#' \itemize{
#'  \item{\code{"eigen"} (Default): Uses a one-shot eigen decomposition strategy based on \code{gmdLA}. It explicitly forms and decomposes a \eqn{p \times p} or \eqn{n \times n} matrix (depending on \code{n} vs \code{p}).}
#'  \item{\code{"spectra"}: Uses a matrix-free iterative approach via the \pkg{RcppSpectra} package to solve the same eigen problem as \code{"eigen"} but without forming the large intermediate matrix. Generally faster and uses less memory for large \code{n} or \code{p}. Requires C++ compiler and \pkg{RcppSpectra}.}
#'  \item{\code{"deflation"}: Uses an iterative power/deflation algorithm. Can be slower but potentially uses less memory than \code{"eigen"} for very large dense problems where \code{ncomp} is small.}
#' }
#'
#' @param X   Numeric matrix n x p.
#' @param A   Column constraint: vector (implies diagonal), dense matrix, or sparse
#'            symmetric p x p PSD matrix. If `NULL`, defaults to identity.
#' @param M   Row constraint: vector (implies diagonal), dense matrix, or sparse
#'            symmetric n x n PSD matrix. If `NULL`, defaults to identity.
#' @param ncomp Number of components to extract. Defaults to `min(dim(X))`. Must be positive.
#' @param method Character string specifying the computation method. One of \code{"eigen"} (default, uses \code{gmdLA}), \code{"spectra"} (uses matrix-free C++/Spectra implementation \code{gmd_fast_cpp}), or \code{"deflation"} (uses \code{gmd_deflationR} or \code{gmd_deflation_cpp}).
#' @param constraints_remedy Character string specifying the remedy for constraints. One of \code{"error"}, \code{"ridge"}, \code{"clip"}, or \code{"identity"}.
#' @param preproc Pre‑processing transformer object from the **multivarious** package
#'                (default `multivarious::pass()`). Use `multivarious::center()` for centered GPCA.
#'                See `?multivarious::prep` for options.
#' @param threshold Convergence tolerance for the \code{"deflation"} method's inner loop. Default `1e-6`.
#' @param use_cpp Logical. If `TRUE` (default) and package was compiled with C++ support,
#'                use faster C++ implementation for \code{method = "deflation"}. Fallback to R otherwise.
#'                (Ignored for \code{method = "eigen"} and \code{method = "spectra"}).
#' @param maxeig Upper bound on subspace dimension for eigen/SVD calculations, primarily for
#'               \code{method = "eigen"}. Affects internal calls to `RSpectra::eigs_sym`. Default `800`.
#' @param maxit_spectra Maximum iterations for the Spectra iterative solver when \code{method = "spectra"}. Default `1000`.
#' @param tol_spectra Tolerance for the Spectra iterative solver when \code{method = "spectra"}. Default `1e-9`.
#' @param verbose Logical. If `TRUE`, print progress messages. Default `FALSE`.
#'
#' @return An object of class `c("genpca", "bi_projector")` inheriting from `multivarious::bi_projector`,
#'   with slots including:
#'   \describe{
#'     \item{u,v}{Left/right singular vectors scaled by the constraint metrics
#'                (MU, AV). These correspond to loadings/components in the original space's geometry.
#'                Use `loadings(fit)` or `components(fit)`.} 
#'     \item{ou,ov}{Orthonormal singular vectors in the constraint metric
#'                  (U, V such that UT M U = I, VT AV = I). These are the core mathematical factors.}
#'     \item{sdev}{Generalised singular values d_k.}
#'     \item{s}{Scores ( X V or equivalently MU D). Represent projection of rows onto components. Use `scores(fit)`.} 
#'     \item{preproc}{The `multivarious` pre‑processing object used.}
#'     \item{A, M}{The constraint matrices used (potentially after coercion to sparse format).}
#'     \item{propv}{Proportion of generalized variance explained by each component.}
#'     \item{cumv}{Cumulative proportion of generalized variance explained.}
#'   }
#'
#' @references
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014).
#' *A Generalized Least‑Squares Matrix Decomposition.* 
#' Journal of the American Statistical Association, 109(505), 145‑159.
#' arXiv:1102.3074.
#'
#' @seealso \code{\link{prep_constraints}}, \code{\link{gmdLA}}, \code{\link{gmd_deflationR}},
#'   \code{\link{gmd_fast_cpp}} (if using method="spectra"), \code{\link{truncate.genpca}},
#'   \code{\link{reconstruct.genpca}},
#'   `multivarious::bi_projector`, `multivarious::project`, `multivarious::scores`,
#'   `multivarious::loadings`, `multivarious::reconstruct`.
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE) &&
#'     requireNamespace("multivarious", quietly = TRUE)) {
#'   set.seed(123)
#'   X <- matrix(stats::rnorm(200 * 100), 200, 100)
#'   rownames(X) <- paste0("R", 1:200)
#'   colnames(X) <- paste0("C", 1:100)
#'
#'   # Standard PCA (A=I, M=I, centered) - using default method="eigen"
#'   gpca_std_eigen <- genpca(X, ncomp = 5, preproc = multivarious::center(), verbose = FALSE)
#'   
#'   # Standard PCA using Spectra method (requires C++ build)
#'   # gpca_std_spectra <- try(genpca(X, ncomp = 5, preproc = multivarious::center(), \n#'   #                              method="spectra", verbose = TRUE))
#'   # if (!inherits(gpca_std_spectra, "try-error")) {
#'   #    print(head(gpca_std_spectra$sdev))
#'   # }
#'
#'   # Compare singular values with prcomp
#'   pr_std <- stats::prcomp(X, center = TRUE, scale. = FALSE)
#'   print("Eigen Method Sdev:")
#'   print(head(gpca_std_eigen$sdev))
#'   print("prcomp Sdev:")
#'   print(head(pr_std$sdev))
#'   print(paste("Total Var Explained (Eigen):"), round(sum(gpca_std_eigen$propv)*100), "%")
#'
#'   # Weighted column PCA (diagonal A, no centering)
#'   col_weights <- stats::runif(100, 0.5, 1.5)
#'   gpca_weighted <- genpca(X, A = col_weights, ncomp = 3, preproc=multivarious::pass(), verbose = FALSE)
#'   print("Weighted GPCA Sdev:")
#'   print(gpca_weighted$sdev)
#'   print(head(loadings(gpca_weighted)))
#' }
#' @useDynLib genpca, .registration = TRUE 
#' @importFrom Rcpp sourceCpp
#' @importFrom multivarious bi_projector init_transform prep pass scores sdev loadings components reconstruct reverse_transform ncomp
#' @importFrom Matrix Matrix isSymmetric isDiagonal diag as t forceSymmetric Diagonal crossprod tcrossprod sweep
#' @importFrom assertthat assert_that
#' @importFrom RSpectra eigs_sym svds
#' @importFrom methods as is
#' @importFrom stats rnorm runif eigen
#' @export
genpca <- function(X, A=NULL, M=NULL, ncomp=NULL,
                   method = c("eigen", "spectra", "deflation"),
                   constraints_remedy = c("error", "ridge", "clip", "identity"),
                   preproc = multivarious::pass(), # Default to pass() for safety
                   threshold = 1e-6, # For deflation
                   use_cpp = TRUE, # For deflation
                   maxeig = 800, # For method="eigen"
                   maxit_spectra = 1000, # For method="spectra"
                   tol_spectra = 1e-9,   # For method="spectra"
                   verbose = FALSE) {

  method <- match.arg(method)
  constraints_remedy <- match.arg(constraints_remedy)

  if (is.null(ncomp)) {
      ncomp <- min(dim(X))
  }
  assert_that(ncomp > 0, msg="ncomp must be positive.")
  ncomp <- min(min(dim(X)), ncomp) # Cannot exceed dimensions

  # Prepare and validate constraints M and A
  if (verbose) message("Preparing constraints...")
  pcon <- prep_constraints(X, A, M, remedy = constraints_remedy, verbose = verbose)
  A <- pcon$A
  M <- pcon$M

  # Prepare pre-processing object
  procres <- multivarious::prep(preproc)
  # Apply pre-processing
  if (verbose) message("Applying pre-processing...")
  Xp <- multivarious::init_transform(procres, X)

  n = nrow(Xp)
  p = ncol(Xp)

  # Check if C++ code is available (specific function name depends on package build)
  # Placeholder check - replace with actual check if package uses compiled code
  cpp_deflation_available <- exists("gmd_deflation_cpp", mode = "function") # Example check
  cpp_spectra_available <- exists("gmd_fast_cpp", mode = "function") # Example check

  if (method == "deflation" && use_cpp && !cpp_deflation_available) {
      if (verbose) message("use_cpp=TRUE but C++ deflation code not found. Falling back to R version.")
      use_cpp <- FALSE # Force R version if C++ not found
  }
  if (method == "spectra" && !cpp_spectra_available) {
      stop("method='spectra' requires the C++ function 'gmd_fast_cpp', which was not found. Ensure the package was compiled correctly with Rcpp/RcppArmadillo/RcppSpectra support.")
  }

  # --- Core Decomposition --- #
  if (method == "deflation") {
    if (verbose) message(paste0("Using iterative deflation (", ifelse(use_cpp, "C++", "R"), ") to extract ", ncomp, " components..."))
    if (n < p) {
        if (use_cpp) {
          svdfit <- gmd_deflation_cpp(Matrix::t(Xp), A, M, ncomp, threshold)
          svdfit$d <- svdfit$d[,1] # Ensure d is vector
          # C++ might not return propv/cumv, need to calculate if necessary
          if (is.null(svdfit$k)) svdfit$k <- length(svdfit$d)
        } else {
          svdfit <- gmd_deflationR(Matrix::t(Xp), A, M, ncomp, threshold, verbose=verbose)
        }
        # Swap u and v back
        svdfit <- list(u=svdfit$v, v=svdfit$u, d=svdfit$d, k=svdfit$k, cumv=svdfit$cumv, propv=svdfit$propv)
    } else {
        if (use_cpp) {
          svdfit <- gmd_deflation_cpp(Xp, M, A, ncomp, threshold)
          svdfit$d <- svdfit$d[,1]
          if (is.null(svdfit$k)) svdfit$k <- length(svdfit$d)
        } else {
          svdfit <- gmd_deflationR(Xp, M, A, ncomp, threshold, verbose=verbose)
        }
    }
    # Deflation methods should return propv/cumv, but double check
    if (is.null(svdfit$propv) || is.null(svdfit$cumv)) {
       if (verbose) message(" Calculating variance explained for deflation method...")
       total_variance <- sum(Matrix::diag(Matrix::crossprod(Xp, M) %*% Xp %*% A))
       if(total_variance > 1e-8) {
          svdfit$propv <- svdfit$d^2 / total_variance
          svdfit$cumv <- cumsum(svdfit$propv)
       } else {
          svdfit$propv <- rep(0, svdfit$k)
          svdfit$cumv <- rep(0, svdfit$k)
       }
    }

  } else if (method == "eigen") { # One-shot eigen-decomposition approach
    if (verbose) message(paste0("Using one-shot eigen decomposition (gmdLA) to extract ", ncomp, " components..."))
    if (n < p) {
        if (verbose) message(" (n < p, using dual formulation)")
        ret <- gmdLA(Matrix::t(Xp), A, M, k=ncomp, n_orig=p, p_orig=n, maxeig=maxeig, use_dual = TRUE, verbose=verbose)
        # Swap u and v back
        svdfit <- list(u=ret$v, v=ret$u, d=ret$d, k=ret$k, cumv=ret$cumv, propv=ret$propv)
    } else {
        svdfit <- gmdLA(Xp, M, A, k=ncomp, n_orig=n, p_orig=p, maxeig=maxeig, use_dual = FALSE, verbose=verbose)
    }

  } else if (method == "spectra") { # Matrix-free C++/Spectra approach
      if (verbose) message(paste0("Using matrix-free Spectra C++ code to extract ", ncomp, " components..."))
      # Ensure Xp is dense matrix for the C++ function
      Xp_dense <- as.matrix(Xp)
      if (any(!is.finite(Xp_dense))) stop("Input matrix X (after preproc) contains non-finite values.")
      
      # Call the C++ function
      spectra_res <- tryCatch(gmd_fast_cpp(Xp_dense, M, A, k=ncomp, tol=tol_spectra, maxit=maxit_spectra),
                              error = function(e) {stop("Call to gmd_fast_cpp failed: ", e$message)}) 
      
      # Calculate variance explained
      if (verbose) message(" Calculating variance explained for Spectra method...")
      total_variance <- sum(Matrix::diag(Matrix::crossprod(Xp, M) %*% Xp %*% A))
      if (total_variance < 1e-8) {
          propv <- rep(0, spectra_res$k)
          warning("Total generalized variance is near zero.")
      } else {
          propv <- (spectra_res$d^2) / total_variance
      }
      cumv <- cumsum(propv)
      
      # Map results to the svdfit structure
      svdfit <- list(d = spectra_res$d,
                     u = spectra_res$u, # ou
                     v = spectra_res$v, # ov
                     k = spectra_res$k,
                     propv = propv,
                     cumv = cumv)
  } else {
      stop("Internal error: Unknown method specified.") # Should not happen due to match.arg
  }

  # Check how many components were actually found
  k_found <- svdfit$k # gmdLA and gmd_deflationR/cpp should return 'k'
  if (is.null(k_found)) k_found <- length(svdfit$d) # Fallback if k not returned

  if (k_found < ncomp) {
      warning("Found only ", k_found, " valid components, less than requested ncomp=", ncomp)
      ncomp <- k_found
      # Trim results if necessary (should be done in helpers, but ensures consistency)
      if (length(svdfit$d) > ncomp) svdfit$d <- svdfit$d[1:ncomp]
      if (ncol(svdfit$u) > ncomp) svdfit$u <- svdfit$u[, 1:ncomp, drop=FALSE]
      if (ncol(svdfit$v) > ncomp) svdfit$v <- svdfit$v[, 1:ncomp, drop=FALSE]
      if (length(svdfit$propv) > ncomp) svdfit$propv <- svdfit$propv[1:ncomp]
      if (length(svdfit$cumv) > ncomp) svdfit$cumv <- svdfit$cumv[1:ncomp]
  }

  if (ncomp == 0) {
      warning("No valid components found.")
      # Return an empty but valid structure
      return(multivarious::bi_projector(v=matrix(0.0, p, 0), s=matrix(0.0, n, 0), sdev=numeric(0),
                                        preproc=procres, ov=matrix(0.0, p, 0), ou=matrix(0.0, n, 0),
                                        u=matrix(0.0, n, 0), classes="genpca", A=A, M=M,
                                        propv=numeric(0), cumv=numeric(0)))
  }

  # --- Construct bi_projector object --- #
  if (verbose) message("Constructing final object...")

  # Scores: F = X V (where V is ov) or F = M U D (where U is ou)
  # Use F = M U D form for consistency with paper's U definition

 

  M_ou <- M %*% svdfit$u
  # Use sweep for element-wise multiplication D onto columns of M_ou
  scores_mat <- sweep(M_ou, 2, svdfit$d, `*`)

  # Assign row/col names if available from original X
  # Get original indices if preproc modified them
  row_indices <- if (!is.null(attr(Xp, "row_indices"))) attr(Xp, "row_indices") else 1:nrow(X)
  col_indices <- if (!is.null(attr(Xp, "col_indices"))) attr(Xp, "col_indices") else 1:ncol(X)
  
  if (!is.null(rownames(X))) {
      rownames(scores_mat) <- rownames(X)[row_indices]
  } else {
      rownames(scores_mat) <- paste0("Obs", 1:n)
  }
  colnames(scores_mat) <- paste0("PC", 1:ncomp)

  # Loadings (components): v = A V (where V is ov)
  loadings_mat <- A %*% svdfit$v
  if (!is.null(colnames(X))) {
      rownames(loadings_mat) <- colnames(X)[col_indices]
  } else {
      rownames(loadings_mat) <- paste0("Var", 1:p)
  }
  colnames(loadings_mat) <- paste0("PC", 1:ncomp)

  # Create the S3 object using the multivarious constructor
  ret <- multivarious::bi_projector(
    v = loadings_mat,     # Loadings = A %*% ov
    s = scores_mat,       # Scores = M %*% ou %*% D
    sdev = svdfit$d,      # Singular values
    preproc = procres,    # Preprocessing object
    ov = svdfit$v,        # Orthonormal V in A metric
    ou = svdfit$u,        # Orthonormal U in M metric
    u = M_ou,             # U scaled by M metric (M %*% ou)
    classes = "genpca",   # Specific class first
    A = A,                # Store constraint matrices
    M = M,
    propv = svdfit$propv, # Proportion of variance
    cumv = svdfit$cumv   # Cumulative variance
  )

  # bi_projector should handle adding "bi_projector", "projector" classes.

  if (verbose) message("GPCA finished.")
  return(ret)
}




#' @rdname genpca
#' @keywords internal
#' @importFrom Matrix Diagonal t crossprod tcrossprod diag solve isDiagonal Matrix
#' @importFrom RSpectra eigs_sym
#' @importFrom methods as is
#' @importFrom stats eigen
#' @details
#' `gmdLA` caches the eigen decomposition of the constraint matrices by
#' storing it as an attribute on the matrix. `compute_sqrtm()` returns this
#' modified matrix so callers can reassign it (e.g. `R <- sqrtm_res$matrix`)
#' to reuse the cached decomposition in subsequent calls.
gmdLA <- function(X, Q, R, k=min(n_orig, p_orig), n_orig, p_orig,
                  maxeig=800, tol=1e-8, use_dual=FALSE, verbose=FALSE) {

  # Caching key based on object ID might be fragile. Attribute caching is used.
  cache_attr_name <- "eigen_decomp_cache"
  eigen_tol <- tol # Tolerance for filtering eigenvalues

  # Helper to compute M^(1/2) and M^(-1/2) or retrieve from cache. The
  # eigen decomposition is cached on the matrix via an attribute so that
  # repeated calls avoid recomputation. The updated matrix with this
  # attribute is returned alongside the decomposition.
  compute_sqrtm <- function(M, cache_attr, mat_name) {
    if (!is.null(attr(M, cache_attr))) {
      if (verbose) message(paste(" Using cached decomposition for matrix", mat_name))
      decomp <- attr(M, cache_attr)
    } else {
      if (verbose) message(paste(" Computing eigen decomposition for matrix", mat_name))
      if (Matrix::isDiagonal(M)) {
        m_diag <- Matrix::diag(M)
        if (any(m_diag < -eigen_tol)) stop(paste(mat_name, "(diagonal) must be PSD (eigenvalues >= -tol)."))
        m_diag_sqrt <- sqrt(pmax(m_diag, 0)) # Ensure non-negative before sqrt
        m_diag_invsqrt <- ifelse(m_diag > eigen_tol, 1 / m_diag_sqrt, 0) # Avoid division by zero
        decomp <- list(values = m_diag, vectors = NULL,
                       sqrtm = Matrix::Diagonal(x=m_diag_sqrt), 
                       invsqrtm = Matrix::Diagonal(x=m_diag_invsqrt))
      } else {
        # Ensure M is symmetric sparse for eigs_sym or dense for eigen
        M_sym <- if(!Matrix::isSymmetric(M)) Matrix::forceSymmetric(M) else M
        if (!is(M_sym, "sparseMatrix")) M_sym <- Matrix::Matrix(M_sym, sparse=TRUE)
        
        decomp_raw <- tryCatch({
            # Determine max rank safely for RSpectra
            safe_k <- min(maxeig, nrow(M_sym)-1)
            if (safe_k < 1) safe_k <- 1 # Ensure k >= 1
            if (nrow(M_sym) > safe_k && safe_k >= 1) { # Use RSpectra if large and maxeig allows k < dim-1
                if (verbose) message(paste(" (", mat_name, "is large, using RSpectra::eigs_sym with k=", safe_k, ")"))
                RSpectra::eigs_sym(M_sym, k=safe_k, which="LM")
            } else {
                if (verbose) message(paste(" (", mat_name, "using base::eigen)"))
                stats::eigen(as.matrix(M_sym), symmetric=TRUE) # Ensure dense for base::eigen
            }
          },
          error = function(e) {stop(paste("Eigen decomposition failed for", mat_name, ":", e$message))} 
        )
        
        valid_idx <- which(decomp_raw$values > eigen_tol)
        if (length(valid_idx) == 0) stop(paste(mat_name, "has no positive eigenvalues > tol."))
        
        vals <- decomp_raw$values[valid_idx]
        vecs <- decomp_raw$vectors[, valid_idx, drop=FALSE]
        k_actual_decomp <- length(vals)
        if (verbose) message(paste("  (Found ", k_actual_decomp, " eigenvalues > tol for ", mat_name, ")"))
        
        vals_sqrt <- sqrt(vals)
        vals_invsqrt <- 1 / vals_sqrt
        
        # Compute sqrtm and invsqrtm using found components
        vecs_mat <- Matrix::Matrix(vecs)
        sqrtm <- vecs_mat %*% Matrix::Diagonal(x=vals_sqrt) %*% Matrix::t(vecs_mat)
        invsqrtm <- vecs_mat %*% Matrix::Diagonal(x=vals_invsqrt) %*% Matrix::t(vecs_mat)
        
        decomp <- list(values = vals, vectors = vecs,
                       sqrtm = sqrtm, invsqrtm = invsqrtm)
      }
      attr(M, cache_attr) <- decomp # Cache the computed decomposition
    }
    # Return matrix with cached decomposition for reuse
    decomp$matrix <- M
    return(decomp)
  }
  
  # --- Main Logic --- #
  if (!use_dual) { # Primal: n_orig >= p_orig
      if (verbose) message(" gmdLA: Using primal approach (n >= p)")
      R_decomp <- compute_sqrtm(R, paste0(cache_attr_name, "_R"), "R")
      R <- R_decomp$matrix
      Rtilde <- R_decomp$sqrtm
      Rtilde.inv <- R_decomp$invsqrtm

      if (verbose) message("  Calculating X'QX...")
      XQX <- Matrix::crossprod(X, Q) %*% X # p x p matrix

      if (verbose) message("  Calculating R(1/2) X'QX R(1/2)...")
      target_mat <- Rtilde %*% XQX %*% Rtilde # p x p symmetric

      if (verbose) message("  Performing eigen decomposition on target matrix (dim: ", nrow(target_mat), ")...")
      k_request <- min(k, p_orig - 1) # Cannot request more than dim-1 for RSpectra
      if (k_request < 1) stop("k_request must be >= 1 in gmdLA (primal)")
      eig_res <- tryCatch(RSpectra::eigs_sym(target_mat, k=k_request, which="LM"),
                          error=function(e){stop("Eigen decomp failed in gmdLA (primal): ", e$message)})

      valid_idx <- which(eig_res$values > eigen_tol)
      if (length(valid_idx) == 0) stop("No positive eigenvalues found in gmdLA (primal).")
      
      eig_vals <- eig_res$values[valid_idx]
      eig_vecs <- eig_res$vectors[, valid_idx, drop = FALSE]
      k_found <- length(eig_vals)
      if (verbose) message(paste("  (Found ", k_found, " eigenvalues > tol)"))
      
      dgmd <- sqrt(eig_vals)
      
      if (verbose) message("  Calculating ov (V = R^(-1/2) * eigenvectors)...")
      vgmd <- Rtilde.inv %*% eig_vecs # ov (p x k_found)
      
      if (verbose) message("  Calculating ou (U)...")
      # Revert to OLD normalization method for ugmd (as per user request)
      XR <- X %*% R
      # RnR = R %*% (X' Q X) %*% R. XQX already calculated above.
      RnR <- R %*% XQX %*% R
      
      ugmd <- matrix(0.0, n_orig, k_found)
      for (i in 1:k_found) {
          # Use the normalization factor from the OLD gmdLA version
          vgmd_i <- vgmd[, i, drop=FALSE]
          normalizing_number_sq <- Matrix::crossprod(vgmd_i, RnR) %*% vgmd_i
          if (as.numeric(normalizing_number_sq) > tol^2) { # Check squared norm > tol^2
             normalizing_number <- sqrt(as.numeric(normalizing_number_sq))
             ugmd[, i] <- as.vector(XR %*% vgmd_i / normalizing_number)
          } else {
             # Handle near-zero norm case (e.g., set column to zero)
             ugmd[, i] <- 0.0 
             warning("Near-zero norm encountered during old ugmd normalization for component ", i)
          }          
      }
  } else { # Dual: n_orig < p_orig
      if (verbose) message(" gmdLA: Using dual approach (n < p)")
      Q_decomp <- compute_sqrtm(Q, paste0(cache_attr_name, "_Q"), "Q")
      Q <- Q_decomp$matrix
      Qtilde <- Q_decomp$sqrtm
      Qtilde.inv <- Q_decomp$invsqrtm

      if (verbose) message("  Calculating X R X'...")
      XRXt <- X %*% R %*% Matrix::t(X) # n x n matrix

      if (verbose) message("  Calculating Q(1/2) X R X' Q(1/2)...")
      target_mat <- Qtilde %*% XRXt %*% Qtilde # n x n symmetric

      if (verbose) message("  Performing eigen decomposition on target matrix (dim: ", nrow(target_mat), ")...")
      k_request <- min(k, n_orig - 1) # Cannot request more than dim-1 for RSpectra
      if (k_request < 1) stop("k_request must be >= 1 in gmdLA (dual)")
      eig_res <- tryCatch(RSpectra::eigs_sym(target_mat, k=k_request, which="LM"),
                          error=function(e){stop("Eigen decomp failed in gmdLA (dual): ", e$message)})

      valid_idx <- which(eig_res$values > eigen_tol)
      if (length(valid_idx) == 0) stop("No positive eigenvalues found in gmdLA (dual).")

      eig_vals <- eig_res$values[valid_idx]
      eig_vecs <- eig_res$vectors[, valid_idx, drop = FALSE]
      k_found <- length(eig_vals)
      if (verbose) message(paste("  (Found ", k_found, " eigenvalues > tol)"))
      
      dgmd <- sqrt(eig_vals)
      
      if (verbose) message("  Calculating ou (U = Q^(-1/2) * eigenvectors)...")
      ugmd <- Qtilde.inv %*% eig_vecs # ou (n x k_found)
      
      if (verbose) message("  Calculating ov (V)...")
      Xt_ugmd <- Matrix::crossprod(X, ugmd) # p x k_found
      vgmd_unnorm <- sweep(Xt_ugmd, 2, dgmd, `/`) # p x k_found
      
      vgmd <- matrix(0.0, p_orig, k_found)
      for (i in 1:k_found) {
          v_i <- vgmd_unnorm[, i, drop=FALSE]
          norm_factor_sq <- as.numeric(Matrix::crossprod(v_i, R %*% v_i))
          if (norm_factor_sq > tol) {
              vgmd[, i] <- as.vector(v_i / sqrt(norm_factor_sq))
          }
      }
  }
  
  # Calculate explained variance
  if (verbose) message(" Calculating explained variance...")
  # Trace(X' Q X R) is invariant under primal/dual formulation (using args passed to gmdLA)
  total_variance <- tryCatch(sum(Matrix::diag(Matrix::crossprod(X, Q) %*% X %*% R)),
                           error = function(e) {
                              warning("Could not compute total variance trace: ", e$message); NA_real_
                           })
  
  if (is.na(total_variance) || total_variance < tol) {
      propv <- rep(0, k_found)
      warning("Total generalized variance is near zero or could not be computed.")
  } else {
      propv <- (dgmd^2) / total_variance
  }
  cumv <- cumsum(propv)
  
  if (k_found < k) {
      warning("gmdLA: Found only ", k_found, " positive eigenvalues > tol, less than requested k=", k)
  }

  list(
    u = ugmd,       # ou (dimension depends on dual/primal)
    v = vgmd,       # ov (dimension depends on dual/primal)
    d = dgmd,       # singular values (k_found)
    k = k_found,    # Number of components found
    cumv = cumv,    # Cumulative variance explained
    propv = propv   # Proportion variance explained
  )
}




#' @rdname genpca
#' @keywords internal
#' @importFrom Matrix diag crossprod t
#' @importFrom stats rnorm
gmd_deflationR <- function(X, Q, R, k, thr = 1e-6, verbose=FALSE) {

  n = nrow(X)
  p = ncol(X)

  ugmd = matrix(0.0, n, k)
  vgmd = matrix(0.0, p, k)
  dgmd = numeric(k)
  propv = numeric(k)
  Xhat <- X

  # Calculate total variance once: trace(X' Q X R)
  qrnorm <- tryCatch(sum(Matrix::diag( Matrix::crossprod(X, Q %*% X) %*% R )),
                    error = function(e) {
                        warning("Could not compute total variance trace: ", e$message)
                        NA_real_
                    })

  if (is.na(qrnorm) || qrnorm < 1e-10) {
      warning("Total generalized variance is near zero or could not be computed.")
      qrnorm <- 1 # Avoid division by zero, propv will be inaccurate
  }

  k_found <- 0
  for(i in 1:k) {
    if (verbose) message(paste(" Deflation component", i))
    # Initialize u, v for power iteration
    u <- matrix(stats::rnorm(n), ncol=1); u <- u / sqrt(sum(u^2))
    v <- matrix(stats::rnorm(p), ncol=1); v <- v / sqrt(sum(v^2))

    iter <- 0
    max_iter_defl <- 500 # Max iterations for inner power method loop
    converged <- FALSE

    while(iter < max_iter_defl) {
      iter <- iter + 1
      oldu <- u
      oldv <- v

      # Update u: u_hat = Xhat R v; normalize u = u_hat / sqrt(u_hat' Q u_hat)
      uhat <- Xhat %*% (R %*% v)
      u_norm_sq <- Matrix::crossprod(uhat, Q) %*% uhat
      if (u_norm_sq < thr^2) { # Check for near zero norm
          if (verbose) message("  u norm near zero, stopping power iteration for component ", i)
          break
      }
      u <- uhat / sqrt(as.numeric(u_norm_sq))

      # Update v: v_hat = Xhat' Q u; normalize v = v_hat / sqrt(v_hat' R v_hat)
      vhat <- Matrix::crossprod(Xhat, (Q %*% u))
      v_norm_sq <- Matrix::crossprod(vhat, R) %*% vhat
       if (v_norm_sq < thr^2) { # Check for near zero norm
          if (verbose) message("  v norm near zero, stopping power iteration for component ", i)
          break
      }
      v <- vhat / sqrt(as.numeric(v_norm_sq))

      # Check convergence (squared norm difference)
      err <- sum((oldu - u)^2) + sum((oldv - v)^2)
      if (err < thr) {
          converged <- TRUE
          if (verbose) message(paste("  Power iteration converged in", iter, "steps."))
          break
      }
    } # End inner while loop

    if (!converged) {
        warning("Power iteration did not converge for component ", i, " within ", max_iter_defl, " iterations.")
        # If not converged, should we stop? Or continue with the current u,v?
        # Let's stop deflation here if power method fails
        warning("Stopping deflation due to non-convergence of power iteration.")
        break # Exit outer for loop
    }

    # Calculate singular value d = u' Q X_hat R v (use Xhat)
    d_i <- Matrix::crossprod(u, Q) %*% Xhat %*% (R %*% v)
    current_d <- as.numeric(d_i)

    # Check for degenerate component
    if (abs(current_d) < thr) {
        warning("Component ", i, " is degenerate (singular value near zero: ", signif(current_d, 3), "). Stopping deflation.")
        break # Exit outer for loop
    }

    # Store results for this valid component
    k_found <- k_found + 1
    dgmd[k_found] <- current_d
    ugmd[, k_found] <- u[,1]
    vgmd[, k_found] <- v[,1]

    # Deflate Xhat = Xhat - d * u v'
    if (k_found < k) { # No need to deflate after the last requested component
        Xhat <- Xhat - dgmd[k_found] * u %*% Matrix::t(v)
    }

    # Calculate proportion of variance for this component
    propv[k_found] <- dgmd[k_found]^2 / qrnorm

  } # End outer for loop (components)

  if (k_found < k) {
      warning("Deflation stopped early. Found ", k_found, " components instead of requested ", k, ".")
      # Trim result arrays
      dgmd <- dgmd[1:k_found]
      ugmd <- ugmd[, 1:k_found, drop=FALSE]
      vgmd <- vgmd[, 1:k_found, drop=FALSE]
      propv <- propv[1:k_found]
  }

  if (k_found == 0) {
      return(list(d=numeric(0), v=matrix(NA_real_, p, 0), u=matrix(NA_real_, n, 0), k=0, cumv=numeric(0), propv=numeric(0)))
  }

  cumv <- cumsum(propv) # Calculate cumulative sum on valid components

  list(d=as.vector(dgmd), v=vgmd, u=ugmd, k=k_found, cumv=cumv, propv=propv)
}


# S3 method for truncate
#' @rdname genpca
#' @importFrom multivarious ncomp scores loadings sdev bi_projector
#' @export
truncate.genpca <- function(x, ncomp) {
  # Check requested ncomp
  current_ncomp <- multivarious::ncomp(x)
  if (missing(ncomp)) stop("Argument 'ncomp' must be provided.")
  if (!is.numeric(ncomp) || length(ncomp) != 1 || ncomp < 1 || ncomp > current_ncomp || ncomp != floor(ncomp)) {
      stop(paste0("Requested ncomp (", ncomp, ") must be a positive integer <= ", current_ncomp, "."))
  }

  if (ncomp == current_ncomp) return(x) # Nothing to do

  # Use the bi_projector constructor to create the truncated object
  # Select the first 'ncomp' components from relevant slots
  ret <- multivarious::bi_projector(
    v = multivarious::loadings(x)[, 1:ncomp, drop=FALSE], # A ov
    s = multivarious::scores(x)[, 1:ncomp, drop=FALSE],   # M ou D
    sdev = multivarious::sdev(x)[1:ncomp],                # d
    preproc = x$preproc,                    # Preprocessing object
    ov = x$ov[, 1:ncomp, drop=FALSE],       # Orthonormal V
    ou = x$ou[, 1:ncomp, drop=FALSE],       # Orthonormal U
    u = x$u[, 1:ncomp, drop=FALSE],         # M ou
    classes = class(x),                     # Keep original classes ("genpca", "bi_projector", ...)
    A = x$A,                                # Constraint matrix A
    M = x$M,                                # Constraint matrix M
    propv = if (!is.null(x$propv)) x$propv[1:ncomp] else NULL, # Proportion variance
    cumv = if (!is.null(x$cumv)) x$cumv[1:ncomp] else NULL    # Cumulative variance
  )
  return(ret) # Explicitly return the new object
}


# S3 method for reconstruct
#' @rdname genpca
#' @importFrom multivarious ncomp sdev scores loadings reverse_transform
#' @importFrom assertthat assert_that
#' @importFrom Matrix Diagonal t %*%
#' @export
reconstruct.genpca <- function(x,
                               comp = 1:multivarious::ncomp(x),
                               rowind = NULL, # Default to all rows
                               colind = NULL) { # Default to all cols

  max_comp <- multivarious::ncomp(x)
  if (max_comp == 0) return(matrix(0,
                                   nrow = length(rowind %||% 1:nrow(x$M)),
                                   ncol = length(colind %||% 1:nrow(x$A)))) # Return empty if no components
  if (length(comp) == 0 || min(comp) < 1 || max(comp) > max_comp) {
      stop("Selected components 'comp' are out of bounds [1, ", max_comp, "].")
  }

  # Reconstruction uses U D V' where U,V are orthonormal in M,A metrics (ou, ov)
  # X_hat_preproc = ou[, comp] %*% D[comp, comp] %*% t(ov[, comp])

  dvals <- multivarious::sdev(x)[comp]
  # Use Matrix::Diagonal for efficiency
  D_comp <- Matrix::Diagonal(n = length(dvals), x = dvals)

  # Determine effective row/col indices for ou/ov matrices
  eff_rowind <- rowind %||% 1:nrow(x$ou)
  eff_colind <- colind %||% 1:nrow(x$ov)
  
  # Helper function for safe indexing
  safe_index <- function(mat, rows, cols) {
     if (is.null(rows) && is.null(cols)) return(mat)
     if (is.null(rows)) rows <- 1:nrow(mat)
     if (is.null(cols)) cols <- 1:ncol(mat)
     mat[rows, cols, drop=FALSE]
  }

  # Select the specified components and indices for ou and ov
  OU_comp <- safe_index(x$ou, eff_rowind, comp)
  OV_comp <- safe_index(x$ov, eff_colind, comp) # ov rows correspond to X columns

  # Perform the core reconstruction: OU %*% D %*% t(OV)
  # Ensure matrix multiplication handles sparse matrices correctly
  reconstructed_data_preproc <- OU_comp %*% D_comp %*% Matrix::t(OV_comp)

  # Apply inverse pre-processing transform
  final_reconstruction <- multivarious::reverse_transform(x$preproc, reconstructed_data_preproc)

  return(final_reconstruction)
}

# Helper for default NULL indexing
`%||%` <- function(a, b) {
  if (!is.null(a)) a else b
}


# S3 method for ncomp
#' @rdname genpca
#' @export
ncomp.genpca <- function(x) {
  # Number of components is determined by the length of singular values
  # or columns in ou/ov/s/v
  length(multivarious::sdev(x))
}

