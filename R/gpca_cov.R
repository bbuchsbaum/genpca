#' Generalized PCA on a covariance matrix
#' 
#' Performs Generalized PCA directly on a pre-computed covariance matrix C with a single 
#' variable-side constraint/metric R. This is useful when you already have C = X'MX 
#' or when X is too large to store but C is manageable. Supports two methods: 
#' "gmd" (Allen et al.'s GMD approach, default) which exactly matches the two-sided 
#' \code{\link{genpca}}, and "geigen" (generalized eigenvalue approach) which solves 
#' C v = lambda R v.
#' 
#' @param C A p x p symmetric positive semi-definite covariance matrix.
#'   Typically C = X'MX where X is the data matrix and M is a row metric.
#' @param R Variable-side constraint/metric. Can be:
#'   \itemize{
#'     \item{NULL: Identity matrix (standard PCA on C)}
#'     \item{A numeric vector of length p: Interpreted as diagonal weights (must be non-negative)}
#'     \item{A p x p symmetric PSD matrix: General metric/smoothing/structure penalties}
#'   }
#' @param ncomp Number of components to return. Default is all positive eigenvalues.
#' @param method Character string specifying the method. One of:
#'   \itemize{
#'     \item{"gmd" (default): Allen et al.'s GMD approach via eigen decomposition of R^{1/2} C R^{1/2}}
#'     \item{"geigen": Generalized eigenvalue approach solving C v = lambda R v}
#'   }
#' @param constraints_remedy How to handle slightly non-PSD inputs (for geigen method). One of:
#'   \itemize{
#'     \item{"error": Stop with an error if constraints are not PSD}
#'     \item{"ridge": Add a small ridge to the diagonal to make PSD}
#'     \item{"clip": Clip negative eigenvalues to zero}
#'     \item{"identity": Replace with identity matrix}
#'   }
#' @param tol Numerical tolerance for PSD checks and filtering small eigenvalues. Default 1e-8.
#' @param verbose Logical. If TRUE, print progress messages. Default FALSE.
#' 
#' @return A list with components:
#'   \describe{
#'     \item{v}{p x k matrix of loadings (R-orthonormal eigenvectors)}
#'     \item{d}{Singular values (square root of eigenvalues lambda)}
#'     \item{lambda}{Eigenvalues (variances under the R-metric)}
#'     \item{k}{Number of components returned}
#'     \item{propv}{Proportion of variance explained by each component}
#'     \item{cumv}{Cumulative proportion of variance explained}
#'     \item{R_rank}{Rank of the constraint matrix R}
#'     \item{method}{The method used ("gmd" or "geigen")}
#'   }
#' 
#' @details
#' \strong{Method Selection Guide:}
#' 
#' Use \code{method = "gmd"} when:
#' \itemize{
#'   \item You need exact equivalence with \code{\link{genpca}(X, M, A)}
#'   \item You're following Allen et al.'s GMD framework
#'   \item You want consistent results with the two-sided decomposition
#' }
#' 
#' Use \code{method = "geigen"} when:
#' \itemize{
#'   \item You specifically need the generalized eigenvalue formulation
#'   \item You're working with legacy code that expects this approach
#'   \item Computational efficiency is critical and R is well-conditioned
#' }
#' 
#' \strong{Method "gmd" (default):}
#' 
#' This method implements Allen et al.'s GMD approach and exactly matches the 
#' two-sided genpca when C = X'MX. It computes the eigendecomposition of 
#' R^{1/2} C R^{1/2} and maps back with V = R^{-1/2} Z, ensuring V'RV = I.
#' The total variance is tr(CR) as in Allen's GPCA (Corollary 5).
#' 
#' \strong{Method "geigen":}
#' 
#' This method solves the generalized eigenproblem C v = lambda R v directly.
#' While mathematically valid, it solves a different optimization than Allen's 
#' GMD and will not, in general, match the two-sided genpca unless R = I or 
#' special commutation conditions hold.
#' 
#' For exact equivalence with genpca(X, M, A), use method="gmd" with C = X'MX and R = A.
#' 
#' @examples
#' # Example 1: Standard PCA on covariance (no constraint)
#' C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
#' fit0 <- genpca_cov(C, R=NULL, ncomp=3)
#' print(fit0$d[1:3])       # first 3 singular values
#' print(fit0$propv[1:3])   # variance explained by first 3 components
#' 
#' # Example 2: Demonstrating equivalence with genpca
#' set.seed(123)
#' X <- matrix(rnorm(50 * 10), 50, 10)
#' M_diag <- runif(50, 0.5, 1.5)  # row weights
#' A_diag <- runif(10, 0.5, 2)    # column weights
#' 
#' # Two-sided GPCA
#' fit_gpca <- genpca(X, M = M_diag, A = A_diag, ncomp = 5,
#'                    preproc = multivarious::pass())
#' 
#' # Equivalent covariance-based GPCA
#' C <- crossprod(X, diag(M_diag) %*% X)  # C = X'MX
#' fit_cov <- genpca_cov(C, R = A_diag, ncomp = 5, method = "gmd")
#' 
#' # These should match exactly
#' all.equal(fit_gpca$sdev, fit_cov$d, tolerance = 1e-10)
#' 
#' # Example 3: Variable weights via a diagonal metric
#' w <- c(1, 1, 0.5, 2)  # emphasize Sepal.Width less, Petal.Width more
#' fitW <- genpca_cov(C, R = w, ncomp=3, method="gmd")
#' print(fitW$d[1:3])
#' 
#' # Example 4: Compare GMD and generalized eigenvalue approaches
#' fit_gmd <- genpca_cov(C, R = w, ncomp=2, method="gmd")
#' fit_geigen <- genpca_cov(C, R = w, ncomp=2, method="geigen")
#' # These will generally differ unless R = I
#' print(paste("GMD singular values:", paste(round(fit_gmd$d, 3), collapse=", ")))
#' print(paste("GEigen singular values:", paste(round(fit_geigen$d, 3), collapse=", ")))
#' 
#' @seealso \code{\link{genpca}} for the standard two-sided GPCA on data matrices,
#'   \code{\link{genpls}} for generalized partial least squares
#' 
#' @references
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014).
#' A Generalized Least-Squares Matrix Decomposition. 
#' Journal of the American Statistical Association, 109(505), 145-159.
#' 
#' @export
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag crossprod
#' @importFrom RSpectra eigs_sym
genpca_cov <- function(C, R = NULL, ncomp = NULL,
                       method = c("gmd", "geigen"),
                       constraints_remedy = c("error","ridge","clip","identity"),
                       tol = 1e-8, verbose = FALSE) {
  
  method <- match.arg(method)
  constraints_remedy <- match.arg(constraints_remedy)
  
  # Dispatch to appropriate implementation
  if (method == "gmd") {
    genpca_cov_gmd(C, R, ncomp, tol, verbose)
  } else {
    genpca_cov_geigen(C, R, ncomp, constraints_remedy, tol, verbose)
  }
}

#' GMD-based covariance GPCA (internal)
#' 
#' Implements Allen et al.'s GMD approach for covariance matrices.
#' Computes eigendecomposition of R^{1/2} C R^{1/2} and maps back.
#' 
#' @keywords internal
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag
genpca_cov_gmd <- function(C, R = NULL, ncomp = NULL, tol = 1e-8, verbose = FALSE) {
  
  # Basic checks & normalization
  stopifnot(is.matrix(C) || inherits(C, "Matrix"))
  p <- nrow(C)
  stopifnot(p == ncol(C))
  if (!inherits(C, "Matrix")) C <- Matrix::Matrix(C, sparse = FALSE)
  if (!Matrix::isSymmetric(C)) C <- Matrix::forceSymmetric(C)
  
  # Column operator R (vector of weights, NULL=I, or PSD matrix)
  if (is.null(R)) {
    R <- Matrix::Diagonal(p)
  } else if (is.vector(R)) {
    stopifnot(length(R) == p)
    if (any(R < -tol)) stop("Negative weights in R.")
    R <- Matrix::Diagonal(p, x = as.numeric(R))
  } else {
    if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)
    if (!Matrix::isSymmetric(R)) R <- Matrix::forceSymmetric(R)
  }
  
  if (verbose) message("Computing eigen factorization of R...")
  
  # Eigen factorization of R to build R^{1/2} and R^{-1/2} on range(R)
  Re <- eigen(as.matrix(R), symmetric = TRUE)
  vals <- pmax(Re$values, 0)
  keep <- which(vals > tol)
  if (length(keep) == 0L) stop("R is (numerically) zero.")
  
  U <- Re$vectors[, keep, drop = FALSE]
  s <- sqrt(vals[keep])
  Rsqrt <- U %*% (s * t(U))                 # R^{1/2}
  Rmhalf <- U %*% ((1/s) * t(U))           # R^{-1/2} on range(R)
  
  if (verbose) message("Computing R^{1/2} C R^{1/2}...")
  
  # Core step per Allen: eigen of B = R^{1/2} C R^{1/2}
  B <- Rsqrt %*% as.matrix(C) %*% Rsqrt
  B <- 0.5 * (B + t(B))  # Ensure symmetry
  
  if (verbose) message("Computing eigendecomposition...")
  
  Ee <- eigen(B, symmetric = TRUE)
  lam_all <- pmax(Ee$values, 0)
  
  # Keep positive eigenvalues
  pos <- which(lam_all > tol)
  if (length(pos) == 0L) stop("No positive eigenvalues in R^{1/2} C R^{1/2}.")
  if (is.null(ncomp)) ncomp <- length(pos)
  ncomp <- min(ncomp, length(pos))
  
  lam <- lam_all[pos][1:ncomp]
  Z <- Ee$vectors[, pos, drop = FALSE][, 1:ncomp, drop = FALSE]
  
  if (verbose) message("Mapping back: V = R^{-1/2} Z...")
  
  # Map back: V = R^{-1/2} Z  (so that V' R V = I)
  V <- Rmhalf %*% Z
  
  # Variance accounting in Allen's GPCA: total = tr(C R)
  total <- sum(Matrix::diag(as.matrix(C) %*% R))
  propv <- as.numeric(lam / ifelse(total > 0, total, 1))
  cumv <- cumsum(propv)
  
  list(
    v      = Matrix::Matrix(V, sparse = FALSE),  # p x k, R-orthonormal
    d      = sqrt(lam),                          # GMD values
    lambda = lam,                                # = d^2
    k      = ncomp,
    propv  = propv,
    cumv   = cumv,
    R_rank = length(keep),
    method = "gmd"
  )
}

#' Generalized eigenvalue-based covariance GPCA (internal)
#' 
#' Solves the generalized eigenproblem C v = lambda R v directly.
#' This is the original implementation that was in gpca.R.
#' 
#' @keywords internal
#' @importFrom Matrix Matrix isSymmetric forceSymmetric Diagonal t diag
#' @importFrom RSpectra eigs_sym
genpca_cov_geigen <- function(C, R = NULL, ncomp = NULL,
                              constraints_remedy = c("error","ridge","clip","identity"),
                              tol = 1e-8, verbose = FALSE) {
  
  constraints_remedy <- match.arg(constraints_remedy)
  
  # --- Basic checks & normalization
  stopifnot(is.matrix(C) || inherits(C, "Matrix"))
  p <- nrow(C)
  stopifnot(p == ncol(C))
  if (!inherits(C, "Matrix")) C <- Matrix::Matrix(C, sparse = FALSE)
  if (!Matrix::isSymmetric(C)) C <- Matrix::forceSymmetric(C)
  
  # Ensure C is (numerically) PSD-ish; clip tiny negatives if needed
  eigC_min <- tryCatch(
    if (p <= 800) min(eigen(as.matrix(C), symmetric=TRUE, only.values=TRUE)$values)
    else RSpectra::eigs_sym(C, k=1, which="SA")$values,
    error = function(e) NA_real_
  )
  if (is.na(eigC_min)) {
    if (verbose) warning("Could not verify PSD of C; proceeding.")
  } else if (eigC_min < -1e-6) {
    warning("C appears non-PSD (min eig: ", signif(eigC_min,4), "). Proceeding but results may be unstable.")
  } else if (eigC_min < 0 && eigC_min > -1e-6) {
    # Symmetrize and small clip
    if (verbose) message("Clipping tiny negative eigenvalues of C.")
    C <- (C + Matrix::t(C)) / 2
  }
  
  # --- Prepare metric R (renamed from G in original)
  if (is.null(R)) {
    R <- Matrix::Diagonal(p)
  } else if (is.vector(R)) {
    if (length(R) != p) stop("Length of weight vector R must match ncol(C).")
    if (any(R < -tol)) stop("Weights in R must be nonnegative (>= -tol).")
    R <- Matrix::Diagonal(p, x = as.numeric(R))
  } else {
    if (!inherits(R, "Matrix")) R <- Matrix::Matrix(R, sparse = FALSE)
    if (!Matrix::isSymmetric(R)) R <- Matrix::forceSymmetric(R)
  }
  
  # --- Remedy for non-PSD R if needed
  ensure_psd <- function(M, name) {
    # quick path for diagonal
    if (Matrix::isDiagonal(M)) {
      d <- Matrix::diag(M)
      minv <- min(d)
      if (minv < -tol) {
        if (verbose) message(name, " diagonal has negatives (min=", signif(minv,4), "). Remedy: ", constraints_remedy)
        if (constraints_remedy == "error") stop(name, " must be PSD.")
        if (constraints_remedy == "ridge") M <- M + Matrix::Diagonal(length(d), x = (tol - minv))
        if (constraints_remedy == "clip")  M <- Matrix::Diagonal(length(d), x = pmax(d, 0))
        if (constraints_remedy == "identity") M <- Matrix::Diagonal(length(d))
      }
      return(M)
    }
    # general case
    minEig <- tryCatch(
      if (nrow(M) <= 800) min(eigen(as.matrix(M), symmetric=TRUE, only.values=TRUE)$values)
      else RSpectra::eigs_sym(M, k=1, which="SA")$values,
      error=function(e) NA_real_
    )
    if (is.na(minEig) || minEig < -tol) {
      if (verbose) message(name, " not PSD (min eig: ", signif(minEig,4), "). Remedy: ", constraints_remedy)
      if (constraints_remedy == "error") stop(name, " must be PSD (>= -tol).")
      if (constraints_remedy == "ridge") {
        ridge <- if (is.na(minEig)) tol else (tol - minEig)
        M <- (M + Matrix::t(M))/2 + Matrix::Diagonal(nrow(M), x = ridge)
      } else if (constraints_remedy == "clip") {
        ed <- eigen(as.matrix(M), symmetric=TRUE)
        ed$values <- pmax(ed$values, 0)
        M <- Matrix::Matrix(ed$vectors %*% (ed$values * t(ed$vectors)), sparse = FALSE)
        M <- (M + Matrix::t(M))/2
      } else if (constraints_remedy == "identity") {
        M <- Matrix::Diagonal(nrow(M))
      }
    } else {
      M <- (M + Matrix::t(M))/2
    }
    M
  }
  R <- ensure_psd(R, "R")
  
  if (verbose) message("Solving generalized eigenproblem C v = lambda R v...")
  
  # --- Solve C v = lambda R v in the range of R (handles semidefinite R)
  # Eigen-decompose R = U diag(gamma) U^T, keep gamma > tol
  ER <- eigen(as.matrix(R), symmetric = TRUE)
  gamma <- pmax(ER$values, 0)
  keep <- which(gamma > tol)
  if (length(keep) == 0L) stop("R is numerically zero; no components can be extracted.")
  U <- ER$vectors[, keep, drop = FALSE]
  gam_sqrt_inv <- 1 / sqrt(gamma[keep])
  Lminushalf <- diag(gam_sqrt_inv)                 # Lambda^{-1/2}
  
  # Work in reduced coordinates: S = Lambda^{-1/2} (U^T C U) Lambda^{-1/2}
  CU <- crossprod(U, as.matrix(C) %*% U)           # U^T C U
  Sred <- Lminushalf %*% CU %*% Lminushalf         # whitened covariance
  ES <- eigen(Sred, symmetric = TRUE)
  
  # Filter positive eigenvalues
  lam_all <- pmax(ES$values, 0)
  pos <- which(lam_all > tol)
  if (length(pos) == 0L) stop("No positive eigenvalues found (within tol).")
  if (is.null(ncomp)) ncomp <- length(pos)
  ncomp <- min(ncomp, length(pos))
  lam <- lam_all[pos][1:ncomp]
  W  <- ES$vectors[, pos, drop = FALSE][, 1:ncomp, drop = FALSE]
  
  # Map back: V = U Lambda^{-1/2} W  (R-orthonormal: V^T R V = I)
  V <- U %*% (Lminushalf %*% W)
  
  # Explained variance: sum(lambda) equals trace(R^{-1/2} C R^{-1/2}) over range(R)
  total_var <- sum(diag(Sred))
  propv <- as.numeric(lam / ifelse(total_var > 0, total_var, 1))
  cumv  <- cumsum(propv)
  
  list(
    v      = Matrix::Matrix(V, sparse = FALSE),  # loadings (p x k), R-orthonormal
    d      = sqrt(lam),                          # singular values (sqrt of variances)
    lambda = lam,                                # variances under the R-metric
    k      = ncomp,                              # number of components returned
    propv  = propv,
    cumv   = cumv,
    R_rank = length(keep),
    method = "geigen"
  )
}