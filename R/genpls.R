#' Generalized PLS via Implicit Operator (PLS-SVD / GPLSSVD)
#'
#' Canonical (two-block) generalized PLS using sparse-friendly implicit
#' matrix–vector products. Solves the SVD of the operator S = t(Xe) %*% Ye
#' without materializing Xe = Mx^{1/2} X Ax^{1/2} or Ye = My^{1/2} Y Ay^{1/2}.
#'
#' This follows the GPLSSVD/PLS-SVD formulation (Beaton, eqs. 10–14):
#' the top `ncomp` singular triplets of S are computed by iterative SVD
#' on the linear maps v -> S v and u -> S^T u, implemented with metric
#' Cholesky multiplies/solves when possible. Works with dense or sparse
#' `Matrix` inputs and constraint metrics.
#'
#' Returns a `multivarious::cross_projector` with X-/Y-weights (vx, vy)
#' chosen to provide natural projection of new data (`X %*% vx`, `Y %*% vy`).
#' Additional GPLSSVD quantities are attached to the object for access:
#' singular values `d`, generalized weights `p`, `q`, variable scores `fi`, `fj`,
#' and row latent variables `lx`, `ly`.
#'
#' @param X Numeric or Matrix, n x p.
#' @param Y Numeric or Matrix, n x q. Must have same n as `X`.
#' @param Ax Column metric for X (W_X): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param Ay Column metric for Y (W_Y): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param Mx Row metric for X (M_X): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param My Row metric for Y (M_Y): vector/diagonal/matrix; `NULL` ⇒ identity.
#' @param ncomp Number of components to extract (rank-k). Default 2.
#' @param preproc_x,preproc_y Optional `multivarious` preprocessors (e.g., `center()`).
#'   Defaults to `multivarious::pass()` (no-op).
#' @param svd_backend Character, one of `"RSpectra"` (default) or `"irlba"` for
#'   iterative SVD. If neither backend is available, a dense fallback is used
#'   for small problems by materializing S.
#' @param svd_opts List of options passed to the SVD backend, e.g., `tol`, `maxitr`.
#' @param verbose Logical; print brief progress messages.
#'
#' @return An object of class `c("genpls", "cross_projector", "projector")` with:
#'   - vx, vy: X- and Y- weights usable with `predict/transfer` (stored in cross_projector)
#'   - d:     singular values (attached field)
#'   - p, q:  generalized weights W_X^{-1/2} u, W_Y^{-1/2} v (attached)
#'   - fi,fj: variable/component scores W_X p diag(d), W_Y q diag(d) (attached)
#'   - lx,ly: row latent variables M_X^{1/2} X W_X p, M_Y^{1/2} Y W_Y q (attached)
#'   - metrics: the supplied metrics (attached)
#'
#' @examples
#' if (requireNamespace("RSpectra", quietly = TRUE) &&
#'     requireNamespace("multivarious", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 100; p <- 40; q <- 30
#'   X <- matrix(rnorm(n*p), n, p)
#'   Y <- matrix(rnorm(n*q), n, q)
#'   w <- runif(n); w <- w/sum(w)
#'   Mx <- My <- Matrix::Diagonal(x = w)
#'   fit <- genpls(X, Y, Mx = Mx, My = My, ncomp = 2,
#'                 preproc_x = multivarious::center(),
#'                 preproc_y = multivarious::center())
#'   fit$d  # singular values
#' }
#'
#' @references
#' Beaton, Dougal. Generalized eigen, singular value, and partial least squares
#' decompositions: The GSVD package. (Eqs. 10–14). 2020.
#'
#' @importFrom Matrix Matrix Diagonal crossprod t forceSymmetric Cholesky solve
#' @importFrom RSpectra svds
#' @importFrom multivarious cross_projector prep init_transform pass
#' @export
genpls <- function(X, Y,
                   Ax = NULL, Ay = NULL,
                   Mx = NULL, My = NULL,
                   ncomp = 2,
                   preproc_x = multivarious::pass(),
                   preproc_y = multivarious::pass(),
                   svd_backend = c("RSpectra", "irlba"),
                   svd_opts = list(tol = 1e-7, maxitr = 1000),
                   verbose = FALSE) {

  svd_backend <- match.arg(svd_backend)
  stopifnot(length(ncomp) == 1L, ncomp >= 1)

  to_Matrix <- function(A) {
    if (inherits(A, "Matrix")) A else Matrix::Matrix(A, sparse = FALSE)
  }

  X <- to_Matrix(X); Y <- to_Matrix(Y)
  n <- nrow(X); px <- ncol(X); py <- ncol(Y)
  if (nrow(Y) != n) stop("X and Y must have the same number of rows.")

  # Preprocess (e.g., centering/scaling) via multivarious
  proc_x <- multivarious::prep(preproc_x)
  proc_y <- multivarious::prep(preproc_y)
  # multivarious expects base matrices
  Xp <- multivarious::init_transform(proc_x, as.matrix(X))
  Yp <- multivarious::init_transform(proc_y, as.matrix(Y))
  # Convert back to Matrix for efficient ops / sparsity
  Xp <- to_Matrix(Xp)
  Yp <- to_Matrix(Yp)

  # Metric builders: identity/diag/symmetric-PSD (sparse-friendly)
  # For general SPD, use sparse Cholesky: if W = R^T R (upper),
  # W^{1/2} x = R^T x; W^{-1/2} x = solve(R, x)
  build_metric <- function(W, dim_expected, name) {
    if (is.null(W)) {
      return(list(mult = function(x) x,
                  mult_sqrt = function(x) x,
                  mult_invsqrt = function(x) x))
    }
    if (is.numeric(W) && length(W) == dim_expected) {
      d <- as.numeric(W); d[d < 0] <- 0
      ds <- sqrt(d); invds <- ifelse(ds > 0, 1/ds, 0)
      D  <- Matrix::Diagonal(x = d)
      Ds <- Matrix::Diagonal(x = ds)
      Dis<- Matrix::Diagonal(x = invds)
      return(list(mult = function(x) D %*% x,
                  mult_sqrt = function(x) Ds %*% x,
                  mult_invsqrt = function(x) Dis %*% x))
    }
    if (inherits(W, "diagonalMatrix")) {
      d <- as.numeric(Matrix::diag(W)); d[d < 0] <- 0
      ds <- sqrt(d); invds <- ifelse(ds > 0, 1/ds, 0)
      Ds <- Matrix::Diagonal(x = ds)
      Dis<- Matrix::Diagonal(x = invds)
      return(list(mult = function(x) W %*% x,
                  mult_sqrt = function(x) Ds %*% x,
                  mult_invsqrt = function(x) Dis %*% x))
    }
    # General symmetric PSD -> use symmetric sqrt via eigen for robustness
    W <- to_Matrix(W)
    W <- Matrix::forceSymmetric(W, uplo = "U")
    if (exists("ensure_spd", mode = "function")) W <- ensure_spd(W)
    Wd <- as.matrix(W)
    es <- eigen(Wd, symmetric = TRUE)
    lam <- pmax(es$values, 0)
    Q   <- es$vectors
    list(
      mult = function(x) W %*% x,
      mult_sqrt = function(x) {
        X <- as.matrix(x)
        alpha <- crossprod(Q, X)
        if (length(lam) > 0) alpha <- diag(sqrt(lam), nrow = length(lam)) %*% alpha
        Matrix::Matrix(Q %*% alpha, sparse = FALSE)
      },
      mult_invsqrt = function(x) {
        X <- as.matrix(x)
        alpha <- crossprod(Q, X)
        invs <- ifelse(lam > 0, 1 / sqrt(lam), 0)
        if (length(invs) > 0) alpha <- diag(invs, nrow = length(invs)) %*% alpha
        Matrix::Matrix(Q %*% alpha, sparse = FALSE)
      }
    )
  }

  WX <- build_metric(Ax, px, "Ax")
  WY <- build_metric(Ay, py, "Ay")
  MX <- build_metric(Mx, n,  "Mx")
  MY <- build_metric(My, n,  "My")

  # Implicit operator for S = t(Xe) %*% Ye, where
  # Xe = Mx^{1/2} X Ax^{1/2},  Ye = My^{1/2} Y Ay^{1/2}
  S_mv <- function(v, args = NULL) {
    v2 <- WY$mult_sqrt(v)             # J x m
    t1 <- Yp %*% v2                   # n x m
    t2 <- MY$mult_sqrt(t1)            # n x m
    t3 <- MX$mult_sqrt(t2)            # n x m
    t4 <- Matrix::crossprod(Xp, t3)   # p x m
    as.matrix(WX$mult_sqrt(t4))       # p x m (base matrix for RSpectra)
  }
  ST_mv <- function(u, args = NULL) {
    u2 <- WX$mult_sqrt(u)             # p x m
    s1 <- Xp %*% u2                   # n x m
    s2 <- MX$mult_sqrt(s1)            # n x m
    s3 <- MY$mult_sqrt(s2)            # n x m
    s4 <- Matrix::crossprod(Yp, s3)   # q x m
    as.matrix(WY$mult_sqrt(s4))       # q x m
  }

  if (verbose) message("Computing top-", ncomp, " GPLSSVD components (", svd_backend, ") ...")

  # Compute top-ncomp SVD of operator
  u <- v <- NULL; d <- rep(NA_real_, ncomp)
  if (svd_backend == "RSpectra" && requireNamespace("RSpectra", quietly = TRUE)) {
    sv <- RSpectra::svds(A = S_mv, k = ncomp, nu = ncomp, nv = ncomp,
                         opts = svd_opts, Atrans = ST_mv, dim = c(px, py))
    u <- Matrix::Matrix(sv$u, sparse = FALSE)
    v <- Matrix::Matrix(sv$v, sparse = FALSE)
    d <- sv$d
  } else if (svd_backend == "irlba" && requireNamespace("irlba", quietly = TRUE)) {
    # Use irlba with a mult() callback and a dummy matrix for dims
    A0 <- matrix(0, px, py)
    mult_fun <- function(x, y) {
      if (is.matrix(x) || inherits(x, "Matrix")) {
        S_mv(y)
      } else {
        ST_mv(x)
      }
    }
    sv <- irlba::irlba(A0,
                       nv = ncomp, nu = ncomp,
                       work = max(3 * ncomp, ncomp + 1),
                       tol = if (!is.null(svd_opts$tol)) svd_opts$tol else 1e-7,
                       maxit = if (!is.null(svd_opts$maxitr)) svd_opts$maxitr else 1000,
                       mult = mult_fun,
                       fastpath = FALSE)
    u <- Matrix::Matrix(sv$u, sparse = FALSE)
    v <- Matrix::Matrix(sv$v, sparse = FALSE)
    d <- sv$d
  } else {
    # Dense fallback for small problems: materialize S
    if (verbose) message("Falling back to dense SVD by forming S (small problems only).")
    # Form Xe and Ye explicitly for fallback
    Xe <- MX$mult_sqrt(Xp %*% WX$mult_sqrt(Matrix::Diagonal(px, x = 1)))
    Ye <- MY$mult_sqrt(Yp %*% WY$mult_sqrt(Matrix::Diagonal(py, x = 1)))
    S <- Matrix::crossprod(Xe, Ye)                          # p x q
    svdS <- svd(as.matrix(S), nu = ncomp, nv = ncomp)
    u <- Matrix::Matrix(svdS$u, sparse = FALSE)
    v <- Matrix::Matrix(svdS$v, sparse = FALSE)
    d <- svdS$d[seq_len(ncomp)]
  }

  # Generalized singular vectors and scores
  p <- WX$mult_invsqrt(u)                         # W_X^{-1/2} u
  q <- WY$mult_invsqrt(v)                         # W_Y^{-1/2} v
  Fi <- WX$mult(p %*% diag(d, nrow = length(d)))  # W_X p diag(d)
  Fj <- WY$mult(q %*% diag(d, nrow = length(d)))  # W_Y q diag(d)
  Lx <- MX$mult_sqrt(Xp %*% WX$mult(p))           # M_X^{1/2} X W_X p
  Ly <- MY$mult_sqrt(Yp %*% WY$mult(q))           # M_Y^{1/2} Y W_Y q

  # Build a cross_projector for multivarious ecosystem
  # Use vx = W_X p and vy = W_Y q so that scores_x = X %*% vx, scores_y = Y %*% vy
  vx <- as.matrix(WX$mult(p))
  vy <- as.matrix(WY$mult(q))

  obj <- multivarious::cross_projector(
    vx = vx,
    vy = vy,
    preproc_x = proc_x,
    preproc_y = proc_y,
    classes = "genpls"
  )

  # Attach GPLSSVD details for direct access
  obj$d  <- as.numeric(d)
  obj$p  <- as.matrix(p)
  obj$q  <- as.matrix(q)
  obj$fi <- as.matrix(Fi)
  obj$fj <- as.matrix(Fj)
  obj$lx <- as.matrix(Lx)
  obj$ly <- as.matrix(Ly)
  obj$metrics <- list(Ax = Ax, Ay = Ay, Mx = Mx, My = My)
  obj$ncomp   <- ncomp
  obj$backend <- svd_backend

  if (verbose) message("genpls finished.")
  obj
}
