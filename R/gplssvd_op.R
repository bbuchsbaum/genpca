#' Generalized PLS-SVD via Implicit Operator (memory-safe)
#'
#' Compute the top-k singular triplets of \eqn{S = Xe' Ye} without
#' materializing the whitened matrices \eqn{Xe = Mx^{1/2} X Wx^{1/2}},
#' \eqn{Ye = My^{1/2} Y Wy^{1/2}}. Works with dense/sparse constraints.
#'
#' @param X n x I matrix (numeric or Matrix)
#' @param Y n x J matrix (numeric or Matrix)
#' @param XLW Row metric for X (M_X): NULL/identity, numeric length-n, diagonalMatrix, or PSD Matrix
#' @param YLW Row metric for Y (M_Y)
#' @param XRW Column metric for X (W_X)
#' @param YRW Column metric for Y (W_Y)
#' @param k  Number of components
#' @param center,scale Logical; pre-center/scale columns of X, Y before metrics
#' @param svd_backend One of "RSpectra" (default) or "irlba"
#' @param svd_opts List of options for the backend (e.g., tol, maxitr)
#' @return list with elements d, u, v, p, q, fi, fj, lx, ly, k, dims, center, scale
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(40 * 6), 40, 6)
#' Y <- matrix(rnorm(40 * 4), 40, 4)
#' op <- gplssvd_op(X, Y, k = 2, center = TRUE)
#' round(op$d, 3)
#' @export
gplssvd_op <- function(X, Y,
                       XLW = NULL, YLW = NULL,
                       XRW = NULL, YRW = NULL,
                       k = 2, center = FALSE, scale = FALSE,
                       svd_backend = c("RSpectra", "irlba"),
                       svd_opts = list(tol = 1e-7, maxitr = 1000)) {

  svd_backend <- match.arg(svd_backend)
  if (!is.numeric(k) || length(k) != 1 || k < 1) {
    stop("k must be a single positive integer >= 1")
  }
  k <- as.integer(k)

  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Matrix package required.")
  if (svd_backend == "RSpectra" && !requireNamespace("RSpectra", quietly = TRUE)) {
    stop("RSpectra package required for svd_backend='RSpectra'.")
  }
  if (svd_backend == "irlba" && !requireNamespace("irlba", quietly = TRUE)) {
    stop("irlba package required for svd_backend='irlba'.")
  }

  to_Matrix <- function(A) {
    if (inherits(A, "Matrix")) A else Matrix::Matrix(A, sparse = FALSE)
  }

  X <- to_Matrix(X)
  Y <- to_Matrix(Y)
  N <- nrow(X)
  I <- ncol(X)
  J <- ncol(Y)
  stopifnot(nrow(Y) == N)

  # Validate k against matrix dimensions
 if (k > min(I, J)) {
    warning("k (", k, ") exceeds min(ncol(X), ncol(Y)) = ", min(I, J),
            "; will return at most ", min(I, J), " components")
    k <- min(I, J)
  }

  # Column center/scale prior to constraints
  cs <- function(A, do_center, do_scale) {
    A <- to_Matrix(A)
    cen <- if (isTRUE(do_center)) Matrix::Matrix(Matrix::colMeans(A), nrow = 1) else NULL
    if (!is.null(cen)) {
      A <- A - Matrix::Matrix(rep(1, nrow(A)), ncol = 1) %*% cen
    }
    scl <- NULL
    if (isTRUE(do_scale)) {
      s <- sqrt(Matrix::colSums(A^2) / pmax(nrow(A) - 1, 1))
      s[s == 0] <- 1
      A <- A %*% Matrix::Diagonal(x = 1 / as.numeric(s))
      scl <- s
    }
    list(A = A, center = if (is.null(cen)) rep(0, ncol(A)) else as.numeric(cen),
         scale = if (is.null(scl)) rep(1, ncol(A)) else as.numeric(scl))
  }

  Xcs <- cs(X, center, scale)
  Ycs <- cs(Y, center, scale)
  X <- Xcs$A
  Y <- Ycs$A

  # Metric operators (shared helper)
  MX <- .metric_operators(XLW, N)
  MY <- .metric_operators(YLW, N)
  WX <- .metric_operators(XRW, I)
  WY <- .metric_operators(YRW, J)

  # Linear operators for S = t(Xe) %*% Ye (shared builder)
  opc <- .build_pls_operator(X, Y, MX, MY, WX, WY)
  S_mv  <- opc$S_mv
  ST_mv <- opc$ST_mv

  # Small-dense fallback for stability on toy sizes
  if (I <= 64 && J <= 64) {
    # Build S explicitly in the same algebra as the implicit operator:
    # S = WX^{1/2} X^T MX^{1/2} MY^{1/2} Y WY^{1/2}
    T1 <- MY$mult_sqrt(Y)              # n x J
    T2 <- MX$mult_sqrt(T1)             # n x J
    T3 <- Matrix::crossprod(X, T2)     # I x J
    T4 <- WX$mult_sqrt(T3)             # I x J
    S  <- t(WY$mult_sqrt(t(T4)))       # I x J (right-mult by WY^{1/2})
    svdS <- svd(as.matrix(S))
    u <- Matrix::Matrix(svdS$u[, seq_len(k), drop = FALSE], sparse = FALSE)
    v <- Matrix::Matrix(svdS$v[, seq_len(k), drop = FALSE], sparse = FALSE)
    d <- svdS$d[seq_len(k)]
  } else if (svd_backend == "RSpectra") {
    sv <- RSpectra::svds(A = S_mv, k = k, nu = k, nv = k,
                         opts = svd_opts, Atrans = ST_mv, dim = c(I, J))
    u <- Matrix::Matrix(sv$u, sparse = FALSE)
    v <- Matrix::Matrix(sv$v, sparse = FALSE)
    d <- sv$d
  } else {
    # irlba operator path via mult() callback. Pass a dummy A for dims.
    A0 <- matrix(0, I, J)
    mult_fun <- function(x, y) {
      # Handle both mult(A, v) and mult(v, A) calling styles
      if (is.matrix(x) || inherits(x, "Matrix")) {
        # x is A, y is vector(s): return A %*% y
        return(S_mv(y))
      } else {
        # x is vector(s), y is A: return t(A) %*% x
        return(ST_mv(x))
      }
    }
    sv <- irlba::irlba(A0,
                       nv = k, nu = k,
                       work = max(3 * k, k + 1),
                       tol = if (!is.null(svd_opts$tol)) svd_opts$tol else 1e-7,
                       maxit = if (!is.null(svd_opts$maxitr)) svd_opts$maxitr else 1000,
                       mult = mult_fun,
                       fastpath = FALSE)
    u <- Matrix::Matrix(sv$u, sparse = FALSE)
    v <- Matrix::Matrix(sv$v, sparse = FALSE)
    d <- sv$d
  }

  # Generalized singular vectors & component scores
  p <- WX$mult_invsqrt(u)
  q <- WY$mult_invsqrt(v)
  # Use Diagonal for efficient scaling by singular values
  Dmat <- Matrix::Diagonal(x = d)
  Fi <- WX$mult(p %*% Dmat)
  Fj <- WY$mult(q %*% Dmat)
  Lx <- MX$mult_sqrt(X %*% WX$mult(p))
  Ly <- MY$mult_sqrt(Y %*% WY$mult(q))

  list(
    d  = as.numeric(d),
    u  = Matrix::Matrix(u, sparse = FALSE),
    v  = Matrix::Matrix(v, sparse = FALSE),
    p  = Matrix::Matrix(p, sparse = FALSE),
    q  = Matrix::Matrix(q, sparse = FALSE),
    fi = Matrix::Matrix(Fi, sparse = FALSE),
    fj = Matrix::Matrix(Fj, sparse = FALSE),
    lx = Matrix::Matrix(Lx, sparse = FALSE),
    ly = Matrix::Matrix(Ly, sparse = FALSE),

    k = k,
    dims = list(N = N, I = I, J = J),
    center = list(X = Xcs$center, Y = Ycs$center),
    scale  = list(X = Xcs$scale,  Y = Ycs$scale)
  )
}
