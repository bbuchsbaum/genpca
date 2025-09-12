#' Generalized PLS-SVD via Implicit Operator (memory-safe)
#'
#' Compute the top-k singular triplets of S = t(Xe) %*% Ye without
#' materializing the whitened matrices Xe = Mx^{1/2} X Wx^{1/2},
#' Ye = My^{1/2} Y Wy^{1/2}. Works with dense/sparse constraints.
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
#' @export
gplssvd_op <- function(X, Y,
                       XLW = NULL, YLW = NULL,
                       XRW = NULL, YRW = NULL,
                       k = 2, center = FALSE, scale = FALSE,
                       svd_backend = c("RSpectra", "irlba"),
                       svd_opts = list(tol = 1e-7, maxitr = 1000)) {

  svd_backend <- match.arg(svd_backend)
  stopifnot(k >= 1)

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

  X <- to_Matrix(X); Y <- to_Matrix(Y)
  N <- nrow(X); I <- ncol(X); J <- ncol(Y)
  stopifnot(nrow(Y) == N)

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
  X <- Xcs$A; Y <- Ycs$A

  # Metric operators
  build_metric <- function(W, n) {
    if (is.null(W)) {
      return(list(
        mult         = function(x) x,
        mult_sqrt    = function(x) x,
        mult_invsqrt = function(x) x
      ))
    }
    if (is.numeric(W) && length(W) == n) {
      d  <- as.numeric(W); d[d < 0] <- 0
      ds <- sqrt(d)
      invds <- ifelse(ds > 0, 1 / ds, 0)
      D  <- Matrix::Diagonal(x = d)
      Ds <- Matrix::Diagonal(x = ds)
      Dis<- Matrix::Diagonal(x = invds)
      return(list(
        mult         = function(x) D  %*% x,
        mult_sqrt    = function(x) Ds %*% x,
        mult_invsqrt = function(x) Dis %*% x
      ))
    }
    if (inherits(W, "diagonalMatrix")) {
      d  <- as.numeric(Matrix::diag(W)); d[d < 0] <- 0
      ds <- sqrt(d)
      invds <- ifelse(ds > 0, 1 / ds, 0)
      Ds <- Matrix::Diagonal(x = ds)
      Dis<- Matrix::Diagonal(x = invds)
      return(list(
        mult         = function(x) W %*% x,
        mult_sqrt    = function(x) Ds %*% x,
        mult_invsqrt = function(x) Dis %*% x
      ))
    }
    # General PSD Matrix -> symmetric sqrt via eigen
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

  MX <- build_metric(XLW, N)
  MY <- build_metric(YLW, N)
  WX <- build_metric(XRW, I)
  WY <- build_metric(YRW, J)

  # Linear operators for S = t(Xe) %*% Ye
  S_mv <- function(v, args = NULL) {
    v2  <- WY$mult_sqrt(v)
    t1  <- Y %*% v2
    t2  <- MY$mult_sqrt(t1)
    t3  <- MX$mult_sqrt(t2)
    t4  <- Matrix::crossprod(X, t3)
    as.matrix(WX$mult_sqrt(t4))
  }
  ST_mv <- function(u, args = NULL) {
    u2  <- WX$mult_sqrt(u)
    s1  <- X %*% u2
    s2  <- MX$mult_sqrt(s1)
    s3  <- MY$mult_sqrt(s2)
    s4  <- Matrix::crossprod(Y, s3)
    as.matrix(WY$mult_sqrt(s4))
  }

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
  Fi <- WX$mult(p %*% diag(d, nrow = length(d)))
  Fj <- WY$mult(q %*% diag(d, nrow = length(d)))
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
