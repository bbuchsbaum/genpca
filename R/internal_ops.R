#' Internal metric and operator helpers
#'
#' These utilities are used across `genpls()` and `gplssvd_op()` to construct
#' metric-aware linear operators for GPLSSVD/PLS-SVD without duplicating logic.
#'
#' @details
#' - `.metric_operators(W, n_expected)` accepts a symmetric positive (semi)definite
#'   matrix or a diagonal specification and returns a list of closures that apply
#'   `W`, `W^{1/2}`, and `W^{-1/2}` to vectors/matrices. Identity and diagonal
#'   cases use fast paths; general SPD uses an eigen-based symmetric square root.
#' - `.build_pls_operator(X, Y, MX, MY, WX, WY)` builds the implicit operator
#'   `v -> WX^{1/2} X^T MX^{1/2} MY^{1/2} Y WY^{1/2} v` and its adjoint.
#'
#' These helpers assume inputs are conformable and (when required) symmetric PSD;
#' callers are responsible for any validation or coercion.
#'
#' @keywords internal
#' @noRd

# Build metric operators for a symmetric (PS) matrix W:
# returns closures: mult (W %*% x), mult_sqrt (W^{1/2} x), mult_invsqrt (W^{-1/2} x)
.metric_operators <- function(W, n_expected = NULL) {
  to_Matrix <- function(A) if (inherits(A, "Matrix")) A else Matrix::Matrix(A, sparse = FALSE)

  # Identity
  if (is.null(W)) {
    return(list(
      mult = function(x) x,
      mult_sqrt = function(x) x,
      mult_invsqrt = function(x) x
    ))
  }

  # Diagonal from numeric
  if (is.numeric(W) && !is.null(n_expected) && length(W) == n_expected) {
    d <- as.numeric(W)
    d[d < 0] <- 0
    ds <- sqrt(d)
    invds <- ifelse(ds > 0, 1 / ds, 0)
    D  <- Matrix::Diagonal(x = d)
    Ds <- Matrix::Diagonal(x = ds)
    Dis <- Matrix::Diagonal(x = invds)
    return(list(
      mult = function(x) D %*% x,
      mult_sqrt = function(x) Ds %*% x,
      mult_invsqrt = function(x) Dis %*% x
    ))
  }

  # DiagonalMatrix
  if (inherits(W, "diagonalMatrix")) {
    d  <- as.numeric(Matrix::diag(W))
    d[d < 0] <- 0
    ds <- sqrt(d)
    invds <- ifelse(ds > 0, 1 / ds, 0)
    Ds <- Matrix::Diagonal(x = ds)
    Dis <- Matrix::Diagonal(x = invds)
    return(list(
      mult = function(x) W %*% x,
      mult_sqrt = function(x) Ds %*% x,
      mult_invsqrt = function(x) Dis %*% x
    ))
  }

  # General symmetric PSD -> symmetric sqrt via eigen
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

#' Build PLS operator closures (internal)
#'
#' @param X Matrix `n x p`
#' @param Y Matrix `n x q`
#' @param MX,MY,WX,WY Lists of metric closures with elements `mult`, `mult_sqrt`, `mult_invsqrt`
#' @return A list with `S_mv`, `ST_mv`, and `dims = c(p, q)`
#' @keywords internal
#' @noRd
#'
# Build PLS operator closures S_mv and ST_mv for given data and metric operators
.build_pls_operator <- function(X, Y, MX, MY, WX, WY) {
  # X: n x p, Y: n x q; MX,MY,WX,WY are metric operator lists (mult/mult_sqrt/...)
  p <- ncol(X)
  q <- ncol(Y)
  S_mv <- function(v, args = NULL) {
    v2 <- WY$mult_sqrt(v)
    t1 <- Y %*% v2
    t2 <- MY$mult_sqrt(t1)
    t3 <- MX$mult_sqrt(t2)
    t4 <- Matrix::crossprod(X, t3)
    as.matrix(WX$mult_sqrt(t4))
  }
  ST_mv <- function(u, args = NULL) {
    u2 <- WX$mult_sqrt(u)
    s1 <- X %*% u2
    s2 <- MX$mult_sqrt(s1)
    s3 <- MY$mult_sqrt(s2)
    s4 <- Matrix::crossprod(Y, s3)
    as.matrix(WY$mult_sqrt(s4))
  }
  list(S_mv = S_mv, ST_mv = ST_mv, dims = c(p, q))
}
