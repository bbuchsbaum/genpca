# Experimental shim: route RSpectra calls through eigencore (dogfood test).
.ec_seed <- 1234L
.ec_with_rng_guard <- function(expr) {
  # eigencore consumes .Random.seed even with seed=; keep caller's stream intact
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (is.null(old)) { if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv()) } else assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  expr
}
.ec_gemm_wrap <- function(f, nrow_out, args) {
  force(f); force(nrow_out); force(args)  # A is rebound by the caller after this closure is made
  function(X, alpha = 1, beta = 0, Y = NULL, ...) {
    X <- as.matrix(X)
    out <- matrix(as.numeric(f(X, args)), nrow = nrow_out, ncol = ncol(X))
    if (!identical(alpha, 1)) out <- alpha * out
    if (!is.null(Y) && !identical(beta, 0)) out <- out + beta * as.matrix(Y)
    out
  }
}
.ec_eigs_sym <- function(A, k, which = "LM", opts = list(), ...) {
  if (methods::is(A, "Matrix") && !methods::is(A, "sparseMatrix")) A <- as.matrix(A)
  a <- list(A, k = k, which = which, seed = .ec_seed)
  if (!is.null(opts$tol)) a$tol <- opts$tol
  if (!is.null(opts$maxitr)) a$maxit <- opts$maxitr
  .ec_with_rng_guard(do.call(eigencore::eigs_sym, c(a, list(...))))
}
.ec_svds <- function(A, k, nu = k, nv = k, opts = list(), Atrans = NULL, dim = NULL, args = NULL, ...) {
  if (is.function(A)) {
    stopifnot(is.function(Atrans), length(dim) == 2L)
    A <- eigencore::linear_operator(dim = dim,
                                    apply = .ec_gemm_wrap(A, dim[1], args),
                                    apply_adjoint = .ec_gemm_wrap(Atrans, dim[2], args))
  } else if (methods::is(A, "Matrix") && !methods::is(A, "sparseMatrix")) {
    A <- as.matrix(A)
  }
  a <- list(A, k = k, nu = nu, nv = nv, seed = .ec_seed)
  if (!is.null(opts$tol)) a$tol <- opts$tol
  .ec_with_rng_guard(do.call(eigencore::svds, c(a, list(...))))
}
