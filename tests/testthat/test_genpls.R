## tests/testthat/test-genpls.R

library(testthat)
library(Matrix)
# suppose your package is "myplspackage"
# library(myplspackage)

test_that("genpls handles basic numeric input and constraints", {
  set.seed(123)
  n <- 10
  px <- 4
  py <- 3
  X <- matrix(rnorm(n*px), n, px)
  Y <- matrix(rnorm(n*py), n, py)
  
  # trivial adjacency => identity
  Ax <- Diagonal(px)
  Mx <- Diagonal(n)
  Ay <- Diagonal(py)
  My <- Diagonal(n)
  
  # fit with 2 comps
  fit <- genpls(X, Y, Mx, Ax, My, Ay, ncomp=2, verbose=FALSE)
  expect_s3_class(fit, c("genpls", "cross_projector", "projector"))
  expect_true(ncol(fit$vx) <= 2)
  expect_true(ncol(fit$vy) <= 2)
  expect_equal(multivarious::ncomp(fit), ncol(fit$vx))
  tol <- 1e-9
  nnz <- sum(colSums(abs(fit$tilde_Px)) > tol)
  expect_equal(multivarious::ncomp(fit), nnz)
})

test_that("genpls catches non-numeric or invalid ncomp", {
  X_bad <- data.frame(a=1:5, b=6:10)
  Y_good <- matrix(rnorm(5*3), 5, 3)
  expect_error(genpls(X_bad, Y_good), "numeric matrix")
  
  X_ok <- matrix(rnorm(5*2), 5, 2)
  expect_error(genpls(X_ok, Y_good, ncomp=0), "positive integer")
})

test_that("genpls warns if ncomp > rank", {
  X <- matrix(rnorm(6*2), 6, 2)
  Y <- matrix(rnorm(6*3), 6, 3)
  # ncomp=5 is more than min(6,2)=2
  expect_warning(
    genpls(X, Y, ncomp=5),
    "exceeds approximate rank"
  )
})

test_that("genpls runs and handles degenerate solutions gracefully", {
  set.seed(234)
  # small example
  X <- matrix(rnorm(8*3), 8, 3)
  Y <- X %*% matrix(rnorm(3*2), 3, 2)  # col-lin combo => rank-limited
  # add small noise
  Y <- Y + rnorm(length(Y), sd=1e-8)
  
  fit <- genpls(X, Y, ncomp=3, verbose=TRUE)
  # might see warnings about degenerate or rank-limited
  expect_true(ncol(fit$vx) <= 3)
  expect_equal(multivarious::ncomp(fit), ncol(fit$vx))
  nnz <- sum(colSums(abs(fit$tilde_Px)) > 1e-9)
  expect_equal(multivarious::ncomp(fit), nnz)
})

test_that("genpls works with adjacency-like constraints (row & col)", {
  skip_on_cran()  # optional if you want to skip for CRAN
  
  set.seed(999)
  n <- 12
  px <- 5
  py <- 4
  
  # 1) Generate random data
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)
  
  # 2) Build adjacency-like constraints:
  #    For row constraints (Mx, My), let's create a random sparse adjacency,
  #    then make it PSD by crossprod. 
  #    For col constraints (Ax, Ay), we'll also create small adjacency or diagonal.
  
  # row adjacency for X
  Mx_tmp <- rsparsematrix(n, n, density=0.2)  
  Mx_tmp <- Matrix::crossprod(Mx_tmp)  # ensures PSD-like
  # ensure it's symmetric:
  Mx_tmp <- forceSymmetric(Mx_tmp)
  # turn into a dsCMatrix if not:
  Mx <- as(Mx_tmp, "dsCMatrix")  
  
  # row adjacency for Y
  My_tmp <- rsparsematrix(n, n, density=0.2)
  My_tmp <- Matrix::crossprod(My_tmp)
  My_tmp <- forceSymmetric(My_tmp)
  My <- as(My_tmp, "dsCMatrix")
  
  # col adjacency for X
  Ax_tmp <- rsparsematrix(px, px, density=0.5)
  Ax_tmp <- Matrix::crossprod(Ax_tmp)
  Ax_tmp <- forceSymmetric(Ax_tmp)
  Ax <- as(Ax_tmp, "dsCMatrix")
  
  # col adjacency for Y
  Ay_tmp <- rsparsematrix(py, py, density=0.5)
  Ay_tmp <- Matrix::crossprod(Ay_tmp)
  Ay_tmp <- forceSymmetric(Ay_tmp)
  Ay <- as(Ay_tmp, "dsCMatrix")
  
  # 3) Fit genpls with partial rank 
  #    (small rank to approximate large adjacency).
  #    ncomp=2 => we want two factors.
  
  fit <- genpls(
    X, Y,
    Mx=Mx, Ax=Ax,
    My=My, Ay=Ay,
    ncomp=2,
    # partial rank for each:
    rank_Mx=3, rank_Ax=2, rank_My=3, rank_Ay=2,
    verbose=FALSE
  )
  
  # Basic checks
  expect_s3_class(fit, c("genpls", "cross_projector", "projector"))
  expect_true(ncol(fit$vx) <= 2)  # loadings in original domain
  expect_true(ncol(fit$vy) <= 2)
  
  # 4) Test that we can project X -> latent:
  Zx <- project(fit, X, source="X")
  # Zx should be n x ncomp
  expect_equal(dim(Zx), c(n, 2))
  
  # 5) Test that we can also project Y -> latent:
  Zy <- project(fit, Y, source="Y")
  # Zy should be n x ncomp
  expect_equal(dim(Zy), c(n, 2))
  
  # 6) Transfer: X->Y 
  #    (the shape should be n x q in the target domain)
  Yhat <- transfer(fit, X, source="X", target="Y")
  expect_equal(dim(Yhat), c(n, py))
})

test_that("genpls works with partial-rank adjacency that is diagonal, used as col constraints", {
  skip_on_cran()
  
  set.seed(234)
  n <- 8
  px <- 4
  py <- 3
  
  X <- matrix(rnorm(n * px), n, px)
  Y <- matrix(rnorm(n * py), n, py)
  
  # row adjacency => random PSD for X:
  Mx_tmp <- rsparsematrix(n, n, density=0.3)
  Mx_tmp <- Matrix::crossprod(Mx_tmp)
  Mx_tmp <- forceSymmetric(Mx_tmp)
  Mx <- as(Mx_tmp, "dsCMatrix")
  
  # col adjacency => diagonal for X => partial-eig skipping
  Ax_diagvals <- c(1, 2, 3, 4)
  Ax <- Diagonal(x=Ax_diagvals)
  
  # row adjacency => random PSD for Y:
  My_tmp <- rsparsematrix(n, n, density=0.3)
  My_tmp <- Matrix::crossprod(My_tmp)
  My_tmp <- forceSymmetric(My_tmp)
  My <- as(My_tmp, "dsCMatrix")
  
  # col adjacency => diagonal for Y => partial-eig skip
  Ay_diagvals <- c(5, 6, 7)
  Ay <- Diagonal(x=Ay_diagvals)
  
  # Fit
  fit <- genpls(
    X, Y,
    Mx=Mx, Ax=Ax,
    My=My, Ay=Ay,
    ncomp=2,
    rank_Mx=3, rank_Ax=2,  # won't matter for diag
    rank_My=3, rank_Ay=2,
    verbose=TRUE
  )
  
  expect_s3_class(fit, c("genpls", "cross_projector", "projector"))
  
  # check dimension for loadings
  expect_equal(dim(fit$vx), c(px, 2))
  expect_equal(dim(fit$vy), c(py, 2))
  
  # project => latent
  Zx <- project(fit, X, source="X")
  expect_equal(dim(Zx), c(n, 2))
  # etc.
})


test_that("rpls and genpls give similar results for identity constraints", {
  set.seed(123)
  
  # 1) Generate synthetic data
  n <- 10
  p <- 5
  q <- 3
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  
  # (Optionally) center columns
  X_centered <- scale(X, center = TRUE, scale = FALSE)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)
  
  # 2) Fit rpls with minimal or zero penalty => effectively classical PLS
  #    If your code allows lambda=0, do that. Otherwise, pick something small like 1e-7.
  fit_rpls_obj <- rpls(
    X_centered, Y_centered,
    K       = 2,                # 2 factors
    lambda  = 1e-7,             # minimal ridge penalty
    penalty = "ridge",
    verbose = FALSE
  )
  
  # 3) Fit genpls with identity constraints => Mx=I, Ax=I, My=I, Ay=I
  #    rank_... is irrelevant for identity, but we specify to show no partial-eig needed
  library(Matrix)
  I_n <- Diagonal(n)  # identity for row constraints
  I_p <- Diagonal(p)  # identity for col constraints
  I_q <- Diagonal(q)
  
  fit_gen_obj <- genpls(
    X_centered, Y_centered,
    Mx=I_n, Ax=I_p, 
    My=I_n, Ay=I_q, 
    ncomp=2,
    rank_Mx=2, rank_Ax=2, rank_My=2, rank_Ay=2,
    verbose=FALSE
  )
  
  tmp = data.frame(X=X_centered, Y=Y_centered)
  fit_plsr <- plsr(Y ~ X, data=tmp, ncomp=2, scale=FALSE)
  
  # 4) Compare loadings:
  rpls_vx  <- fit_rpls_obj$vx  # p x K
  genpls_vx<- fit_gen_obj$vx   # p x K
  
  # Typically loadings can differ by a sign or scalar factor. Let's do 
  # a correlation approach or a small difference approach. We'll flatten:
  rpls_vx_vec   <- as.numeric(rpls_vx)
  genpls_vx_vec <- as.numeric(genpls_vx)
  
  # If all zero => skip. Otherwise check correlation:
  if (sd(rpls_vx_vec) < 1e-15 || sd(genpls_vx_vec) < 1e-15) {
    diff_abs <- sum(abs(rpls_vx_vec - genpls_vx_vec))
    expect_lt(diff_abs, 1e-6)
  } else {
    cval <- cor(rpls_vx_vec, genpls_vx_vec)
    # We expect a fairly high correlation if they're effectively the same PLS 
    expect_gt(cval, 0.9) 
  }
  
  # 5) Compare latent scores by projecting X->latent
  #    For rpls: we can do X_centered %*% fit_rpls_obj$vx
  #    For genpls: use project(fit_gen_obj, X_centered, source="X")
  Z_rpls <- X_centered %*% fit_rpls_obj$vx
  Z_gen  <- project(fit_gen_obj, X_centered, source="X")
  
  Z_rpls_vec <- as.numeric(Z_rpls)
  Z_gen_vec  <- as.numeric(Z_gen)
  if (sd(Z_rpls_vec)<1e-15 || sd(Z_gen_vec)<1e-15) {
    diff_scores <- sum(abs(Z_rpls_vec - Z_gen_vec))
    expect_lt(diff_scores, 1e-6)
  } else {
    cval_scores <- cor(Z_rpls_vec, Z_gen_vec)
    expect_gt(cval_scores, 0.9)
  }
})