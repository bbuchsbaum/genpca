context("genpca")

mat_10_10 <- matrix(rnorm(10*10), 10, 10)

test_that("ncomp must be integer", {
  expect_error(genpca(mat_10_10, ncomp = 2.5), "single positive integer")
  expect_error(genpca(mat_10_10, ncomp = c(1, 2)), "single positive integer")
})

test_that("pca and genpca have same results with identity matrix for row and column constraints", {
  res1 <- genpca(mat_10_10, preproc=multivarious::center())
  res2 <- multivarious::pca(mat_10_10,ncomp=ncomp(res1), preproc=multivarious::center())
  
  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(sum(diffscores) < 1e-5)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2))
  
  expect_equal(unname(apply(multivarious::components(res1), 2, function(x) sum(x^2))), 
               rep(1, ncomp(res1)))
})

test_that("gen_pca with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=multivarious::center())
  res2 <- multivarious::pca(mat_10_10, preproc=multivarious::standardize())
  
  # Compare absolute scores element-wise with tolerance
  # Scores should match up to sign flips
  expect_equal(abs(multivarious::scores(res1)), abs(multivarious::scores(res2)), tolerance = 1e-6)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), check.attributes=FALSE)
  
})

test_that("gen_pca with use_cpp with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=multivarious::center(), use_cpp=TRUE)
  res2 <- multivarious::pca(mat_10_10, preproc=multivarious::standardize())
  
  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(abs(sum(diffscores)) < 1e-5)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), check.attributes=FALSE)
  
})

test_that("gen_pca with use_cpp (+ deflation) with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=multivarious::center(), ncomp=9, 
  method="deflation", use_cpp=TRUE, threshold=1e-7)
  res2 <- multivarious::pca(mat_10_10, preproc=multivarious::standardize(), ncomp=9)
  
  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), tolerance=.01)
  
})

test_that("gen_pca with use_cpp (+ deflation and n < p) with column variances is equivalent to a scaled pca", {
  mat_10_20 <- matrix(rnorm(10*20), 10, 20)
  wts <- 1/apply(mat_10_20, 2, var)
  res1 <- genpca(mat_10_20, A=wts, preproc=multivarious::center(), ncomp=9, method="deflation", 
  use_cpp=FALSE, threshold=1e-8)
  res2 <- multivarious::pca(mat_10_20, preproc=multivarious::standardize(), ncomp=9)
  
  diffscores <- abs(multivarious::scores(res1)) - abs(multivarious::scores(res2))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(multivarious::sdev(res1), multivarious::sdev(res2), tolerance=.01)
  
})

test_that("gen_pca with dense column and row constraints works", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  res1 <- genpca(mat_10_10, A=A, M=M, preproc=multivarious::center())
  expect_equal(ncomp(res1),length(multivarious::sdev(res1)))
})

test_that("gen_pca with sparse column and row constraints works", {
  A <- neighborweights:::adjacency.neighbor_graph(neighborweights::graph_weights(mat_10_10, k=8))
  Matrix::diag(A) <- 1
  M <- neighborweights:::adjacency.neighbor_graph(neighborweights::graph_weights(t(mat_10_10), k=3))
  Matrix::diag(M) <- 1.5
  res1 <- genpca(mat_10_10, A=A, M=M, preproc=multivarious::center())
})


test_that("can reconstruct a genpca model with component selection", {
  A <- cov(matrix(rnorm(20*10), 20,10))
  M <- cov(matrix(rnorm(20*10), 20,10))
  res1 <- genpca(mat_10_10, preproc=multivarious::center())
  recon1 <- reconstruct(res1)
  expect_equal(as.matrix(recon1), mat_10_10, check.attributes=FALSE)
  
  res1 <- genpca(mat_10_10, A=A, M=M, ncomp=10, preproc=multivarious::center())
  res2 <- multivarious::pca(mat_10_10,ncomp=10, preproc=multivarious::center())
  recon2 <- reconstruct(res2)
  
   
  
})

test_that("can project a row vector", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  
  res1 <- genpca(mat_10_10, A=A, M=M)
  p <- multivarious::project(res1, mat_10_10[1,])
  expect_equal(dim(p), c(1,ncomp(res1)))
})

#test_that("can extract residuals", {
#  res1 <- genpca(mat_10_10)
#  resid <- residuals(res1, ncomp=2, mat_10_10)
#  d <- multivarious::sdev(res1)
#  expect_equal(sum(d[3:length(d)] ^2), sum(resid^2))
#})

test_that("can run genpca with deflation", {
  X <- matrix(rnorm(100),10,10)
  res1 <- genpca(X, preproc=multivarious::center(), ncomp=5,deflation=TRUE)
  res2 <- genpca(X, preproc=multivarious::center(), ncomp=5)
  expect_true(sum(abs(res1$u) - abs(res2$u)) < 1)
})

test_that("can run genpca with sparse weighting matrix", {
  X <- matrix(rnorm(10000*20),10000,20)
  A <- neighborweights::temporal_adjacency(1:20)
  A <- cov(as.matrix(A))
  M <- neighborweights::temporal_adjacency(1:10000)
  res1 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), M=M, preproc=multivarious::center(), ncomp=5,deflation=TRUE)
  res2 <- genpca(X, A=A, M=M, preproc=multivarious::center(), ncomp=5)
  expect_true(!is.null(res1))
})

test_that("can run genpca on a largeish matrix with deflation", {
  nr <- 1000
  nc <- 500
  X <- matrix(rnorm(nr*nc),nr,nc)
  A <- neighborweights::temporal_adjacency(1:nc)
  A <- t(A) %*% A
  
  M <- neighborweights::temporal_adjacency(1:nr)
  M <- t(M) %*% M
  
  res1 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=multivarious::center(), ncomp=5,deflation=TRUE, threshold=1e-5)
  res2 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=multivarious::center(), ncomp=5,deflation=TRUE, 
                 threshold=1e-5, use_cpp=FALSE)
  
  res3 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=multivarious::center(), ncomp=20,deflation=FALSE)
  
  expect_true(!is.null(res1))
})

## tests/testthat/test-genpca-spatial.R

test_that("genpca with spatial adjacency recovers a smooth temporal blob better than no constraints", {
  
  skip_if_not_installed("Matrix")
  library(Matrix)
  
  set.seed(1234)
  
  ## 1) Generate small 2D grid, e.g. 8x8
  nr <- 8
  nc <- 8
  P  <- nr * nc          # total number of pixels
  T  <- 20               # number of time points
  
  ## 2) Construct a ground-truth smooth blob
  ##    We'll define a circular-ish blob in the center
  grid_x <- matrix(rep(1:nr, each=nc), nrow=nr, ncol=nc)
  grid_y <- matrix(rep(1:nc, nr),      nrow=nr, ncol=nc, byrow=TRUE)
  
  center_r <- floor(nr / 2)
  center_c <- floor(nc / 2)
  radius   <- 2.5
  blob     <- exp(-((grid_x - center_r)^2 + (grid_y - center_c)^2) / (2*radius^2))
  
  ## Flatten to length=P
  blob_vec <- as.vector(blob)  # shape (P)
  
  ## 3) Create a time-varying amplitude with a sinusoid
  tt <- seq(0, 2*pi, length.out=T)
  amp <- 2 + sin(tt)   # shape (T)
  
  ## Expand to a T x P matrix
  ## Each row t is the blob * amp[t]
  signal_mat <- outer(amp, blob_vec)  # shape (T x P)
  
  ## 4) Add noise
  noise_mat <- matrix(rnorm(T * P, sd=0.5), T, P)
  X <- signal_mat + noise_mat   # final observed data
  
  rownames(X) <- paste0("Time", seq_len(T))
  
  ## 5) Build a "spatial adjacency" or Laplacian matrix A for the 8x8 grid
  ##    - We'll do adjacency for up/down/left/right neighbors
  ##    - Then we might create a Laplacian from it (D - A, etc.)
  ## Here, let's do adjacency directly: A_ij = 1 if i & j are neighbors
  
  adj_list <- list()
  index_2d_to_1d <- function(r, c) (r-1)*nc + c
  for (r in seq_len(nr)) {
    for (c in seq_len(nc)) {
      cur_id <- index_2d_to_1d(r, c)
      neighbors <- c()
      if (r > 1)         neighbors <- c(neighbors, index_2d_to_1d(r-1, c))
      if (r < nr)        neighbors <- c(neighbors, index_2d_to_1d(r+1, c))
      if (c > 1)         neighbors <- c(neighbors, index_2d_to_1d(r, c-1))
      if (c < nc)        neighbors <- c(neighbors, index_2d_to_1d(r, c+1))
      for (ngb in neighbors) {
        adj_list[[length(adj_list)+1]] <- c(cur_id, ngb)
      }
    }
  }
  
  ## Build a sparse adjacency matrix from these edges
  row_inds <- sapply(adj_list, `[[`, 1)
  col_inds <- sapply(adj_list, `[[`, 2)
  ones     <- rep(1, length(adj_list))
  
  # A_sp is shape P x P
  A_sp <- sparseMatrix(i=row_inds, j=col_inds, x=ones, dims=c(P, P))
  
  ## (Optional) Laplacian = Diag(rowSums(A_sp)) - A_sp
  d_vec <- rowSums(A_sp)
  Lap_sp <- sparseMatrix(i=1:P, j=1:P, x=d_vec) - A_sp
  
  
  ## 6) Run genpca with no constraint  (A=I, M=I)
  # We'll do just 1 component for demonstration
  A_id <- Diagonal(x=rep(1, P))  # identity
  M_id <- Diagonal(x=rep(1, T))  # identity
  # Possibly center or no preproc
  fit_no_constraint <- genpca(X, A=A_id, M=M_id, ncomp=1,
                              preproc=multivarious::center(),
                              deflation=FALSE, use_cpp=FALSE)
  
  ## 7) Run genpca with adjacency or Laplacian
  ## Let's do adjacency or Laplacian. Here we do adjacency:
  fit_adj <- genpca(X, A=A_sp, M=M_id, ncomp=1,
                    preproc=multivarious::center(),
                    deflation=FALSE, use_cpp=FALSE)
  
  ## 8) Compare reconstruction MSE for rank-1 approximation
  #   We'll do Xhat = reconstruct(..., comp=1)
  #   Then measure mean((X - Xhat)^2)
  Xhat_no_constr <- reconstruct(fit_no_constraint, comp=1)
  Xhat_adj       <- reconstruct(fit_adj, comp=1)
  
  mse_no_constr <- mean((X - Xhat_no_constr)^2)
  mse_adj       <- mean((X - Xhat_adj)^2)
  
  ## 9) Expect adjacency to do better => MSE should be smaller
  cat("MSE no constraint:", mse_no_constr, "\n")
  cat("MSE adjacency:    ", mse_adj, "\n")
  
  expect_lt(mse_adj, mse_no_constr * 0.9,
            info="Adjacency-constrained genpca should yield lower MSE than no-constraint for a smooth blob.")
})

############################################################
## Direct tests for gmd_fast_cpp implementation
############################################################

context("gmd_fast_cpp direct tests")

# Helper function to compare subspaces (ignoring column signs)
compare_subspaces <- function(U1, U2, tol = 1e-6) {
  expect_equal(ncol(U1), ncol(U2))
  if (ncol(U1) == 0) return(TRUE)
  
  # Normalize columns just in case
  U1 <- apply(U1, 2, function(x) x / sqrt(sum(x^2)))
  U2 <- apply(U2, 2, function(x) x / sqrt(sum(x^2)))
  
  # Correlation matrix between columns
  corr_mat <- abs(t(U1) %*% U2)
  
  # Check if it's close to a permutation matrix
  # Each row and column should have exactly one element close to 1
  row_max_close_to_1 <- all(abs(apply(corr_mat, 1, max) - 1) < tol)
  col_max_close_to_1 <- all(abs(apply(corr_mat, 2, max) - 1) < tol)
  
  # Check if the sum of squares of correlations is close to the number of components
  sum_sq_corr_close_to_k <- abs(sum(corr_mat^2) - ncol(U1)) < tol * ncol(U1)
  
  expect_true(row_max_close_to_1, info = "Max correlation per row not close to 1")
  expect_true(col_max_close_to_1, info = "Max correlation per col not close to 1")
  expect_true(sum_sq_corr_close_to_k, info = "Sum of squared correlations not close to k")
}

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p <= n, dense constraints", {
  set.seed(1)
  n <- 20
  p <- 15
  k <- 5
  X <- matrix(rnorm(n*p), n, p)
  Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1 # Dense SPD
  R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1 # Dense SPD
  
  # Center data as gmd_fast_cpp assumes centered data
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, deflation=FALSE, preproc=NULL)
  # Directly call C++ function (make sure it's exported properly)
  res_cpp <- gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)
  
  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  compare_subspaces(res_cpp$u, multivarious::scores(res_r), tol = 1e-5)
  compare_subspaces(res_cpp$v, multivarious::components(res_r), tol = 1e-5)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p > n, dense constraints", {
  set.seed(2)
  n <- 15
  p <- 20
  k <- 5
  X <- matrix(rnorm(n*p), n, p)
  Q <- crossprod(matrix(rnorm(n*n), n, n)) + diag(n)*0.1 # Dense SPD
  R <- crossprod(matrix(rnorm(p*p), p, p)) + diag(p)*0.1 # Dense SPD
  
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, deflation=FALSE, preproc=NULL)
  res_cpp <- gmd_fast_cpp(X_centered, Q=Matrix(Q), R=Matrix(R), k=k)
  
  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  compare_subspaces(res_cpp$u, multivarious::scores(res_r), tol = 1e-5)
  compare_subspaces(res_cpp$v, multivarious::components(res_r), tol = 1e-5)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p <= n, sparse constraints", {
  set.seed(3)
  n <- 25
  p <- 20
  k <- 4
  X <- matrix(rnorm(n*p), n, p)
  # Create sparse constraints (e.g., diagonal)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))
  
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, deflation=FALSE, preproc=NULL)
  res_cpp <- gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)
  
  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  compare_subspaces(res_cpp$u, multivarious::scores(res_r), tol = 1e-5)
  compare_subspaces(res_cpp$v, multivarious::components(res_r), tol = 1e-5)
})

test_that("gmd_fast_cpp matches genpca (use_cpp=TRUE) for p > n, sparse constraints", {
  set.seed(4)
  n <- 20
  p <- 25
  k <- 4
  X <- matrix(rnorm(n*p), n, p)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))
  
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, deflation=FALSE, preproc=NULL)
  res_cpp <- gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)
  
  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  compare_subspaces(res_cpp$u, multivarious::scores(res_r), tol = 1e-5)
  compare_subspaces(res_cpp$v, multivarious::components(res_r), tol = 1e-5)
})

test_that("gmd_fast_cpp handles k=1 correctly", {
  set.seed(5)
  n <- 10
  p <- 8
  k <- 1
  X <- matrix(rnorm(n*p), n, p)
  Q <- Diagonal(n, x = runif(n, 0.5, 1.5))
  R <- Diagonal(p, x = runif(p, 0.5, 1.5))
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  res_r <- genpca(X_centered, M=Q, A=R, ncomp=k, use_cpp=TRUE, deflation=FALSE, preproc=NULL)
  res_cpp <- gmd_fast_cpp(X_centered, Q=Q, R=R, k=k)
  
  expect_equal(res_cpp$d, multivarious::sdev(res_r), tolerance = 1e-6)
  expect_equal(res_cpp$k, k)
  expect_equal(ncol(res_cpp$u), k)
  expect_equal(ncol(res_cpp$v), k)
  compare_subspaces(res_cpp$u, multivarious::scores(res_r), tol = 1e-5)
  compare_subspaces(res_cpp$v, multivarious::components(res_r), tol = 1e-5)
})

test_that("gmd_fast_cpp returns fewer components if necessary", {
  # Test case where numerical rank might be less than k requested
  set.seed(6)
  n <- 10
  p <- 10
  k <- 8
  # Create a low-rank matrix
  rank_true <- 5
  X <- matrix(rnorm(n * rank_true), n, rank_true) %*% matrix(rnorm(rank_true * p), rank_true, p)
  Q <- Diagonal(n)
  R <- Diagonal(p)
  X_centered <- scale(X, center=TRUE, scale=FALSE)
  
  # Expect warning when k > rank? Maybe not from C++ directly
  # The C++ code filters eigenvalues close to zero
  res_cpp <- gmd_fast_cpp(X_centered, Q=Q, R=R, k=k, tol = 1e-9)
  
  expect_lte(res_cpp$k, k)
  expect_lte(res_cpp$k, min(n, p))
  # Actual rank might depend on numerical tolerance
  expect_true(res_cpp$k >= rank_true - 1 && res_cpp$k <= rank_true + 1) 
  expect_equal(length(res_cpp$d), res_cpp$k)
  expect_equal(ncol(res_cpp$u), res_cpp$k)
  expect_equal(ncol(res_cpp$v), res_cpp$k)
})

## ────────────────────────────────────────────────────────────────────────────────
##  tests/testthat/test-genpca-advanced.R
## ────────────────────────────────────────────────────────────────────────────────
context("genpca – advanced properties")

set.seed(42)
Xsmall <- matrix(rnorm(50 * 30), 50, 30)      # n > p       (spectra: right‑side)
Xwide  <- matrix(rnorm(40 * 120), 40, 120)    # n < p       (spectra: left‑side)

## -------------------------------------------------------------------------------
test_that("Spectra method matches eigen method on modest problems", {

  skip_if_not_installed("RSpectra")
  skip_if_not_installed("genpca")             # skip if C++ code was not built
  skip_on_cran()                              # heavy-ish

  ## 1) n >= p   (right‑side operator)
  fit_eig  <- genpca(Xsmall, ncomp = 10, method = "eigen",
                     preproc = multivarious::center())
  fit_spc  <- genpca(Xsmall, ncomp = 10, method = "spectra",
                     preproc = multivarious::center(), tol_spectra = 1e-10)

  expect_equal(fit_eig$sdev,           fit_spc$sdev,           tolerance = 1e-6)
  expect_equal(abs(fit_eig$multivarious::scores()),  abs(fit_spc$multivarious::scores()),  tolerance = 1e-5)

  ## 2) n < p    (left‑side operator)
  fit_eig_w <- genpca(Xwide,  ncomp = 15, method = "eigen",
                      preproc = multivarious::center())
  fit_spc_w <- genpca(Xwide,  ncomp = 15, method = "spectra",
                      preproc = multivarious::center(), tol_spectra = 1e-10)

  expect_equal(fit_eig_w$sdev,          fit_spc_w$sdev,          tolerance = 1e-6)
  expect_equal(abs(fit_eig_w$multivarious::scores()), abs(fit_spc_w$multivarious::scores()), tolerance = 1e-5)
})

## -------------------------------------------------------------------------------
test_that("Orthonormality holds in (M,A) metrics", {

  Mrow <- neighborweights::temporal_adjacency(1:nrow(Xsmall))
  Mrow <- t(Mrow) %*% Mrow                 # PSD & dense
  Acol <- cov(matrix(rnorm(ncol(Xsmall)^2), ncol(Xsmall)))

  fit <- genpca(Xsmall, M = Mrow, A = Acol,
                ncomp = 12, preproc = multivarious::center())

  QtU <- t(fit$ou) %*% Mrow %*% fit$ou     # should be ~ I
  RtV <- t(fit$ov) %*% Acol %*% fit$ov     # should be ~ I

  expect_true(max(abs(QtU - diag(ncol(QtU)))) < 1e-8)
  expect_true(max(abs(RtV - diag(ncol(RtV)))) < 1e-8)
})

## -------------------------------------------------------------------------------
test_that("Reconstruction error decreases monotonically with ncomp", {

  X <- matrix(rnorm(120 * 45), 120, 45)
  fit_all <- genpca(X, ncomp = 20, preproc = multivarious::center())

  rss <- vapply(1:20, function(k) {
           recon <- reconstruct(fit_all, comp = 1:k)
           sum((X - recon)^2)
         }, numeric(1))

  ## RSS should strictly decrease until floating‑point noise kicks in
  expect_true(all(diff(rss) <= 1e-12))
})
