## tests/testthat/test_genplscor.R

library(testthat)
library(Matrix)
library(multivarious)

test_that("genplscor works with identity constraints (basic usage)", {
  set.seed(123)
  n <- 10
  p <- 5
  q <- 3
  
  # Generate random data X, Y
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  
  # Typically, we center for PLS
  Xc <- scale(X, center=TRUE, scale=FALSE)
  Yc <- scale(Y, center=TRUE, scale=FALSE)
  
  # Identity constraints => Mx=I_n, Ax=I_p, My=I_n, Ay=I_q
  # We'll do no constraints => default is NULL => identity
  fit_id <- genplsc(Xc, Yc, ncomp=2, verbose=FALSE)
  
  expect_s3_class(fit_id, c("genplscorr","cross_projector","projector"))
  # check dimension of loadings
  expect_equal(dim(fit_id$vx), c(p, 2))
  expect_equal(dim(fit_id$vy), c(q, 2))
  
  # project X -> latent => dimension n x ncomp
  Zx <- project(fit_id, Xc, source="X")
  expect_equal(dim(Zx), c(n, 2))
  
  # project Y -> latent => dimension n x ncomp
  Zy <- project(fit_id, Yc, source="Y")
  expect_equal(dim(Zy), c(n, 2))
  
  # transfer X->Y => dimension n x q
  Yhat <- transfer(fit_id, Xc, source="X", target="Y")
  expect_equal(dim(Yhat), c(n, q))
})


test_that("genplscor works with diagonal col constraints for X, Y", {
  set.seed(999)
  n <- 12
  p <- 4
  q <- 3
  
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*q), n, q)
  
  # Let's define col constraints: Ax=diag(1, p) => identity? Let's do a random diag
  Ax_vals <- c(1.0, 2.5, 0.5, 4.0) # just random positive
  Ax <- Diagonal(x=Ax_vals)
  # Similarly for Y
  Ay_vals <- c(2.0, 1.5, 3.2)
  Ay <- Diagonal(x=Ay_vals)
  
  # row constraints => identity
  # pass them as NULL => identity
  fit_diag <- genplsc(X, Y, Ax=Ax, Ay=Ay, ncomp=2, verbose=FALSE)
  
  expect_s3_class(fit_diag, c("genplscorr","cross_projector","projector"))
  # dimension checks
  expect_equal(dim(fit_diag$vx), c(p, 2))
  expect_equal(dim(fit_diag$vy), c(q, 2))
  
  # project
  Zx <- project(fit_diag, X, source="X")
  expect_equal(dim(Zx), c(n, 2))
  
  Zy <- project(fit_diag, Y, source="Y")
  expect_equal(dim(Zy), c(n, 2))
})

test_that("genplscor with random PSD row adjacency constraints for X, Y", {
  skip_on_cran()  # skip for CRAN if bigger
  
  set.seed(234)
  n <- 15
  p <- 6
  q <- 5
  
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*q), n, q)
  
  # row constraints => random adjacency => crossprod => PSD
  Mx_tmp <- rsparsematrix(n, n, density=0.2)
  Mx_tmp <- crossprod(Mx_tmp)
  Mx_tmp <- forceSymmetric(Mx_tmp)
  Mx <- as(Mx_tmp, "dsCMatrix")
  
  My_tmp <- rsparsematrix(n, n, density=0.2)
  My_tmp <- crossprod(My_tmp)
  My_tmp <- forceSymmetric(My_tmp)
  My <- as(My_tmp, "dsCMatrix")
  
  # col constraints => identity => pass as NULL
  fit_adj <- genplsc(X, Y, Mx=Mx, My=My, ncomp=3, verbose=FALSE)
  
  expect_s3_class(fit_adj, c("genplscorr","cross_projector","projector"))
  expect_equal(dim(fit_adj$vx), c(p, 3))
  expect_equal(dim(fit_adj$vy), c(q, 3))
  
  # project => see if we can do it
  Zx <- project(fit_adj, X, source="X")
  Zy <- project(fit_adj, Y, source="Y")
  expect_equal(dim(Zx), c(n, 3))
  expect_equal(dim(Zy), c(n, 3))
})

test_that("genplscor handles partial rank (ncomp < min(p,q))", {
  skip_on_cran()
  
  set.seed(345)
  n <- 10
  p <- 6
  q <- 7
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*q), n, q)
  # # no constraints => identity
  
  # ncomp=2 < min(p,q)=6 or 7
  fit_prank <- genplsc(X, Y, ncomp=2, verbose=TRUE)
  
  expect_s3_class(fit_prank, c("genplscorr","cross_projector","projector"))
  expect_true(ncol(fit_prank$vx) <= 2)
  expect_true(ncol(fit_prank$vy) <= 2)
  Zx <- project(fit_prank, X, source="X")
  Zy <- project(fit_prank, Y, source="Y")
  expect_equal(dim(Zx), c(n, 2))
  expect_equal(dim(Zy), c(n, 2))
})

test_that("genplscor internal: correlation of loadings with direct SVD approach", {
  skip_on_cran()
  
  set.seed(456)
  n <- 8
  p <- 4
  q <- 3
  X <- matrix(rnorm(n*p), n, p)
  Y <- matrix(rnorm(n*q), n, q)
  
  # let us define a quick "classical PLS correlation" approach:
  # 1) center => Xc, Yc
  # 2) crossprod => do SVD => top ncomp
  Xc <- scale(X, center=TRUE, scale=FALSE)
  Yc <- scale(Y, center=TRUE, scale=FALSE)
  # direct "pls correlation" => Mx=I,Ax=I, My=I,Ay=I
  # => genplscor
  fit_gcor <- genplsc(Xc, Yc, ncomp=2, verbose=FALSE)
  
  # as a naive check, let's do crossprod(Xc, Yc) => SVD => top2
  C_mat <- crossprod(Xc, Yc) # shape (p x q)
  sres <- svd(C_mat)
  # top 2
  d_naive <- sres$d[1:2]
  u_naive <- sres$u[,1:2, drop=FALSE]
  v_naive <- sres$v[,1:2, drop=FALSE]
  
  # compare fit_gcor$vx => p x 2 vs u_naive in correlation sense
  # because they can differ by sign or scale
  # flatten & cor
  vx_vec <- as.numeric(fit_gcor$vx)
  u_vec  <- as.numeric(u_naive)
  if(sd(vx_vec)<1e-15 || sd(u_vec)<1e-15) {
    # check difference
    dval <- sum(abs(vx_vec - u_vec))
    expect_lt(dval, 1e-6)
  } else {
    cor_val <- cor(vx_vec, u_vec)
    # we expect correlation > 0.9 or so
    expect_gt(cor_val, 0.9)
  }
})


test_that("genplscorr dual=FALSE vs dual=TRUE yield similar results on small data", {
  set.seed(123)
  n <- 8
  p <- 5
  q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  
  # (1) No constraints => identity
  fit_no_dual <- genplsc(X, Y, ncomp=3, dual=FALSE, verbose=FALSE)
  fit_yes_dual<- genplsc(X, Y, ncomp=3, dual=TRUE,  verbose=FALSE)
  
  # We'll compare some basic items:
  #  - singular values (these are d_full in each fit)
  d_nodual <- fit_no_dual$d_full
  d_dual   <- fit_yes_dual$d_full
  
  # Because there's an inherent "sign or orientation" ambiguity in SVD/PLS loadings,
  # we can't expect exact equality. We'll check that they're not drastically different.
  # We'll do a length check first:
  expect_true(abs(length(d_nodual) - length(d_dual)) <= 1)
  
  # Then compare first 2 or 3 singular values if they exist:
  K <- min(length(d_nodual), length(d_dual))
  expect_true(K >= 2)
  
  ratio_sv <- d_dual[1:K] / d_nodual[1:K]
  # ratio should be close to 1 if they're truly aligned
  # We'll allow a tolerance:
  expect_true(all(ratio_sv > 0.7 & ratio_sv < 1.4)) 
  
  # Compare factor scores in embedded domain: Tx, Ty
  # We'll do correlation check between the columns from each approach.
  # Because sign flips can happen, we do an absolute correlation check.
  # We'll define a function to compare columns via max corr:
  col_cors <- function(mat1, mat2, tol=0.7) {
    # For each col in mat1, pick best correlation in mat2
    # Return min of best cor across columns
    stopifnot(ncol(mat1) == ncol(mat2))
    out <- numeric(ncol(mat1))
    for(j in seq_len(ncol(mat1))) {
      corr_vals <- sapply(seq_len(ncol(mat2)), function(k) {
        cor(mat1[, j], mat2[, k])
      })
      out[j] <- max(abs(corr_vals), na.rm=TRUE)
    }
    min(out)
  }
  
  score_corX <- col_cors(fit_no_dual$Tx, fit_yes_dual$Tx)
  score_corY <- col_cors(fit_no_dual$Ty, fit_yes_dual$Ty)
  expect_true(score_corX > 0.7)  # somewhat aligned
  expect_true(score_corY > 0.7)
  
  # Compare loadings in original space
  loading_corX <- col_cors(fit_no_dual$vx, fit_yes_dual$vx)
  loading_corY <- col_cors(fit_no_dual$vy, fit_yes_dual$vy)
  expect_true(loading_corX > 0.7)
  expect_true(loading_corY > 0.7)
})

test_that("genplscorr dual=FALSE vs dual=TRUE with constraints yield similar results", {
  set.seed(456)
  n <- 10
  p <- 6
  q <- 5
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  
  # Let's define some basic PSD constraints, e.g. random cov for columns:
  # We'll do it for Ax only, for demonstration
  Ax_0 <- crossprod(matrix(rnorm(p * p, sd=0.2), p, p))
  # reduce rank => set small values below threshold
  Ax_0[Ax_0 < 0.15] <- 0
  # ensure symmetry
  Ax_0 <- 0.5 * (Ax_0 + t(Ax_0))
  
  # For Mx, My => identity, Ay => identity
  # to keep test simpler:
  Mx_0 <- diag(n)
  My_0 <- diag(n)
  Ay_0 <- diag(q)
  
  fit_no_dual <- genplsc(X, Y,
                            Mx=Mx_0, Ax=Ax_0,
                            My=My_0, Ay=Ay_0,
                            ncomp=3, dual=FALSE, verbose=TRUE)
  fit_yes_dual<- genplsc(X, Y,
                            Mx=Mx_0, Ax=Ax_0,
                            My=My_0, Ay=Ay_0,
                            ncomp=3, dual=TRUE,  verbose=TRUE)
  
  # do the same checks: singular values, factor scores, loadings
  d1 <- fit_no_dual$d_full
  d2 <- fit_yes_dual$d_full
  K  <- min(length(d1), length(d2))
  ratio_sv <- d2[1:K] / pmax(d1[1:K], 1e-9)
  # loosen tolerance as constraints can shift results more
  expect_true(all(ratio_sv > 0.5 & ratio_sv < 2.0))
  
  # factor scores correlation
  col_cors <- function(mat1, mat2, tol=0.5) {
    stopifnot(ncol(mat1) == ncol(mat2))
    out <- numeric(ncol(mat1))
    for(j in seq_len(ncol(mat1))) {
      corr_vals <- sapply(seq_len(ncol(mat2)), function(k) {
        cor(mat1[, j], mat2[, k])
      })
      out[j] <- max(abs(corr_vals), na.rm=TRUE)
    }
    min(out)
  }
  
  scX_cor <- col_cors(fit_no_dual$Tx, fit_yes_dual$Tx)
  scY_cor <- col_cors(fit_no_dual$Ty, fit_yes_dual$Ty)
  expect_true(scX_cor > 0.5)
  expect_true(scY_cor > 0.5)
  
  ldX_cor <- col_cors(fit_no_dual$vx, fit_yes_dual$vx)
  ldY_cor <- col_cors(fit_no_dual$vy, fit_yes_dual$vy)
  expect_true(ldX_cor > 0.5)
  expect_true(ldY_cor > 0.5)
})

test_that("genplscorr handles ncomp > rank gracefully", {
  set.seed(999)
  n <- 8
  p <- 5
  q <- 4
  X <- matrix(rnorm(n * p), n, p)
  Y <- matrix(rnorm(n * q), n, q)
  
  # Basic constraints (identity)
  # We'll choose ncomp= larger than min(n,p,q)
  # min(n,p,q) = 4, so let's pick ncomp=6
  ncomp_too_big <- 6
  
  fit_nodual <- genplsc(
    X, Y,
    ncomp = ncomp_too_big,
    dual = FALSE,
    verbose = FALSE
  )
  fit_dual <- genplsc(
    X, Y,
    ncomp = ncomp_too_big,
    dual = TRUE,
    verbose = FALSE
  )
  
  # Because ncomp > rank, each method might produce fewer actual components 
  # or produce degenerate columns in factor scores.
  # We'll check that the code didn't crash and that the final number of 
  # components is <= the requested ncomp.
  expect_true(fit_nodual$ncomp <= ncomp_too_big)
  expect_true(fit_dual$ncomp <= ncomp_too_big)
  
  # Also check that the singular values are not huge or NA
  expect_true(all(!is.na(fit_nodual$d_full)))
  expect_true(all(!is.na(fit_dual$d_full)))
  
  # Quick check that at least something non-degenerate emerged
  # We'll ensure they produce at least 1 component if the data isn't rank-0
  expect_true(fit_nodual$ncomp >= 1)
  expect_true(fit_dual$ncomp >= 1)
  
})

test_that("genplscorr dual=TRUE works with random sparse PSD constraints for X and Y", {
  set.seed(1001)
  n <- 12
  p <- 10
  q <- 8
  
  X <- matrix(rnorm(n * p, sd=0.5), n, p)
  Y <- matrix(rnorm(n * q, sd=0.5), n, q)
  
  # Build random sparse PSD constraints for Ax, Ay, Mx, My.
  # We'll do something like:
  # Ax0: a random p x p, then threshold to sparse; symmetrize => PSD approx
  # Then refine by multiplying with t(...) if needed.
  # A real PSD approach might do crossprod( ) then zero out small entries.
  
  random_psd_sparse <- function(dimsize=10, threshold=0.3) {
    mat <- matrix(rnorm(dimsize^2, sd=0.2), dimsize, dimsize)
    mat[abs(mat) < threshold] <- 0
    # symmetrize
    mat <- 0.5 * (mat + t(mat))
    # crossprod to ensure PSD
    # for a small example, crossprod is enough
    psd <- crossprod(mat)
    # keep it sparse
    psd <- Matrix::Matrix(psd, sparse=TRUE)
    psd
  }
  
  Ax_sp <- random_psd_sparse(dimsize=p, threshold=0.25)
  Ay_sp <- random_psd_sparse(dimsize=q, threshold=0.2)
  Mx_sp <- random_psd_sparse(dimsize=n, threshold=0.1)
  My_sp <- random_psd_sparse(dimsize=n, threshold=0.1)
  
  # Just do a small number of components
  K <- 3
  
  # We'll try dual=TRUE 
  fit_sp_dual <- genplsc(
    X, Y,
    Mx=Mx_sp, Ax=Ax_sp,
    My=My_sp, Ay=Ay_sp,
    ncomp=K,
    dual=TRUE,
    rank_Mx=NULL, rank_Ax=NULL, rank_My=NULL, rank_Ay=NULL, # => adaptive expansions
    var_threshold=0.9,  # let's capture 90% for speed
    max_k=5,           # small partial expansions
    verbose=FALSE
  )
  
  # If all is well, we won't crash or produce obviously degenerate results.
  expect_true(fit_sp_dual$ncomp <= K)
  expect_true(all(!is.na(fit_sp_dual$d_full)))
  expect_true(length(fit_sp_dual$d_full) >= 1)
  
  # Check that the embedding factor scores, etc., are of correct dimension
  # Factor scores Tx => n x (some #), loadings vx => p x (some #)
  expect_equal(dim(fit_sp_dual$Tx), c(n, fit_sp_dual$ncomp))
  expect_equal(dim(fit_sp_dual$vx), c(p, fit_sp_dual$ncomp))
  
  # Not easily verifiable if the decomposition is "right," but 
  # at least we confirm the code handles random sparse PSD constraints and 
  # the partial expansions didn't blow up. 
  # In real usage, you'd do deeper validations or compare with a known solution on smaller data.
  
  # Additional sanity checks:
  # The d_full typically won't be NA or infinite
  expect_true(all(is.finite(fit_sp_dual$d_full)))
  # Some measure that the loadings aren't all zero
  loadnormsX <- apply(fit_sp_dual$vx, 2, function(cc) sqrt(sum(cc^2)))
  expect_true(all(loadnormsX > 1e-8))
})


test_that("genplscorr: dual=TRUE vs dual=FALSE produce closely equivalent results on moderate data", {
  set.seed(222)
  n <- 12
  p <- 7
  q <- 6
  
  # Generate moderately sized data
  X <- matrix(rnorm(n * p, mean=0, sd=0.6), n, p)
  Y <- matrix(rnorm(n * q, mean=1, sd=0.8), n, q)
  
  # Let's define modest PSD constraints for Ax, Ay, Mx, My:
  # We'll do a crossprod approach for each, plus a small threshold to ensure some sparsity.
  
  make_psd <- function(dimsize, threshold=0.1, sdev=0.2) {
    mat <- matrix(rnorm(dimsize^2, sd=sdev), dimsize, dimsize)
    # zero out small absolute entries
    mat[abs(mat) < threshold] <- 0
    # symmetrize
    mat <- 0.5*(mat + t(mat))
    # crossprod => ensure PSD
    out <- crossprod(mat)
    # small dimension => we can keep as is for demonstration
    # optionally convert to sparse
    # out <- Matrix::Matrix(out, sparse=TRUE)
    out
  }
  
  # Let n=12 => Mx, My are 12x12. p=7 => Ax is 7x7. q=6 => Ay is 6x6
  Mx_0 <- make_psd(n, threshold=0.08)
  My_0 <- make_psd(n, threshold=0.08)
  Ax_0 <- make_psd(p, threshold=0.07)
  Ay_0 <- make_psd(q, threshold=0.07)
  
  # We choose a small-ish rank expansions but let them be adaptive:
  # We'll do an example ncomp=4
  K <- 4
  
  # dual=FALSE
  fit_nodual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp = K,
    dual = FALSE,
    rank_Mx=NULL, rank_Ax=NULL, rank_My=NULL, rank_Ay=NULL,
    var_threshold=0.95,
    max_k=10,  # partial expansions up to 10 
    verbose=FALSE
  )
  
  # dual=TRUE
  fit_dual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp = K,
    dual = TRUE,
    rank_Mx=NULL, rank_Ax=NULL, rank_My=NULL, rank_Ay=NULL,
    var_threshold=0.95,
    max_k=10,
    verbose=FALSE
  )
  
  expect_true(fit_nodual$ncomp <= K)
  expect_true(fit_dual$ncomp <= K)
  
  # 1) Compare singular values
  d1 <- fit_nodual$d_full
  d2 <- fit_dual$d_full
  expect_true(length(d1) > 0 && length(d2) > 0)
  
  # We'll compare up to min of lengths
  L <- min(length(d1), length(d2))
  ratio_sv <- d2[1:L] / pmax(d1[1:L], 1e-12)
  # Because constraints can shift them somewhat, we just expect no wild mismatch:
  expect_true(all(ratio_sv > 0.5 & ratio_sv < 2.0))
  
  # 2) Compare factor scores or loadings
  # We'll define a helper again for correlation across columns:
  col_cors <- function(matA, matB) {
    stopifnot(ncol(matA) == ncol(matB))
    out <- numeric(ncol(matA))
    for (j in seq_len(ncol(matA))) {
      # best correlation with any col in matB
      cvals <- sapply(seq_len(ncol(matB)), function(k) {
        cor(matA[, j], matB[, k])
      })
      out[j] <- max(abs(cvals))
    }
    min(out)
  }
  
  # Factor scores Tx
  if (!is.null(fit_nodual$Tx) && !is.null(fit_dual$Tx)) {
    scX_cor <- col_cors(fit_nodual$Tx, fit_dual$Tx)
    expect_true(scX_cor > 0.5)
  }
  
  # Factor scores Ty
  if (!is.null(fit_nodual$Ty) && !is.null(fit_dual$Ty)) {
    scY_cor <- col_cors(fit_nodual$Ty, fit_dual$Ty)
    expect_true(scY_cor > 0.5)
  }
  
  # Loadings vx, vy in original domain
  ldX_cor <- col_cors(fit_nodual$vx, fit_dual$vx)
  ldY_cor <- col_cors(fit_nodual$vy, fit_dual$vy)
  # We'll want these to be decently correlated
  expect_true(ldX_cor > 0.5)
  expect_true(ldY_cor > 0.5)
})



test_that("genplscorr2 dual=TRUE vs dual=FALSE produce similar results for moderate data and PSD constraints", {
  set.seed(227)
  n <- 34
  p <- 16
  q <- 15
  
  # Generate moderately sized data
  X <- matrix(rnorm(n * p, mean=0, sd=0.6), n, p)
  Y <- matrix(rnorm(n * q, mean=1, sd=0.8), n, q)
  
  # Let's define modest PSD constraints for Ax, Ay, Mx, My:
  # We'll do a crossprod approach for each, plus a small threshold to ensure some sparsity.
  make_psd <- function(dimsize, threshold=0.1, sdev=0.2) {
    mat <- matrix(rnorm(dimsize^2, sd=sdev), dimsize, dimsize)
    # zero out small entries
    mat[abs(mat) < threshold] <- 0
    # symmetrize
    mat <- 0.5 * (mat + t(mat))
    # ensure PSD by crossprod
    out <- crossprod(mat)
    out
  }
  
  # Let n=12 => Mx, My are 12x12. p=7 => Ax=7x7, q=6 => Ay=6x6
  Mx_0 <- make_psd(n, threshold=0.08)
  My_0 <- make_psd(n, threshold=0.08)
  Ax_0 <- make_psd(p, threshold=0.07)
  Ay_0 <- make_psd(q, threshold=0.07)
  
  # We'll do an example with ncomp=4
  K <- 4
  
  # 1) dual=FALSE
  fit_nodual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp = K,
    dual = FALSE,
    rank_Mx=3L, rank_Ax=3, rank_My=3, rank_Ay=3,
    var_threshold=1,
    max_k=100,  # partial expansions up to 10 
    verbose=FALSE
  )
  
  # 2) dual=TRUE
  fit_dual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp = K,
    dual = TRUE,
    rank_Mx=3, rank_Ax=3, rank_My=3, rank_Ay=3,
    var_threshold=1,
    max_k=100,
    verbose=FALSE
  )
  
  # Basic checks
  expect_true(fit_nodual$ncomp <= K)
  expect_true(fit_dual$ncomp <= K)
  expect_s3_class(fit_nodual, "cross_projector")
  expect_s3_class(fit_dual, "cross_projector")
  
  # Compare singular values (d_full)
  d1 <- fit_nodual$d_full
  d2 <- fit_dual$d_full
  expect_true(length(d1) > 0 && length(d2) > 0)
  L <- min(length(d1), length(d2))
  ratio_sv <- d2[1:L] / pmax(1e-12, d1[1:L])
  # We just want them within a modest ratio
  expect_true(all(ratio_sv > 0.3 & ratio_sv < 3.0))
  
  # A helper to measure cross-column correlation
  # allowing sign flips or different order
  col_cors <- function(matA, matB) {
    # returns min of max correlation per column in matA
    # ignoring sign or col permutations
    stopifnot(ncol(matA) > 0, ncol(matB) > 0)
    out <- numeric(ncol(matA))
    for (j in seq_len(ncol(matA))) {
      # correlate matA col j with *all* columns in matB, take max of absolute
      cors <- sapply(seq_len(ncol(matB)), function(k) {
        cc <- cor(matA[, j], matB[, k])
        if (is.na(cc)) 0 else abs(cc)
      })
      out[j] <- max(cors)
    }
    mean(out)
  }
  
  # Compare factor scores if available
  if (!is.null(fit_nodual$Tx) && !is.null(fit_dual$Tx)) {
    score_cor_x <- col_cors(fit_nodual$Tx, fit_dual$Tx)
    expect_true(score_cor_x > 0.5) 
  }
  if (!is.null(fit_nodual$Ty) && !is.null(fit_dual$Ty)) {
    score_cor_y <- col_cors(fit_nodual$Ty, fit_dual$Ty)
    expect_true(score_cor_y > 0.5) 
  }
  
  # Compare loadings in original domain
  expect_equal(nrow(fit_nodual$vx), p)
  expect_equal(nrow(fit_dual$vx), p)
  if (ncol(fit_nodual$vx) > 0 && ncol(fit_dual$vx) > 0) {
    cor_vx <- col_cors(fit_nodual$vx, fit_dual$vx)
    expect_true(cor_vx > 0.5)
  }
  
  expect_equal(nrow(fit_nodual$vy), q)
  expect_equal(nrow(fit_dual$vy), q)
  if (ncol(fit_nodual$vy) > 0 && ncol(fit_dual$vy) > 0) {
    cor_vy <- col_cors(fit_nodual$vy, fit_dual$vy)
    expect_true(cor_vy > 0.5)
  }
})


test_that("genplscorr2 handles low-rank data + constraints consistently", {
  set.seed(999)  # for reproducibility
  
  #### 1) Generate low-rank data (X, Y) + noise ####
  n <- 30
  p <- 20
  q <- 15
  
  # We assume rank ~ 4 for the "true" latent structure.
  rank_true <- 4
  
  # (a) Build a latent factor matrix for X: shape n x rank_true
  #     plus a separate factor for Y.
  #     Then create X, Y by (FactorX * LoadingsX) + noise, etc.
  
  # "true" latent factor scores for X, dimension n x rank_true
  Fx <- matrix(rnorm(n * rank_true, sd=1.0), n, rank_true)
  # loadings for X in p-d space
  Lx <- matrix(rnorm(p * rank_true, sd=1.0), p, rank_true)
  
  # "true" latent factor for Y, dimension n x rank_true
  # optionally correlated with X's factor: e.g. let's take Fx plus small random
  # or you can do a separate factor. We'll do correlated for demonstration:
  Fy <- Fx + matrix(rnorm(n * rank_true, sd=0.2), n, rank_true)
  Ly <- matrix(rnorm(q * rank_true, sd=1.0), q, rank_true)
  
  # Construct X, Y (rank ~ rank_true) + random noise
  Xbase <- Fx %*% t(Lx)   # shape n x p
  Ybase <- Fy %*% t(Ly)   # shape n x q
  
  noiseX <- matrix(rnorm(n * p, sd=0.5), n, p)
  noiseY <- matrix(rnorm(n * q, sd=0.5), n, q)
  
  X <- Xbase + noiseX
  Y <- Ybase + noiseY
  
  # Optionally center columns (typical for PLS):
  # We'll rely on genplscorr2's own 'preproc_x=pass()' or we can do:
  # X <- scale(X, center=TRUE, scale=FALSE)
  # Y <- scale(Y, center=TRUE, scale=FALSE)
  
  #### 2) Build low-rank constraints for Mx, Ax, My, Ay ####
  # Let's produce, e.g., rank 5 expansions for each, plus small noise => shape n x n or p x p, etc.
  # We'll do something like U diag(...) V^T, with rank=5. Then we can add a dash of noise.
  
  make_lowrank_psd <- function(dimsize, r=5, noise_sd=0.01) {
    # generate a random U, V => shape dimsize x r
    U <- matrix(rnorm(dimsize*r), dimsize, r)
    V <- matrix(rnorm(dimsize*r), dimsize, r)
    # build rank-r mat => U diag(...) V^T
    # simplest: let diag(...) be positives (random?)
    diagvals <- runif(r, min=0.5, max=2.0)
    lowrank <- U %*% (diag(diagvals)) %*% t(V)  # shape dimsize x dimsize
    # symmetrize
    lowrank <- 0.5 * (lowrank + t(lowrank))
    # add small noise
    noise <- matrix(rnorm(dimsize*dimsize, sd=noise_sd), dimsize, dimsize)
    out <- lowrank + noise
    # make sure PSD => do a quick shift or crossprod approach
    # to ensure non-negative eigenvalues
    # simplest: crossprod => or we can do an eigen-cleaning
    out <- 0.5*(out + t(out))  # sym
    # We'll do a small "eigen-truncation" approach:
    eout <- eigen(out, symmetric=TRUE)
    lampos <- pmax(eout$values, 1e-5)
    out_psd <- eout$vectors %*% diag(lampos) %*% t(eout$vectors)
    out_psd
  }
  
  # Mx, My => n x n
  # Ax => p x p,   Ay => q x q
  Mx_0 <- make_lowrank_psd(n, r=5)
  My_0 <- make_lowrank_psd(n, r=5)
  Ax_0 <- make_lowrank_psd(p, r=5)
  Ay_0 <- make_lowrank_psd(q, r=5)
  
  #### 3) Run genplscorr2 => dual=FALSE vs. dual=TRUE ####
  # We'll request ncomp=6, see how many we actually get
  K <- 6
  
  fit_nodual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp=K,
    dual=FALSE,
    rank_Mx=NULL, rank_Ax=NULL, rank_My=NULL, rank_Ay=NULL,
    var_threshold=.9,
    max_k=100,
    verbose=FALSE
  )
  
  fit_dual <- genplsc(
    X, Y,
    Mx=Mx_0, Ax=Ax_0,
    My=My_0, Ay=Ay_0,
    ncomp=K,
    dual=TRUE,
    rank_Mx=NULL, rank_Ax=NULL, rank_My=NULL, rank_Ay=NULL,
    var_threshold=.9,
    max_k=100,
    verbose=FALSE
  )
  
  #### 4) Basic checks ####
  # Check that each fit has # of components <= K
  expect_true(fit_nodual$ncomp <= K)
  expect_true(fit_dual$ncomp   <= K)
  
  # Confirm we have loadings of appropriate dimension
  expect_equal(dim(fit_nodual$vx)[1], p)  # p x ncomp
  expect_equal(dim(fit_dual$vx)[1],   p)
  
  # 5) Compare first few directions 
  # They may differ due to partial expansions, sign flips, etc.
  # We can do a rough check on correlation among the first 2-3 columns of vx or vy
  # e.g., the absolute correlation => should be fairly large if they're capturing similar directions.
  # but won't be identical => we just check if there's some overlap. 
  library(stats)
  ncomp_compare <- min(fit_nodual$ncomp, fit_dual$ncomp, 3)
  cor_thresh <- 0.5  # somewhat arbitrary
  for(j in seq_len(ncomp_compare)) {
    # correlation of loadings on X
    cval_x <- cor(fit_nodual$vx[,j], fit_dual$vx[,j])
    # correlation of loadings on Y
    cval_y <- cor(fit_nodual$vy[,j], fit_dual$vy[,j])
    # Just expect they're not trivial => we won't do expect_* but we can check 
    # they aren't negligible
    testthat::expect_gt(abs(cval_x), 0.2, 
                        info=paste("component", j, "X-loadings correlation is non-trivial"))
    testthat::expect_gt(abs(cval_y), 0.2, 
                        info=paste("component", j, "Y-loadings correlation is non-trivial"))
  }
  
  # You might also check factor scores or singular values, etc.
  # The point is: in a low-rank scenario, we'd see that both approaches 
  # pick up somewhat similar directions, though not identical.
})