lib <- Sys.getenv("EC_LIB", "devlib"); .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(eigencore)})
cat("### eigencore", lib, as.character(packageVersion("eigencore")), "\n")
set.seed(3)
# certificate notes on the native block matrix-free hermitian path
S <- crossprod(matrix(rnorm(400 * 300), 400)) / 400; refS <- eigen(S, symmetric = TRUE)$values
ops <- linear_operator(c(300, 300), function(X_, alpha = 1, beta = 0, Y = NULL, ...) S %*% X_, function(X_, alpha = 1, beta = 0, Y = NULL, ...) S %*% X_, structure = hermitian())
r <- eigencore::eigs_sym(ops, 3, "LA", method = lanczos(block = 4))
cat("-- mf block4: passed", isTRUE(r$certificate$passed), "| norm_bound", r$certificate$norm_bound_type, "| scale_is_estimate", r$certificate$scale_is_estimate, "| notes:", paste(r$certificate$notes, collapse = "; "), "\n")
# SVD-formulated GMD via eigencore vs genpca C++ spectra kernel
suppressMessages(pkgload::load_all("~/code/genpca", quiet = TRUE))
n <- 4000; p <- 600; k <- 10
Xg <- matrix(rnorm(n * p), n); Rg <- crossprod(matrix(rnorm(p * p), p)) / p + diag(p) * 0.1; Qg <- runif(n, 0.5, 2)
tm <- function(expr) { t0 <- proc.time()[[3]]; force(expr); proc.time()[[3]] - t0 }
ts <- tm(fs <- genpca(Xg, A = Rg, M = Qg, ncomp = k, method = "spectra", preproc = multivarious::pass()))
LR <- t(chol(Rg)); sq <- sqrt(Qg)
opg <- eigencore::linear_operator(c(n, p), apply = function(V, alpha = 1, beta = 0, Y = NULL, ...) sq * (Xg %*% (LR %*% V)),
                                  apply_adjoint = function(U, alpha = 1, beta = 0, Y = NULL, ...) crossprod(LR, crossprod(Xg, sq * U)))
te <- tm(fe <- eigencore::svds(opg, k)); te0 <- tm(fe0 <- eigencore::svds(opg, k, certify = FALSE))
ove <- solve(t(LR), fe$v)
sinang <- function(U, V) { s <- svd(crossprod(qr.Q(qr(as.matrix(U))), qr.Q(qr(as.matrix(V)))))$d; sqrt(max(0, 1 - min(s)^2)) }
cat(sprintf("-- GMD 4000x600 k=10: C++ Spectra %.2fs | eigencore svd-operator %.2fs (certify=FALSE %.2fs; passed %s, %s) | sdev max rel diff %.1e | V subspace sin %.1e\n",
            ts, te, te0, isTRUE(fe$certificate$passed), fe$certificate$norm_bound_type, max(abs(fs$sdev - fe$d) / fe$d), sinang(fs$ov, ove)))
# dense-metric diagonal case: genpca spectra (RSpectra::svds on Xw) vs eigencore svds on Xw
Xw <- sq * Xg
cat(sprintf("-- diag-metric svds 4000x600 k=10: RSpectra %.2fs | eigencore %.2fs\n", tm(RSpectra::svds(Xw, k)), tm(eigencore::svds(Xw, k))))
# fair timings (installed lib)
Sd <- crossprod(matrix(rnorm(1600 * 1500), 1600)) / 1600
cat(sprintf("-- dense sym 1500 k=10: RSpectra %.2fs | eigencore %.2fs | certify=FALSE %.2fs\n", tm(RSpectra::eigs_sym(Sd, 10, "LA")), tm(eigencore::eigs_sym(Sd, 10, "LA")), tm(eigencore::eigs_sym(Sd, 10, "LA", certify = FALSE))))
nn <- 20000; Sp <- rsparsematrix(nn, nn, density = 5e-4); Sp <- as(as(forceSymmetric(Sp + Diagonal(nn)), "generalMatrix"), "CsparseMatrix")
cat(sprintf("-- sparse 20000 k=10: RSpectra %.2fs | eigencore %.2fs\n", tm(RSpectra::eigs_sym(Sp, 10, "LM")), tm(eigencore::eigs_sym(Sp, 10, "LM"))))
Xb <- matrix(rnorm(5000 * 300), 5000)
cat(sprintf("-- svds 5000x300 k=10: RSpectra %.2fs | eigencore %.2fs\n", tm(RSpectra::svds(Xb, 10)), tm(eigencore::svds(Xb, 10))))
Xe <- matrix(rnorm(3000 * 400), 3000); Ye <- matrix(rnorm(3000 * 300), 3000)
cat(sprintf("-- op svds 400x300 k=5: RSpectra %.2fs | eigencore %.2fs\n",
  tm(RSpectra::svds(function(x, args) crossprod(Xe, Ye %*% x), 5, Atrans = function(x, args) crossprod(Ye, Xe %*% x), dim = c(400, 300))),
  tm(eigencore::svds(eigencore::linear_operator(c(400, 300), function(V, alpha = 1, beta = 0, Y = NULL, ...) crossprod(Xe, Ye %*% V), function(U, alpha = 1, beta = 0, Y = NULL, ...) crossprod(Ye, Xe %*% U)), 5))))
# small-problem overhead (typical test sizes)
Ss <- crossprod(matrix(rnorm(200 * 120), 200)); Xs <- matrix(rnorm(150 * 40), 150)
cat(sprintf("-- small: eigs_sym 120 k=3 RSpectra %.3fs | eigencore %.3fs ; svds 150x40 k=3 RSpectra %.3fs | eigencore %.3fs\n",
  tm(for (i in 1:20) RSpectra::eigs_sym(Ss, 3, "LA")) / 20, tm(for (i in 1:20) eigencore::eigs_sym(Ss, 3, "LA")) / 20,
  tm(for (i in 1:20) RSpectra::svds(Xs, 3)) / 20, tm(for (i in 1:20) eigencore::svds(Xs, 3)) / 20))
