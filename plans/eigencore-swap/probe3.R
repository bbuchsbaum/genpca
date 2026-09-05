lib <- Sys.getenv("EC_LIB", "devlib"); .libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({library(Matrix); library(eigencore)})
cat("### eigencore", lib, as.character(packageVersion("eigencore")), "\n")
set.seed(3); X <- matrix(rnorm(400 * 60), 400); refX <- svd(X)
# (i) alpha/beta/Y values actually passed
ab <- list()
op <- linear_operator(dim(X), apply = function(X_, alpha = 1, beta = 0, Y = NULL, ...) { ab[[length(ab) + 1]] <<- c(alpha, beta, is.null(Y)); X %*% X_ },
                      apply_adjoint = function(X_, alpha = 1, beta = 0, Y = NULL, ...) crossprod(X, X_))
r <- svds(op, 3); cat("-- op svds d err:", max(abs(r$d - refX$d[1:3])), " unique (alpha,beta,Ynull):", paste(unique(sapply(ab, paste, collapse = "/")), collapse = "  "), "\n")
# (ii) matrix-free hermitian with block lanczos
S <- crossprod(matrix(rnorm(400 * 300), 400)) / 400; refS <- eigen(S, symmetric = TRUE)$values
ops <- linear_operator(c(300, 300), function(X_, alpha = 1, beta = 0, Y = NULL, ...) S %*% X_, function(X_, alpha = 1, beta = 0, Y = NULL, ...) S %*% X_, structure = hermitian())
for (b in c(1, 2, 4)) { r <- tryCatch(eigs_sym(ops, 3, "LA", method = lanczos(block = b)), error = function(e) e)
  if (inherits(r, "error")) cat(sprintf("-- mf hermitian lanczos(block=%d): ERROR %s\n", b, conditionMessage(r))) else
  cat(sprintf("-- mf hermitian lanczos(block=%d): err %.1e passed %s nconv %d method '%s'\n", b, max(abs(r$values - refS[1:3])), isTRUE(r$certificate$passed), r$nconv, r$diagnostics$method)) }
# (iv) seed determinism magnitude
v1 <- eigs_sym(S, 3, "LA", seed = 1L)$vectors; v2 <- eigs_sym(S, 3, "LA", seed = 1L)$vectors; v3 <- eigs_sym(S, 3, "LA", seed = 2L)$vectors
cat(sprintf("-- same seed vec maxdiff %.1e ; different seed %.1e\n", max(abs(abs(v1) - abs(v2))), max(abs(abs(v1) - abs(v3)))))
# (iii) SVD-formulated GMD operator vs genpca C++ spectra kernel
suppressMessages(pkgload::load_all("~/code/genpca", quiet = TRUE))
n <- 4000; p <- 600; k <- 10
Xg <- matrix(rnorm(n * p), n); Rg <- crossprod(matrix(rnorm(p * p), p)) / p + diag(p) * 0.1; Qg <- runif(n, 0.5, 2)
t0 <- proc.time()[[3]]; fs <- genpca(Xg, A = Rg, M = Qg, ncomp = k, method = "spectra", preproc = multivarious::pass()); ts <- proc.time()[[3]] - t0
LR <- t(chol(Rg)); sq <- sqrt(Qg)
opg <- linear_operator(c(n, p), apply = function(V, alpha = 1, beta = 0, Y = NULL, ...) sq * (Xg %*% (LR %*% V)),
                       apply_adjoint = function(U, alpha = 1, beta = 0, Y = NULL, ...) crossprod(LR, crossprod(Xg, sq * U)))
t0 <- proc.time()[[3]]; fe <- svds(opg, k); te <- proc.time()[[3]] - t0
ovs <- fs$ov; ove <- solve(t(LR), fe$v)   # V = L^{-T} Z
cat(sprintf("-- GMD 4000x600 k=10: spectra C++ %.2fs, eigencore svd-operator %.2fs (passed %s); sdev maxreldiff %.1e; V subspace sin %.1e\n", ts, te, isTRUE(fe$certificate$passed), max(abs(fs$sdev - fe$d) / fe$d),
            { s <- svd(crossprod(qr.Q(qr(as.matrix(ovs))), qr.Q(qr(ove))))$d; sqrt(max(0, 1 - min(s)^2)) }))
t0 <- proc.time()[[3]]; fe2 <- svds(opg, k, certify = FALSE); cat(sprintf("-- same, certify=FALSE: %.2fs\n", proc.time()[[3]] - t0))
# (v) fair timings with installed lib
Sd <- crossprod(matrix(rnorm(1600 * 1500), 1600)) / 1600
tm <- function(expr) { t0 <- proc.time()[[3]]; force(expr); proc.time()[[3]] - t0 }
cat(sprintf("-- dense sym 1500 k=10: RSpectra %.2fs | eigencore %.2fs | certify=FALSE %.2fs\n", tm(RSpectra::eigs_sym(Sd, 10, "LA")), tm(eigs_sym(Sd, 10, "LA")), tm(eigs_sym(Sd, 10, "LA", certify = FALSE))))
nn <- 20000; Sp <- rsparsematrix(nn, nn, density = 5e-4); Sp <- as(as(forceSymmetric(Sp + Diagonal(nn)), "generalMatrix"), "CsparseMatrix")
cat(sprintf("-- sparse 20000 k=10: RSpectra %.2fs | eigencore %.2fs\n", tm(RSpectra::eigs_sym(Sp, 10, "LM")), tm(eigs_sym(Sp, 10, "LM"))))
Xb <- matrix(rnorm(5000 * 300), 5000)
cat(sprintf("-- svds 5000x300 k=10: RSpectra %.2fs | eigencore %.2fs\n", tm(RSpectra::svds(Xb, 10)), tm(svds(Xb, 10))))
Xe <- matrix(rnorm(3000 * 400), 3000); Ye <- matrix(rnorm(3000 * 300), 3000)
cat(sprintf("-- op svds 400x300 k=5: RSpectra %.2fs | eigencore %.2fs\n",
  tm(RSpectra::svds(function(x, args) crossprod(Xe, Ye %*% x), 5, Atrans = function(x, args) crossprod(Ye, Xe %*% x), dim = c(400, 300))),
  tm(svds(linear_operator(c(400, 300), function(V, alpha = 1, beta = 0, Y = NULL, ...) crossprod(Xe, Ye %*% V), function(U, alpha = 1, beta = 0, Y = NULL, ...) crossprod(Ye, Xe %*% U)), 5))))
