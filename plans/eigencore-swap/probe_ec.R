mode <- Sys.getenv("EC_MODE", "dev")
suppressPackageStartupMessages(library(Matrix))
if (mode == "dev") {
  suppressMessages(pkgload::load_all("~/code/eigencore", quiet = TRUE))
  ver <- read.dcf("~/code/eigencore/DESCRIPTION")[, "Version"]
} else {
  suppressPackageStartupMessages(library(eigencore, lib.loc = "cranlib"))
  ver <- as.character(packageVersion("eigencore", lib.loc = "cranlib"))
}
cat("### eigencore", mode, ver, "  RSpectra", as.character(packageVersion("RSpectra")), "\n")
P <- function(label, expr) {
  nw <- 0L; t0 <- proc.time()[["elapsed"]]
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "perr")),
    warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  dt <- proc.time()[["elapsed"]] - t0
  if (inherits(r, "perr")) cat(sprintf("[%-30s] ERROR: %s\n", label, substr(gsub("\n", " ", r$msg), 1, 150)))
  else cat(sprintf("[%-30s] %s  (%.2fs%s)\n", label, paste(format(r, digits = 4), collapse = " "), dt, if (nw) paste0(", ", nw, " warn") else ""))
  invisible(r)
}
sinang <- function(U, V) { s <- svd(crossprod(qr.Q(qr(as.matrix(U))), qr.Q(qr(as.matrix(V)))))$d; sqrt(max(0, 1 - min(s)^2)) }
set.seed(42)
# --- 1. matrix-class acceptance
S0 <- crossprod(matrix(rnorm(300 * 60), 300)) / 300; S0 <- S0 - 0.2 * diag(60)   # indefinite-ish, 60x60
ref <- eigen(S0, symmetric = TRUE)
Ssp <- as(Matrix(S0 * (abs(S0) > 0.05)), "CsparseMatrix"); Ssp <- forceSymmetric(Ssp)
classes <- list(base = S0, dge = Matrix(S0, sparse = FALSE), dsy = forceSymmetric(Matrix(S0, sparse = FALSE)),
                dgC = as(as(Ssp, "generalMatrix"), "CsparseMatrix"), dsC = Ssp)
for (nm in names(classes)) P(paste0("class ", nm, " (", class(classes[[nm]])[1], ")"), {
  r <- eigs_sym(classes[[nm]], k = 3, which = "LM"); refc <- eigen(as.matrix(classes[[nm]]), symmetric = TRUE)
  c(maxabs_err = max(abs(sort(abs(r$values), TRUE) - sort(abs(refc$values), TRUE)[1:3]))) })
# --- 2. accuracy vs eigen and vs RSpectra, dense 300x300
S <- crossprod(matrix(rnorm(400 * 300), 400)) / 400; S <- S - 0.5 * diag(300); refS <- eigen(S, symmetric = TRUE)
P("LM k=5 values err", { r <- eigs_sym(S, 5, "LM"); c(ec = max(abs(r$values - refS$values[1:5])), rs = max(abs(RSpectra::eigs_sym(S, 5, "LM")$values - refS$values[1:5]))) })
P("LM k=5 subspace sin", { r <- eigs_sym(S, 5, "LM"); c(ec = sinang(r$vectors, refS$vectors[, 1:5]), rs = sinang(RSpectra::eigs_sym(S, 5, "LM")$vectors, refS$vectors[, 1:5])) })
P("LA k=5 values err", { r <- eigs_sym(S, 5, "LA"); c(ec = max(abs(r$values - refS$values[1:5])), rs = max(abs(RSpectra::eigs_sym(S, 5, "LA")$values - refS$values[1:5]))) })
P("SA k=1 min eigenvalue", { r <- eigs_sym(S, 1, "SA"); c(ec = r$values, rs = RSpectra::eigs_sym(S, 1, "SA")$values, true = min(refS$values)) })
P("order of values LM", { r <- eigs_sym(S, 5, "LM"); c(decreasing = !is.unsorted(rev(r$values)), first = r$values[1]) })
# --- 3. svds dense
X <- matrix(rnorm(400 * 60), 400); refX <- svd(X)
P("svds dense k=5 d err", { r <- svds(X, k = 5); c(ec = max(abs(r$d - refX$d[1:5])), rs = max(abs(RSpectra::svds(X, 5)$d - refX$d[1:5]))) })
P("svds dense u/v sin", { r <- svds(X, k = 5); c(u = sinang(r$u, refX$u[, 1:5]), v = sinang(r$v, refX$v[, 1:5]), dims = paste(dim(r$u), collapse = "x")) })
P("svds nu=0 nv=5", { r <- svds(X, k = 5, nu = 0, nv = 5); c(u_null = is.null(r$u), v_cols = ncol(r$v)) })
P("svds dgeMatrix input", { r <- svds(Matrix(X, sparse = FALSE), k = 3); max(abs(r$d - refX$d[1:3])) })
# --- 4. passing tol/maxit via ... and opts
P("eigs_sym ... tol=1e-12", { r <- eigs_sym(S, 3, "LA", tol = 1e-12); max(abs(r$values - refS$values[1:3])) })
P("eigs_sym ... maxit=200", { r <- eigs_sym(S, 3, "LA", maxit = 200); r$niter })
P("eigs_sym opts=list(tol,maxitr)", { r <- eigs_sym(S, 3, "LA", opts = list(tol = 1e-6, maxitr = 100)); "accepted" })
P("svds ... tol=1e-6", { r <- svds(X, 3, tol = 1e-6); "accepted" })
P("svds Atrans/dim args (RSpectra style)", { r <- svds(function(x, args) X %*% x, k = 2, Atrans = function(x, args) crossprod(X, x), dim = dim(X)); "accepted" })
# --- 5. operator interface: block semantics
seen <- list()
op <- linear_operator(dim = dim(X),
  apply = function(x, ...) { seen[[length(seen) + 1]] <<- c(class(x)[1], paste(dim(x), collapse = "x")); X %*% x },
  apply_adjoint = function(x, ...) crossprod(X, x))
P("operator svds k=2 d err", { r <- svds(op, k = 2); max(abs(r$d - refX$d[1:2])) })
P("operator apply() input shapes", unique(vapply(seen, function(z) paste(z, collapse = " "), "")))
opv <- linear_operator(dim = dim(X), apply = function(x, ...) as.numeric(X %*% x), apply_adjoint = function(x, ...) as.numeric(crossprod(X, x)))
P("operator returning numeric vec", { r <- svds(opv, k = 1); max(abs(r$d - refX$d[1])) })
ops <- linear_operator(dim = c(300, 300), apply = function(x, ...) S %*% x, structure = hermitian())
P("sym operator eigs_sym k=3", { r <- eigs_sym(ops, 3, "LA"); max(abs(r$values - refS$values[1:3])) })
P("symmetric_operator() k=3", { r <- eigs_sym(symmetric_operator(linear_operator(c(300,300), function(x, ...) S %*% x, function(x, ...) S %*% x)), 3, "LA"); max(abs(r$values - refS$values[1:3])) })
# --- 6. RNG and determinism
P("consumes .Random.seed?", { set.seed(7); a <- runif(1); set.seed(7); invisible(eigs_sym(S, 3, "LA")); invisible(svds(X, 3)); b <- runif(1); c(rng_consumed = !identical(a, b)) })
P("deterministic repeat (values)", { a <- eigs_sym(S, 3, "LA"); b <- eigs_sym(S, 3, "LA"); c(identical_values = identical(a$values, b$values), vec_maxdiff = max(abs(abs(a$vectors) - abs(b$vectors)))) })
P("RSpectra deterministic repeat", { a <- RSpectra::eigs_sym(S, 3, "LA"); b <- RSpectra::eigs_sym(S, 3, "LA"); c(identical_values = identical(a$values, b$values)) })
# --- 7. k edge cases
S50 <- S[1:50, 1:50]; ref50 <- eigen(S50, symmetric = TRUE)
P("eigs_sym k = n-1", { r <- eigs_sym(S50, 49, "LA"); c(n = length(r$values), err = max(abs(r$values - ref50$values[1:49]))) })
P("eigs_sym k = n", { r <- eigs_sym(S50, 50, "LA"); c(n = length(r$values), err = max(abs(r$values - ref50$values))) })
P("RSpectra eigs_sym k = n", { r <- RSpectra::eigs_sym(S50, 50, "LA"); length(r$values) })
X30 <- X[1:30, 1:10]
P("svds k = min(n,p)", { r <- svds(X30, 10); c(n = length(r$d), err = max(abs(r$d - svd(X30)$d))) })
P("svds k = 1 on 1-col", { r <- svds(X[, 1, drop = FALSE], 1); c(d = r$d, true = sqrt(sum(X[, 1]^2))) })
# --- 8. certificate behaviour on a hard case: repeated eigenvalue split by k
Sc <- diag(c(rep(2, 5), 1.999, 1, 0.5, rep(0.1, 42))); Sc <- Sc + matrix(rnorm(2500, sd = 1e-9), 50); Sc <- (Sc + t(Sc)) / 2
P("clustered k=3 (cluster of 5)", { r <- eigs_sym(Sc, 3, "LA"); c(vals = r$values, passed = isTRUE(r$certificate$passed)) })
P("clustered k=6 (splits 5+1)", { r <- eigs_sym(Sc, 6, "LA"); c(passed = isTRUE(r$certificate$passed), err = max(abs(r$values - sort(diag(Sc), TRUE)[1:6]))) })
P("what does a failed cert do", { r <- eigs_sym(Sc, 6, "LA", tol = 1e-15); c(passed = isTRUE(r$certificate$passed), nconv = r$nconv, class = class(r)[1]) })
# --- 9. timing
Sd <- crossprod(matrix(rnorm(1600 * 1500), 1600)) / 1600
P("time dense sym 1500 k=10 RSpectra", { t0 <- proc.time()[[3]]; r <- RSpectra::eigs_sym(Sd, 10, "LA"); c(sec = proc.time()[[3]] - t0, v1 = r$values[1]) })
P("time dense sym 1500 k=10 eigencore", { t0 <- proc.time()[[3]]; r <- eigs_sym(Sd, 10, "LA"); c(sec = proc.time()[[3]] - t0, v1 = r$values[1], passed = isTRUE(r$certificate$passed)) })
n <- 20000; Sp <- rsparsematrix(n, n, density = 5e-4, rand.x = function(m) rnorm(m)); Sp <- forceSymmetric(Sp + Diagonal(n, 1)); SpC <- as(as(Sp, "generalMatrix"), "CsparseMatrix")
P("time sparse 20000 k=10 RSpectra", { t0 <- proc.time()[[3]]; r <- RSpectra::eigs_sym(SpC, 10, "LM"); c(sec = proc.time()[[3]] - t0, v1 = r$values[1]) })
P("time sparse 20000 k=10 eigencore", { t0 <- proc.time()[[3]]; r <- eigs_sym(SpC, 10, "LM"); c(sec = proc.time()[[3]] - t0, v1 = r$values[1], passed = isTRUE(r$certificate$passed)) })
Xb <- matrix(rnorm(5000 * 300), 5000)
P("time svds 5000x300 k=10 RSpectra", { t0 <- proc.time()[[3]]; r <- RSpectra::svds(Xb, 10); c(sec = proc.time()[[3]] - t0, d1 = r$d[1]) })
P("time svds 5000x300 k=10 eigencore", { t0 <- proc.time()[[3]]; r <- svds(Xb, 10); c(sec = proc.time()[[3]] - t0, d1 = r$d[1], passed = isTRUE(r$certificate$passed)) })
Xe <- matrix(rnorm(3000 * 400), 3000); Ye <- matrix(rnorm(3000 * 300), 3000)
Sm <- function(v) crossprod(Xe, Ye %*% v); STm <- function(u) crossprod(Ye, Xe %*% u)
P("time op svds 400x300 k=5 RSpectra", { t0 <- proc.time()[[3]]; r <- RSpectra::svds(function(x, args) Sm(x), 5, Atrans = function(x, args) STm(x), dim = c(400, 300)); c(sec = proc.time()[[3]] - t0, d1 = r$d[1]) })
P("time op svds 400x300 k=5 eigencore", { t0 <- proc.time()[[3]]; r <- svds(linear_operator(c(400, 300), function(x, ...) Sm(x), function(x, ...) STm(x)), 5); c(sec = proc.time()[[3]] - t0, d1 = r$d[1], passed = isTRUE(r$certificate$passed)) })
