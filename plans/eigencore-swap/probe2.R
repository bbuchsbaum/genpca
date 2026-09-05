mode <- Sys.getenv("EC_MODE", "dev"); suppressPackageStartupMessages(library(Matrix))
if (mode == "dev") suppressMessages(pkgload::load_all("~/code/eigencore", quiet = TRUE)) else suppressPackageStartupMessages(library(eigencore, lib.loc = "cranlib"))
cat("### mode", mode, "\n")
set.seed(42); X <- matrix(rnorm(400 * 60), 400); S <- crossprod(matrix(rnorm(400 * 300), 400)) / 400; refS <- eigen(S, symmetric = TRUE)$values
# callback convention
cat("-- callback args seen:\n")
op <- try(linear_operator(dim = dim(X), apply = function(...) { a <- list(...); cat("   apply called with names:", paste(sprintf("%s<%s %s>", names(a), sapply(a, function(z) class(z)[1]), sapply(a, function(z) paste(dim(z), collapse="x"))), collapse=", "), "\n"); X %*% a[[1]] },
                          apply_adjoint = function(...) { a <- list(...); crossprod(X, a[[1]]) }), silent = TRUE)
r <- try(svds(op, k = 2), silent = TRUE); if (inherits(r, "try-error")) cat("   svds(op) error:", conditionMessage(attr(r, "condition")), "\n") else cat("   svds(op) d:", r$d, "\n")
# matrix-free hermitian accuracy
ops <- linear_operator(c(300, 300), function(x, ...) S %*% x, function(x, ...) S %*% x, structure = hermitian())
for (tl in c(1e-8, 1e-12)) { r <- eigs_sym(ops, 3, "LA", tol = tl); cat(sprintf("-- matrix-free hermitian tol=%g: max value err %.2e, cert passed %s, max residual %s, method %s\n", tl, max(abs(r$values - refS[1:3])), isTRUE(r$certificate$passed), format(r$certificate$max_residual, digits=3), r$diagnostics$method %||% NA)) }
r <- eigs_sym(ops, 3, "LA"); cat("-- diagnostics names:", paste(names(r$diagnostics), collapse=", "), "\n"); cat("-- certificate names:", paste(names(r$certificate), collapse=", "), "\n")
str(r$diagnostics, max.level = 1, give.attr = FALSE)
r2 <- eigs_sym(ops, 3, "LA", method = lanczos()); cat(sprintf("-- matrix-free lanczos(): err %.2e\n", max(abs(r2$values - refS[1:3]))))
# seed= and RNG
set.seed(7); a <- runif(1); set.seed(7); invisible(eigs_sym(S, 3, "LA", seed = 1L)); invisible(svds(X, 3, seed = 1L)); b <- runif(1)
cat("-- with seed=1: .Random.seed consumed?", !identical(a, b), "\n")
v1 <- eigs_sym(S, 3, "LA", seed = 1L)$vectors; v2 <- eigs_sym(S, 3, "LA", seed = 1L)$vectors; cat("-- seed=1 repeat vectors identical?", identical(v1, v2), "\n")
# timing dense with/without certify
Sd <- crossprod(matrix(rnorm(1600 * 1500), 1600)) / 1600
for (cf in c(TRUE, FALSE)) { t0 <- proc.time()[[3]]; r <- eigs_sym(Sd, 10, "LA", certify = cf); cat(sprintf("-- dense 1500 k=10 certify=%s: %.2fs\n", cf, proc.time()[[3]] - t0)) }
t0 <- proc.time()[[3]]; r <- eigs_sym(Sd, 10, "LA", method = lanczos()); cat(sprintf("-- dense 1500 k=10 method=lanczos(): %.2fs err %.1e\n", proc.time()[[3]] - t0, max(abs(r$values - eigen(Sd, TRUE, TRUE)$values[1:10]))))
cat("-- plan for dense 1500:\n"); print(plan_solver(eigen_problem(Sd, structure = hermitian(), target = largest()), k = 10))
