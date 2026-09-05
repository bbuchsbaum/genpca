suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/genpca", quiet = TRUE))
set.seed(1); X <- matrix(rnorm(40*2), 40)
A_ind <- matrix(c(1, 2, 2, 1), 2)   # eigenvalues -1, 3: indefinite
msgs <- character()
f <- withCallingHandlers(genpca(X, A = A_ind, ncomp = 1),
  warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("indefinite A, default remedy: warnings/messages =", length(msgs), "\n")
cat("A actually used:\n"); print(as.matrix(f$A))
cat("min eig of A used:", min(eigen(as.matrix(f$A))$values), "\n")
# randomized backend at small scale: is the absolute jitter the cause?
set.seed(1); X <- matrix(rnorm(60*8), 60); Am <- crossprod(matrix(rnorm(8*8), 8))/8
ref <- genpca(X, A = Am, ncomp = 4, method = "eigen", preproc = multivarious::pass())$sdev
for (jit in c(1e-10, 1e-20)) {
  f <- suppressWarnings(genpca(X * 1e-5, A = Am, ncomp = 4, method = "randomized", preproc = multivarious::pass(),
                               seed_randomized = 1L, jitter_metric = jit))
  cat(sprintf("randomized c=1e-5 jitter_metric=%g: sdev/c = %s  (eigen ref: %s)\n", jit,
              paste(signif(f$sdev / 1e-5, 4), collapse = ","), paste(signif(ref, 4), collapse = ",")))
}
