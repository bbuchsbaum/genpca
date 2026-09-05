suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/genpca", quiet = TRUE))
run <- function(expr) { ws <- character()
  val <- withCallingHandlers(tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(val = val, ws = ws) }
Z <- matrix(0, 3, 3)
cat("is_spd(zero 3x3):", is_spd(Z), "\n")
r <- ensure_spd(Z, tol = 1e-6)
cat("ensure_spd(zero) min eig:", min(eigen(as.matrix(r))$values), "\n")
M <- diag(c(1, -5e-7))
cat("ensure_spd(diag(1,-5e-7), tol=1e-12) min eig:", min(eigen(as.matrix(ensure_spd(M, tol = 1e-12)))$values), "\n")
cat("clip_psd(diag(1,-5e-7)) min eig:", min(eigen(as.matrix(clip_psd(M)))$values), "\n")
set.seed(1); X <- matrix(rnorm(40*5), 40)
w <- c(1, 1, 1, 1, -5e-7)
for (m in c("eigen", "spectra", "deflation", "randomized")) {
  r <- run({ f <- genpca(X, A = w, ncomp = 2, method = m); paste0("ok  k=", length(f$sdev), " sdev=", paste(signif(f$sdev, 4), collapse = ",")) })
  cat(sprintf("A=vector(-5e-7) %-10s %s | warnings: %s\n", m, r$val, paste(unique(r$ws), collapse = " / ")))
}
for (m in c("eigen", "spectra", "deflation")) {
  r <- run({ f <- genpca(X, A = w, ncomp = 2, method = m, use_cpp = FALSE); paste0("ok  k=", length(f$sdev)) })
  cat(sprintf("A=vector(-5e-7) %-10s use_cpp=FALSE %s | warnings: %s\n", m, r$val, paste(unique(r$ws), collapse = " / ")))
}
