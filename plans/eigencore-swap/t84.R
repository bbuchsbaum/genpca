suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/genpca", quiet = TRUE))
run <- function(expr) { ws <- character()
  val <- withCallingHandlers(tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e))),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") })
  list(val = val, ws = ws) }
set.seed(1); X <- matrix(rnorm(60*8), 60)
Am <- crossprod(matrix(rnorm(8*8), 8))/8; Mm <- diag(runif(60, 0.5, 1.5))
for (c in c(1, 1e-3, 1e-5, 1e-7)) for (m in c("eigen", "spectra", "deflation", "randomized")) {
  r <- run({ f <- genpca(X * c, A = Am, M = Mm, ncomp = 4, method = m, preproc = multivarious::pass(), seed_randomized = 1L)
             paste0("k=", length(f$sdev), " sdev/c=", paste(signif(f$sdev / c, 4), collapse = ",")) })
  cat(sprintf("c=%-6g %-10s %s | warn: %s\n", c, m, r$val, paste(unique(substr(r$ws, 1, 60)), collapse = " / ")))
}
