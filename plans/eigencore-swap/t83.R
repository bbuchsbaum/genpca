suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/genpca", quiet = TRUE))
S <- matrix(c(2, 0.9, 0.1, 2), 2)   # S[1,2]=0.1 (upper), S[2,1]=0.9 (lower)
cat("ensure_spd(asym) ->\n"); print(as.matrix(ensure_spd(S)))
set.seed(1); X <- matrix(rnorm(40*2), 40)
msgs <- character()
f <- withCallingHandlers(genpca(X, A = S, ncomp = 1),
  warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
  message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
cat("genpca(X, A = asym) default remedy: warnings/messages =", length(msgs), "\n")
fu <- genpca(X, A = matrix(c(2,0.1,0.1,2),2), ncomp = 1); fl <- genpca(X, A = matrix(c(2,0.9,0.9,2),2), ncomp = 1)
fa <- genpca(X, A = matrix(c(2,0.5,0.5,2),2), ncomp = 1)
cat("sdev asym:", f$sdev, " upper-sym:", fu$sdev, " lower-sym:", fl$sdev, " averaged:", fa$sdev, "\n")
C <- crossprod(X); Ca <- C; Ca[1,2] <- C[1,2] + 1.0   # perturb upper only
cat("genpca_cov(asym C) d:", genpca_cov(Ca, ncomp = 1)$d, "  from upper:", genpca_cov(Matrix::forceSymmetric(Ca, "U"), ncomp = 1)$d,
    "  from lower:", genpca_cov(Matrix::forceSymmetric(Ca, "L"), ncomp = 1)$d, "\n")
cat("genpca_cov(asym C, geigen) d:", genpca_cov(Ca, ncomp = 1, method = "geigen")$d, "\n")
