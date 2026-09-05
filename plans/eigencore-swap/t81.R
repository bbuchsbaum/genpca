suppressMessages(pkgload::load_all("/Users/bbuchsbaum/code/genpca", quiet = TRUE))
set.seed(1)
n <- 50; p <- 6; eps <- 1e-2
Qo <- qr.Q(qr(matrix(rnorm(p*p), p)))
lam <- c(1, 0.9, 0.8, 0.7, 0.6, eps)
A <- Qo %*% diag(lam) %*% t(Qo); A <- (A + t(A))/2
Z <- scale(matrix(rnorm(n*p), n), scale = FALSE)
X <- Z %*% diag(c(1,1,1,1,1, 1/eps)) %*% t(Qo)   # X'X ~ n*Qo diag(1,..,eps^-2) Qo'
Ah <- Qo %*% diag(sqrt(lam)) %*% t(Qo)
ev <- eigen(Ah %*% crossprod(X) %*% Ah, symmetric = TRUE)$values
cat("true sqrt(eigen) top3 :", signif(sqrt(ev[1:3]), 5), "\n")
f_full  <- genpca(X, A = A, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = 800)
f_trunc <- suppressWarnings(genpca(X, A = A, ncomp = 3, method = "eigen", preproc = multivarious::pass(), maxeig = 4))
cat("sdev maxeig=800       :", signif(f_full$sdev, 5), "\n")
cat("sdev maxeig=4 (trunc) :", signif(f_trunc$sdev, 5), "\n")
v1 <- multivarious::components(f_full)[,1]; v1t <- multivarious::components(f_trunc)[,1]
cat("|cos| angle of PC1 loadings full vs trunc:", abs(sum(v1*v1t))/sqrt(sum(v1^2)*sum(v1t^2)), "\n")
