pkgload::load_all('.', quiet = TRUE)

cat('Explicit clip retains a negative eigenvalue\n')
B <- repair_metric(diag(c(1, -1e-10)), method = 'clip')
print(diag(B))
print(attr(B, 'repair_report'))

cat('\nSingular generalized eigenproblem residual\n')
C <- matrix(c(2, 1, 1, 2), 2)
R <- diag(c(1, 0))
f <- geigen_cov(C, R, ncomp = 1)
print(list(v = f$v, lambda = f$lambda,
           residual = C %*% f$v - R %*% f$v * f$lambda))

cat('\nCommuting C and R can select different leading directions\n')
C <- diag(c(4, 1)); R <- diag(c(9, 1))
print(list(gmd = genpca_cov(C, R, ncomp = 1)$v,
           geigen = geigen_cov(C, R, ncomp = 1)$v))
