pkgload::load_all('/Users/bbuchsbaum/code/genpca', quiet=TRUE)
cat('EIGENCORE', as.character(packageVersion('eigencore')), '\n')
cat('\nRANK CUTOFF\n')
X <- diag(c(10,1,0.1))
for (m in c('eigen','spectra','randomized','deflation')) for (cpp in if(m=='deflation') c(TRUE,FALSE) else TRUE) {
 set.seed(42)
 f <- suppressWarnings(genpca(X, ncomp=3, method=m, rank_rtol=.2, use_cpp=cpp))
 cat(m, 'cpp=',cpp, 'd=', f$sdev, '\n')
}
cat('\nSMALL POSITIVE METRIC DIRECTION\n')
X <- diag(c(1, 1e6)); A <- diag(c(1, 1e-10))
for (m in c('eigen','spectra','randomized','deflation')) {
 f <- suppressWarnings(genpca(X,A=A,ncomp=1,method=m))
 cat(m,'d=',f$sdev,'ov=', f$ov, 'ou=',f$ou, 'VAV=',as.numeric(crossprod(f$ov,A%*%f$ov)), '\n')
 cat('project-score error=',max(abs(multivarious::project(f,X)-multivarious::scores(f))), '\n')
}
f <- genpca_cov(crossprod(X), R=A, ncomp=1)
cat('cov d=',f$d, 'v=',as.numeric(f$v),'VAV=',as.numeric(t(f$v)%*%A%*%f$v),'\n')
cat('\nMLE DELTA WITHOUT RESCALE\n')
set.seed(92); X <- matrix(rnorm(120),20)
f <- gpca_mle(X,ncomp=2,max_iter=1)
cat('delta=',f$loglik_rescale_delta,'loglik=',f$loglik,'path=',f$loglik_path,'\n')
