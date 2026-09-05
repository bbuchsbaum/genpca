pkgload::load_all('/Users/bbuchsbaum/code/genpca',quiet=TRUE)
set.seed(55);X<-matrix(rnorm(110*80),110);Y<-matrix(rnorm(110*70),110);S<-crossprod(X,Y)
for (tol in c(1e-12,1e-16,1e-20)) {
 cat('tol=',tol,'\n'); sv <- genpca:::.top_svd(function(z,args)S%*%z,2,adjoint=function(z,args)crossprod(S,z),dim=dim(S),tol=tol)
 print(list(nconv=sv$nconv,converged=sv$converged,d=sv$d))
 w<-character();f<-withCallingHandlers(gplssvd_op(X,Y,k=2,svd_opts=list(tol=tol)),warning=function(e){w<<-c(w,conditionMessage(e));invokeRestart('muffleWarning')})
 cat('gplssvd warnings=',paste(w,collapse=' | '),'d=',f$d,'\n')
}
