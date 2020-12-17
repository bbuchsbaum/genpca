
#' @keywords internal
prep_constraints <- function(X, A, M) {
  if (is.null(A)) {
    A <- sparseMatrix(i=1:ncol(X), j=1:ncol(X),x=rep(1, ncol(X)))
  }
  
  if (is.null(M)) {
    M <- sparseMatrix(i=1:nrow(X), j=1:nrow(X),x=rep(1,nrow(X)))
  }
  
  if (is.vector(A)) {
    assert_that(length(A) == ncol(X))
    A <- sparseMatrix(i=1:length(A), j=1:length(A),x=A)
  } else {
    assert_that(nrow(A) == ncol(X), msg=paste("nrow(A) != ncol(X) -- ", nrow(A), " != ", ncol(X)))
    assert_that(ncol(A) == ncol(X), msg=paste("ncol(A) != ncol(X) -- ", ncol(A), " != ", ncol(X)))
  }
  
  if (is.vector(M)) {
    assert_that(length(M) == nrow(X))
    M <- sparseMatrix(i=1:length(M), j=1:length(M),x=M)
  } else {
    assert_that(nrow(M) == nrow(X))
    assert_that(ncol(M) == nrow(X))
  }
  
  list(A=A, M=M)
  
}

#' Generalized Principal Components Analysis
#' 
#' Compute a PCA in a inner-product space defined by row and column constraint matrices.
#' 
#' @param X the data matrix
#' @param A the column constraints. Can be a \code{vector}, symmetric \code{matrix}, or symmetric sparse matrix with \code{ncol(X)} rows and columns.
#' @param M the row constraints. Can be a \code{vector}, symmetric \code{matrix}, or symmetric sparse matrix with \code{nrow(X)} rows and columns.
#' @param ncomp the number of components to return
#' @param preproc a pre-processing function (e.g. `center()` or `standardize()` from `multivarious` package, see Details)
#' @param deflation use the iterative deflation method
#' @param threshold convergence threshold for deflation method
#' @param use_cpp use the C++ implementation
#' @importFrom assertthat assert_that
#' @importFrom Matrix sparseMatrix t isDiagonal
#' @importFrom multivarious bi_projector init_transform prep
#' 
#' @return an instance of type `genpca`, extending `bi_projector`
#' @export
#' 
#' @details Pre-processing such as scaling and centering is carried out using special functions from the `multivarious` package. 
#' Basic options are `multivarious::center()`, `multivarious::standardize()`, and for no pre-processing at all, `multivarious::pass()`.
#' 
#' 
#' @references 
#' 
#' Abdi, H. (2007). Singular value decomposition (SVD) and generalized singular value decomposition. \emph{Encyclopedia of measurement and statistics}, 907-912.
#' 
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014). A generalized least-square matrix decomposition. \emph{Journal of the American Statistical Association}, 109(505), 145-159.
#' 
#' 
#' @examples 
#' X <- matrix(rnorm(100*100), 100,100)
#' A <- cov(X)
#' M <- cov(t(X))
#' gp1 <- genpca(X, A=A, M=M, ncomp=100)
#' 
#' Xrecon <- reconstruct(gp1)
#' 
#' @import multivarious
#' @importFrom assertthat assert_that
genpca <- function(X, A=NULL, M=NULL, ncomp=min(dim(X)), 
                   preproc=multivarious::center(), deflation=FALSE, 
                   threshold=1e-06, 
                   use_cpp=TRUE, maxeig=800) {
  
  
  pcon <- prep_constraints(X, A, M)
 
  A <- pcon$A
  M <- pcon$M
  
  procres <- multivarious::prep(preproc)
  Xp <- multivarious::init_transform(procres, X)
  
  assert_that(ncomp > 0)
  ncomp <- min(min(dim(Xp)), ncomp)
  
  n = nrow(Xp)
  p = ncol(Xp)
  
  
  if (deflation) {
    if (n < p) {
      
      if (use_cpp) {
        svdfit <- gmd_deflation_cpp(t(Xp), A, M, ncomp, threshold)
        #svdfit <- gmd_deflationR(t(Xp), A, M, ncomp, threshold,svd_init)
        svdfit$d <- svdfit$d[,1]
      } else {
        svdfit <- gmd_deflationR(t(Xp), A, M, ncomp, threshold)
      }
      ## swap v and u ...
      svdfit <- list(u=svdfit$v, v=svdfit$u,d=svdfit$d, cumv=svdfit$cumv,propv=svdfit$propv)
    } else {
      if (use_cpp) {
        svdfit <- gmd_deflation_cpp(Xp, M, A, ncomp, threshold)
        svdfit$d <- svdfit$d[,1]
      } else {
        svdfit <- gmd_deflationR(Xp, M, A, ncomp, threshold)
      }
    }
  } else { 
    if(n < p){
      ret = gmdLA(t(Xp), A,M, ncomp,p,n, maxeig=maxeig)
      ## swap v and u ...
      svdfit <- list(u=ret$v, v=ret$u,d=ret$d, cumv=ret$cumv,propv=ret$propv)
    } else{
      svdfit <- gmdLA(Xp, M,A,ncomp,n,p, maxeig=maxeig)
    }
  }
  

  scores <- t(t(as.matrix(M %*% svdfit$u)) * svdfit$d)
  #col_scores <- t(t(as.matrix(A %*% svdfit$v)) * svdfit$d)
  
  #scores <- t(t(as.matrix(svdfit$u)) * svdfit$d)
  row.names(scores) <- row.names(X)
  #norm_loadings <- t(t(as.matrix(svdfit$v)) * svdfit$d)
  
  ret <- bi_projector(
    v = A %*% svdfit$v,
    s = scores,
    sdev=svdfit$d, 
    preproc=procres,
    ncomp=length(svdfit$d),
    ## generalized singular vector v
    ov=svdfit$v, 
    ## generalized singular vector u
    ou=svdfit$u, 
    u=M %*% svdfit$u,
    #col_scores=col_scores,
    classes=c("genpca"),
    A=A,
    M=M)
  
  ret
}




#' @keywords internal
gmdLA <- function(X, Q, R, k=min(n,p), n, p, maxeig=800) {
  ##computation
  if (isDiagonal(R)) {
    Rtilde <- Matrix::Diagonal(x=sqrt(Matrix::diag(R)))
    Rtilde.inv = Matrix::Diagonal(x=1/sqrt(Matrix::diag(R)))
  } else {
    decomp <- if (!is.null(attr(R, "decomp"))) {
      ## cached decomposition (should validate)
      attr(R, "decomp")
    } else {
      if (nrow(R) > maxeig) {
        ## large matrix, truncate for speed
        RSpectra::eigs_sym(R, k=maxeig)
      } else {
        eigen(R, symmetric=TRUE)
      }
    }
    
    if (length(decomp$values) > 1) {
      v <- decomp$values
      if (sum(v < 0) > 1) {
        warning(paste("genpca: removing ", sum(v<0), 
                      " negative eigenvalues when computing inverse of constraint matrix."))
      }
    }
    
    keep <- which(decomp$values > 0 & (abs(decomp$values) > 1e-6))
    
    decomp$vectors <- decomp$vectors[, keep]
    decomp$values <- decomp$values[keep]
    
    ## R^(1/2)
    Rtilde <- decomp$vectors %*% diag(sqrt(decomp$values)) %*% t(decomp$vectors)
    
    inv.values = 1 / sqrt(decomp$values)
    ## R^(-1/2)
    Rtilde.inv = decomp$vectors %*% diag(inv.values) %*% t(decomp$vectors)
  }
  
  ## X' %*% Q %*% X
  inmat <- Matrix::crossprod(X, Q) %*% X
  
  ## R^(1/2) %*% X' %*% Q %*% X %*% R^(1/2) 
  RtinRt <- Rtilde %*% inmat %*% Rtilde
  
  XR <- X %*% R
  ## nXn * n*n * pXp
  ## R %*% X' %*% Q %*% X %*% R
  RnR <- R %*% inmat %*% R
  
  xtilde.decomp <- if (k == min(n,p)) {
    eigen(RtinRt)
  } else {
    #RSpectra::eigs_sym(RtinRt, k=k)
    ret <- RSpectra::eigs_sym(RtinRt, k=k+1)
    list(vectors=ret$vectors, values=ret$values)
  }
  
  keep <- which(abs(xtilde.decomp$values) > 1e-6 & (xtilde.decomp$values > 0 ))
  k <- min(k, length(keep))
  xtilde.decomp$vectors <- xtilde.decomp$vectors[, 1:k]
  xtilde.decomp$values <- xtilde.decomp$values[1:k]
  
  #Rtilde.inv %*% xtilde.decomp$vectors
  
  ## R^(-1/2) %*% v
  vgmd <- Rtilde.inv %*% xtilde.decomp$vectors
  dgmd <- sqrt(xtilde.decomp$values[1:k])
  ugmd <- matrix(nrow = n, ncol = k)
  cumv <- rep(0, k)
  propv <- dgmd ^ 2 / sum(diag(as.matrix(inmat %*% R)))
  normalizing.number <- 1
  
  
  for (i in 1:k) {
    normalizing.number = sqrt(vgmd[, i] %*% RnR %*% vgmd[, i])
    ugmd[, i] = as.matrix(XR %*% vgmd[, i]) / as.double(normalizing.number)
    cumv[i] = sum(propv[1:i])
  }
  
  list(
    u = ugmd[, 1:k, drop=FALSE],
    v = vgmd[, 1:k, drop=FALSE],
    d = dgmd[1:k],
    k = k,
    cumv = cumv,
    propv = propv
  )
  
}

#' @importFrom multivarious reconstruct
reconstruct.genpca <- function(x, 
                               comp=1:ncomp(x), 
                               rowind=1:nrow(scores(x)), colind=1:nrow(components(x))) {
  
  ## X = FV
  ## F = M(XAV)
  ## F = MUD
  ## X = M(AV')(F')
  
  #Xr = x$M %*% (Fp %*% t(x$A %*% x$ov))
  
  ## Xr = x$ou %*% diag(x$d) %*% t(x$ov)
  ## does not work when X is uncentered
  
  assert_that(max(comp) <= ncomp(x))
  
  out <- x$ou[rowind,comp,drop=FALSE] %*% diag(sdev(x)[comp]) %*% t(x$ov[,comp,drop=FALSE])[,colind]
  reverse_transform(x$preproc, out)
  
}



gmd_deflationR <- function(X, Q, R, k, thr = 1e-6) {
  
  n = nrow(X)
  p = ncol(X)
  #arma::mat ugmd(n, k, fill::zeros);
  #arma::mat vgmd(p, k, fill::zeros);
  
  ugmd = matrix(nrow=n, ncol = k)
  vgmd = matrix(nrow=p, ncol = k)
  dgmd = rep(0,k)
  propv = rep(0,k)
  Xhat <- X
  

  u <- cbind(rnorm(n))
  v <- cbind(rnorm(p))
  
  
  #browser()
  print("qnorm")
  
  #if (p > n) {
  #  qrnorm = trace(X * R * X.t() * Q);
  #} else {
  #  qrnorm = trace(X.t() * Q * X * R);
  #}
  
  
  qrnorm = sum(Matrix::diag(Matrix::crossprod(X,Q) %*% X %*% R))
  cumv = rep(0,k)
  
  #print("begin loop")
  for(i in 1:k){
    print(i)
    err <- 1
    
    #if (svd_init && i > 1) {
    #  init <- rsvd::rsvd(Xhat,k=1)
    #  v <- init$v[,1,drop=FALSE]
    #} 
    
    while(err > thr){
      oldu = u
      oldv = v
      
      uhat = Xhat %*% (R %*% v)
      #u = uhat/as.double(sqrt(t(uhat)%*% Q %*% uhat))
      u = uhat/as.double(sqrt(Matrix::crossprod(uhat, Q) %*% uhat))
      #vhat = t(Xhat) %*% Q %*% u
      #vhat2 <- t(Xhat) %*% (Q %*% u)
      vhat <- Matrix::crossprod(Xhat, (Q %*% u)) 
      #v = vhat/as.double(sqrt(t(vhat) %*% R %*% vhat))
      v <- vhat / as.double(sqrt(Matrix::crossprod(vhat, R) %*% vhat))
      err = as.numeric(t(oldu - u) %*% (oldu - u) + t(oldv -v ) %*% (oldv - v))
      
      #print(paste("err: ", err))
      if (is.na(err)) {
        message("error is NA, convergence failed")
        message(paste("t(oldu - u) %*% (oldu - u)", t(oldu - u) %*% (oldu - u)))
        message(paste("t(oldv -v ) %*% (oldv - v)", t(oldv -v ) %*% (oldv - v)))
        stop()
      }
    }
    
    
    dgmd[i] <- Matrix::crossprod(u, Q) %*% X %*% (R %*% v)
    ugmd[,i] <- u[,1]
    vgmd[,i] <- v[,1]
    Xhat <- Xhat - dgmd[i] *  u %*% t(v)
    propv[i] <- dgmd[i]^2/as.double(qrnorm)
    cumv[i] <- sum(propv[1:i])
  }
  
  list(d=as.vector(dgmd), v=vgmd, u=ugmd, cumv=cumv, propv=propv)
  
}

