#' @keywords internal
nipals_deflation_gs <- function(Xmat, Ymat, ncomp, maxiter=200, tol=1e-9, verbose=FALSE) {
  n <- nrow(Xmat)
  p <- ncol(Xmat)
  q <- ncol(Ymat)
  
  loadX <- matrix(0, p, ncomp)  # P
  loadY <- matrix(0, q, ncomp)  # Q
  scoreX<- matrix(0, n, ncomp)  # T
  scoreY<- matrix(0, n, ncomp)  # U
  
  Xres <- Xmat
  Yres <- Ymat
  
  for(k in seq_len(ncomp)) {
    if(verbose) cat(sprintf("Extracting comp %d/%d...\n", k, ncomp))
    
    # initialize u => col j in Yres with max var
    var_ycols <- apply(Yres, 2, var)
    jmax <- which.max(var_ycols)
    u_n <- Yres[, jmax, drop=FALSE]
    
    for(iter in seq_len(maxiter)) {
      # (1) p_k
      denom_u <- sum(u_n^2)
      if(denom_u<1e-15) break
      p_k <- crossprod(Xres, u_n)/ denom_u
      
      # Gram–Schmidt p_k
      for(j in seq_len(k-1)) {
        dot_pp<- sum(p_k * loadX[, j])
        p_k <- p_k - dot_pp* loadX[, j]
      }
      norm_p<- sqrt(sum(p_k^2))
      if(norm_p>1e-15) p_k<- p_k/norm_p else p_k[]=0
      
      # (2) t_k
      denom_p<- sum(p_k^2)
      if(denom_p<1e-15) break
      t_n <- Xres %*% p_k / denom_p
      
      # Gram–Schmidt t_n
      for(j in seq_len(k-1)) {
        dot_tt<- sum(t_n * scoreX[, j])
        t_n <- t_n - dot_tt* scoreX[, j]
      }
      norm_t<- sqrt(sum(t_n^2))
      if(norm_t>1e-15) t_n<- t_n/norm_t else t_n[]=0
      
      # (3) q_k
      denom_t<- sum(t_n^2)
      if(denom_t<1e-15) break
      q_k <- crossprod(Yres, t_n)/ denom_t
      
      # Gram–Schmidt q_k
      for(j in seq_len(k-1)) {
        dot_qq<- sum(q_k * loadY[, j])
        q_k <- q_k - dot_qq* loadY[, j]
      }
      norm_q<- sqrt(sum(q_k^2))
      if(norm_q>1e-15) q_k<- q_k/norm_q else q_k[]=0
      
      # (4) u_n new
      denom_q<- sum(q_k^2)
      if(denom_q<1e-15) break
      u_n_new <- Yres %*% q_k / denom_q
      
      diff_val<- sqrt(sum((u_n_new - u_n)^2))
      u_n <- u_n_new
      if(diff_val<tol) break
    }
    
    loadX[,k] <- p_k
    loadY[,k] <- q_k
    scoreX[,k]<- t_n
    scoreY[,k]<- u_n
    
    # deflate
    Xres <- Xres - t_n %*% t(p_k)
    Yres <- Yres - t_n %*% t(q_k)
  }
  
  list(P=loadX, Q=loadY, T=scoreX, U=scoreY)
}


#' Generalized Partial Least Squares with Deflation + Gram-Schmidt, Adaptive Rank by Default
#'
#' This function forms the embedded data \eqn{\tilde{X} = \sqrt{M_x}\,X_{\mathrm{proc}}\,\sqrt{A_x}}
#' and \eqn{\tilde{Y} = \sqrt{M_y}\,Y_{\mathrm{proc}}\,\sqrt{A_y}} explicitly in memory, then
#' performs a multi-component NIPALS with rank-1 \emph{deflation} \emph{and} a Gram–Schmidt
#' (GS) step to keep each new component orthonormal to previously extracted components.
#'
#' By default, we do \emph{adaptive} partial eigen expansions for each constraint matrix:
#' \itemize{
#'   \item If \code{rank_Mx} (etc.) is a positive integer, we use that rank directly.
#'   \item Otherwise (if \code{NULL}, \code{NA}, or \code{0}), we do an initial partial eigen
#'         up to \code{max_k}, then keep as many eigenvalues as needed to capture
#'         \code{var_threshold} fraction of the total variance. This approach is repeated
#'         for \code{Mx}, \code{Ax}, \code{My}, \code{Ay}.
#' }
#'
#' The multi-factor solution avoids the repeated-first-factor problem and yields more stable
#' loadings for each component.
#'
#' @param X numeric matrix \eqn{(n \times p)}.
#' @param Y numeric matrix \eqn{(n \times q)}.
#' @param Mx,Ax row/column constraints for \eqn{X}, each either diagonal/identity
#'        or PSD of size \eqn{(n \times n)/(p \times p)}. If \code{Mx} is \code{NULL},
#'        we default to identity. If \code{rank_Mx} is \code{NULL}, \code{NA}, or \code{0},
#'        the adaptive approach is used.
#' @param My,Ay row/column constraints for \eqn{Y}, each either diagonal/identity
#'        or PSD of size \eqn{(n \times n)/(q \times q)}. Same adaptive logic as above.
#' @param ncomp integer, number of factors (components) to extract.
#'   After the dimensions `p` and `q` are determined, oversized
#'   values are reduced to `min(ncomp, p, q)`.
#'
#' @param preproc_x, preproc_y \code{\link[multivarious]{pre_processor}} objects or
#'        \code{\link[multivarious]{pass}} for no preprocessing.
#'
#' @param rank_Mx, rank_Ax, rank_My, rank_Ay integer or \code{NULL}, \code{NA}, \code{0}
#'        controlling the partial eigen expansions for each constraint matrix.
#'        If a positive integer, that rank is used. Otherwise, we use an adaptive approach
#'        capturing \code{var_threshold} fraction of variance, capped by \code{max_k}.
#'
#' @param var_threshold numeric fraction of total variance to capture if ranks are \code{NULL}/\code{NA}/\code{0}.
#'        Defaults to 0.99 (99\%).
#' @param max_k integer maximum cap for the partial eigen expansions. Defaults to 200.
#'
#' @param maxiter, tol numeric parameters for the internal NIPALS iteration.
#' @param verbose logical; if TRUE, prints progress messages.
#' @param ... additional arguments stored in the returned object.
#'
#' @return A \code{cross_projector} object of class \code{c("genpls","cross_projector","projector")},
#' whose loadings \code{vx, vy} are distinct for each factor and orthonormal
#' across components (within numerical tolerance).
#'
#' @importFrom PRIMME eigs_sym
#' @importFrom Matrix diag forceSymmetric crossprod isDiagonal
#' @importFrom multivarious prep pass init_transform
#' @export
genpls <- function(X, Y,
                   Mx=NULL, Ax=NULL,
                   My=NULL, Ay=NULL,
                   ncomp=2,
                   preproc_x=pass(),
                   preproc_y=pass(),
                   rank_Mx=NULL, rank_Ax=NULL,
                   rank_My=NULL, rank_Ay=NULL,
                   var_threshold=0.99,
                   max_k=200,
                   maxiter=200, tol=1e-9,
                   verbose=FALSE,
                   ...)
{
  # 1) Basic checks, same as before
  if(!is.matrix(X) || !is.numeric(X))
    stop("X must be a numeric matrix.")
  if(!is.matrix(Y) || !is.numeric(Y))
    stop("Y must be a numeric matrix.")
  if(nrow(X)!=nrow(Y))
    stop("X,Y must have the same number of rows.")
  if(ncomp<1 || floor(ncomp)!=ncomp)
    stop("ncomp must be a positive integer.")
  
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)
  approx_rank <- min(n,p,q)
  if(ncomp > approx_rank && verbose) {
    warning(sprintf("Requested ncomp=%d > min(n,p,q)=%d => degenerate?", ncomp, approx_rank))
  }
  ncomp <- min(ncomp, p, q)  # cannot exceed dimensions
  
  # 2) Default constraints => identity if missing
  if(is.null(Mx)) Mx <- Matrix::Diagonal(n)
  if(is.null(Ax)) Ax <- Matrix::Diagonal(p)
  if(is.null(My)) My <- Matrix::Diagonal(n)
  if(is.null(Ay)) Ay <- Matrix::Diagonal(q)
  
  # 3) Preprocess
  px_obj <- prep(preproc_x)
  py_obj <- prep(preproc_y)
  Xp <- init_transform(px_obj, X)
  Yp <- init_transform(py_obj, Y)
  
  # 4) Build approximate sqrt transforms => from plsutils
  Mx_sqrt <- build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol)
  Ax_sqrt <- build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
  My_sqrt <- build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol)
  Ay_sqrt <- build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
  
  # 5) Build Xtilde, Ytilde in memory => from plsutils row_transform, col_transform
  Xrw <- row_transform(Xp, Mx_sqrt)
  Xtilde <- col_transform(Xrw, Ax_sqrt)
  Yrw <- row_transform(Yp, My_sqrt)
  Ytilde <- col_transform(Yrw, Ay_sqrt)
  
  nipres <- nipals_deflation_gs(Xtilde, Ytilde, ncomp, maxiter, tol, verbose)
  
  Ptilde <- nipres$P
  Qtilde <- nipres$Q
  Ttilde <- nipres$T
  Utilde <- nipres$U
  
  # 7) build_invsqrt_mult => from plsutils
  Ax_invsqrt <- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
  Ay_invsqrt <- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
  
  vx <- sapply(seq_len(ncomp), function(j) Ax_invsqrt(Ptilde[,j]))
  vy <- sapply(seq_len(ncomp), function(j) Ay_invsqrt(Qtilde[,j]))
  
  out <- cross_projector(
    vx=vx,
    vy=vy,
    preproc_x=px_obj,
    preproc_y=py_obj,
    classes=c("genpls"),
    tilde_Px=Ptilde,
    tilde_Py=Qtilde,
    tilde_Tx=Ttilde,
    tilde_Ty=Utilde,
    ncomp=ncomp,
    rank_Mx=rank_Mx,
    rank_Ax=rank_Ax,
    rank_My=rank_My,
    rank_Ay=rank_Ay,
    var_threshold=var_threshold,
    max_k=max_k,
    maxiter=maxiter,
    tol=tol,
    verbose=verbose,
    ...
  )
  out
}