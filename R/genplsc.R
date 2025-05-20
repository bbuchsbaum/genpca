#' Generalised Partial‑Least‑Squares **Correlation** (GPLS‑C)
#' 
#' `genplsc()` extracts latent **correlation directions** between two data
#' sets \eqn{X\in\mathbb{R}^{n\times p}} and \eqn{Y\in\mathbb{R}^{n\times q}}
#' while respecting *row* and *column* metric constraints
#' \eqn{M_x,\,A_x,\,M_y,\,A_y}.  
#' It unifies the “direct cross‑product” formulation and the
#' “row/column‑likeness operator” trick in a single function:
#' 
#' \itemize{
#'   \item **`dual = FALSE`** – explicitly forms the embedded cross‑product
#'         \eqn{C = (\sqrt{M_x}XA_x^{1/2})^\top\,(\sqrt{M_y}Y A_y^{1/2})}
#'         and does (partial) SVD on the \eqn{p\times q} matrix.
#'   \item **`dual = TRUE`**  – avoids the large \eqn{p\times q} object by
#'         projecting onto row‑ or column‑space operators (à la *Allen
#'         et al., 2014*), chosen automatically or via `force_row_likeness`
#'         / `force_col_likeness`.
#' }
#' 
#' For every constraint matrix you may either
#' 
#' * give a **fixed rank** (`rank_Mx = 50`) or  
#' * leave it `NULL` / `NA` / `0` to let an **adaptive** scheme keep the
#'   smallest eigen‑subspace that captures `var_threshold` of its total
#'   variance (capped by `max_k`).
#' 
#' A truncated SVD is used when `svd_method = "RSpectra"` or `"irlba"` and the
#' package is installed; otherwise base `svd()` is used.
#' 
#' @section Returned object:
#' A [`cross_projector`][multivarious::cross_projector] with slots
#' 
#' \describe{
#'   \item{`vx`, `vy`}{original‑space loadings \eqn{p\times k}, \eqn{q\times k}}
#'   \item{`Tx`, `Ty`}{embedded factor scores \eqn{n\times k}}
#'   \item{`d_full`}{all singular values (length \eqn{\le k})}
#' }
#' and a class vector `c("genplsc", "cross_projector", "projector")`.
#' 
#' @inheritParams build_sqrt_mult
#' @inheritParams build_invsqrt_mult
#' 
#' @param X,Y              Numeric matrices with identical numbers of rows.
#' @param Mx,Ax,My,Ay      Row/column constraint matrices (NULL ⇒ identity).
#' @param ncomp            Number of latent factors.  If \eqn{\le 0},
#'                         uses \eqn{\min(p,q)}.
#' @param preproc_x,preproc_y  [`multivarious`] preprocessing objects
#'                         (e.g. `center()`, `standardize()`, `pass()`).
#' @param rank_Mx,rank_Ax,rank_My,rank_Ay
#'                         Fixed ranks; `NULL`/`NA`/`0` ⇒ adaptive.
#' @param var_threshold    Fraction of variance to keep in adaptive mode.
#' @param max_k            Hard cap for adaptive sub‑space dimension.
#' @param dual             `FALSE` for direct cross‑product; `TRUE` for
#'                         operator trick.
#' @param force_row_likeness,force_col_likeness
#'                         Override automatic choice when `dual = TRUE`.
#' @param svd_method       `"base"`, `"RSpectra"`, or `"irlba"`.
#' @param tol              Numerical tolerance (e.g. tiny eigenvalues).
#' @param verbose          Logical; print progress if `TRUE`.
#' @param ...              Ignored; added for forward compatibility.
#' 
#' @references
#' Allen, G. I., Grosenick, L., & Taylor, J. (2014).  
#' *A Generalized Least‑Squares Matrix Decomposition*. **JASA**, 109(505), 145–159.
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 50; p <- 40; q <- 30
#' X <- matrix(rnorm(n*p), n, p)
#' Y <- matrix(rnorm(n*q), n, q)
#' 
#' # unconstrained, direct
#' fit1 <- genplsc(X, Y, ncomp = 5)
#' 
#' # unconstrained, operator trick (row likeness here)
#' fit2 <- genplsc(X, Y, ncomp = 5, dual = TRUE)
#' 
#' # diagonal row constraints and adaptive rank
#' Mx <- diag(abs(rnorm(n)))
#' My <- diag(abs(rnorm(n)))
#' fit3 <- genplsc(X, Y, Mx = Mx, My = My, dual = TRUE,
#'                 var_threshold = 0.95, ncomp = 4,
#'                 svd_method = "RSpectra")
#' }
#' 
#' @aliases genplscorr genplscorr2
#' @export
genplsc <- function(X, Y,
                             Mx=NULL, Ax=NULL, My=NULL, Ay=NULL,
                             ncomp = 2,
                             preproc_x = pass(), preproc_y = pass(),
                             rank_Mx=NULL, rank_Ax=NULL,
                             rank_My=NULL, rank_Ay=NULL,
                             var_threshold = .99, max_k = 200,
                             dual = FALSE,
                             force_row_likeness = NULL,
                             force_col_likeness = NULL,
                             svd_method = c("base","RSpectra","irlba"),
                             tol = 1e-9,
                             verbose = FALSE)
{
  ## ---- 0. sanity ---------------------------------------------------
  if(!is.numeric(X) || !is.matrix(X)) stop("X must be a numeric matrix")
  if(!is.numeric(Y) || !is.matrix(Y)) stop("Y must be a numeric matrix")
  if(nrow(X) != nrow(Y)) stop("X and Y must have the same number of rows")
  if(ncomp <= 0) ncomp <- min(ncol(X), ncol(Y))

  svd_method <- match.arg(svd_method)
  n <- nrow(X); p <- ncol(X); q <- ncol(Y)

  if(!is.null(force_row_likeness) && !is.null(force_col_likeness))
    stop("Specify only one of force_row_likeness or force_col_likeness")

  ## ---- 1. defaults for constraints --------------------------------
  if(is.null(Mx)) Mx <- Matrix::Diagonal(n)
  if(is.null(My)) My <- Matrix::Diagonal(n)
  if(is.null(Ax)) Ax <- Matrix::Diagonal(p)
  if(is.null(Ay)) Ay <- Matrix::Diagonal(q)

  ## ---- 2. preprocess ----------------------------------------------
  px <- prep(preproc_x);  py <- prep(preproc_y)
  Xp <- init_transform(px, X)
  Yp <- init_transform(py, Y)

  ## ---- 3. helper closures & local SVD function ---------------------
  Mx_sqrt <- build_sqrt_mult( Mx, rank_Mx, var_threshold, max_k, tol=tol)
  My_sqrt <- build_sqrt_mult( My, rank_My, var_threshold, max_k, tol=tol)
  Ax_sqrt <- build_sqrt_mult( Ax, rank_Ax, var_threshold, max_k, tol=tol)
  Ay_sqrt <- build_sqrt_mult( Ay, rank_Ay, var_threshold, max_k, tol=tol)

  Ax_isqrt<- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
  Ay_isqrt<- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)

  partial_svd <- function(M, k) {
    k_eff <- min(k, nrow(M), ncol(M))
    sv <- base::svd(M, nu = k_eff, nv = k_eff)
    list(u = sv$u[, seq_len(k_eff), drop = FALSE],
         d = sv$d[seq_len(k_eff)],
         v = sv$v[, seq_len(k_eff), drop = FALSE])
  }

  ## ---- 4. decide row/col likeness (for dual path) ------------------
  want_row_like <- if(!dual) NA else
    if(!is.null(force_row_likeness)) TRUE else
    if(!is.null(force_col_likeness)) FALSE else n < (p+q)/2

  ## ---- 5. DIRECT path (dual = FALSE) -------------------------------
  if(!dual){
      Xtil <- col_transform(row_transform(Xp,Mx_sqrt), Ax_sqrt)   # n × p
      Ytil <- col_transform(row_transform(Yp,My_sqrt), Ay_sqrt)   # n × q

      C    <- crossprod(Xtil, Ytil)                               # p × q
      sv   <- partial_svd(C, ncomp)

      Tt   <- Xtil %*% sv$u
      Ut   <- Ytil %*% sv$v

      vx   <- Ax_isqrt(sv$u)      # original-space loadings
      vy   <- Ay_isqrt(sv$v)

      return(cross_projector(
               vx=vx, vy=vy, Tx=Tt, Ty=Ut,
               preproc_x=px, preproc_y=py,
               d_full=sv$d, ncomp=ncol(vx),
               classes=c("genplscorr"),
               dual=FALSE, rank_Mx=rank_Mx, rank_Ax=rank_Ax,
               rank_My=rank_My, rank_Ay=rank_Ay,
               var_threshold=var_threshold, max_k=max_k,
               tol=tol, verbose=verbose))
  }

  ## ---- 6. OPERATOR path (dual = TRUE) ------------------------------
  if(want_row_like){                       # ------ 6A. row likeness -----
      Xr  <- row_transform(Xp, Mx_sqrt)                   # n × p
      Yr  <- row_transform(Yp, My_sqrt)                   # n × q

      Gx  <- tcrossprod( Ax_sqrt(t(Xr)) )                 # n × n
      Gy  <- tcrossprod( Ay_sqrt(t(Yr)) )                 # n × n

      Qx  <- partial_eig_approx(Gx, rank_Mx,
                                var_threshold,max_k,tol=tol)$Q
      Qy  <- partial_eig_approx(Gy, rank_My,
                                var_threshold,max_k,tol=tol)$Q

      XYR <- tcrossprod(Xr, Yr)                           # n × n
      sv  <- partial_svd( crossprod(Qx, XYR %*% Qy), ncomp)

      Tt  <- Qx %*% sv$u
      Ut  <- Qy %*% sv$v

      Xstar<- col_transform(Xr, Ax_sqrt)                  # n × p
      Ystar<- col_transform(Yr, Ay_sqrt)                  # n × q

      vx  <- Ax_isqrt( crossprod(Xstar,  Tt) )
      vy  <- Ay_isqrt( crossprod(Ystar,  Ut) )

  }else{                                   # ------ 6B. column likeness ---
      Xc  <- col_transform(Xp, Ax_sqrt)                  # n × p
      Yc  <- col_transform(Yp, Ay_sqrt)                  # n × q

      Hx  <- crossprod( Mx_sqrt(Xc) )                    # p × p
      Hy  <- crossprod( My_sqrt(Yc) )                    # q × q

      Qx  <- partial_eig_approx(Hx, rank_Ax,
                                var_threshold,max_k,tol=tol)$Q
      Qy  <- partial_eig_approx(Hy, rank_Ay,
                                var_threshold,max_k,tol=tol)$Q

      XYc <- crossprod( row_transform(Xc,Mx_sqrt),
                        row_transform(Yc,My_sqrt) )      # p × q

      sv  <- partial_svd( crossprod(Qx, XYc %*% Qy), ncomp)

      vx_emb <- Qx %*% sv$u
      vy_emb <- Qy %*% sv$v

      Tt  <- row_transform(col_transform(Xp,Ax_sqrt),Mx_sqrt) %*% vx_emb
      Ut  <- row_transform(col_transform(Yp,Ay_sqrt),My_sqrt) %*% vy_emb

      vx <- Ax_isqrt(vx_emb)
      vy <- Ay_isqrt(vy_emb)
  }

  ## ---- 7. return ---------------------------------------------------
  cross_projector(
        vx=vx, vy=vy, Tx=Tt, Ty=Ut,
        preproc_x=px, preproc_y=py,
        d_full=sv$d, ncomp=ncol(vx),
        classes=c("genplscorr"),
        dual=TRUE,  want_row_like=want_row_like,
        rank_Mx=rank_Mx, rank_Ax=rank_Ax,
        rank_My=rank_My, rank_Ay=rank_Ay,
        var_threshold=var_threshold, max_k=max_k,
        tol=tol, verbose=verbose)
}

