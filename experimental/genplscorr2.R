genplscorr2 <- function(X, Y,
                        Mx=NULL, Ax=NULL,
                        My=NULL, Ay=NULL,
                        ncomp=2,
                        preproc_x=pass(),
                        preproc_y=pass(),
                        rank_Mx=NULL, rank_Ax=NULL,
                        rank_My=NULL, rank_Ay=NULL,
                        var_threshold=0.99,
                        max_k=200,
                        dual=FALSE,
                        # optional param to refine row_likeness heuristic:
                        force_row_likeness = NULL,
                        force_col_likeness = NULL,
                        # choose partial or full svd method
                        svd_method = c("base", "RSpectra", "irlba"),
                        tol=1e-9,
                        maxiter=200,  # not used
                        verbose=FALSE)
{
  #### 1) Basic checks
  if(!is.matrix(X) || !is.numeric(X)) {
    stop("genplscorr2: X must be a numeric matrix.")
  }
  if(!is.matrix(Y) || !is.numeric(Y)) {
    stop("genplscorr2: Y must be a numeric matrix.")
  }
  if(nrow(X) != nrow(Y)) {
    stop("genplscorr2: X,Y must have the same number of rows.")
  }
  if(ncomp <= 0) {
    ncomp <- min(ncol(X), ncol(Y))
    if(verbose) message("genplscorr2: ncomp <= 0 => using all up to min(p,q).")
  }
  
  # handle user-chosen SVD method
  svd_method <- match.arg(svd_method)
  
  library(multivarious)  # for cross_projector, etc. (ideally at top of file)
  
  #### 2) Preprocessing
  px_obj <- prep(preproc_x)
  py_obj <- prep(preproc_y)
  Xp <- init_transform(px_obj, X)
  Yp <- init_transform(py_obj, Y)
  n <- nrow(Xp)
  p <- ncol(Xp)
  q <- ncol(Yp)
  
  #### 3) Identity defaults for constraints
  if(is.null(Mx)) Mx <- Matrix::Diagonal(n)
  if(is.null(Ax)) Ax <- Matrix::Diagonal(p)
  if(is.null(My)) My <- Matrix::Diagonal(n)
  if(is.null(Ay)) Ay <- Matrix::Diagonal(q)
  
  #### Helper for partial SVD on a (m x n) matrix 'mat' => return (u, d, v)
  partial_svd <- function(mat, k, method = "base") {
    # k: how many vectors we want
    # mat can be crossOp or crossprod(...) etc.
    # returns list(u=..., d=..., v=...) each truncated to k
    if(k <= 0) k <- min(nrow(mat), ncol(mat))  # safety
    if(method == "base") {
      # do full base::svd and truncate
      sfull <- base::svd(mat, nu=min(nrow(mat), k), nv=min(ncol(mat), k))
      sfull$u <- sfull$u[, seq_len(min(k, ncol(sfull$u))), drop=FALSE]
      sfull$v <- sfull$v[, seq_len(min(k, ncol(sfull$v))), drop=FALSE]
      sfull$d <- sfull$d[seq_len(min(k, length(sfull$d)))]
      return(sfull)
    } else if(method == "RSpectra") {
      if(!requireNamespace("RSpectra", quietly=TRUE)) {
        warning("RSpectra not installed; falling back to base::svd")
        return(partial_svd(mat, k, method="base"))
      }
      # We want top k largest singular values => compute partial SVD
      # RSpectra::svds => we have to choose 'nu' and 'nv'
      sv <- RSpectra::svds(mat, k=k)
      # note that RSpectra::svds returns a list(u, v, d), each dimension truncated to k
      return(list(u=sv$u, d=sv$d, v=sv$v))
    } else if(method == "irlba") {
      if(!requireNamespace("irlba", quietly=TRUE)) {
        warning("irlba not installed; falling back to base::svd")
        return(partial_svd(mat, k, method="base"))
      }
      # use irlba::irlba
      sv <- irlba::irlba(mat, nv=k, nu=k)
      return(list(u=sv$u, d=sv$d, v=sv$v))
    } else {
      stop("Unrecognized svd_method: ", method)
    }
  }
  
  #### Decide row-likeness or column-likeness if dual=TRUE
  # user can override with force_row_likeness or force_col_likeness
  decide_row_col <- function(n, p, q) {
    if(!is.null(force_row_likeness) && force_row_likeness) {
      return(TRUE)
    }
    if(!is.null(force_col_likeness) && force_col_likeness) {
      return(FALSE)
    }
    # fallback to the old heuristic
    (n < (p + q)/2)
  }
  
  #### 4) if not dual => direct expansions
  if(!dual) {
    if(verbose) message("[genplscorr2] dual=FALSE => direct expansions + crossprod + partial SVD approach.")
    Xrw  <- row_transform(Xp, build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol))
    Xstar<- col_transform(Xrw, build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol))
    Yrw  <- row_transform(Yp, build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol))
    Ystar<- col_transform(Yrw, build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol))
    
    if(verbose) {
      message(sprintf(" shape(Xstar)=(%d,%d), shape(Ystar)=(%d,%d). forming crossprod => p x q",
                      nrow(Xstar), ncol(Xstar), nrow(Ystar), ncol(Ystar)))
    }
    # crossprod => p x q
    C <- crossprod(Xstar, Ystar)
    
    # partial SVD => keep top ncomp
    sres <- partial_svd(C, k=ncomp, method=svd_method)
    d_full <- sres$d
    r0     <- length(d_full)  # actual # singular values returned
    u_trunc<- sres$u  # p x r0
    v_trunc<- sres$v  # q x r0
    
    # factor scores in embedded domain
    Tstar <- Xstar %*% u_trunc  # n x r0
    Ustar <- Ystar %*% v_trunc  # n x r0
    
    # original loadings => multiply by Ax^-1/2, Ay^-1/2
    Ax_inv <- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
    Ay_inv <- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
    
    # instead of apply(..., 2, Ax_inv), we do a matrix approach:
    vx_orig <- Ax_inv(u_trunc)  # p x r0
    vy_orig <- Ay_inv(v_trunc)  # q x r0
    
    out <- cross_projector(
      vx = vx_orig,  vy = vy_orig,
      preproc_x = px_obj, preproc_y = py_obj,
      classes=c("genplscorr"),
      tilde_vx = u_trunc, tilde_vy = v_trunc,
      Tx = Tstar,    Ty = Ustar,
      ncomp = r0,
      d_full=d_full,
      rank_Mx=rank_Mx, rank_Ax=rank_Ax,
      rank_My=rank_My, rank_Ay=rank_Ay,
      var_threshold=var_threshold, max_k=max_k,
      dual=FALSE, tol=tol, verbose=verbose
    )
    return(out)
  }
  
  #### 5) Otherwise, if dual=TRUE => row-likeness or column-likeness approach
  if(verbose) message("[genplscorr2] dual=TRUE => operator-based expansions (like gmdLA).")
  
  row_likeness <- decide_row_col(n, p, q)
  if(row_likeness) {
    if(verbose) message(" row-likeness: building operators Gx, Gy in n x n space.")
    ## step 1: partial expansions for row-likeness
    #   Gx = sqrt(Mx) X Ax X^T sqrt(Mx) => n x n
    Xrow    <- row_transform(Xp, build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol))
    Ax_sqrt <- build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
    midX    <- Ax_sqrt(t(Xrow)) # shape p x n
    Gx      <- crossprod(midX)  # n x n
    outGx   <- partial_eig_approx(Gx, user_rank=rank_Mx,  # respect rank_Mx
                                  var_threshold, max_k, tol=tol)
    Qx      <- outGx$Q
    lamx    <- outGx$lams
    
    #   Gy = sqrt(My) Y Ay Y^T sqrt(My) => n x n
    Yrow    <- row_transform(Yp, build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol))
    Ay_sqrt <- build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
    midY    <- Ay_sqrt(t(Yrow)) # shape q x n
    Gy      <- crossprod(midY)  # n x n
    outGy   <- partial_eig_approx(Gy, user_rank=rank_My,
                                  var_threshold, max_k, tol=tol)
    Qy      <- outGy$Q
    lamy    <- outGy$lams
    
    # step 2: define cross-operator => XY_row = Xrow Yrow^T => n x n
    XY_row  <- crossprod(Xrow, Yrow)  # n x n
    crossOp <- crossprod(Qx, XY_row %*% Qy)
    
    # step 3: partial SVD => correlation directions alpha_x, alpha_y
    sres    <- partial_svd(crossOp, k=ncomp, method=svd_method)
    d_full  <- sres$d
    r0      <- length(d_full)
    alpha_x <- sres$u  # rx x r0
    alpha_y <- sres$v  # ry x r0
    
    # step 4: factor scores => T^* = Qx alpha_x, U^* = Qy alpha_y
    Tstar <- Qx %*% alpha_x
    Ustar <- Qy %*% alpha_y
    
    # step 5: original loadings
    #   Xstar2 = (n x p) after row_transform(Xp, Mx^(1/2)) then col_transform(..., Ax^(1/2))
    Xstar2 <- col_transform(Xrow, Ax_sqrt)  # shape n x p
    p_k    <- crossprod(Xstar2, Tstar)      # => p x r0
    
    Ax_inv <- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
    # do a matrix approach if p_k is p x r0
    vx_orig <- Ax_inv(p_k)
    
    Ystar2 <- col_transform(Yrow, Ay_sqrt)  # shape n x q
    q_k    <- crossprod(Ystar2, Ustar)      # => q x r0
    Ay_inv <- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
    vy_orig<- Ay_inv(q_k)
    
  } else {
    if(verbose) message(" column-likeness: building Hx, Hy in p x p and q x q space.")
    ## step 1: partial expansions for column-likeness
    Xcol <- col_transform(Xp, build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol))
    Mx_sqrt <- build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol)
    midX  <- Mx_sqrt(Xcol)
    Hx    <- crossprod(midX)  # p x p
    outHx <- partial_eig_approx(Hx, user_rank=rank_Ax,
                                var_threshold, max_k, tol=tol)
    Qx    <- outHx$Q
    lamx  <- outHx$lams
    
    Ycol   <- col_transform(Yp, build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol))
    My_sqrt<- build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol)
    midY   <- My_sqrt(Ycol)
    Hy     <- crossprod(midY)  # q x q
    outHy  <- partial_eig_approx(Hy, user_rank=rank_Ay,
                                 var_threshold, max_k, tol=tol)
    Qy     <- outHy$Q
    lamy   <- outHy$lams
    
    # step 2: cross-operator in col-likeness
    Xcol2  <- row_transform(Xcol, Mx_sqrt)
    Ycol2  <- row_transform(Ycol, My_sqrt)
    XY_col <- crossprod(Xcol2, Ycol2)  # p x q
    crossOp<- crossprod(Qx, XY_col %*% Qy)
    
    # step 3: partial SVD => correlation directions
    sres   <- partial_svd(crossOp, k=ncomp, method=svd_method)
    d_full <- sres$d
    r0     <- length(d_full)
    alpha_x<- sres$u
    alpha_y<- sres$v
    
    # step 4: factor scores
    vx_emb <- Qx %*% alpha_x  # p x r0
    vy_emb <- Qy %*% alpha_y  # q x r0
    
    # define Tstar, Ustar similarly
    #   Tstar= Xcol4 %*% vx_emb
    Xcol3  <- col_transform(Xp, build_sqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol))
    Xcol4  <- row_transform(Xcol3, build_sqrt_mult(Mx, rank_Mx, var_threshold, max_k, tol=tol))
    Tstar  <- Xcol4 %*% vx_emb
    
    Ycol3  <- col_transform(Yp, build_sqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol))
    Ycol4  <- row_transform(Ycol3, build_sqrt_mult(My, rank_My, var_threshold, max_k, tol=tol))
    Ustar  <- Ycol4 %*% vy_emb
    
    # step 5: map loadings back
    Ax_inv <- build_invsqrt_mult(Ax, rank_Ax, var_threshold, max_k, tol=tol)
    Ay_inv <- build_invsqrt_mult(Ay, rank_Ay, var_threshold, max_k, tol=tol)
    vx_orig<- Ax_inv(vx_emb)
    vy_orig<- Ay_inv(vy_emb)
  }
  
  #### 6) final cross_projector
  out <- cross_projector(
    vx = vx_orig,
    vy = vy_orig,
    preproc_x = px_obj,
    preproc_y = py_obj,
    classes   = c("genplscorr"),
    Tx = Tstar,
    Ty = Ustar,
    ncomp = r0,
    d_full= d_full,
    rank_Mx=rank_Mx,
    rank_Ax=rank_Ax,
    rank_My=rank_My,
    rank_Ay=rank_Ay,
    var_threshold=var_threshold,
    max_k=max_k,
    dual=dual,
    tol=tol,
    verbose=verbose
  )
  out
}