#' @export
#'
#' @title Generalized partial least squares "regression decomposition" (GPLSREG)
#'
#' @description Computes generalized partial least squares "regression decomposition" between two data matrices.
#' GPLSREG allows for the use of left (row) and right (column) weights for each data matrix.
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param XLW An \emph{I} by \emph{I} matrix of row weights for \code{X}. Default is \code{diag(nrow(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YLW An \emph{I} by \emph{I} matrix of row weights for \code{Y}. Default is \code{diag(nrow(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param XRW A \emph{J} by \emph{J} matrix of row weights for \code{X}. Default is \code{diag(ncol(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YRW A \emph{K} by \emph{K} matrix of row weights for \code{Y}. Default is \code{diag(ncol(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#'
#'
#' @return A list of outputs
#' \item{d}{A vector containing the singular values from each iteration.}
#' \item{u}{Left (rows) singular vectors.}
#' \item{v}{Right (columns) singular vectors. In GPLSREG sometimes called "weight matrix".}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{fj}{Right (columns) component scores.}
#' \item{tx}{X "Latent vectors": A normed version of \code{lx} for use in rebuilding \code{X} data}
#' \item{u_hat}{X "Loading matrix": A "predicted" version of \code{u} for use in rebuilding \code{X} data}
#' \item{betas}{"Regression weights": Akin to betas for use in rebuilding \code{Y}}
#' \item{X_reconstructeds}{A version of \code{X} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{Y_reconstructeds}{A version of \code{Y} reconstructed for each iteration (i.e., latent variable/component)}
#' \item{X_residuals}{The residualized (i.e., \code{X - X_reconstructeds}) version of \code{X} for each iteration (i.e., latent variable/component)}
#' \item{Y_residuals}{The residualized (i.e., \code{Y - Y_reconstructeds}) version of \code{Y} for each iteration (i.e., latent variable/component)}
#' \item{r2_x}{Proportion of explained variance from \code{X} to each latent variable/component.}
#' \item{r2_y}{Proportion of explained variance from \code{Y} to each latent variable/component.}
#' \item{Y_reconstructed}{A version of \code{Y} reconstructed from all iterations (i.e., latent variables/components); see \code{components}.}
#' \item{Y_residual}{The residualized (i.e., \code{Y - Y_reconstructed} from all iterations (i.e., latent variables/components); see \code{components}.}
#'
#' @seealso \code{\link{gpls_can}} \code{\link{gpls_cor}} \code{\link{pls_reg}} \code{\link{plsca_reg}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#'
#' @examples
#'
#'  \dontrun{
#'  library(GSVD)
#'  data("wine", package = "GSVD")
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'
#'  ## standard partial least squares "regression decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  gplsreg_pls_optimization <- gpls_reg(X, Y)
#'
#'  ## partial least squares "regression decomposition"
#'  ### but with the optimization per latent variable of CCA
#'  #### because of optimization, this ends up identical to
#'  #### cca(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplsreg_cca_optimization <- gpls_reg( t(MASS::ginv(X)), t(MASS::ginv(Y)),
#'       XRW = crossprod(X), YRW = crossprod(Y))
#'
#'  ## partial least squares "regression decomposition"
#'  ### but with the optimization per latent variable of RRR/RDA
#'  #### because of optimization, this ends up identical to
#'  #### rrr(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  #### or rda(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplsreg_rrr_optimization <- gpls_reg( t(MASS::ginv(X)), Y, XRW = crossprod(X))
#'
#'  rm(X)
#'  rm(Y)
#'
#'  ## partial least squares-correspondence analysis "regression decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  X_ca_preproc <- ca_preproc(X)
#'  Y_ca_preproc <- ca_preproc(Y)
#'
#'  gplsreg_plsca <- gpls_reg( X = X_ca_preproc$Z, Y = Y_ca_preproc$Z,
#'      XLW = diag(1/X_ca_preproc$m), YLW = diag(1/Y_ca_preproc$m),
#'      XRW = diag(1/X_ca_preproc$w), YRW = diag(1/Y_ca_preproc$w)
#'  )
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares
#' @importFrom MASS ginv
#' @importFrom GSVD invsqrt_psd_matrix sqrt_psd_matrix

gpls_reg <- function(X, Y,
                     XLW = diag(nrow(X)), YLW = diag(nrow(Y)),
                     XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                     components = 0, tol = .Machine$double.eps){


  X_gsvd <- gsvd(X, XLW, XRW)
    X_rank <- length(X_gsvd$d)
    X_trace <- sum(X_gsvd$d^2)
  rm(X_gsvd)

  Y_gsvd <- gsvd(Y, YLW, YRW)
    Y_rank <- length(Y_gsvd$d)
    Y_trace <- sum(Y_gsvd$d^2)
  rm(Y_gsvd)

  stopped_early <- F

  if((components > X_rank) | (components < 1)){
    components <- X_rank
  }

  lx <- tx <- matrix(NA,nrow(X),components)
  ly <- matrix(NA,nrow(Y),components)

  u_hat <- fi <- p <- u <- matrix(NA,ncol(X),components)
  fj <- q <- v <- matrix(NA,ncol(Y),components)

  r2_x_cumulative <- r2_y_cumulative <- d <- betas <- rep(NA, components)

  X_reconstructeds <- X_residuals <- array(NA,dim=c(nrow(X),ncol(X),components))
  Y_reconstructeds <- Y_residuals <- array(NA,dim=c(nrow(Y),ncol(Y),components))

  X_deflate <- X
  Y_deflate <- Y

  for(i in 1:components){

    gplssvd_results <- gpls_cor(X_deflate, Y_deflate,
                                XRW = XRW, YRW = YRW,
                                XLW = XLW, YLW = YLW,
                                components = 1, tol = tol)

    u[,i] <- gplssvd_results$u
    p[,i] <- gplssvd_results$p
    fi[,i] <- gplssvd_results$fi
    v[,i] <- gplssvd_results$v
    q[,i] <- gplssvd_results$q
    fj[,i] <- gplssvd_results$fj
    d[i] <- gplssvd_results$d
    lx[,i] <- gplssvd_results$lx
    ly[,i] <- gplssvd_results$ly

    tx[,i] <- lx[,i] / sqrt(sum(lx[,i]^2))
    betas[i] <- t(ly[,i]) %*% tx[,i]
    u_hat[,i] <- t(tx[,i]) %*% ( GSVD::sqrt_psd_matrix(XLW) %*% X_deflate %*% GSVD::sqrt_psd_matrix(XRW) )


    X_reconstructeds[,,i] <- GSVD::invsqrt_psd_matrix(XLW) %*% (tx[,i] %o% u_hat[,i]) %*% GSVD::invsqrt_psd_matrix(XRW)
      X_reconstructeds[abs(X_reconstructeds) < tol] <- 0
    Y_reconstructeds[,,i] <- GSVD::invsqrt_psd_matrix(YLW) %*% ((tx[,i] * betas[i]) %o% v[,i]) %*% GSVD::invsqrt_psd_matrix(YRW)
      Y_reconstructeds[abs(Y_reconstructeds) < tol] <- 0


    X_residuals[,,i] <- (X_deflate - X_reconstructeds[,,i])
      X_residuals[abs(X_residuals) < tol] <- 0
    Y_residuals[,,i] <- (Y_deflate - Y_reconstructeds[,,i])
      Y_residuals[abs(Y_residuals) < tol] <- 0

    X_deflate <- X_residuals[,,i]
    Y_deflate <- Y_residuals[,,i]


    r2_x_cumulative[i] <- (X_trace-sum( ( GSVD::sqrt_psd_matrix(XLW) %*%  X_deflate %*% GSVD::sqrt_psd_matrix(XRW) ) ^2)) / X_trace
    r2_y_cumulative[i] <- (Y_trace-sum( ( GSVD::sqrt_psd_matrix(YLW) %*%  Y_deflate %*% GSVD::sqrt_psd_matrix(YRW) ) ^2)) / Y_trace


    if( (sum(Y_deflate^2) < tol) & (i < components) ){

      stopped_early <- T
      warning("gpls_reg: Y is fully deflated. Stopping early.")

    }
    if( (sum(X_deflate^2) < tol) & (i < components) ){

      stopped_early <- T
      warning("gpls_reg: X is fully deflated. Stopping early.")

    }

    if(stopped_early){
      break
    }

  }

  if(stopped_early){
    u <- u[,1:i]
    p <- p[,1:i]
    fi <- fi[,1:i]
    v <- v[,1:i]
    q <- q[,1:i]
    fj <- fj[,1:i]
    d <- d[1:i]
    lx <- lx[,1:i]
    ly <- ly[,1:i]

    tx <- tx[,1:i]
    betas <- betas[1:i]
    u_hat <- u_hat[,1:i]

    X_reconstructeds <- X_reconstructeds[,,1:i]
    Y_reconstructeds <- Y_reconstructeds[,,1:i]
    X_residuals <- X_residuals[,,1:i]
    Y_residuals <- Y_residuals[,,1:i]

    r2_x_cumulative <- r2_x_cumulative[1:i]
    r2_y_cumulative <- r2_y_cumulative[1:i]
  }


  Y_reconstructed <- GSVD::invsqrt_psd_matrix(YLW) %*% (tx %*% diag(betas,length(betas),length(betas)) %*% t(v)) %*% GSVD::invsqrt_psd_matrix(YRW)
    Y_reconstructed[abs(Y_reconstructed) < tol] <- 0
  Y_residual <- Y - Y_reconstructed



  rownames(tx) <- rownames(lx) <- rownames(X_reconstructeds) <- rownames(X_residuals) <- rownames(X)
  rownames(fi) <- rownames(u) <- rownames(p) <- rownames(u_hat) <- colnames(X_reconstructeds) <- colnames(X_residuals) <- colnames(X)
  rownames(Y_reconstructed) <- rownames(Y_residual) <- rownames(Y_reconstructeds) <- rownames(Y_residuals) <- rownames(Y)
  rownames(fj) <- rownames(v) <- rownames(q) <- colnames(Y_reconstructed) <- colnames(Y_residual) <- colnames(Y_reconstructeds) <- colnames(Y_residuals) <- colnames(Y)


  return( list(
    d = d, u = u, v = v, lx = lx, ly = ly,
    p = p, q = q, fi = fi, fj = fj,
    tx = tx, u_hat = u_hat, betas = betas,
    X_reconstructeds = X_reconstructeds, X_residuals = X_residuals,
    Y_reconstructeds = Y_reconstructeds, Y_residuals = Y_residuals,
    r2_x = diff(c(0,r2_x_cumulative)), r2_y = diff(c(0,r2_y_cumulative)),
    Y_reconstructed = Y_reconstructed, Y_residual = Y_residual
  ) )


}

#' @export
#'
#' @title Generalized partial least squares "correlation decomposition" (GPLSCOR)
#'
#' @description Computes generalized partial least squares "correlation decomposition" between two data matrices.
#' GPLSCOR allows for the use of left (row) and right (column) weights for each data matrix.
#'
#' @param X Data matrix with \emph{I} rows and \emph{J} columns
#' @param Y Data matrix with \emph{I} rows and \emph{K} columns
#' @param XLW An \emph{I} by \emph{I} matrix of row weights for \code{X}. Default is \code{diag(nrow(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YLW An \emph{I} by \emph{I} matrix of row weights for \code{Y}. Default is \code{diag(nrow(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param XRW A \emph{J} by \emph{J} matrix of row weights for \code{X}. Default is \code{diag(ncol(X))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param YRW A \emph{K} by \emph{K} matrix of row weights for \code{Y}. Default is \code{diag(ncol(Y))} (i.e., all ones on the diagonal; zeros off-diagonal).
#' @param components The number of components to return. If < 1 then the maximum components will be returned. Default = 0.
#' @param tol default is .Machine$double.eps. A parameter to pass through to \code{\link[GSVD]{gplssvd}}; eliminates singular values that are effectively zero (and thus drops null components).
#'
#'
#' @return a list of outputs; see also \code{\link[GSVD]{gplssvd}}
#' \item{d.orig}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l.orig}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{tau}{A vector that contains the (original) explained variance per component (via eigenvalues: \code{$l.orig}.}
#' \item{d}{A vector of length \code{min(length(d.orig), k)} containing the retained singular values}
#' \item{l}{A vector of length \code{min(length(l.orig), k)} containing the retained eigen values}
#' \item{u}{Left (rows) singular vectors.}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{v}{Right (columns) singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fj}{Right (columns) component scores.}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#'
#' @seealso \code{\link{gpls_reg}} \code{\link{gpls_can}} \code{\link{pls_cor}} \code{\link{plsca_cor}} \code{\link{cca}} \code{\link{rrr}} \code{\link{rda}} \code{\link[GSVD]{gplssvd}}
#'
#' @references
#' Beaton, D., ADNI, Saporta, G., Abdi, H. (2019). A generalization of partial least squares regression and correspondence analysis for categorical and mixed data: An application with the ADNI data. \emph{bioRxiv}, 598888.
#' Beaton, D., Dunlop, J., & Abdi, H. (2016). Partial least squares correspondence analysis: A framework to simultaneously analyze behavioral and genetic data. \emph{Psychological methods}, \bold{21} (4), 621.
#'
#' @examples
#'
#'  \dontrun{
#'  library(GSVD)
#'  data("wine", package = "GSVD")
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'  ## partial least squares "correlation decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  #### this is pls_cor(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_pls_optimization <- gpls_cor(X, Y)
#'
#'  ## partial least squares "correlation decomposition"
#'  ### but with the optimization per latent variable of CCA
#'  #### this is cca(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_cca_optimization <- gpls_cor( t(MASS::ginv(X)), t(MASS::ginv(Y)),
#'       XRW = crossprod(X), YRW = crossprod(Y))
#'
#'  ## partial least squares "correlation decomposition"
#'  ### but with the optimization per latent variable of RRR/RDA
#'  #### this is rrr(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F) or
#'  #### rda(X, Y, center_X = F, center_Y = F, scale_X = F, scale_Y = F)
#'  gplscor_rrr_optimization <- gpls_cor( t(MASS::ginv(X)), Y, XRW = crossprod(X))
#'
#'
#'  rm(X)
#'  rm(Y)
#'
#'  ## partial least squares-correspondence analysis "correlation decomposition"
#'  #### the first latent variable from reg & cor & can are identical in all PLSs.
#'  data("snps.druguse", package = "GSVD")
#'  X <- make_data_disjunctive(snps.druguse$DATA1)
#'  Y <- make_data_disjunctive(snps.druguse$DATA2)
#'
#'  X_ca_preproc <- ca_preproc(X)
#'  Y_ca_preproc <- ca_preproc(Y)
#'
#'  #### this is plsca_cor(X, Y)
#'  gplscor_plsca <- gpls_cor( X = X_ca_preproc$Z, Y = Y_ca_preproc$Z,
#'      XLW = diag(1/X_ca_preproc$m), YLW = diag(1/Y_ca_preproc$m),
#'      XRW = diag(1/X_ca_preproc$w), YRW = diag(1/Y_ca_preproc$w)
#'  )
#'
#'  }
#'
#' @keywords multivariate, diagonalization, partial least squares
#' @importFrom MASS ginv

gpls_cor <- function(X, Y,
                  XLW = diag(nrow(X)), YLW = diag(nrow(Y)),
                  XRW = diag(ncol(X)), YRW = diag(ncol(Y)),
                  components = 0, tol = .Machine$double.eps){

    ## that's it.
    gplssvd(X, Y, XLW = XLW, YLW = YLW, XRW = XRW, YRW = YRW, k = components, tol = tol)

}

' @export
#'
#' @title Generalized partial least squares singular value decomposition
#'
#' @description
#' \code{gplssvd} takes in left (\code{XLW}, \code{YLW}) and right (\code{XRW}, \code{YRW}) constraints (usually diagonal matrices, but any positive semi-definite matrix is fine) that are applied to each data matrix (\code{X} and \code{Y}), respectively.
#'   Right constraints for each data matrix are used for the orthogonality conditions.
#'
#' @param X a data matrix
#' @param Y a data matrix
#' @param XLW \bold{X}'s \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the \code{X} matrix and thus left singular vectors.
#' @param YLW \bold{Y}'s \bold{L}eft \bold{W}eights -- the constraints applied to the left side (rows) of the \code{Y} matrix and thus left singular vectors.
#' @param XRW \bold{X}'s \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the \code{X} matrix and thus right singular vectors.
#' @param YRW \bold{Y}'s \bold{R}ight \bold{W}eights -- the constraints applied to the right side (columns) of the \code{Y} matrix and thus right singular vectors.
#' @param k total number of components to return though the full variance will still be returned (see \code{d_full}). If 0, the full set of components are returned.
#' @param tol default is .Machine$double.eps. A tolerance level for eliminating effectively zero (small variance), negative, imaginary eigen/singular value components.
#'
#' @return A list with thirteen elements:
#' \item{d_full}{A vector containing the singular values of DAT above the tolerance threshold (based on eigenvalues).}
#' \item{l_full}{A vector containing the eigen values of DAT above the tolerance threshold (\code{tol}).}
#' \item{d}{A vector of length \code{min(length(d_full), k)} containing the retained singular values}
#' \item{l}{A vector of length \code{min(length(l_full), k)} containing the retained eigen values}
#' \item{u}{Left (rows) singular vectors. Dimensions are \code{nrow(DAT)} by k.}
#' \item{p}{Left (rows) generalized singular vectors.}
#' \item{fi}{Left (rows) component scores.}
#' \item{lx}{Latent variable scores for rows of \code{X}}
#' \item{v}{Right (columns) singular vectors.}
#' \item{q}{Right (columns) generalized singular vectors.}
#' \item{fj}{Right (columns) component scores.}
#' \item{ly}{Latent variable scores for rows of \code{Y}}
#'
#' @seealso \code{\link{tolerance_svd}}, \code{\link{geigen}} and \code{\link{gsvd}}
#'
#' @examples
#'
#'  # Three "two-table" technique examples
#'  data(wine)
#'  X <- scale(wine$objective)
#'  Y <- scale(wine$subjective)
#'
#'  ## Partial least squares (correlation)
#'  pls.res <- gplssvd(X, Y)
#'
#'
#'  ## Canonical correlation analysis (CCA)
#'  ### NOTE:
#'  #### This is not "traditional" CCA because of the generalized inverse.
#'  #### However the results are the same as standard CCA when data are not rank deficient.
#'  #### and this particular version uses tricks to minimize memory & computation
#'  cca.res <- gplssvd(
#'      X = MASS::ginv(t(X)),
#'      Y = MASS::ginv(t(Y)),
#'      XRW=crossprod(X),
#'      YRW=crossprod(Y)
#'  )
#'
#'
#'  ## Reduced rank regression (RRR) a.k.a. redundancy analysis (RDA)
#'  ### NOTE:
#'  #### This is not "traditional" RRR because of the generalized inverse.
#'  #### However the results are the same as standard RRR when data are not rank deficient.
#'  rrr.res <- gplssvd(
#'      X = MASS::ginv(t(X)),
#'      Y = Y,
#'      XRW=crossprod(X)
#'  )
#'
#' @author Derek Beaton
#' @keywords multivariate

gplssvd <- function(X, Y, XLW, YLW, XRW, YRW, k = 0, tol = .Machine$double.eps){

  # preliminaries
  X_dimensions <- dim(X)
  ## stolen from MASS::ginv()
  if (length(X_dimensions) > 2 || !(is.numeric(X) || is.complex(X))){
    stop("gplssvd: 'X' must be a numeric or complex matrix")
  }
  if (!is.matrix(X)){
    X <- as.matrix(X)
  }

  # preliminaries
  Y_dimensions <- dim(Y)
  ## stolen from MASS::ginv()
  if (length(Y_dimensions) > 2 || !(is.numeric(Y) || is.complex(Y))){
    stop("gplssvd: 'Y' must be a numeric or complex matrix")
  }
  if (!is.matrix(Y)){
    Y <- as.matrix(Y)
  }

  # check that row dimensions match
  if(X_dimensions[1] != Y_dimensions[1]){
    stop("gplssvd: X and Y must have the same number of rows")
  }

  # a few things about XLW for stopping conditions
  XLW_is_missing <- missing(XLW)
  if(!XLW_is_missing){

    XLW_is_vector <- is.vector(XLW)

    if(!XLW_is_vector){

      if( nrow(XLW) != ncol(XLW) | nrow(XLW) != X_dimensions[1] ){
        stop("gplssvd: nrow(XLW) does not equal ncol(XLW) or nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(XLW)){
        stop("gplssvd: XLW is empty (i.e., all 0s")
      }
    }

    if(XLW_is_vector){
      if(length(XLW)!=X_dimensions[1]){
        stop("gplssvd: length(XLW) does not equal nrow(X)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(XLW)<=tol)){
      if(!are_all_values_positive(XLW)){
        stop("gplssvd: XLW is not strictly positive values")
      }
    }
  }

  # a few things about XRW for stopping conditions
  XRW_is_missing <- missing(XRW)
  if(!XRW_is_missing){

    XRW_is_vector <- is.vector(XRW)

    if(!XRW_is_vector){

      if( nrow(XRW) != ncol(XRW) | nrow(XRW) != X_dimensions[2] ){
        stop("gplssvd: nrow(XRW) does not equal ncol(XRW) or ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(XRW)){
        stop("gplssvd: XRW is empty (i.e., all 0s")
      }
    }

    if(XRW_is_vector){
      if(length(XRW)!=X_dimensions[2]){
        stop("gplssvd: length(XRW) does not equal ncol(X)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(XRW)<=tol)){
      if(!are_all_values_positive(XRW)){
        stop("gplssvd: XRW is not strictly positive values")
      }
    }
  }

  # a few things about YLW for stopping conditions
  YLW_is_missing <- missing(YLW)
  if(!YLW_is_missing){

    YLW_is_vector <- is.vector(YLW)

    if(!YLW_is_vector){

      if( nrow(YLW) != ncol(YLW) | nrow(YLW) != Y_dimensions[1] ){
        stop("gplssvd: nrow(YLW) does not equal ncol(YLW) or nrow(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(YLW)){
        stop("gplssvd: YLW is empty (i.e., all 0s")
      }
    }

    if(YLW_is_vector){
      if(length(YLW)!=Y_dimensions[1]){
        stop("gplssvd: length(YLW) does not equal nrow(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(YLW)<=tol)){
      if(!are_all_values_positive(YLW)){
        stop("gplssvd: YLW is not strictly positive values")
      }
    }
  }

  # a few things about YRW for stopping conditions
  YRW_is_missing <- missing(YRW)
  if(!YRW_is_missing){

    YRW_is_vector <- is.vector(YRW)

    if(!YRW_is_vector){

      if( nrow(YRW) != ncol(YRW) | nrow(YRW) != Y_dimensions[2] ){
        stop("gplssvd: nrow(YRW) does not equal ncol(YRW) or ncol(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      if(is_empty_matrix(YRW)){
        stop("gplssvd: YRW is empty (i.e., all 0s")
      }
    }

    if(YRW_is_vector){
      if(length(YRW)!=Y_dimensions[2]){
        stop("gplssvd: length(YRW) does not equal ncol(Y)")
      }

      # if you gave me all zeros, I'm stopping.
      # if(all(abs(YRW)<=tol)){
      if(!are_all_values_positive(YRW)){
        stop("gplssvd: YRW is not strictly positive values")
      }
    }
  }

  ## convenience checks *could* be removed* if problematic
  # convenience checks & conversions; these are meant to minimize XLW's memory footprint
  if(!XLW_is_missing){

    ## this is a matrix call; THIS MAKES IT GO MISSING AND NEXT CHECKS FAIL.
      ### technically, the following two checks actually handle this.
      ### it is a diagonal, and then it's a vector that's all 1s
    if( !XLW_is_vector ){

      # if(is_identity_matrix(XLW) ){
      #   XLW_is_missing <- T
      #   XLW <- substitute() # neat! this makes it go missing
      # }

      ## this is a matrix call
      if( is_diagonal_matrix(XLW) ){
        XLW <- diag(XLW)
        XLW_is_vector <- T  # now it's a vector
      }
    }

    ## this is a vector call
    if( XLW_is_vector & all(XLW==1) ){
      XLW_is_missing <- T
      XLW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize XRW's memory footprint
  if(!XRW_is_missing){
    if( !XRW_is_vector ){

      # if( is_identity_matrix(XRW) ){
      #   XRW_is_missing <- T
      #   XRW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(XRW) ){
        XRW <- diag(XRW)
        XRW_is_vector <- T  # now it's a vector
      }
    }

    if( XRW_is_vector & all(XRW==1) ){
      XRW_is_missing <- T
      XRW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize YLW's memory footprint
  if(!YLW_is_missing){
    if( !YLW_is_vector ){

      # if( is_identity_matrix(YLW) ){
      #   YLW_is_missing <- T
      #   YLW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(YLW) ){
        YLW <- diag(YLW)
        YLW_is_vector <- T  # now it's a vector
      }

    }

    if( YLW_is_vector & all(YLW==1) ){
      YLW_is_missing <- T
      YLW <- substitute() # neat! this makes it go missing
    }
  }

  # convenience checks & conversions; these are meant to minimize YRW's memory footprint
  if(!YRW_is_missing){
    if( !YRW_is_vector ){

      # if( is_identity_matrix(YRW) ){
      #   YRW_is_missing <- T
      #   YRW <- substitute() # neat! this makes it go missing
      # }

      if( is_diagonal_matrix(YRW) ){
        YRW <- diag(YRW)
        YRW_is_vector <- T  # now it's a vector
      }

    }

    if( YRW_is_vector & all(YRW==1) ){
      YRW_is_missing <- T
      YRW <- substitute() # neat! this makes it go missing
    }
  }


  # this manipulates X as needed based on XLW
  if(!XLW_is_missing){

    if( XLW_is_vector ){
      # X <- sweep(X,1,sqrt(XLW),"*") ## replace the sweep with * & t()
      X <- X * sqrt(XLW)
    }else{
      XLW <- as.matrix(XLW)
      # X <- (XLW %^% (1/2)) %*% X
      X <- sqrt_psd_matrix(XLW) %*% X
    }

  }
  # this manipulates X as needed based on XRW
  if(!XRW_is_missing){

    if( XRW_is_vector ){
      sqrt_XRW <- sqrt(XRW)
      # X <- sweep(X,2, sqrt_XRW,"*") ## replace the sweep with * & t()
      X <- t(t(X) * sqrt_XRW)
    }else{
      XRW <- as.matrix(XRW)
      # X <- X %*% (XRW %^% (1/2))
      X <- X %*% sqrt_psd_matrix(XRW)
    }

  }
  # this manipulates Y as needed based on YLW
  if(!YLW_is_missing){

    if( YLW_is_vector ){
      # Y <- sweep(Y,1,sqrt(YLW),"*")  ## replace the sweep with * & t()
      Y <- Y * sqrt(YLW)
    }else{
      YLW <- as.matrix(YLW)
      Y <- sqrt_psd_matrix(YLW) %*% Y
    }

  }
  # this manipulates Y as needed based on YRW
  if(!YRW_is_missing){

    if( YRW_is_vector ){
      sqrt_YRW <- sqrt(YRW)
      # Y <- sweep(Y,2,sqrt_YRW,"*")  ## replace the sweep with * & t()
      Y <- t(t(Y) * sqrt_YRW)
    }else{
      YRW <- as.matrix(YRW)
      # Y <- Y %*% (YRW %^% (1/2))
      Y <- Y %*% sqrt_psd_matrix(YRW)
    }

  }


  # all the decomposition things
  if(k<=0){
    k <- min(X_dimensions, Y_dimensions)
  }

  res <- tolerance_svd( t(X) %*% Y, nu=k, nv=k, tol=tol)
  res$d_full <- res$d
  res$l_full <- res$d_full^2
  # res$tau <- (res$l_full/sum(res$l_full)) * 100
  components.to.return <- min(length(res$d_full),k) #a safety check
  res$d <- res$d_full[1:components.to.return]
  res$l <- res$d^2
  res$u <- res$u[,1:components.to.return, drop = FALSE]
  res$v <- res$v[,1:components.to.return, drop = FALSE]

  res$lx <- X %*% res$u
  res$ly <- Y %*% res$v


  # make scores according to weights
  if(!XRW_is_missing){
    if(XRW_is_vector){

      # res$p <- sweep(res$u,1,1/sqrt_XRW,"*")
      res$p <- res$u / sqrt_XRW
      # res$fi <- sweep(sweep(res$p,1,XRW,"*"),2,res$d,"*")
      res$fi <- t(t(res$p * XRW) * res$d)

    }else{

      # res$p <- (XRW %^% (-1/2)) %*% res$u
      res$p <- invsqrt_psd_matrix(XRW) %*% res$u
      # res$fi <- sweep((XRW %*% res$p),2,res$d,"*")
      res$fi <- t(t(XRW %*% res$p) * res$d)

    }
  }else{

    res$p <- res$u
    # res$fi <- sweep(res$p,2,res$d,"*")
    res$fi <- t(t(res$p) * res$d)

  }

  if(!YRW_is_missing){
    if(YRW_is_vector){

      # res$q <- sweep(res$v,1,1/sqrt_YRW,"*")
      res$q <- res$v /sqrt_YRW
      # res$fj <- sweep(sweep(res$q,1,YRW,"*"),2,res$d,"*")
      res$fj <- t(t(res$q * YRW) * res$d)

    }else{

      # res$q <- (YRW %^% (-1/2)) %*% res$v
      res$q <- invsqrt_psd_matrix(YRW) %*% res$v
      # res$fj <- sweep((YRW %*% res$q),2,res$d,"*")
      res$fj <- t(t(YRW %*% res$q) * res$d)

    }
  }else{

    res$q <- res$v
    # res$fj <- sweep(res$q,2,res$d,"*")
    res$fj <- t(t(res$q) * res$d)

  }

  rownames(res$fi) <- rownames(res$u) <- rownames(res$p) <- colnames(X)
  rownames(res$fj) <- rownames(res$v) <- rownames(res$q) <- colnames(Y)

  rownames(res$lx) <- rownames(X)
  rownames(res$ly) <- rownames(Y)

  class(res) <- c("gplssvd", "GSVD", "list")
  return(res)

}

