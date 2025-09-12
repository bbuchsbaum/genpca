#' Canonical Generalized PLS (alias)
#'
#' Convenience alias for `genpls()`; computes canonical generalized PLS
#' (PLS-SVD/GPLSSVD). See `?genpls` for full documentation.
#'
#' @inheritParams genpls
#' @return See `genpls()`
#' @export
genplsc <- function(X, Y,
                    Ax = NULL, Ay = NULL,
                    Mx = NULL, My = NULL,
                    ncomp = 2,
                    preproc_x = multivarious::pass(),
                    preproc_y = multivarious::pass(),
                    svd_backend = c("RSpectra", "irlba"),
                    svd_opts = list(tol = 1e-7, maxitr = 1000),
                    verbose = FALSE) {
  genpls(X = X, Y = Y,
         Ax = Ax, Ay = Ay,
         Mx = Mx, My = My,
         ncomp = ncomp,
         preproc_x = preproc_x,
         preproc_y = preproc_y,
         svd_backend = svd_backend,
         svd_opts = svd_opts,
         verbose = verbose)
}

