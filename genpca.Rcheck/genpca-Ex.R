pkgname <- "genpca"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('genpca')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("genpca")
### * genpca

flush(stderr()); flush(stdout())

### Name: genpca
### Title: Generalised Principal Components Analysis (GPCA)
### Aliases: genpca gmdLA gmd_deflationR truncate.genpca reconstruct.genpca
###   ncomp.genpca
### Keywords: internal

### ** Examples

if (requireNamespace("RSpectra", quietly = TRUE) &&
    requireNamespace("multivarious", quietly = TRUE)) {
  set.seed(123)
  X <- matrix(stats::rnorm(200 * 100), 200, 100)
  rownames(X) <- paste0("R", 1:200)
  colnames(X) <- paste0("C", 1:100)

  # Standard PCA (A=I, M=I, centered) - using default method="eigen"
  gpca_std_eigen <- genpca(X, ncomp = 5, preproc = multivarious::center(), verbose = FALSE)
  
  # Standard PCA using Spectra method (requires C++ build)
  # gpca_std_spectra <- try(genpca(X, ncomp = 5, preproc = multivarious::center(), \n#'   #                              method="spectra", verbose = TRUE))
  # if (!inherits(gpca_std_spectra, "try-error")) {
  #    print(head(gpca_std_spectra$sdev))
  # }

  # Compare singular values with prcomp
  pr_std <- stats::prcomp(X, center = TRUE, scale. = FALSE)
  print("Eigen Method Sdev:")
  print(head(gpca_std_eigen$sdev))
  print("prcomp Sdev:")
  print(head(pr_std$sdev))
  print(paste("Total Var Explained (Eigen):", round(sum(gpca_std_eigen$propv)*100), "%"))

  # Weighted column PCA (diagonal A, no centering)
  col_weights <- stats::runif(100, 0.5, 1.5)
  gpca_weighted <- genpca(X, A = col_weights, ncomp = 3, preproc=multivarious::pass(), verbose = FALSE)
  print("Weighted GPCA Sdev:")
  print(gpca_weighted$sdev)
  print(head(loadings(gpca_weighted)))
}



cleanEx()
nameEx("genpca_cov")
### * genpca_cov

flush(stderr()); flush(stdout())

### Name: genpca_cov
### Title: Generalized PCA on a covariance matrix
### Aliases: genpca_cov

### ** Examples

# Standard PCA on covariance (no constraint)
C <- cov(scale(iris[,1:4], center=TRUE, scale=FALSE))
fit0 <- genpca_cov(C, R=NULL, ncomp=3)
print(fit0$d[1:3])       # first 3 singular values
print(fit0$propv[1:3])   # variance explained by first 3 components

# Variable weights via a diagonal metric
w <- c(1, 1, 0.5, 2)  # emphasize Sepal.Width less, Petal.Width more
fitW <- genpca_cov(C, R = w, ncomp=3, method="gmd")
print(fitW$d[1:3])

# Compare GMD and generalized eigenvalue approaches
fit_gmd <- genpca_cov(C, R = w, ncomp=2, method="gmd")
fit_geigen <- genpca_cov(C, R = w, ncomp=2, method="geigen")
# These will generally differ unless R = I
print(paste("GMD singular values:", paste(round(fit_gmd$d, 3), collapse=", ")))
print(paste("GEigen singular values:", paste(round(fit_geigen$d, 3), collapse=", ")))




cleanEx()
nameEx("sfpca")
### * sfpca

flush(stderr()); flush(stdout())

### Name: second_diff_matrix
### Title: Sparse and Functional Principal Components Analysis (SFPCA) with
###   Spatial Coordinates
### Aliases: second_diff_matrix sfpca construct_spatial_penalty
### Keywords: internal

### ** Examples

library(Matrix)
set.seed(123)
# Simulate a small example due to resource constraints
n <- 100  # Number of time points
p <- 50   # Number of spatial locations
X <- matrix(rnorm(n * p), n, p)
# Simulate spatial coordinates
spat_cds <- matrix(runif(p * 3), nrow = 3, ncol = p)  # 3D coordinates
result <- sfpca(X, K = 2, spat_cds = spat_cds)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
