context("genpca")

mat_10_10 <- matrix(rnorm(10*10), 10, 10)
library(multivarious)

test_that("pca and genpca have same results with identity matrix for row and column constraints", {
  res1 <- genpca(mat_10_10)
  res2 <- pca(mat_10_10,ncomp=ncomp(res1))
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(sum(diffscores) < 1e-5)
  expect_equal(sdev(res1), sdev(res2))
  
  expect_equal(apply(components(res1), 2, function(x) sum(x^2)), rep(1, ncomp(res1)))
})

test_that("gen_pca with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=center())
  res2 <- pca(mat_10_10, preproc=standardize())
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(abs(sum(diffscores)) < 1e-5)
  expect_equal(sdev(res1), sdev(res2))
  
})

test_that("gen_pca with use_cpp with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=center(), use_cpp=TRUE)
  res2 <- pca(mat_10_10, preproc=standardize())
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(abs(sum(diffscores)) < 1e-5)
  expect_equal(sdev(res1), sdev(res2))
  
})

test_that("gen_pca with use_cpp (+ deflation) with column variances is equivalent to a scaled pca", {
  wts <- 1/apply(mat_10_10, 2, var)
  res1 <- genpca(mat_10_10, A=wts, preproc=center(), ncomp=9, deflation=TRUE, use_cpp=TRUE, threshold=1e-7)
  res2 <- pca(mat_10_10, preproc=standardize(), ncomp=9)
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(sdev(res1), sdev(res2), tolerance=.01)
  
})

test_that("gen_pca with use_cpp (+ deflation and n < p) with column variances is equivalent to a scaled pca", {
  mat_10_20 <- matrix(rnorm(10*20), 10, 20)
  wts <- 1/apply(mat_10_20, 2, var)
  res1 <- genpca(mat_10_20, A=wts, preproc=center(), ncomp=9, deflation=TRUE, use_cpp=FALSE, threshold=1e-8)
  res2 <- pca(mat_10_20, preproc=standardize(), ncomp=9)
  
  diffscores <- abs(scores(res1)) - abs(scores(res2))
  expect_true(mean(abs(diffscores)) < .01)
  expect_equal(sdev(res1), sdev(res2), tolerance=.01)
  
})

test_that("gen_pca with dense column and row constraints works", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  res1 <- genpca(mat_10_10, A=A, M=M, preproc=center())
  expect_equal(ncomp(res1),length(sdev(res1)))
})

test_that("gen_pca with sparse column and row constraints works", {
  A <- neighborweights:::adjacency.neighbor_graph(neighborweights::graph_weights(mat_10_10, k=8))
  Matrix::diag(A) <- 1
  M <- neighborweights:::adjacency.neighbor_graph(neighborweights::graph_weights(t(mat_10_10), k=3))
  Matrix::diag(M) <- 1.5
  res1 <- genpca(mat_10_10, A=A, M=M, preproc=center())
})


test_that("can reconstruct a genpca model with component selection", {
  A <- cov(matrix(rnorm(20*10), 20,10))
  M <- cov(matrix(rnorm(20*10), 20,10))
  res1 <- genpca(mat_10_10, preproc=center())
  recon1 <- reconstruct(res1)
  expect_equal(as.matrix(recon1), mat_10_10, check.attributes=FALSE)
  
  res1 <- genpca(mat_10_10, A=A, M=M, ncomp=10, preproc=center())
  res2 <- pca(mat_10_10,ncomp=10, preproc=center())
  recon2 <- reconstruct(res2)
  
   
  
})

test_that("can project a row vector", {
  A <- cov(matrix(rnorm(10*10),10,10))
  M <- cov(matrix(rnorm(10*10),10,10))
  
  res1 <- genpca(mat_10_10, A=A, M=M)
  p <- project(res1, mat_10_10[1,])
  expect_equal(dim(p), c(1,ncomp(res1)))
})

#test_that("can extract residuals", {
#  res1 <- genpca(mat_10_10)
#  resid <- residuals(res1, ncomp=2, mat_10_10)
#  d <- sdev(res1)
#  expect_equal(sum(d[3:length(d)] ^2), sum(resid^2))
#})

test_that("can run genpca with deflation", {
  X <- matrix(rnorm(100),10,10)
  res1 <- genpca(X, preproc=center(), ncomp=5,deflation=TRUE)
  res2 <- genpca(X, preproc=center(), ncomp=5)
  expect_true(sum(abs(res1$u) - abs(res2$u)) < 1)
})

test_that("can run genpca with sparse weighting matrix", {
  X <- matrix(rnorm(10000*20),10000,20)
  A <- neighborweights::temporal_adjacency(1:20)
  A <- cov(as.matrix(A))
  M <- neighborweights::temporal_adjacency(1:10000)
  res1 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), M=M, preproc=center(), ncomp=5,deflation=TRUE)
  res2 <- genpca(X, A=A, M=M, preproc=center(), ncomp=5)
  expect_true(!is.null(res1))
})

test_that("can run genpca on a largeish matrix with deflation", {
  nr <- 1000
  nc <- 500
  X <- matrix(rnorm(nr*nc),nr,nc)
  A <- neighborweights::temporal_adjacency(1:nc)
  A <- t(A) %*% A
  
  M <- neighborweights::temporal_adjacency(1:nr)
  M <- t(M) %*% M
  
  res1 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=center(), ncomp=5,deflation=TRUE, threshold=1e-5)
  res2 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=center(), ncomp=5,deflation=TRUE, 
                 threshold=1e-5, use_cpp=FALSE)
  
  res3 <- genpca(X, A=Matrix::Matrix(A, sparse=TRUE), 
                 M=M, preproc=center(), ncomp=20,deflation=FALSE)
  
  expect_true(!is.null(res1))
})
