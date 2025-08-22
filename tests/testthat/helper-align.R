# Align columns of A to B up to permutation and sign; return correlation diag.
align_perm_sign <- function(A, B) {
  stopifnot(ncol(A) == ncol(B))
  
  # Check if clue package is available
  if (!requireNamespace("clue", quietly = TRUE)) {
    warning("Package 'clue' not available for optimal alignment. Using simpler matching.")
    # Simple greedy matching based on absolute correlations
    C <- abs(cor(A, B))
    perm <- integer(ncol(A))
    for (i in seq_len(ncol(A))) {
      j <- which.max(C[i, ])
      perm[i] <- j
      C[, j] <- 0  # Mark as used
    }
    return(list(corr = diag(abs(cor(A, B[, perm]))), perm = perm))
  }
  
  # Optimal matching using Hungarian algorithm
  C <- abs(cor(A, B))
  perm <- clue::solve_LSAP(C, maximum = TRUE)
  list(corr = C[cbind(seq_len(ncol(A)), perm)],
       perm = perm)
}