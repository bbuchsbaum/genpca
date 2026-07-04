# Coordinate descent for the SFPCA penalized quadratic subproblem

Internal solver for \`min_x 0.5 x'Sx - b'x + P(x; lambda)\` with sparse
SPD \`S\` and an L1 or SCAD penalty.

## Usage

``` r
sfpca_cd_solve_cpp(S, b, x0, lambda, penalty, scad_a, max_sweeps, tol)
```

## Arguments

- S:

  sparse SPD matrix (\`dgCMatrix\`)

- b:

  numeric vector, linear term

- x0:

  numeric vector, warm start

- lambda:

  penalty level (must be \>= 0)

- penalty:

  0 for L1, 1 for SCAD

- scad_a:

  SCAD shape parameter (\> 2)

- max_sweeps:

  maximum number of full-equivalent sweeps

- tol:

  convergence tolerance on the KKT residual (gradient units)

## Value

list with \`x\`, \`sweeps\`, and \`kkt\` (max KKT residual)
