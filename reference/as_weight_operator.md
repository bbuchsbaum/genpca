# Create Weight Operator Function

Returns a closure that applies a weight matrix W or its transformations
(square root, inverse, or combinations thereof) to vectors/matrices.

## Usage

``` r
as_weight_operator(W, transpose = FALSE, sqrt = FALSE, inverse = FALSE)
```

## Arguments

- W:

  A weight matrix (SPD) or NULL for identity

- transpose:

  Logical, whether to transpose W before applying

- sqrt:

  Logical, whether to use square root of W

- inverse:

  Logical, whether to use inverse of W

## Value

A function that applies the requested transformation
