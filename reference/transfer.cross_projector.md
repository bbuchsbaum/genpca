# Compatibility wrapper for multivarious::transfer

Provides a transfer method that accepts \`source\`/\`target\` arguments
used in the tests. The method forwards to \`multivarious\`'s transfer
method which expects \`from\`/\`to\`.

## Usage

``` r
# S3 method for class 'cross_projector'
transfer(x, new_data, from = NULL, to = NULL, opts = list(), ...)
```

## Arguments

- x:

  A cross_projector object.

- new_data:

  New data to transfer.

- from:

  Source space ("X" or "Y").

- to:

  Target space ("X" or "Y").

- opts:

  Options list passed to multivarious::transfer.

- ...:

  Additional arguments. Legacy parameters \`source\` and \`target\` are
  accepted here for backwards compatibility but are deprecated; use
  \`from\` and \`to\` instead.

## Value

Matrix with transferred data.

## Examples

``` r
# \donttest{
# Generate sample data
set.seed(123)
n <- 50
X <- matrix(rnorm(n * 10), n, 10)
Y <- matrix(rnorm(n * 8), n, 8)

# Create a cross projector using rpls
fit <- rpls(X, Y, K = 2)

# Transfer new X data to Y space
new_X <- matrix(rnorm(10 * 10), 10, 10)
transferred <- transfer(fit, new_X, from = "X", to = "Y")
#> Error in transfer(fit, new_X, from = "X", to = "Y"): could not find function "transfer"
# }
```
