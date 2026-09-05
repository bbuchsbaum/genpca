# Symmetrize a nearly symmetric matrix or stop

Measures `||A - A'||_F / ||A||_F`. Below `rtol` the two triangles are
averaged and the result is marked symmetric; above it the function
stops. Asymmetry is an input error, never something a PSD remedy
repairs.

## Usage

``` r
symmetrize_or_stop(A, rtol = 1e-10, name = "A")
```

## Arguments

- A:

  square numeric matrix or Matrix::Matrix

- rtol:

  relative asymmetry allowed (default 1e-10)

- name:

  label used in the error message

## Value

a symmetric Matrix (dense or sparse as supplied)
