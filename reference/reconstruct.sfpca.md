# Reconstruct data from an sfpca fit

Reconstructs the rank-`K`
[`sfpca`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md) model
as \\U D V'\\, using the stored (non-orthogonal) factors directly.

## Usage

``` r
# S3 method for class 'sfpca'
reconstruct(
  x,
  comp = seq_len(multivarious::ncomp(x)),
  rowind = NULL,
  colind = NULL,
  ...
)
```

## Arguments

- x:

  An `sfpca` object.

- comp:

  Integer vector of components to use (default: all).

- rowind:

  Optional integer vector of rows to reconstruct (default: all).

- colind:

  Optional integer vector of columns to reconstruct (default: all).

- ...:

  Ignored.

## Value

A numeric matrix of dimension `length(rowind) x length(colind)`, the
rank-`length(comp)` reconstruction \\U D V'\\ using sfpca's stored
non-orthogonal factors.

## Details

sfpca components are Euclidean unit vectors but are **not** mutually
orthogonal, so \\V'V \ne I\\. The inherited `reconstruct.bi_projector()`
method reconstructs through the Moore-Penrose pseudoinverse of the
loadings (`scores \%*\% pinv(V)`), which for non-orthogonal `V` does
**not** return the rank-`comp` model \\U D V'\\ that
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
actually fits and deflates with. This method instead computes
`scores(x)[rowind, comp] \%*\% t(components(x)[colind, comp])`, i.e. \\U
D V'\\ restricted to the requested rows/columns/components.
[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md) does
no preprocessing, so no inverse transform is applied.

## See also

[`sfpca()`](https://bbuchsbaum.github.io/genpca/reference/sfpca.md)
