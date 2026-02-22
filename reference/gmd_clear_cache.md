# Clear internal cache for matrix decompositions

Clears the internal cache used by generalized matrix decomposition
functions. This can be useful to free up memory or when working with
different datasets.

## Usage

``` r
gmd_clear_cache()
```

## Value

Invisibly returns TRUE after clearing the cache.

## Examples

``` r
# Clear the internal cache
gmd_clear_cache()
```
