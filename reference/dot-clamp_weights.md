# Validate diagonal weights

Weights must be finite and non-negative. Entries in
`[-rtol * max(abs(w)), 0)` are set to exactly zero with a message;
anything below is an error. Every backend sees the same cleaned vector.

## Usage

``` r
.clamp_weights(w, rtol = .metric_rtol_default(), name = "A")
```
