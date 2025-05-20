# ──────────────────────────────────────────────────────────────
#  benchmark‑genpca.R   –  runtime comparison of algorithm paths
# ──────────────────────────────────────────────────────────────

library(genpca)        # your dev package
library(Matrix)
library(multivarious)
library(bench)
library(dplyr)
library(tidyr)
library(ggplot2)

## helper ──────────────────────────────────────────────────────
make_constraint <- function(type, dim, guarantee_psd = TRUE) {
  switch(
    type,
    "identity" = Diagonal(dim, 1),
    "dense"    = {
      M <- matrix(rnorm(dim * dim, sd = 0.2), dim, dim)
      if (guarantee_psd) {
        tcrossprod(M) / dim  # Guaranteed PSD
      } else {
        # Generate a symmetric matrix, likely indefinite
        (M + t(M)) / 2 
      }
    },
    "sparse"   = {
      # simple 1‑D chain Laplacian (tridiagonal)  dim × dim
      i <- rep(1:dim, each = 3L)
      j <- c(i, i - 1, i + 1)
      keep <- j >= 1 & j <= dim
      i <- i[keep]; j <- j[keep]
      x <- ifelse(i == j, 1, -0.5) # Base Laplacian weights
      S <- sparseMatrix(i, j, x = x, dims = c(dim, dim)) # PSD
      if (!guarantee_psd) {
         # Perturb diagonal slightly to make it potentially indefinite
         diag_perturb <- rnorm(dim, mean = 0, sd = 0.1) 
         S <- S + Diagonal(dim, x = diag_perturb)
         S <- (S + t(S)) / 2 # Ensure symmetry after perturbation
      }
      S
    }
  )
}

## benchmark grid ──────────────────────────────────────────────
grid  <- expand_grid(
  n        = c(200, 500, 1000),
  p        = c(100, 300, 600),
  constraint_type = c("identity", "dense", "sparse"), # Renamed for clarity
  psd_status      = c("psd", "non-psd"),
  method          = c("eigen", "spectra", "deflation")
) |>
  # Exclude identity for non-psd case as identity is always PSD
  filter(!(constraint_type == "identity" & psd_status == "non-psd")) |>
  arrange(n, p, constraint_type, psd_status, method)

## sanity checks for C++ availability so we don't hard‑fail
has_cpp_defl   <- exists("gmd_deflation_cpp", mode = "function")
has_cpp_spectr <- exists("gmd_fast_cpp",      mode = "function")

## run! ────────────────────────────────────────────────────────
run_case <- function(n, p, constraint_type, psd_status, method) {

  if (method == "spectra"   && !has_cpp_spectr) return(NULL)
  if (method == "deflation" && !has_cpp_defl)   return(NULL)

  browser()
  X <- matrix(rnorm(n * p), n, p)
  
  # Generate constraints based on psd_status
  guarantee_psd_flag <- (psd_status == "psd")
  A <- make_constraint(constraint_type, p, guarantee_psd = guarantee_psd_flag)
  M <- make_constraint(constraint_type, n, guarantee_psd = guarantee_psd_flag)

  ## 10 components or as many as possible
  nc <- min(10L, n, p)
  
  # Set remedy if testing non-PSD case
  remedy_arg <- if (psd_status == "non-psd") "ridge" else "error"

  ## keep everything inside bench::mark so timing is clean
  bench::mark(
    genpca(
      X,
      A = A,
      M = M,
      ncomp  = nc,
      method = method,
      constraints_remedy = remedy_arg,
      preproc = center(),
      use_cpp = TRUE,          # safe; ignored unless needed
      verbose = FALSE
    ),
    iterations = 1,
    check = FALSE
  ) |>
    mutate(n = n, p = p, constraint_type = constraint_type, 
           psd_status = psd_status, method = method)
}

bench_tbl <- pmap_dfr(grid, ~run_case(..1, ..2, ..3, ..4, ..5)) |>
             filter(!is.null(expression))      # drop skipped specs

## ──────────────────────────────────────────────────────────────
##  visualise
## ──────────────────────────────────────────────────────────────
plt <- bench_tbl |>
  mutate(dim_pair = paste0("n=", n, "\np=", p)) |>
  # Add psd_status to interaction for distinct bars
  mutate(method_psd = interaction(method, psd_status, sep = "\n(")) |>
  mutate(method_psd = gsub("\\)$", ")", method_psd)) |> # Clean up label
  ggplot(aes(x = method, # Keep original method on x-axis label
             y = median, # Using median time
             fill = psd_status)) + # Color by PSD status
  geom_col(aes(group = method_psd), # Group bars by method & psd status
           position = position_dodge(width = 0.8), width = 0.7) + 
  scale_y_continuous(trans = "log10",
                     name = "runtime  (median,  log₁₀ seconds)") +
  scale_fill_brewer(palette = "Paired", name = "Constraint Status") + # Legend for color
  facet_grid(constraint_type ~ dim_pair, labeller = label_parsed) +
  labs(
    title    = "genpca: runtime by algorithm path and constraint PSD status",
    subtitle = "Bars = 1 run per scenario (median wall-clock); log-scale. Remedy='ridge' used for non-PSD.",
    x        = "Algorithm Method"
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1),
        panel.spacing = unit(0.6, "lines"),
        legend.position = "top") # Place legend at top

print(plt)

## save objects if you like
# readr::write_rds(bench_tbl, "benchmarks/genpca-benchmark.rds"

ggsave("bench/genpca-runtime.png", plt, width = 11, height = 6, dpi = 300)