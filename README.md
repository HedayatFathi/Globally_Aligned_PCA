# Globally Aligned PCA for Multi-Group Data

R code for the paper *"Globally aligned Principal Component Analysis for
multi-group data"* by Hedayat Fathi, Marzia A. Cremona, and Federico Severino.

Globally aligned PCA interpolates between group-wise and pooled PCA by penalizing
each group's covariance toward the global principal subspace through a single
alignment-strength parameter τ.

## Contents

- `Global_PCA.R` — the method (`Global_PCA`) and the performance criteria:
  average within-group variance (W), proportion of variance explained (PVE),
  alignment index (A), and stability index (S).
- `simulation.Rmd` — the Monte Carlo simulation study, which reproduces the
  simulation tables in the paper.

## Requirements

R (≥ 4.2) with the following packages:

```r
install.packages(c("ggplot2", "MASS", "dplyr", "purrr",
                   "tidyr", "tibble", "knitr"))
```

## Usage

```r
# Run the simulation study
rmarkdown::render("simulation.Rmd")
```

The core function can also be used on your own data:

```r
source("Global_PCA.R")

fit <- Global_PCA(X, groups, r = 3, rho = tau * lambda_global[1:3])
# X       : numeric data matrix
# groups  : grouping factor (one entry per row of X)
# r       : number of global directions to align with
# rho     : per-direction alignment strengths
```

## Citation

> Fathi, H., Cremona, M. A., and Severino, F. (2026).
> Globally aligned Principal Component Analysis for multi-group data.
> submitted.

## License

MIT License — see `LICENSE`.
