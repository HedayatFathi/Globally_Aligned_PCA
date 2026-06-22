# Globally Aligned PCA for Multi-Group Data

R code for the paper *"Globally aligned Principal Component Analysis for
multi-group data"* by Hedayat Fathi, Marzia A. Cremona, and Federico Severino.

Globally aligned PCA interpolates between group-wise and pooled PCA by penalizing
each group's covariance toward the global principal subspace through a single
alignment-strength parameter τ.

## Contents

- `penalized_pca.R` — the method (`penalized_pca`) and the performance criteria:
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
source("penalized_pca.R")

fit <- penalized_pca(X, groups, r = 3, rho = tau * lambda_global[1:3])
# X       : numeric data matrix
# groups  : grouping factor (one entry per row of X)
# r       : number of global directions to align with
# rho     : per-direction alignment strengths
```

## Citation

> Fathi, H., Cremona, M. A., and Severino, F. (2026).
> Globally aligned Principal Component Analysis for multi-group data.
> *Journal of Computational and Graphical Statistics* (under review).

## License

MIT License — see `LICENSE`.
