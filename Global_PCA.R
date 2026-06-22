

# =============================================================================
# penalized_pca.R
#
# Globally aligned PCA for multi-group data -- core method and criteria.
#
# Accompanies the paper
#   "Globally aligned Principal Component Analysis for multi-group data"
#   Hedayat Fathi, Marzia A. Cremona, Federico Severino.
#
# This file defines:
#   * penalized_pca()          -- the globally aligned PCA estimator
#   * the subspace (rank-r) performance criteria used in the paper:
#       compute_W_r()          -- average within-group variance (W)
#       compute_PVE_r()        -- average proportion of variance explained (PVE)
#       compute_A_r()          -- alignment index (A)
#       compute_S_r()          -- stability index (S)
#   * angle_error() / compute_angle_error_r()
#                              -- eigenvector-recovery angle (simulation only,
#                                 where the true directions are known)
#
# These functions are sourced by simulation.Rmd.
# Base R only; no packages required to define them.
# =============================================================================


# -----------------------------------------------------------------------------
# Core method: globally aligned PCA
# -----------------------------------------------------------------------------
# For each group g, the globally aligned covariance is
#       Sigma_g^(rho) = Sigma_g + V_r diag(rho) V_r^T,
# where V_r holds the first r global (pooled) principal directions. The leading
# eigenvectors of Sigma_g^(rho) are the globally aligned principal components of
# group g.
#
# Arguments:
#   X       numeric matrix (n x p): observations stacked across all groups.
#   groups  factor/vector of length n giving each row's group membership.
#   r       number of global directions to align with (r <= p).
#   rho     alignment strengths; length r, or a single value recycled to r.
#           In the paper, rho_m = tau * lambda_global_m.
#   center  logical; center variables using the pooled sample.
#   scale.  logical; scale variables using the pooled sample.
#
# Returns a list:
#   $global     pooled Sigma, its eigenvalues (eigvals) and eigenvectors (eigvecs)
#   $per_group  per group: Sigma_k, aligned Sigma_k_rho, and its eigendecomposition
#   $params     the r and rho actually used

penalized_pca <- function(X, groups, r = 1, rho = rep(1, r),
                          center = TRUE, scale. = TRUE) {
  X      <- as.matrix(X)
  groups <- as.factor(groups)
  p      <- ncol(X)

  if (center || scale.) {
    X <- scale(X, center = center, scale = scale.)
  }

  # Pooled (global) covariance and its eigendecomposition
  Sigma_global <- cov(X)
  eig_global   <- eigen(Sigma_global, symmetric = TRUE)
  Vg <- eig_global$vectors
  Lg <- eig_global$values

  if (r > p) stop("r cannot exceed number of variables p.")
  if (length(rho) == 1) rho <- rep(rho, r)

  # Low-rank penalty V_r diag(rho) V_r^T, shared by all groups
  Vr          <- Vg[, 1:r, drop = FALSE]
  penalty_mat <- Vr %*% diag(rho, nrow = r) %*% t(Vr)

  # Per-group aligned covariance and its eigendecomposition
  per_group <- lapply(levels(groups), function(g) {
    idx         <- which(groups == g)
    Xg          <- X[idx, , drop = FALSE]
    Sigma_k     <- cov(Xg)
    Sigma_k_rho <- Sigma_k + penalty_mat
    eig_k       <- eigen(Sigma_k_rho, symmetric = TRUE)
    list(group       = g,
         Sigma_k     = Sigma_k,
         Sigma_k_rho = Sigma_k_rho,
         eig         = eig_k)
  })
  names(per_group) <- levels(groups)

  list(
    global    = list(Sigma = Sigma_global, eigvals = Lg, eigvecs = Vg),
    per_group = per_group,
    params    = list(r = r, rho = rho)
  )
}


# -----------------------------------------------------------------------------
# Subspace (rank-r) performance criteria
# -----------------------------------------------------------------------------
# Each criterion takes:
#   V_list          list of K matrices, each p x r (one r-dim subspace per group)
#   Sigma_hat_list  list of K sample group covariances Sigma_k
#   V_global_r      p x r matrix of the leading r global directions (for A)

# W: average over groups of the variance captured by the r-dim subspace.
compute_W_r <- function(V_list, Sigma_hat_list) {
  K  <- length(V_list)
  Wk <- numeric(K)
  for (k in seq_len(K))
    Wk[k] <- sum(diag(t(V_list[[k]]) %*% Sigma_hat_list[[k]] %*% V_list[[k]]))
  mean(Wk)
}

# PVE: W normalised by the total variance in each group, averaged over groups.
compute_PVE_r <- function(V_list, Sigma_hat_list) {
  K <- length(V_list); pve_k <- numeric(K)
  for (k in seq_len(K))
    pve_k[k] <- sum(diag(t(V_list[[k]]) %*% Sigma_hat_list[[k]] %*% V_list[[k]])) /
                sum(diag(Sigma_hat_list[[k]]))
  mean(pve_k)
}

# A (alignment index): average mean squared cosine of the principal angles
# between each group subspace and the global r-dim subspace. A in [0, 1].
compute_A_r <- function(V_list, V_global_r) {
  K <- length(V_list); r <- ncol(V_global_r); Ak <- numeric(K)
  for (k in seq_len(K)) {
    M     <- t(V_global_r) %*% qr.Q(qr(V_list[[k]]))
    Ak[k] <- (1 / r) * sum(M * M)
  }
  mean(Ak)
}

# S (stability index): average over group pairs of the mean squared sine of the
# principal angles between the two group subspaces. S in [0, 1]; 0 iff identical.
compute_S_r <- function(V_list) {
  r <- ncol(V_list[[1]])
  K <- length(V_list); if (K <= 1L) return(0)
  total <- 0; count <- 0
  for (k in 1:(K - 1)) {
    for (l in (k + 1):K) {
      M     <- t(V_list[[k]]) %*% V_list[[l]]
      total <- total + (r - sum(M * M)); count <- count + 1
    }
  }
  (total / count) / r
}


# -----------------------------------------------------------------------------
# Eigenvector-recovery angle (simulation only)
# -----------------------------------------------------------------------------
# Used to measure how well an estimated direction recovers a known true
# direction; meaningful only in the simulation, where u_true is available.

# Sign-invariant angle (radians) between two vectors.
angle_error <- function(u_true, v_hat) {
  u_true  <- as.numeric(u_true); v_hat <- as.numeric(v_hat)
  u_true  <- u_true / sqrt(sum(u_true^2))
  v_hat   <- v_hat  / sqrt(sum(v_hat^2))
  cos_val <- min(1, max(-1, abs(sum(u_true * v_hat))))
  acos(cos_val)
}

# Mean recovery angle over all groups and the r leading directions.
#   V_list         list of K matrices, each p x r (estimated directions)
#   u_true_array   p x r x K array of true group directions
compute_angle_error_r <- function(V_list, u_true_array) {
  K  <- length(V_list); r <- ncol(V_list[[1]])
  angles <- numeric(K * r); idx <- 1
  for (k in seq_len(K))
    for (j in seq_len(r)) {
      angles[idx] <- angle_error(u_true_array[, j, k], V_list[[k]][, j])
      idx <- idx + 1
    }
  mean(angles)
}
