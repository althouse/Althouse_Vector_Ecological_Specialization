################################################################################
# compute_bifurcation_coefficients.R
#
# Reproducible numerical computation of the center-manifold bifurcation
# coefficients a_W, a_M, and b for the 2-host, 2-vector SIR-SI model.
#
# Method:  Castillo-Chavez & Song (2004) Theorem 4.1
#          "Dynamical models of tuberculosis and their applications"
#          Math. Biosci. Eng. 1(2):361-404
#
# Applied to: Althouse et al. (2012) multi-host, multi-vector dengue model
#   PLOS Negl. Trop. Dis. 6(11): e1928
#
# Parameters match Figure 2 of the main text.
#
# Author: Benjamin M. Althouse
# Date:   2026-02-25
################################################################################

cat("=========================================================================\n")
cat("  Center-manifold bifurcation coefficients\n")
cat("  Castillo-Chavez & Song (2004), Math. Biosci. Eng. 1(2):361-404\n")
cat("=========================================================================\n\n")

# ═══════════════════════════════════════════════════════════════════════════════
# 1. PARAMETERS (from Table 1 / Figure 2 caption of the main text)
# ═══════════════════════════════════════════════════════════════════════════════
H <- 2   # number of host species
V <- 2   # number of vector species

# Biting preference weights r[j,i]: vector j on host i  (day^{-1} vector^{-1})
r <- matrix(c(
  0.5,   0.001,   # vector 1: specialist on host 1
  0.001, 0.5      # vector 2: specialist on host 2
), nrow = 2, byrow = TRUE)

N_h <- c(1000, 1000)      # host populations (large primate, small primate)
N_v <- c(25000, 25000)    # vector populations

gamma_h <- c(1/4, 1/4)                    # host recovery rates (day^{-1})
mu_h    <- c(1/(60*365), 1/(15*365))      # host natural mortality (day^{-1})
alpha_h <- c(0, 0)                         # disease-induced mortality
mu_v    <- c(1/7, 1/7)                    # vector mortality (day^{-1})

# Derived quantities
Nbar <- sum(N_h)                            # total host population
alpha_norm <- r / rowSums(r)                # normalized preference weights
D <- as.numeric(alpha_norm %*% N_h)         # weighted effective host pop D_j

cat("Parameters:\n")
cat(sprintf("  Hosts:   N_h = (%d, %d)\n", N_h[1], N_h[2]))
cat(sprintf("  Vectors: N_v = (%d, %d)\n", N_v[1], N_v[2]))
cat(sprintf("  Biting prefs (diag):  r_11=%.3f, r_22=%.3f\n", r[1,1], r[2,2]))
cat(sprintf("  Biting prefs (off):   r_12=%.3f, r_21=%.3f\n", r[1,2], r[2,1]))
cat(sprintf("  Recovery:    gamma = (%.4f, %.4f)\n", gamma_h[1], gamma_h[2]))
cat(sprintf("  Host mort:   mu_h  = (%.6e, %.6e)\n", mu_h[1], mu_h[2]))
cat(sprintf("  Vector mort: mu_v  = (%.4f, %.4f)\n", mu_v[1], mu_v[2]))
cat(sprintf("  Nbar = %d,  D = (%.4f, %.4f)\n\n", Nbar, D[1], D[2]))

# ═══════════════════════════════════════════════════════════════════════════════
# 2. ODE RIGHT-HAND SIDE
# ═══════════════════════════════════════════════════════════════════════════════
# State vector: x = (I_h1, I_h2, R_h1, R_h2, I_v1, I_v2)
# With S_hi = N_hi - I_hi - R_hi  and  S_vj = N_vj - I_vj
#
# The model from Supplement S1 (with c_j = 0, intro = 0, sigma = kappa = 1):
#   dI_hi/dt = lambda_hi * S_hi - (gamma_i + mu_hi) * I_hi
#   dR_hi/dt = gamma_i * I_hi - mu_hi * R_hi
#   dI_vj/dt = lambda_vj * S_vj - mu_vj * I_vj
#
# where the FOIs differ ONLY in denominator:
#   Weighted:   denom_j = D_j = sum_k alpha_kj * N_hk
#   Unweighted: denom_j = Nbar = sum_k N_hk
#
# lambda_hi = sum_j  r[j,i] * beta * I_vj / denom_j
# lambda_vj = sum_k  r[j,k] * beta * I_hk / denom_j

f_rhs <- function(x, beta, model = "weighted") {
  I_h <- x[1:2]
  R_h <- x[3:4]
  I_v <- x[5:6]

  S_h <- N_h - I_h - R_h
  S_v <- N_v - I_v

  denom <- if (model == "weighted") D else rep(Nbar, V)

  # Host FOI: lambda_h[i] = sum_j r[j,i] * beta * I_v[j] / denom[j]
  lambda_h <- numeric(H)
  for (i in 1:H)
    lambda_h[i] <- sum(r[, i] * beta * I_v / denom)

  # Vector FOI: lambda_v[j] = sum_k r[j,k] * beta * I_h[k] / denom[j]
  lambda_v <- numeric(V)
  for (j in 1:V)
    lambda_v[j] <- sum(r[j, ] * beta * I_h / denom[j])

  dI_h <- lambda_h * S_h - (gamma_h + mu_h + alpha_h) * I_h
  dR_h <- gamma_h * I_h - mu_h * R_h
  dI_v <- lambda_v * S_v - mu_v * I_v

  c(dI_h, dR_h, dI_v)
}

n <- 6  # dimension of state space

# ═══════════════════════════════════════════════════════════════════════════════
# 3. NEXT-GENERATION MATRIX → R0  AND  CRITICAL BETA
# ═══════════════════════════════════════════════════════════════════════════════
# Following van den Driessche & Watmough (2002).
# For the infected subsystem (I_h1, I_h2, I_v1, I_v2):
#   F = new-infection matrix, V = transition matrix
#   K = F V^{-1},  R0 = rho(K)

compute_R0 <- function(beta, model = "weighted") {
  denom <- if (model == "weighted") D else rep(Nbar, V)

  # F: new infections at DFE (S_hi = N_hi, S_vj = N_vj)
  # F is 4x4 for the infection subsystem (I_h1, I_h2, I_v1, I_v2)
  F_mat <- matrix(0, 4, 4)

  # Host rows (i=1,2): new host infections come from vectors
  for (i in 1:H) for (j in 1:V)
    F_mat[i, H + j] <- r[j, i] * beta * N_h[i] / denom[j]

  # Vector rows (j=1,2): new vector infections come from hosts
  for (j in 1:V) for (i in 1:H)
    F_mat[H + j, i] <- r[j, i] * beta * N_v[j] / denom[j]

  # V: transitions OUT of infected classes
  V_diag <- c(gamma_h + mu_h + alpha_h, mu_v)
  V_mat  <- diag(V_diag)

  K <- F_mat %*% solve(V_mat)
  max(Mod(eigen(K, only.values = TRUE)$values))
}

# Find critical beta where R0 = 1
find_beta_crit <- function(model) {
  f <- function(beta) compute_R0(beta, model) - 1
  uniroot(f, interval = c(1e-8, 2), tol = 1e-14)$root
}

beta_crit_W <- find_beta_crit("weighted")
beta_crit_M <- find_beta_crit("unweighted")

cat("Next-generation matrix analysis:\n")
cat(sprintf("  R0(beta=0.1, weighted)   = %.6f\n", compute_R0(0.1, "weighted")))
cat(sprintf("  R0(beta=0.1, unweighted) = %.6f\n", compute_R0(0.1, "unweighted")))
cat(sprintf("  Critical beta (weighted):   beta* = %.10f  (R0 = %.10f)\n",
            beta_crit_W, compute_R0(beta_crit_W, "weighted")))
cat(sprintf("  Critical beta (unweighted): beta* = %.10f  (R0 = %.10f)\n\n",
            beta_crit_M, compute_R0(beta_crit_M, "unweighted")))

# ═══════════════════════════════════════════════════════════════════════════════
# 4. JACOBIAN AT DFE (numerical, via central differences)
# ═══════════════════════════════════════════════════════════════════════════════

jacobian_num <- function(beta, model, eps = 1e-7) {
  x0 <- rep(0, n)
  J <- matrix(0, n, n)
  for (i in 1:n) {
    xp <- xm <- x0
    xp[i] <- eps
    xm[i] <- -eps
    J[, i] <- (f_rhs(xp, beta, model) - f_rhs(xm, beta, model)) / (2 * eps)
  }
  J
}

# Also: analytical Jacobian for verification
jacobian_analytic <- function(beta, model) {
  denom <- if (model == "weighted") D else rep(Nbar, V)
  J <- matrix(0, n, n)

  gamma_tot <- gamma_h + mu_h + alpha_h   # total removal from I_h

  # Row 1 (I_h1), Row 2 (I_h2)
  for (i in 1:H) {
    J[i, i] <- -(gamma_tot[i])        # d(dI_hi)/d(I_hi)
    for (j in 1:V)
      J[i, H + H + j] <- r[j, i] * beta * N_h[i] / denom[j]  # d(dI_hi)/d(I_vj)
  }

  # Row 3 (R_h1), Row 4 (R_h2)
  for (i in 1:H) {
    J[H + i, i]     <- gamma_h[i]     # d(dR_hi)/d(I_hi)
    J[H + i, H + i] <- -mu_h[i]       # d(dR_hi)/d(R_hi)
  }

  # Row 5 (I_v1), Row 6 (I_v2)
  for (j in 1:V) {
    for (i in 1:H)
      J[2*H + j, i] <- r[j, i] * beta * N_v[j] / denom[j]  # d(dI_vj)/d(I_hi)
    J[2*H + j, 2*H + j] <- -mu_v[j]   # d(dI_vj)/d(I_vj)
  }

  J
}

# ═══════════════════════════════════════════════════════════════════════════════
# 5. EIGENVECTORS FOR THE ZERO EIGENVALUE
# ═══════════════════════════════════════════════════════════════════════════════

get_null_eigenvectors <- function(J) {
  ev <- eigen(J)

  # Find eigenvalue closest to zero
  idx <- which.min(abs(Re(ev$values)))
  cat(sprintf("    Eigenvalue closest to 0: %+.4e + %.4ei\n",
              Re(ev$values[idx]), Im(ev$values[idx])))

  # Right eigenvector
  v <- Re(ev$vectors[, idx])

  # Left eigenvector (= right eigenvector of J^T)
  ev_left <- eigen(t(J))
  idx_left <- which.min(abs(Re(ev_left$values)))
  w <- Re(ev_left$vectors[, idx_left])

  # Normalize: w . v = 1
  wv <- sum(w * v)
  w <- w / wv

  # Convention: make dominant v-component positive
  if (sum(v[c(1,2,5,6)]) < 0) {  # infection compartments
    v <- -v
    w <- -w / sum(w * v)     # re-normalize
    w <- w / sum(w * v)
  }

  list(v = v, w = w, eigenvalue = ev$values[idx])
}

# ═══════════════════════════════════════════════════════════════════════════════
# 6. COMPUTE BIFURCATION COEFFICIENTS  a  AND  b
#    a = sum_{k,i,j} w_k v_i v_j  d²f_k / (dx_i dx_j)   evaluated at DFE
#    b = sum_{k,i}   w_k v_i      d²f_k / (dx_i dbeta)   evaluated at DFE
# ═══════════════════════════════════════════════════════════════════════════════

compute_ab <- function(model) {
  beta_c <- if (model == "weighted") beta_crit_W else beta_crit_M

  cat(sprintf("\n── %s model ─────────────────────────────────\n", toupper(model)))
  cat(sprintf("  beta* = %.10f\n", beta_c))

  # Jacobian (use analytic for accuracy; verify with numerical)
  J_a  <- jacobian_analytic(beta_c, model)
  J_n  <- jacobian_num(beta_c, model)
  cat(sprintf("  ||J_analytic - J_numeric|| = %.2e\n", max(abs(J_a - J_n))))

  # Eigenvectors
  eig <- get_null_eigenvectors(J_a)
  v <- eig$v
  w <- eig$w

  cat(sprintf("  w . v = %.10f  (should be 1.0)\n", sum(w * v)))
  cat("  v = ["); cat(sprintf(" %+.6e", v)); cat(" ]\n")
  cat("  w = ["); cat(sprintf(" %+.6e", w)); cat(" ]\n")

  # ── Coefficient 'a': second partial derivatives of f w.r.t. state ──
  x0 <- rep(0, n)
  h  <- 1e-5

  a_val <- 0
  for (i in 1:n) {
    for (j in 1:n) {
      # Central difference for d²f_k / (dx_i dx_j)
      xpp <- x0; xpp[i] <- xpp[i] + h; xpp[j] <- xpp[j] + h
      xpm <- x0; xpm[i] <- xpm[i] + h; xpm[j] <- xpm[j] - h
      xmp <- x0; xmp[i] <- xmp[i] - h; xmp[j] <- xmp[j] + h
      xmm <- x0; xmm[i] <- xmm[i] - h; xmm[j] <- xmm[j] - h

      d2f <- (f_rhs(xpp, beta_c, model) - f_rhs(xpm, beta_c, model) -
              f_rhs(xmp, beta_c, model) + f_rhs(xmm, beta_c, model)) / (4*h^2)

      a_val <- a_val + sum(w * d2f) * v[i] * v[j]
    }
  }

  # ── Coefficient 'b': mixed partial d²f_k / (dx_i  dbeta) ──
  db <- 1e-7

  b_val <- 0
  for (i in 1:n) {
    xp <- xm <- x0
    xp[i] <- h
    xm[i] <- -h

    d2f_beta <- (f_rhs(xp, beta_c + db, model) - f_rhs(xp, beta_c - db, model) -
                 f_rhs(xm, beta_c + db, model) + f_rhs(xm, beta_c - db, model)) /
                (4 * h * db)

    b_val <- b_val + sum(w * d2f_beta) * v[i]
  }

  cat(sprintf("\n  a = %+.6e\n", a_val))
  cat(sprintf("  b = %+.6e\n", b_val))

  if (a_val > 0 && b_val > 0)
    cat("  ==> BACKWARD (subcritical) transcritical bifurcation\n")
  else if (a_val < 0 && b_val > 0)
    cat("  ==> FORWARD (supercritical) transcritical bifurcation\n")
  else
    cat(sprintf("  ==> Other case (a %s 0, b %s 0)\n",
                ifelse(a_val > 0, ">", "<"), ifelse(b_val > 0, ">", "<")))

  list(a = a_val, b = b_val, beta_crit = beta_c, v = v, w = w)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 7. ANALYTICAL VERIFICATION
# ═══════════════════════════════════════════════════════════════════════════════
# For this model the only nonzero second derivatives at the DFE come from
# the bilinear terms  lambda * S  (susceptible depletion), since D_j and
# Nbar are both constant.
#
# Specifically, in the reduced system x = (I_h1, I_h2, R_h1, R_h2, I_v1, I_v2):
#
#   d²f_{I_hi} / (dI_vj  dI_hi) = -r[j,i]*beta / denom_j
#   d²f_{I_hi} / (dI_vj  dR_hi) = -r[j,i]*beta / denom_j
#   d²f_{I_vj} / (dI_hk  dI_vj) = -r[j,k]*beta / denom_j
#
# (and the symmetric partners; all other d²f = 0).

compute_ab_analytic <- function(model) {
  beta_c <- if (model == "weighted") beta_crit_W else beta_crit_M
  denom <- if (model == "weighted") D else rep(Nbar, V)

  J <- jacobian_analytic(beta_c, model)
  eig <- get_null_eigenvectors(J)
  v <- eig$v
  w <- eig$w

  # Compute 'a' analytically from the known nonzero second derivatives
  a_val <- 0

  for (i in 1:H) {    # host index
    for (j in 1:V) {   # vector index
      coeff <- -r[j, i] * beta_c / denom[j]

      # f_{I_hi}: d²/dI_vj dI_hi  and  d²/dI_vj dR_hi
      vj_idx <- 2*H + j   # I_vj position in state vector
      ih_idx <- i          # I_hi position
      rh_idx <- H + i      # R_hi position

      # d²f_{ih}/d(I_vj)d(I_hi): contributes 2*w[ih]*v[vj]*v[ih]*coeff
      a_val <- a_val + 2 * w[ih_idx] * v[vj_idx] * v[ih_idx] * coeff
      # d²f_{ih}/d(I_vj)d(R_hi): contributes 2*w[ih]*v[vj]*v[rh]*coeff
      a_val <- a_val + 2 * w[ih_idx] * v[vj_idx] * v[rh_idx] * coeff

      # f_{I_vj}: d²/dI_hi dI_vj
      a_val <- a_val + 2 * w[vj_idx] * v[ih_idx] * v[vj_idx] * coeff
    }
  }

  # Compute 'b' analytically
  # d²f_{I_hi}/(dI_vj  dbeta) = r[j,i]*N_hi/denom_j   (at DFE)
  # d²f_{I_vj}/(dI_hk  dbeta) = r[j,k]*N_vj/denom_j   (at DFE)
  b_val <- 0

  for (i in 1:H) for (j in 1:V) {
    ih_idx <- i
    vj_idx <- 2*H + j

    # Host direction
    b_val <- b_val + w[ih_idx] * v[vj_idx] * r[j,i] * N_h[i] / denom[j]
    # Vector direction
    b_val <- b_val + w[vj_idx] * v[ih_idx] * r[j,i] * N_v[j] / denom[j]
  }

  list(a = a_val, b = b_val)
}

# ═══════════════════════════════════════════════════════════════════════════════
# 8. RUN THE COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════════

cat("── Verify Jacobians ──\n")
J_W <- jacobian_analytic(beta_crit_W, "weighted")
J_M <- jacobian_analytic(beta_crit_M, "unweighted")
cat("  Eigenvalues of J (weighted at beta*):\n    ")
cat(sprintf("%+.4e  ", sort(Re(eigen(J_W)$values), decreasing=TRUE)), "\n")
cat("  Eigenvalues of J (unweighted at beta*):\n    ")
cat(sprintf("%+.4e  ", sort(Re(eigen(J_M)$values), decreasing=TRUE)), "\n")

results_W <- compute_ab("weighted")
results_M <- compute_ab("unweighted")

# Cross-check with analytic
cat("\n── Analytic cross-check ──\n")
cat("  (Using closed-form second derivatives)\n")
suppressMessages({
  sink("/dev/null")  # suppress eigenvector printout
  ana_W <- compute_ab_analytic("weighted")
  ana_M <- compute_ab_analytic("unweighted")
  sink()
})
cat(sprintf("  a_W (numeric) = %+.6e    a_W (analytic) = %+.6e    diff = %.2e\n",
            results_W$a, ana_W$a, abs(results_W$a - ana_W$a)))
cat(sprintf("  a_M (numeric) = %+.6e    a_M (analytic) = %+.6e    diff = %.2e\n",
            results_M$a, ana_M$a, abs(results_M$a - ana_M$a)))
cat(sprintf("  b_W (numeric) = %+.6e    b_W (analytic) = %+.6e    diff = %.2e\n",
            results_W$b, ana_W$b, abs(results_W$b - ana_W$b)))
cat(sprintf("  b_M (numeric) = %+.6e    b_M (analytic) = %+.6e    diff = %.2e\n",
            results_M$b, ana_M$b, abs(results_M$b - ana_M$b)))

# ═══════════════════════════════════════════════════════════════════════════════
# 9. REPORT RESULTS
# ═══════════════════════════════════════════════════════════════════════════════

cat("\n=========================================================================\n")
cat("  RESULTS\n")
cat("=========================================================================\n")
cat(sprintf("  a_W = %+.4e    (weighted / specialist)\n", results_W$a))
cat(sprintf("  a_M = %+.4e    (unweighted / opportunistic)\n", results_M$a))
cat(sprintf("  b_W = %+.4e\n", results_W$b))
cat(sprintf("  b_M = %+.4e\n", results_M$b))
cat("=========================================================================\n")
