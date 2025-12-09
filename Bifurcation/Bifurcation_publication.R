################################################################################
# bifurcation_publication.R
# 
# Produce bifurcation analysis figure for the manuscript:
#   Althouse BM, *et al.* (2025) — FOI denominators and multi-host dengue.
# 
# Generates a two-panel bifurcation diagram showing endemic equilibria as a
# function of R0, contrasting the weighted (Eq. 3) and unweighted (Eq. 4)
# force-of-infection formulations.
#
# Simulation results from:
#   Althouse BM, Lessler J, Sall AA, Diallo M, Hanley KA, et al. (2012)
#   *Synchrony of Sylvatic Dengue Isolations: A Multi-Host, Multi-Vector
#   SIR Model of Dengue Virus Transmission in Senegal.
#   PLOS Neglected Tropical Diseases 6(11): e1928.
#   10.1371/journal.pntd.0001928
#
# The script:
#   1. Creates a date-stamped output folder.
#   2. Loads cached bifurcation data (or reads CSV files if cache absent).
#   3. Calculates analytical R0 for each beta value.
#   4. Generates a 6.8 × 5-inch two-panel PDF showing prevalence vs. R0.
#
# Source this file or run it with `Rscript bifurcation_publication.R`.
#
# © 2025 Benjamin M. Althouse – CC BY 4.0
################################################################################

## -------- 1  housekeeping ----------------------------------------------------
rm(list = ls())

# resolve script location even outside RStudio -------------------------------
this.file <- if (interactive() && requireNamespace("rstudioapi", quietly = TRUE)) {
  tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
} else {
  parent.frame(2)$ofile %||% ""
}
if (nzchar(this.file)) setwd(dirname(this.file))

# output folder yyyy-mm-dd ----------------------------------------------------
out.dir <- file.path(sprintf("output_%s", format(Sys.Date(), "%Y-%m-%d")))
if (!dir.exists(out.dir)) {
  dir.create(out.dir, recursive = TRUE)
  cat("Created output directory:", out.dir, "\n")
}

# ---- packages ----------------------------------------------------------------
library(RColorBrewer)

## -------- 2  source analytic R0 functions ------------------------------------
source(file.path("..", "..", "R0Compare", "r0numerical2specIntro.r"))
source(file.path("..", "..", "R0Compare", "r0numerical2specIntroOldDenom.r"))

## -------- 3  parameters & data paths -----------------------------------------
new.dir  <- file.path("..", "NewForm-FixSmallPrimate")
old.dir  <- file.path("..", "OldForm-FixSmallPrimate")

n.sim   <- 100                       # number of beta lines (1…100)
disc    <- 50                        # matches original beta step
beta.x  <- ((disc / n.sim) * (1:n.sim)) / 100 + .025 - disc / n.sim / 100

# Fixed epidemiological parameters -------------------------------------------
rrmh  <- 0.5                         # biting preference: vector 1 → large host
rrmp  <- 0.001                       # biting preference: vector 1 → small host
rrm2h <- 0.001                       # biting preference: vector 2 → large host
rrm2p <- 0.5                         # biting preference: vector 2 → small host

vh  <- 1 / (60 * 365)                # mortality rate: large host
vp  <- 1 / (15 * 365)                # mortality rate: small host
vm  <- 1 / 7                         # mortality rate: vector 1
vm2 <- 1 / 7                         # mortality rate: vector 2

gh <- 1 / 4                          # recovery rate: large host
gp <- 1 / 4                          # recovery rate: small host

Nh  <- 1000                          # population size: large host
Np  <- 1000                          # population size: small host
Nm1 <- 25000                         # population size: vector 1
Nm2 <- 25000                         # population size: vector 2
NN  <- Nh + Np                       # total host population

## -------- 4  calculate R0 for each beta value --------------------------------
message("Calculating R0 values …")

R0.new <- sapply(beta.x, function(bb) {
  r0new(bm1h = bb, bm1p1 = bb, bm2h = bb, bm2p1 = bb,
        bhm1 = bb, bp1m1 = bb, bhm2 = bb, bp1m2 = bb,
        rm1h = rrmh, rm1p1 = rrmp, rm2h = rrm2h, rm2p1 = rrm2p,
        Vm = vm, Vm2 = vm2, gh = gh, mh = vh, gp1 = gp, mp1 = vp,
        Nh = Nh, Np1 = Np, Nm1 = Nm1, Nm2 = Nm2, NN = NN)
})

R0.old <- sapply(beta.x, function(bb) {
  r0old(bm1h = bb, bm1p1 = bb, bm2h = bb, bm2p1 = bb,
        bhm1 = bb, bp1m1 = bb, bhm2 = bb, bp1m2 = bb,
        rm1h = rrmh, rm1p1 = rrmp, rm2h = rrm2h, rm2p1 = rrm2p,
        Vm = vm, Vm2 = vm2, gh = gh, mh = vh, gp1 = gp, mp1 = vp,
        Nh = Nh, Np1 = Np, Nm1 = Nm1, Nm2 = Nm2, NN = NN)
})

## -------- 5  load bifurcation data (cached or from CSV) ----------------------
rds.file <- "bifurcation_data_cache.rds"

if (file.exists(rds.file)) {
  message("Loading cached bifurcation data …")
  cached <- readRDS(rds.file)
  bifH.new <- cached$bifH_new
  bifP.new <- cached$bifP_new
  bifH.old <- cached$bifH_old
  bifP.old <- cached$bifP_old
  
} else {
  message("Cache not found. Loading CSV files …")
  
  # determine sampling indices (annual snapshots from 50-year runs) ----------
  temp <- read.csv(file.path(old.dir, "Human", "H-bb1.csv"), header = FALSE)
  idx  <- seq(1, nrow(temp), by = 365)
  
  # read large-host prevalence: weighted FOI ----------------------------------
  bifH.new <- matrix(NA, nrow = n.sim, ncol = length(idx))
  for (i in seq_len(n.sim)) {
    temp <- read.csv(file.path(new.dir, "Human", paste0("H-bb", i, ".csv")),
                     header = FALSE)
    bifH.new[i, ] <- temp[idx, 1]
  }
  
  # read small-host prevalence: weighted FOI ----------------------------------
  bifP.new <- matrix(NA, nrow = n.sim, ncol = length(idx))
  for (i in seq_len(n.sim)) {
    temp <- read.csv(file.path(new.dir, "Primate", paste0("P-bb", i, ".csv")),
                     header = FALSE)
    bifP.new[i, ] <- temp[idx, 1]
  }
  
  # read large-host prevalence: unweighted FOI --------------------------------
  bifH.old <- matrix(NA, nrow = n.sim, ncol = length(idx))
  for (i in seq_len(n.sim)) {
    temp <- read.csv(file.path(old.dir, "Human", paste0("H-bb", i, ".csv")),
                     header = FALSE)
    bifH.old[i, ] <- temp[idx, 1]
  }
  
  # read small-host prevalence: unweighted FOI --------------------------------
  bifP.old <- matrix(NA, nrow = n.sim, ncol = length(idx))
  for (i in seq_len(n.sim)) {
    temp <- read.csv(file.path(old.dir, "Primate", paste0("P-bb", i, ".csv")),
                     header = FALSE)
    bifP.old[i, ] <- temp[idx, 1]
  }
  
  # cache for future runs -----------------------------------------------------
  saveRDS(list(bifH_new = bifH.new, bifP_new = bifP.new, 
               bifH_old = bifH.old, bifP_old = bifP.old), 
          file = rds.file)
  message("Saved cache: ", rds.file)
}

message("Data loading complete.")

## -------- 6  plotting helper --------------------------------------------------
#' Overlay bifurcation branches for weighted vs. unweighted FOI
#'
#' Plots prevalence (per 1000) as a function of R0, showing endemic equilibria
#' for both the weighted-sum (Eq. 3) and simple-sum (Eq. 4) formulations.
#' Duplicate equilibria are removed for clarity.
#'
#' @param R0.new   numeric vector of R0 values (weighted FOI)
#' @param R0.old   numeric vector of R0 values (unweighted FOI)
#' @param bif.new  matrix of prevalence values (weighted FOI)
#' @param bif.old  matrix of prevalence values (unweighted FOI)
#' @param host.nm  character string for the panel title
#'
#' @keywords internal
plot_overlay <- function(R0.new, R0.old, bif.new, bif.old, host.nm) {
  # colour palette ------------------------------------------------------------
  cols      <- brewer.pal(11, name = "Paired")
  color.new <- cols[2]                      # light blue for weighted
  color.old <- cols[6]                      # orange for unweighted
  
  # scale prevalence to per 1000 ----------------------------------------------
  bif.new.1k <- bif.new * 1000
  bif.old.1k <- bif.old * 1000
  
  ylim <- range(bif.new.1k, bif.old.1k, na.rm = TRUE)
  xlim <- range(c(R0.new, R0.old), na.rm = TRUE)
  
  # round to remove near-duplicate equilibria ---------------------------------
  bif.new.rnd <- round(bif.new, 5)
  bif.old.rnd <- round(bif.old, 5)
  
  # empty frame ---------------------------------------------------------------
  plot(NULL, xlim = xlim, ylim = ylim, xlab = "", ylab = "", main = "",
       cex.lab = 1.2, cex.axis = 0.85, cex.main = 1.1)
  
  # plot weighted FOI (squares) -----------------------------------------------
  for (j in seq_len(nrow(bif.new))) {
    unique.pts <- !duplicated(bif.new.rnd[j, ])
    points(rep(R0.new[j], sum(unique.pts)), bif.new.1k[j, unique.pts], 
           pch = 0, cex = 1, col = paste0(color.new, "85"))
  }
  
  # plot unweighted FOI (circles) ---------------------------------------------
  for (j in seq_len(nrow(bif.old))) {
    unique.pts <- !duplicated(bif.old.rnd[j, ])
    points(rep(R0.old[j], sum(unique.pts)), bif.old.1k[j, unique.pts], 
           pch = 20, cex = 1, col = paste0(color.old, "85"))
  }
  
  # reference line at R0 = 1 --------------------------------------------------
  abline(v = 1, lty = 2, col = "gray50", lwd = 2)
  
  # panel title ---------------------------------------------------------------
  legend("topright", legend = paste("Comparison:", host.nm), 
         cex = 1, bty = "n")
  
  # main legend ---------------------------------------------------------------
  legend("topleft", 
         legend = c("Weighted Sum (Eq. 3)", "Simple Sum (Eq. 4)", 
                    expression(R[0] == 1)),
         pch = c(0, 20, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
         col = c(color.new, color.old, "gray50"), 
         bty = "n", cex = 0.9, pt.cex = 1.2)
}

## -------- 7  generate two-panel figure ----------------------------------------
message("\nGenerating bifurcation plot …")

output.file <- file.path(out.dir, "bifurcation_R0_two_panel_overlay.pdf")
pdf(output.file, width = 6.8, height = 5, useDingbats = FALSE)

par(mfrow = c(2, 1),
    oma  = c(1.5, 1.5, 0, 0),
    mar  = c(0.5, 0.5, 0.5, 0.5),
    mgp  = c(2, 0.05, 0),
    tck  = 0.02)

# panel 1: large primate (human) --------------------------------------------
plot_overlay(R0.new, R0.old, bifH.new, bifH.old, "Large Primate")

# panel 2: small primate ----------------------------------------------------
plot_overlay(R0.new, R0.old, bifP.new, bifP.old, "Small Primate")

# shared axis labels --------------------------------------------------------
mtext(expression(R[0]), side = 1, outer = TRUE, line = 0.5, cex = 1)
mtext("Prevalence (per 1000)", side = 2, outer = TRUE, line = 0.5, cex = 1)

dev.off()

message("PDF written to: ", output.file)
message("\nDone!")
