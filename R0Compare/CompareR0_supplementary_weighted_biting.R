################################################################################
# CompareR0.R
#
# Produce Figure 3 (“CompareR0”) for the manuscript:
#   Althouse BM, *et al.* (2025) — FOI denominators and multi-host dengue.
#
# Requires:
#   • r0numerical2specIntro.r  (weighted FOI, Eq 5) — exported function r0new()
#   • r0numerical2specIntroOldDenom.r  (unweighted FOI, Eq 6) — exported function r0old()
#
# © 2025 Benjamin M. Althouse  •  CC BY 4.0
################################################################################

rm(list = ls())

# Set working directory
# For RStudio: uses active document location
# For command line: assumes script is run from its own directory
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# Get today's date
today <- format(Sys.Date(), "%Y-%m-%d")
output.dir <- paste0("output_", today, "_supp_weighted/")
# Check if the directory exists
if (!dir.exists(output.dir)) {
  # Create the directory if it doesn't exist
  dir.create(output.dir)
}


# ---- packages ----------------------------------------------------------------
library(fields)          # contour & colour support
library(RColorBrewer)
# library(rgl)  # Not needed for 2D contour plots

# ---- source analytic R0 functions --------------------------------------------
source('r0numerical2specIntro.r', chdir = F)
source('r0numerical2specIntroOldDenom.r', chdir = F)


# ---- helper: shaded “R0 ≈ 1” band -------------------------------------------
#' Shade the contour ribbon where R0 falls in the interval \[1 − ε, 1\]
#'
#' @param R0mat  numeric matrix of R0 values
#' @param epsilon half-width of the shading band (default 0.005)
#' @param div     divisor to rescale the x/y axes for plotting
#' @param flip    logical, reverse `y` ribbon for the second panel
#' @keywords internal
makeshade <- function(R0mat, epsilon = 0.005, div = 1000, flip = FALSE)
{
  idx <- which(R0mat < 1 & R0mat >= (1 - epsilon), arr.ind = TRUE)
  if (!length(idx)) return(invisible(NULL))
  
  xidx <- sort(idx[, 1])
  xx   <- c(x[xidx] / div, rev(x[xidx] / div))
  
  yy.idx <- idx[order(idx[, 1]), 2]
  mx <- if (flip) 0 else max(y) / div
  yy <- c(y[yy.idx] / div, rev(rep(mx, length(yy.idx))))
  
  polygon(xx, yy, col = grey(0), border = FALSE, density = 10)
}

# ---- helper: twin-panel contour with empirical ellipse -----------------------
#' Draw side-by-side contours of R0 (weighted vs. unweighted) and add an
#' empirical ellipse showing field uncertainty.
#'
#' @inheritParams makeshade
#' @param newR0,oldR0 matrices of R0 values
#' @param circ.x,circ.y ellipse centre (native axis units)
#' @param circ.rx,circ.ry horizontal / vertical semi-axes
#' @param tit   y-axis title (per two-panel row)
#' @param lab   two-element vector of subplot labels, e.g. c("a.","b.")
#' @export
shadedcontour <- function(newR0, oldR0,
                          epsilon1 = 0.005, epsilon2 = 0.5,
                          flip = FALSE, ax = NULL, div = 1000,
                          tit = NULL, CX.lab = 0.6, lab.line = 1,
                          lab = NULL, box.size = 0.85, cont.cx = 0.75,
                          circ.x = NULL, circ.y = NULL,
                          circ.rx = NULL, circ.ry = NULL,
                          circ.col = "red", circ.lwd = 2, circ.n = 60)
{
  # internal ellipse drawer
  ellipse <- function(cx, cy, rx, ry) {
    th <- seq(0, 2 * pi, length.out = circ.n)
    lines(cx + rx * cos(th), cy + ry * sin(th),
          col = circ.col, lwd = circ.lwd)
  }
  drawEllipse <- function() {
    if (is.null(circ.x) || is.null(circ.y)) return()
    if (is.null(circ.rx) || is.null(circ.ry)) {
      points(circ.x, circ.y, pch = 21, bg = NA,
             col = circ.col, lwd = circ.lwd)
    } else ellipse(circ.x, circ.y, circ.rx, circ.ry)
  }
  
  levels <- c(1, 2, 5)
  
  # ---- Panel 1: weighted FOI -------------------------------------------------
  contour(x / div, y / div, newR0, levels = levels,
          col = grey(0), labcex = cont.cx, axes = FALSE, lwd = 2)
  if (is.null(ax)) { axis(1, lwd = 0); axis(2, lwd = 0) } else {
    axis(1, at = ax, lwd = 0)
  }
  makeshade(newR0, epsilon1, div, flip)
  drawEllipse()
  
  if (!is.null(lab))
    text(par("usr")[2], par("usr")[4] * .93, lab[1],
         font = 2, cex = 1, pos = 2, xpd = NA)
  if (!is.null(tit))
    mtext(tit, side = 2, line = lab.line, cex = CX.lab)
  box()
  
  # ---- Panel 2: unweighted FOI ----------------------------------------------
  contour(x / div, y / div, oldR0, levels = levels,
          col = grey(0), labcex = cont.cx, axes = FALSE, lwd = 2)
  if (is.null(ax)) { axis(1, lwd = 0); axis(2, lwd = 0, labels = FALSE) } else {
    axis(1, at = ax, lwd = 0)
  }
  makeshade(oldR0, epsilon2, div, flip)
  drawEllipse()
  
  if (!is.null(lab))
    text(par("usr")[2], par("usr")[4] * .93, lab[2],
         font = 2, cex = 1, pos = 2, xpd = NA)
  box()
}




################################################################################
# Build the three R0 contour surfaces used in Figure 3
#
#   (i)   small- vs. large-host abundances          →  newR0hosts, oldR0hosts
#   (ii)  large-host abundance vs. vector abundance →  newR0hostvec, oldR0hostvec
#   (iii) host-to-vector transmission probabilities →  newR0trans,  oldR0trans
#
# All grids are 100 × 100 and are filled with the byte-compiled versions of
# r0new() and r0old() via outer() — 10× faster than explicit loops.
#
# © 2025 Benjamin M. Althouse • CC BY 4.0
################################################################################

library(compiler)

# -----------------------------------------------------------------------------#
# 1  Byte-compile the analytic R0 functions
# -----------------------------------------------------------------------------#
r0new <- cmpfun(r0new)
r0old <- cmpfun(r0old)


# -------------------------------------------------------------------
#  After:  source('r0numerical2specIntro.r')      #  r0new()
#          source('r0numerical2specIntroOldDenom.r')  #  r0old()
#
#  Patch:  create a wrapper that enforces "equal‐bite" normalization
# -------------------------------------------------------------------

equalise_bites <- function(f_old) {
  
  function(...) {
    # -------- pull arguments into a named list -------------
    args <- list(...)
    
    # -------- convenience vectors --------------------------
    Nh_vec <- c(args$Nh,  args$Np1)      # host abundances (extend if >2)
    # biting constants for vector 1 (primary = h, secondary = p1)
    r_v1   <- c(args$rm1h, args$rm1p1)
    # biting constants for vector 2
    r_v2   <- c(args$rm2h, args$rm2p1)
    
    # -------- helper to rescale one vector -----------------
    rescale <- function(r_vec) {
      # Reviewer #2's normalization: c_j = N_j^M / N_j^W
      # where N_j^M = total host abundance (unweighted)
      #       N_j^W = weighted denominator = sum(r_vj_hi * N_hi)
      c_j <- sum(Nh_vec) / sum(r_vec * Nh_vec)
      r_vec * c_j                                    # multiply by c_j
    }
    
    # enforce equal‐bite totals
    r_v1 <- rescale(r_v1)
    r_v2 <- rescale(r_v2)
    
    # write the scaled rates back into the arg list
    args$rm1h  <- r_v1[1];  args$rm1p1 <- r_v1[2]
    args$rm2h  <- r_v2[1];  args$rm2p1 <- r_v2[2]
    
    # -------- call the original r0old() with updated args --
    do.call(f_old, args)
  }
}

# Build the "equal‑bite" version of the un‑weighted FOI
r0old.eq <- equalise_bites(r0old)
# NOTE: r0new already has the correct weighted denominator, no normalization needed



# -----------------------------------------------------------------------------#
# 2  Global constants that never change across any surface
# -----------------------------------------------------------------------------#
BB     <- 0.15
cross  <- 0.0005

core.args <- list(
  bm1h = BB, bhm1 = BB, bm1p1 = BB,
  bm2h = BB, bhm2 = BB, bm2p1 = BB,
  bp1m2 = BB, bp1m1 = BB,
  bm2p2 = BB, bp2m1 = BB, bm1p2 = BB, bp2m2 = BB,
  rm1h = 0.5, rm1p1 = cross, rm2h = cross,
  rm2p1 = 0.5, rm2p2 = 0.5, rm1p2 = 0.5,
  Vm  = 1/7, Vm2 = 1/7,
  gh  = 1/4,  mh  = 1/(60 * 365),
  gp1 = 1/4,  mp1 = 1/(15 * 365),
  mp2 = 1/(30 * 365), gp2 = 1/5,
  Np2 = 0, Nm2 = 25000
)

# -----------------------------------------------------------------------------#
# 3  Surface (i):  small- vs. large-host abundances
# -----------------------------------------------------------------------------#
xs.host <- seq(1, 50 * 100, length.out = 100)   # Nh
ys.host <- xs.host                              # Np1

f_new_host <- function(Nh, Np1)
  do.call(r0new, modifyList(core.args,
                            list(Nh = Nh, Np1 = Np1, Nm1 = 25000, NN = Nh + Np1)))

f_old_host <- function(Nh, Np1)
  do.call(r0old.eq, modifyList(core.args,
                            list(Nh = Nh, Np1 = Np1, Nm1 = 25000, NN = Nh + Np1)))

newR0hosts <- outer(xs.host, ys.host, Vectorize(f_new_host))
oldR0hosts <- outer(xs.host, ys.host, Vectorize(f_old_host))

# # NEW – point f_old_host at r0old.eq instead of r0old
# f_old_host <- function(Nh, Np1)
#   do.call(r0old.eq, modifyList(core.args,
#                                list(Nh = Nh, Np1 = Np1,
#                                     Nm1 = 25000, NN = Nh + Np1)))
# oldR0hosts <- outer(xs.host, ys.host, Vectorize(f_old_host))
# 
# f_new_host <- function(Nh, Np1)
#   do.call(r0new.eq, modifyList(core.args,
#                                list(Nh = Nh, Np1 = Np1,
#                                     Nm1 = 25000, NN = Nh + Np1)))
# newR0hosts <- outer(xs.host, ys.host, Vectorize(f_new_host))


# -----------------------------------------------------------------------------#
# 4  Surface (ii):  large-host abundance vs. vector abundance
# -----------------------------------------------------------------------------#
xs.hv <- xs.host                               # Nh: 1–5 ×10^3
ys.hv <- seq(1, 500 * 100, length.out = 100)  # Nm1: 1–5 ×10^4

f_new_hv <- function(Nh, Nm1)
  do.call(r0new, modifyList(core.args,
                            list(Nh = Nh, Np1 = 1000, Nm1 = Nm1, NN = Nh)))

f_old_hv <- function(Nh, Nm1)
  do.call(r0old.eq, modifyList(core.args,
                            list(Nh = Nh, Np1 = 1000, Nm1 = Nm1, NN = Nh)))

newR0hostvec <- outer(xs.hv, ys.hv, Vectorize(f_new_hv))
oldR0hostvec <- outer(xs.hv, ys.hv, Vectorize(f_old_hv))



# -----------------------------------------------------------------------------#
# 5  Surface (iii):  transmission-probability grid
# -----------------------------------------------------------------------------#
xs.tr <- seq(0, 0.5, length.out = 100)   # β_small
ys.tr <- xs.tr                           # β_large

trans.args <- modifyList(core.args,
                         list(rm1p1 = 0.001, rm2h = 0.001,
                              Nh = 1000, Np1 = 1000, Nm1 = 25000, NN = 2000))


f_new_tr <- function(betaSmall, betaLarge)
  do.call(r0new, modifyList(trans.args,
                            list(bm1h = betaSmall, bhm1 = betaSmall,
                                 bm2h = betaSmall, bhm2 = betaSmall,
                                 bm1p1 = betaLarge, bm2p1 = betaLarge,
                                 bp1m1 = betaLarge, bp1m2 = betaLarge)))

f_old_tr <- function(betaSmall, betaLarge)
  do.call(r0old.eq, modifyList(trans.args,
                            list(bm1h = betaSmall, bhm1 = betaSmall,
                                 bm2h = betaSmall, bhm2 = betaSmall,
                                 bm1p1 = betaLarge, bm2p1 = betaLarge,
                                 bp1m1 = betaLarge, bp1m2 = betaLarge)))

newR0trans <- outer(xs.tr, ys.tr, Vectorize(f_new_tr))
oldR0trans <- outer(xs.tr, ys.tr, Vectorize(f_old_tr))

# -----------------------------------------------------------------------------#
# Done — the six matrices are ready for shadedcontour()
#   newR0hosts, oldR0hosts,
#   newR0hostvec, oldR0hostvec,
#   newR0trans,  oldR0trans
# -----------------------------------------------------------------------------#





###############################################################################
## Figure 3: side-by-side comparison of weighted vs. unweighted FOI contours
##
## Produces a 4 × 6-inch PDF with six panels:
##   a,b – small- vs. large-host numbers
##   c,d – large-host vs. vector numbers
##   e,f – host-to-vector transmission probabilities
##
## Assumes that the six matrices
##   newR0hosts, oldR0hosts,
##   newR0hostvec, oldR0hostvec,
##   newR0trans,  oldR0trans
## are already in the workspace, plus the helper `shadedcontour()`.
###############################################################################

output.file <- file.path(output.dir, "Combined-WithData.pdf")
pdf(output.file, width = 4, height = 6)

par(mfrow = c(3, 2),
    mar  = c(2, 1, 1, 0),
    oma  = c(1, 1.5, 0, 2),
    mgp  = c(2, 0.15, 0),
    tck  = 0,
    cex.axis = 0.85)

CX  <- 0.7    # axis-title size
lin <- 1.5    # axis-title offset

## -------------------------------------------------------------------------
## Panel a / b  –  small- vs. large-host abundances
## -------------------------------------------------------------------------
x<-seq(1,50*100,length.out=100)
y<-seq(1,50*100,length.out=100)

## Panel a/b: small vs large hosts
shadedcontour(newR0hosts, oldR0hosts,
              epsilon1 = .05, epsilon2 = .25,
              tit = "Number Small Host (1000s)",
              CX.lab = CX, lab.line = lin,
              box.size = .875, lab = c("a.","b."), cont.cx = .65,
              circ.x = 2.5, circ.y = 1.25,
              circ.rx = 1.25, circ.ry = 0.3)   # oval spans ±1.25k × ±0.3k


mtext("Number Large Host (1000s)", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)

## -------------------------------------------------------------------------
## Panel c / d  –  large-host abundance vs. vector abundance
## -------------------------------------------------------------------------
x <- seq(1, 50*100, length.out = 100)
y <- seq(1, 500*100, length.out = 100)


shadedcontour(newR0hostvec, oldR0hostvec,
              epsilon1 = .005, epsilon2 = .25, flip = TRUE,
              tit = "Number Vectors (1000s)",
              CX.lab = CX, lab.line = lin,
              lab = c("c.","d."), cont.cx = .65,
              circ.x = 1.5, circ.y = 25,
              circ.rx = 0.3, circ.ry = 9)

mtext("Number Large Host (1000s)", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)


## -------------------------------------------------------------------------
## Panel e / f  –  transmission-probability grid
## -------------------------------------------------------------------------
x <- seq(0, .5, length.out = 100)
y <- seq(0, .5, length.out = 100)


shadedcontour(newR0trans, oldR0trans,
              epsilon1 = .25, epsilon2 = .25, flip = TRUE, div = 1,
              tit = "Small Host Transmission Probability",
              CX.lab = CX, lab.line = lin,
              lab = c("e.","f."), cont.cx = .65,
              circ.x = 0.12, circ.y = 0.05,
              circ.rx = 0.03, circ.ry = 0.02)


mtext("Large Host Transmission Probability", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)

dev.off()
message("Figure written to: ", output.file)






# ----------- add 1-1 line
# ---- helper: twin‑panel contour with empirical ellipse ----------------------
#' Draw side‑by‑side contours of R0 (weighted vs. unweighted) and add an
#' empirical ellipse showing field uncertainty.
#'
#' @inheritParams makeshade
#' @param newR0,oldR0 matrices of R0 values
#' @param circ.x,circ.y ellipse centre (native axis units)
#' @param circ.rx,circ.ry horizontal / vertical semi‑axes
#' @param tit   y‑axis title (per two‑panel row)
#' @param lab   two‑element vector of subplot labels, e.g. c("a.","b.")
#' @param refline logical; draw a 1–1 line in both panels?
#' @param ref.col,ref.lwd color and width for the 1–1 line
#' @export
shadedcontour <- function(newR0, oldR0,
                          epsilon1 = 0.005, epsilon2 = 0.5,
                          flip = FALSE, ax = NULL, div = 1000,
                          tit = NULL, CX.lab = 0.6, lab.line = 1,
                          lab = NULL, box.size = 0.85, cont.cx = 0.75,
                          circ.x = NULL, circ.y = NULL,
                          circ.rx = NULL, circ.ry = NULL,
                          circ.col = "red", circ.lwd = 2, circ.n = 60,
                          refline = FALSE, ref.col = "blue", ref.lwd = 3)
{
  # internal ellipse drawer ----------------------------------------------------
  ellipse <- function(cx, cy, rx, ry) {
    th <- seq(0, 2 * pi, length.out = circ.n)
    lines(cx + rx * cos(th), cy + ry * sin(th),
          col = circ.col, lwd = circ.lwd)
  }
  drawEllipse <- function() {
    if (is.null(circ.x) || is.null(circ.y)) return()
    if (is.null(circ.rx) || is.null(circ.ry)) {
      points(circ.x, circ.y, pch = 21, bg = NA,
             col = circ.col, lwd = circ.lwd)
    } else ellipse(circ.x, circ.y, circ.rx, circ.ry)
  }
  
  levels <- c(1, 2, 5)
  
  # ---- Panel 1: weighted FOI -------------------------------------------------
  contour(x / div, y / div, newR0, levels = levels,
          col = grey(0), labcex = cont.cx, axes = FALSE, lwd = 2)
  if (is.null(ax)) { axis(1, lwd = 0); axis(2, lwd = 0) } else {
    axis(1, at = ax, lwd = 0)
  }
  makeshade(newR0, epsilon1, div, flip)
  if (refline) abline(0, 1, col = ref.col, lwd = ref.lwd)
  drawEllipse()
  
  if (!is.null(lab))
    text(par("usr")[2], par("usr")[4] * .93, lab[1],
         font = 2, cex = 1, pos = 2, xpd = NA)
  if (!is.null(tit))
    mtext(tit, side = 2, line = lab.line, cex = CX.lab)
  box()
  
  # ---- Panel 2: unweighted FOI ----------------------------------------------
  contour(x / div, y / div, oldR0, levels = levels,
          col = grey(0), labcex = cont.cx, axes = FALSE, lwd = 2)
  if (is.null(ax)) { axis(1, lwd = 0); axis(2, lwd = 0, labels = FALSE) } else {
    axis(1, at = ax, lwd = 0)
  }
  makeshade(oldR0, epsilon2, div, flip)
  if (refline) abline(0, 1, col = ref.col, lwd = ref.lwd)
  drawEllipse()
  
  if (!is.null(lab))
    text(par("usr")[2], par("usr")[4] * .93, lab[2],
         font = 2, cex = 1, pos = 2, xpd = NA)
  box()
}



###############################################################################
## Figure 3: side-by-side comparison of weighted vs. unweighted FOI contours
##
## Produces a 4 × 6-inch PDF with six panels:
##   a,b – small- vs. large-host numbers
##   c,d – large-host vs. vector numbers
##   e,f – host-to-vector transmission probabilities
##
## Assumes that the six matrices
##   newR0hosts, oldR0hosts,
##   newR0hostvec, oldR0hostvec,
##   newR0trans,  oldR0trans
## are already in the workspace, plus the helper `shadedcontour()`.
###############################################################################

output.file <- file.path(output.dir, "Combined-WithData-1-1.pdf")
pdf(output.file, width = 4, height = 6)

par(mfrow = c(3, 2),
    mar  = c(2, 1, 1, 0),
    oma  = c(1, 1.5, 0, 2),
    mgp  = c(2, 0.15, 0),
    tck  = 0,
    cex.axis = 0.85)

CX  <- 0.7    # axis-title size
lin <- 1.5    # axis-title offset

## -------------------------------------------------------------------------
## Panel a / b  –  small‑ vs. large‑host abundances
## -------------------------------------------------------------------------
x <- seq(1, 50 * 100, length.out = 100)
y <- seq(1, 50 * 100, length.out = 100)

shadedcontour(newR0hosts, oldR0hosts,
              epsilon1 = .05, epsilon2 = .25,
              tit = "Number Small Host (1000s)",
              CX.lab = CX, lab.line = lin,
              box.size = .875, lab = c("a.","b."), cont.cx = .65,
              circ.x = 2.5, circ.y = 1.25,
              circ.rx = 1.25, circ.ry = 0.3,
              refline = TRUE)                     # <‑‑ 1–1 line ON

mtext("Number Large Host (1000s)", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)

## -------------------------------------------------------------------------
## Panel c / d  –  large‑host abundance vs. vector abundance
## -------------------------------------------------------------------------
x <- seq(1, 50 * 100,  length.out = 100)
y <- seq(1, 500 * 100, length.out = 100)

shadedcontour(newR0hostvec, oldR0hostvec,
              epsilon1 = .005, epsilon2 = .25, flip = TRUE,
              tit = "Number Vectors (1000s)",
              CX.lab = CX, lab.line = lin,
              lab = c("c.","d."), cont.cx = .65,
              circ.x = 1.5, circ.y = 25,
              circ.rx = 0.3, circ.ry = 9,
              refline = TRUE)                     # <‑‑ 1–1 line ON

mtext("Number Large Host (1000s)", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)

## -------------------------------------------------------------------------
## Panel e / f  –  transmission‑probability grid
## -------------------------------------------------------------------------
x <- seq(0, .5, length.out = 100)
y <- seq(0, .5, length.out = 100)

shadedcontour(newR0trans, oldR0trans,
              epsilon1 = .25, epsilon2 = .25, flip = TRUE, div = 1,
              tit = "Small Host Transmission Probability",
              CX.lab = CX, lab.line = lin,
              lab = c("e.","f."), cont.cx = .65,
              circ.x = 0.12, circ.y = 0.05,
              circ.rx = 0.03, circ.ry = 0.02,
              refline = FALSE)                    # <‑‑ keep OFF here



mtext("Large Host Transmission Probability", side = 1,
      line = lin, cex = CX, outer = FALSE, at = 0)

dev.off()
message("Figure written to: ", output.file)







