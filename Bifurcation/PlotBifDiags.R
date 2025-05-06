
################################################################################
# compare_bifurcation.R
# 
# Produce Figure 2 (“CompareBifurcation”) for the manuscript:
#   Althouse BM, *et al.* (2025) — FOI denominators and multi-host dengue.
# 
# 
#   Simulation results from:
#   Althouse BM, Lessler J, Sall AA, Diallo M, Hanley KA, et al. (2012)
#   *Synchrony of Sylvatic Dengue Isolations: A Multi-Host, Multi-Vector
#   SIR Model of Dengue Virus Transmission in Senegal.
#   PLOS Neglected Tropical Diseases 6(11): e1928.
#   10.1371/journal.pntd.0001928
# 
# Build the human-prevalence bifurcation diagram that contrasts the weighted
# (Eq. 3) and unweighted (Eq. 4) FOI formulations.  The script:
#
#   1. Creates a date-stamped output folder.
#   2. Reads every quarterly point from 100 deterministic simulations stored as
#        NewForm-FixSmallPrimate/{Human,Primate}/H-bb<i>.csv
#        OldForm-FixSmallPrimate/{Human,Primate}/H-bb<i>.csv
#   3. Generates a 6.8 × 3-inch PDF (CompareBifurcation.pdf) that overlays the
#      two bifurcation branches.
#
# Source this file or run it with `Rscript compare_bifurcation.R`.
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
if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)

## -------- 2  parameters & helpers -------------------------------------------
n.sim   <- 100                       # number of beta lines (1…100)
disc    <- 50                        # matches original beta step
beta.x  <- ((disc / n.sim) * (1:n.sim)) / 100 + .025 - disc / n.sim / 100

quarter.idx <- seq(1, 365 * 50, by = 365)  # every 365 days ≈ “3 months” in code

read_matrix <- function(folder) {
  # returns n.sim × length(idx) matrix for Human prevalence
  file.vec <- sprintf("%s/Human/H-bb%d.csv", folder, 1:n.sim)
  vapply(file.vec, function(f) {
    scan(f, what = numeric(), quiet = TRUE)[quarter.idx]
  }, numeric(length(quarter.idx)))
}

## -------- 3  load data -------------------------------------------------------
new.dir  <- "NewForm-FixSmallPrimate"
old.dir  <- "OldForm-FixSmallPrimate"

message("Reading CSVs …")
bifH.new  <- read_matrix(new.dir)
bifH.old  <- read_matrix(old.dir)

## -------- 4  plot ------------------------------------------------------------
color1 <- grey(0.65)
color2 <- grey(0.15)
col.new  <- paste(color1, "25", sep="")   # light grey, 85 % opacity
col.old  <- paste(color2, "25", sep="")


yrange   <- c(0, 0.0015)

pdf(file.path(out.dir, "CompareBifurcation.pdf"),
    width = 6.8, height = 3, useDingbats = FALSE)

par(mar = c(1, 1, 1, 1), oma = c(2, 2, 0, 0),
    mgp = c(2, 0.25, 0),  tck = -0.02,
    cex.axis = 0.75)

plot(beta.x[1], bifH.new[1, 1],                    # empty frame via first point
     type = "n", xlim = range(beta.x), ylim = yrange,
     xlab = "", ylab = "")

# draw weighted (squares) then unweighted (dots) ----------------------------
for (i in seq_len(n.sim)) {
  points(rep(beta.x[i], length(unique(bifH.new[, i]))),
         unique(bifH.new[, i]),
         pch = 0, cex = 1, col = col.new)
  points(rep(beta.x[i], length(unique(bifH.old[, i]))),
         unique(bifH.old[, i]),
         pch = 20, cex = 0.9, col = col.old)
}

mtext(expression(beta), side = 1, line = 1.5)
mtext("Prevalence",      side = 2, line = 1.5, cex = 0.85)

legend("topleft",
       legend = c("Weighted Sum (Eq. 3)", "Sum (Eq. 4)"),
       pch = c(0, 20), col = c(color1, color2),
       bty = "n", cex = 0.75)

dev.off()
message("PDF written to: ", file.path(out.dir, "CompareBifurcation.pdf"))
