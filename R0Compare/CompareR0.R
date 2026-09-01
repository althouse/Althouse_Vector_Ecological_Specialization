#!/usr/bin/env Rscript

# Reproduce main-text Figure 3: R0 landscapes under the weighted and
# unweighted force-of-infection formulations.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

source("figure_helpers.R", chdir = FALSE)

message("Computing main-text Figure 3 surfaces ...")
surfaces <- compute_surfaces(direction_multiplier = 1, equal_bites = FALSE)
output_file <- file.path("..", "outputs", "Figure3_R0_comparison.pdf")
write_comparison_figure(surfaces, output_file)
message("Wrote ", normalizePath(output_file))
