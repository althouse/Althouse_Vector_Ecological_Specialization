#!/usr/bin/env Rscript

# Reproduce Supplementary Figure S2 with asymmetric transmission:
# host-to-vector transmission is one quarter of vector-to-host transmission.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

source("figure_helpers.R", chdir = FALSE)

message("Computing Supplementary Figure S2 surfaces ...")
surfaces <- compute_surfaces(direction_multiplier = 0.25, equal_bites = FALSE)
output_file <- file.path("..", "outputs", "FigureS2_asymmetric_transmission.pdf")
write_comparison_figure(surfaces, output_file)
message("Wrote ", normalizePath(output_file))
