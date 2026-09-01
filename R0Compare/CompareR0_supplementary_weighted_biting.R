#!/usr/bin/env Rscript

# Reproduce Supplementary Figure S1. The unweighted biting rates are rescaled
# point-by-point so the weighted and unweighted formulations deliver equal
# total bites; their R0 surfaces must therefore coincide.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

source("figure_helpers.R", chdir = FALSE)

message("Computing Supplementary Figure S1 equal-bite surfaces ...")
surfaces <- compute_surfaces(direction_multiplier = 1, equal_bites = TRUE)

differences <- c(
  max(abs(surfaces$weighted_hosts - surfaces$unweighted_hosts)),
  max(abs(surfaces$weighted_vectors - surfaces$unweighted_vectors)),
  max(abs(surfaces$weighted_transmission - surfaces$unweighted_transmission))
)
if (max(differences) > 1e-8) {
  stop("Equal-bite surfaces do not coincide; maximum difference = ",
       format(max(differences), scientific = TRUE))
}

output_file <- file.path("..", "outputs", "FigureS1_equal_bite_R0_comparison.pdf")
write_comparison_figure(surfaces, output_file)
message("Maximum weighted/unweighted difference: ",
        format(max(differences), scientific = TRUE))
message("Wrote ", normalizePath(output_file))
