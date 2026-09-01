#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_dir <- if (length(file_arg)) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  normalizePath(getwd())
}

scripts <- c(
  "Bifurcation/compute_bifurcation_coefficients.R",
  "R0Compare/test_normalization.R",
  "Bifurcation/Bifurcation_publication.R",
  "R0Compare/CompareR0.R",
  "R0Compare/CompareR0_supplementary_weighted_biting.R",
  "R0Compare/CompareR0_supplementary.R"
)

rscript <- file.path(R.home("bin"), "Rscript")
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(repo_dir)
for (script in scripts) {
  cat("\n=== Running", script, "===\n")
  status <- system2(rscript, script)
  if (status != 0) stop("Failed: ", script, call. = FALSE)
}

cat("\nAll analyses completed successfully. Outputs are in outputs/.\n")
