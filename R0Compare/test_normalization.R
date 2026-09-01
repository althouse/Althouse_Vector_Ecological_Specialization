#!/usr/bin/env Rscript

# Numerical regression test for the equal-bite identity.

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg)) {
  setwd(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
}

source("figure_helpers.R", chdir = FALSE)

scenarios <- list(
  equal = c(Nh = 1000, Np1 = 1000),
  unequal = c(Nh = 2500, Np1 = 500),
  highly_asymmetric = c(Nh = 4000, Np1 = 200)
)

differences <- vapply(names(scenarios), function(name) {
  populations <- scenarios[[name]]
  model_args <- modifyList(
    base_arguments(),
    list(Nh = populations[["Nh"]], Np1 = populations[["Np1"]],
         NN = sum(populations))
  )
  weighted <- evaluate_r0(r0new, model_args)
  rescaled_unweighted <- evaluate_r0(r0old, equal_bite_arguments(model_args))
  difference <- abs(weighted - rescaled_unweighted)
  message(sprintf("%-19s weighted=%0.12f  rescaled-unweighted=%0.12f  diff=%0.3e",
                  name, weighted, rescaled_unweighted, difference))
  difference
}, numeric(1))

tolerance <- 1e-10
if (max(differences) > tolerance) {
  stop("Equal-bite normalization regression failed; maximum difference = ",
       format(max(differences), scientific = TRUE))
}
message("PASS: equal-bite identity holds in all scenarios.")
