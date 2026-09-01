# Shared helpers for Figures 3, S1, and S2.

source("r0numerical2specIntro.r", chdir = FALSE)
source("r0numerical2specIntroOldDenom.r", chdir = FALSE)

r0new <- compiler::cmpfun(r0new)
r0old <- compiler::cmpfun(r0old)

base_arguments <- function(direction_multiplier = 1, cross_biting = 0.001) {
  beta <- 0.15
  list(
    bm1h = beta, bhm1 = beta * direction_multiplier,
    bm1p1 = beta, bp1m1 = beta * direction_multiplier,
    bm2h = beta, bhm2 = beta * direction_multiplier,
    bm2p1 = beta, bp1m2 = beta * direction_multiplier,
    bm2p2 = beta, bp2m1 = beta,
    bm1p2 = beta, bp2m2 = beta,
    rm1h = 0.5, rm1p1 = cross_biting,
    rm2h = cross_biting, rm2p1 = 0.5,
    rm2p2 = 0.5, rm1p2 = 0.5,
    Vm = 1 / 7, Vm2 = 1 / 7,
    gh = 1 / 4, mh = 1 / (60 * 365),
    gp1 = 1 / 4, mp1 = 1 / (15 * 365),
    gp2 = 1 / 5, mp2 = 1 / (30 * 365),
    Nh = 1000, Np1 = 1000, Np2 = 0,
    Nm1 = 25000, Nm2 = 25000,
    NN = 2000, intro = 0
  )
}

# Rescale the unweighted formulation so that each vector species delivers the
# same total number of bites as in the weighted formulation:
# c_j = Nbar / D_j, D_j = sum_i[(r_ji / sum_i r_ji) N_i].
equal_bite_arguments <- function(args) {
  host_abundance <- c(args$Nh, args$Np1)

  rescale <- function(rates) {
    weighted_denominator <- sum((rates / sum(rates)) * host_abundance)
    rates * (sum(host_abundance) / weighted_denominator)
  }

  vector_1 <- rescale(c(args$rm1h, args$rm1p1))
  vector_2 <- rescale(c(args$rm2h, args$rm2p1))
  args$rm1h <- vector_1[1]
  args$rm1p1 <- vector_1[2]
  args$rm2h <- vector_2[1]
  args$rm2p1 <- vector_2[2]
  args
}

evaluate_r0 <- function(fun, args) do.call(fun, args)

compute_surfaces <- function(direction_multiplier = 1,
                             equal_bites = FALSE,
                             grid_points = 100,
                             cross_biting = 0.001) {
  core <- base_arguments(direction_multiplier, cross_biting)

  old_value <- function(args) {
    if (equal_bites) args <- equal_bite_arguments(args)
    evaluate_r0(r0old, args)
  }

  host_axis <- seq(1, 5000, length.out = grid_points)
  vector_axis <- seq(1, 50000, length.out = grid_points)
  transmission_axis <- seq(0, 0.5, length.out = grid_points)

  weighted_hosts <- outer(host_axis, host_axis, Vectorize(function(Nh, Np1) {
    args <- modifyList(core, list(Nh = Nh, Np1 = Np1, NN = Nh + Np1))
    evaluate_r0(r0new, args)
  }))
  unweighted_hosts <- outer(host_axis, host_axis, Vectorize(function(Nh, Np1) {
    args <- modifyList(core, list(Nh = Nh, Np1 = Np1, NN = Nh + Np1))
    old_value(args)
  }))

  weighted_vectors <- outer(host_axis, vector_axis, Vectorize(function(Nh, Nm1) {
    args <- modifyList(core, list(Nh = Nh, Np1 = 1000, Nm1 = Nm1,
                                  NN = Nh + 1000))
    evaluate_r0(r0new, args)
  }))
  unweighted_vectors <- outer(host_axis, vector_axis, Vectorize(function(Nh, Nm1) {
    args <- modifyList(core, list(Nh = Nh, Np1 = 1000, Nm1 = Nm1,
                                  NN = Nh + 1000))
    old_value(args)
  }))

  transmission_args <- modifyList(
    core,
    list(Nh = 1000, Np1 = 1000, NN = 2000, Nm1 = 25000, Nm2 = 25000)
  )
  transmission_arguments <- function(beta_large, beta_small) {
    modifyList(
      transmission_args,
      list(
        bm1h = beta_large, bm2h = beta_large,
        bhm1 = beta_large * direction_multiplier,
        bhm2 = beta_large * direction_multiplier,
        bm1p1 = beta_small, bm2p1 = beta_small,
        bp1m1 = beta_small * direction_multiplier,
        bp1m2 = beta_small * direction_multiplier
      )
    )
  }
  weighted_transmission <- outer(
    transmission_axis, transmission_axis,
    Vectorize(function(beta_large, beta_small) {
      evaluate_r0(r0new, transmission_arguments(beta_large, beta_small))
    })
  )
  unweighted_transmission <- outer(
    transmission_axis, transmission_axis,
    Vectorize(function(beta_large, beta_small) {
      old_value(transmission_arguments(beta_large, beta_small))
    })
  )

  matrices <- list(
    weighted_hosts, unweighted_hosts,
    weighted_vectors, unweighted_vectors,
    weighted_transmission, unweighted_transmission
  )
  stopifnot(all(vapply(matrices, function(x) all(is.finite(x)), logical(1))))

  list(
    host_axis = host_axis,
    vector_axis = vector_axis,
    transmission_axis = transmission_axis,
    weighted_hosts = weighted_hosts,
    unweighted_hosts = unweighted_hosts,
    weighted_vectors = weighted_vectors,
    unweighted_vectors = unweighted_vectors,
    weighted_transmission = weighted_transmission,
    unweighted_transmission = unweighted_transmission
  )
}

shade_subthreshold <- function(x, y, z, epsilon, divisor, flip) {
  index <- which(z < 1 & z >= 1 - epsilon, arr.ind = TRUE)
  if (!nrow(index)) return(invisible(NULL))

  groups <- split(index[, 2], index[, 1])
  x_index <- as.integer(names(groups))
  y_index <- vapply(groups, function(values) as.integer(round(mean(values))),
                    integer(1))
  xx <- x[x_index] / divisor
  yy <- y[y_index] / divisor
  edge <- if (flip) min(y) / divisor else max(y) / divisor
  polygon(c(xx, rev(xx)), c(yy, rep(edge, length(yy))),
          col = "black", border = FALSE, density = 10)
}

draw_ellipse <- function(center_x, center_y, radius_x, radius_y) {
  theta <- seq(0, 2 * pi, length.out = 100)
  lines(center_x + radius_x * cos(theta), center_y + radius_y * sin(theta),
        col = "red", lwd = 2)
}

draw_pair <- function(x, y, weighted, unweighted, divisor,
                      epsilon_weighted, epsilon_unweighted, flip,
                      y_title, labels, ellipse, contour_size = 0.65) {
  draw_panel <- function(z, epsilon, show_y, label) {
    contour(x / divisor, y / divisor, z, levels = c(1, 2, 5),
            col = "black", labcex = contour_size, axes = FALSE, lwd = 2)
    axis(1, lwd = 0)
    axis(2, lwd = 0, labels = show_y)
    shade_subthreshold(x, y, z, epsilon, divisor, flip)
    do.call(draw_ellipse, ellipse)
    text(par("usr")[2], par("usr")[4] * 0.93, label,
         font = 2, cex = 1, pos = 2, xpd = NA)
    box()
  }

  draw_panel(weighted, epsilon_weighted, TRUE, labels[1])
  mtext(y_title, side = 2, line = 1.5, cex = 0.7)
  draw_panel(unweighted, epsilon_unweighted, FALSE, labels[2])
}

write_comparison_figure <- function(surfaces, output_file) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  pdf(output_file, width = 4, height = 6, useDingbats = FALSE)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(3, 2), mar = c(2, 1, 1, 0), oma = c(1, 1.5, 0, 2),
      mgp = c(2, 0.15, 0), tck = 0, cex.axis = 0.85)

  draw_pair(
    surfaces$host_axis, surfaces$host_axis,
    surfaces$weighted_hosts, surfaces$unweighted_hosts,
    divisor = 1000, epsilon_weighted = 0.05, epsilon_unweighted = 0.25,
    flip = FALSE, y_title = "Number Small Host (1000s)", labels = c("a.", "b."),
    ellipse = list(center_x = 1.5, center_y = 2.5,
                   radius_x = 0.3, radius_y = 1.25)
  )
  mtext("Number Large Host (1000s)", side = 1, line = 1.5,
        cex = 0.7, outer = FALSE, at = 0)

  draw_pair(
    surfaces$host_axis, surfaces$vector_axis,
    surfaces$weighted_vectors, surfaces$unweighted_vectors,
    divisor = 1000, epsilon_weighted = 0.005, epsilon_unweighted = 0.25,
    flip = TRUE, y_title = "Number Vectors (1000s)", labels = c("c.", "d."),
    ellipse = list(center_x = 1.5, center_y = 25,
                   radius_x = 0.3, radius_y = 9)
  )
  mtext("Number Large Host (1000s)", side = 1, line = 1.5,
        cex = 0.7, outer = FALSE, at = 0)

  draw_pair(
    surfaces$transmission_axis, surfaces$transmission_axis,
    surfaces$weighted_transmission, surfaces$unweighted_transmission,
    divisor = 1, epsilon_weighted = 0.25, epsilon_unweighted = 0.25,
    flip = TRUE, y_title = "Small Host Transmission Probability",
    labels = c("e.", "f."),
    ellipse = list(center_x = 0.05, center_y = 0.12,
                   radius_x = 0.02, radius_y = 0.03)
  )
  mtext("Large Host Transmission Probability", side = 1, line = 1.5,
        cex = 0.7, outer = FALSE, at = 0)

  invisible(output_file)
}
