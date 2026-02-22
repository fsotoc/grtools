# Individual SvA curve helper used by sva_curve.grt_wind_fit
sva_curve_ind <- function(fitted_model, bootstrap_samples, sub, curve = c("SvA", "SvM"),
                          stim_level = c(1, 2), alpha = 0.05, log_scale = FALSE,
                          range_fact = 3) {
  curve <- match.arg(curve)
  if (any(fitted_model$restrictions == "r3x2")) {
    model_type <- 2
  } else {
    model_type <- 1
  }

  ymin <- min(fitted_model$means[, 2]) - range_fact
  ymax <- max(fitted_model$means[, 2]) + range_fact
  y_steps <- seq(ymax, ymin, by = -0.01)
  n_steps <- length(y_steps)

  rlna_ind_vector <- rep(NA_real_, n_steps)
  d_vector <- rep(NA_real_, n_steps)
  dlow <- rep(NA_real_, n_steps)
  dhigh <- rep(NA_real_, n_steps)

  for (i in seq_len(n_steps)) {
    if (curve == "SvA") {
      rlna_ind_vector[i] <- sva_rlna_ind(y_steps[i], fitted_model$fullpars, sub, model_type, stim_level = stim_level)
    } else {
      rlna_ind_vector[i] <- sva_rlm_ind(y_steps[i], fitted_model$fullpars, sub, model_type, stim_level = stim_level)
    }
    d_vector[i] <- sva_d_prime_cond_ind(y_steps[i], fitted_model$fullpars, sub, model_type)

    rlna_ind_edf <- rep(NA_real_, nrow(bootstrap_samples$samples))
    for (s in seq_len(nrow(bootstrap_samples$samples))) {
      rlna_ind_edf[s] <- sva_d_prime_cond_ind(y_steps[i], bootstrap_samples$samples[s, ], sub, model_type)
    }
    ci <- stats::quantile(rlna_ind_edf, c(alpha / 2, 1 - alpha / 2))
    dlow[i] <- ci[1]
    dhigh[i] <- ci[2]
  }

  dmin <- min(dlow)
  dmax <- max(dhigh)
  if (dmin > 0) {
    dmin <- 0
  }

  oc <- 1
  x_label <- "Relative Likelihood of No Awareness"
  ind_c <- sva_rlna_ind(-fitted_model$indpar$a2[sub], fitted_model$fullpars, sub, model_type, stim_level = stim_level)

  if (curve == "SvM") {
    x_label <- "Relative Likelihood of Metacognition"
    rlna_ind_vector <- sva_safe_log2(rlna_ind_vector, log_scale)
    ind_c <- sva_safe_log2(sva_rlm_ind(-fitted_model$indpar$a2[sub], fitted_model$fullpars, sub, model_type, stim_level = stim_level), log_scale)
    if (isTRUE(log_scale)) {
      x_label <- "Relative Likelihood of Metacognition (Log2 Scale)"
      oc <- 0
    }
  }

  graphics::plot(
    rlna_ind_vector,
    d_vector,
    type = "l",
    xlab = x_label,
    ylab = "Sensitivity (d')",
    main = paste("Participant", sub),
    ylim = c(dmin, dmax),
    lwd = 2.5,
    cex.lab = 1.5,
    col = grDevices::rgb(0.8, 0.25, 0.33, 1)
  )
  graphics::polygon(
    c(rlna_ind_vector, rev(rlna_ind_vector)),
    c(dlow, rev(dhigh)),
    col = grDevices::rgb(1, 0, 0, 0.3),
    border = NA
  )
  graphics::abline(v = if (curve == "SvM") oc else 1, col = grDevices::rgb(0, 0, 0.8, 1), lty = 2, lwd = 1.5)
  graphics::abline(
    v = ind_c,
    col = grDevices::rgb(0, 0.6, 0, 1),
    lty = 3
  )
  graphics::abline(h = 0, lty = 2, lwd = 1.5)

  invisible(NULL)
}
