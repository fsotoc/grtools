#' Plot Sensitivity Versus Awareness (SvA) Curves from a Fitted GRTapas Model
#'
#' Plots a sensitivity-versus-awareness (SvA) curve from a fitted GRTapas model,
#' including a bootstrap confidence band for conditional sensitivity.
#'
#' @param x A fitted model returned by \code{grt_wind_fit} or
#'   \code{grt_wind_fit_parallel} under GRTapas restrictions.
#' @param bootstrap_samples An object returned by \code{bootstrap_sample} for
#'   the same fitted model.
#' @param sub A flag indicating whether to plot the group-level SvA curve or an
#'   individual curve. Use \code{"group"} (default) for the group model. If a
#'   number is provided, the function plots the SvA curve for that participant
#'   using \code{sva_curve_ind}.
#' @param curve Type of curve to plot. Defaults to \code{"SvA"}. Use
#'   \code{"SvM"} to follow the sensitivity-versus-metacognition logic from
#'   \code{svm_curve}.
#' @param stim_level Levels for which RLNA is computed and averaged. Defaults to
#'   \code{c(1, 2)}.
#' @param range_fact Controls the range of values in the y-axis of the model
#'   used to obtain the RLNA axis in the SvA curve.
#' @param alpha Error probability used to compute two-sided confidence intervals.
#' @param log_scale Logical flag used only for \code{curve = "SvM"}. If
#'   \code{TRUE}, plots the x-axis on a log2 scale.
#' @param ... Additional arguments passed to methods.
#'
#' @details
#' Reference: Pournaghdali, A., Schwartz, B. L., Hays, J., \& Soto, F. A.
#' (2023). Sensitivity vs. awareness curve: A novel model-based analysis to
#' uncover the processes underlying nonconscious perception. \emph{Psychonomic
#' Bulletin \& Review, 30}(2), 553-563.
#'
#' @return Produces a base R plot. Invisibly returns \code{NULL}.
#'
#' @examples
#' # Example workflow:
#' # fitted_model <- grt_wind_fit(cmats, model = "GRTapas")
#' # boot <- bootstrap_sample(fitted_model, N_rows = lapply(cmats, rowSums), N = 50, n_cores = 1)
#' # sva_curve(fitted_model, boot)       # group-level curve
#' # sva_curve(fitted_model, boot, sub=1)  # participant-level curve
#'
#' @export
sva_curve <- function(x, ...) UseMethod("sva_curve")

#' @export
sva_curve.grt_wind_fit <- function(x, bootstrap_samples, sub = "group", curve = c("SvA", "SvM"),
                                   stim_level = c(1, 2), range_fact = 3, alpha = 0.05,
                                   log_scale = FALSE, ...) {
  fitted_model <- x
  curve <- match.arg(curve)

  if (!any(fitted_model$restrictions == "GRTapas")) {
    stop("Function sva_curve requires a GRTapas model")
  }

  if (!identical(sub, "group")) {
    if (!is.numeric(sub) || length(sub) != 1 || is.na(sub)) {
      stop("'sub' must be 'group' or a single participant number")
    }
    sub <- as.integer(sub)
    if (sub < 1 || sub > fitted_model$N) {
      stop("Participant index in 'sub' is out of range")
    }
    sva_curve_ind(
      fitted_model = fitted_model,
      bootstrap_samples = bootstrap_samples,
      sub = sub,
      curve = curve,
      stim_level = stim_level,
      alpha = alpha,
      log_scale = log_scale,
      range_fact = range_fact
    )
    return(invisible(NULL))
  }

  ymin <- min(fitted_model$means[, 2]) - range_fact
  ymax <- max(fitted_model$means[, 2]) + range_fact
  y_steps <- seq(ymax, ymin, by = -0.01)
  n_steps <- length(y_steps)

  rlna_vector <- rep(NA_real_, n_steps)
  d_vector <- rep(NA_real_, n_steps)
  dlow <- rep(NA_real_, n_steps)
  dhigh <- rep(NA_real_, n_steps)

  for (i in seq_len(n_steps)) {
    if (curve == "SvA") {
      rlna_vector[i] <- sva_rlna(y_steps[i], fitted_model$fullpars, stim_level = stim_level)
    } else {
      rlna_vector[i] <- sva_rlm(y_steps[i], fitted_model$fullpars, stim_level = stim_level)
    }
    d_vector[i] <- sva_d_prime_cond(y_steps[i], fitted_model$fullpars)

    rlna_edf <- rep(NA_real_, nrow(bootstrap_samples$samples))
    for (s in seq_len(nrow(bootstrap_samples$samples))) {
      rlna_edf[s] <- sva_d_prime_cond(y_steps[i], bootstrap_samples$samples[s, ])
    }

    ci <- stats::quantile(rlna_edf, c(alpha / 2, 1 - alpha / 2))
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
  ind_c <- sva_rlna(-fitted_model$indpar$a2, fitted_model$fullpars, stim_level = stim_level)

  if (curve == "SvM") {
    x_label <- "Relative Likelihood of Metacognition"
    rlna_vector <- sva_safe_log2(rlna_vector, log_scale)
    ind_c <- sva_safe_log2(sva_rlm(-fitted_model$indpar$a2, fitted_model$fullpars, stim_level = stim_level), log_scale)
    if (isTRUE(log_scale)) {
      x_label <- "Relative Likelihood of Metacognition (Log2 Scale)"
      oc <- 0
    }
  }

  graphics::plot(
    rlna_vector,
    d_vector,
    type = "l",
    xlab = x_label,
    ylab = "Sensitivity (d')",
    ylim = c(dmin, dmax),
    lwd = 2.5,
    cex.lab = 1.5,
    col = grDevices::rgb(0.8, 0.25, 0.33, 1)
  )
  graphics::polygon(
    c(rlna_vector, rev(rlna_vector)),
    c(dlow, rev(dhigh)),
    col = grDevices::rgb(1, 0, 0, 0.3),
    border = NA
  )
  graphics::abline(v = if (curve == "SvM") oc else 1, col = grDevices::rgb(0, 0, 0.8, 1), lty = 2, lwd = 1.5)

  if (any(fitted_model$restrictions == "DS(B)")) {
    graphics::abline(
      v = ind_c,
      col = grDevices::rgb(0, 0.6, 0, 1),
      lty = 3
    )
  }
  graphics::abline(h = 0, lty = 2, lwd = 1.5)

  invisible(NULL)
}
