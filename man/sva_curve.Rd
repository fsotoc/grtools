\name{sva_curve}
\alias{sva_curve}
\alias{sva_curve.grt_wind_fit}
\title{Plot Sensitivity Versus Awareness (SvA) Curves}
\usage{
sva_curve(x, ...)
}
\arguments{
\item{x}{A fitted model returned by \code{\link{grt_wind_fit}} or
\code{\link{grt_wind_fit_parallel}} under GRTapas restrictions.}

\item{\dots}{Additional arguments passed to methods. For
\code{sva_curve.grt_wind_fit}, these include:
\itemize{
\item \code{bootstrap_samples}: output from \code{\link{bootstrap_sample}}
for the same model.
\item \code{sub}: flag indicating whether to plot the group-level curve or an
individual participant curve. Use \code{"group"} (default) for the group
model; provide a participant number for an individual model curve.
\item \code{curve}: curve type. Use \code{"SvA"} (default) for the original
sensitivity-versus-awareness curve, or \code{"SvM"} for
sensitivity-versus-metacognition logic.
\item \code{stim_level}: stimulus levels used to compute RLNA
(default \code{c(1, 2)}).
\item \code{range_fact}: range factor for the y-axis values used to compute
the SvA curve (default \code{3}).
\item \code{alpha}: error probability used for two-sided confidence intervals
(default \code{0.05}).
\item \code{log_scale}: logical flag used only when \code{curve = "SvM"}; if
\code{TRUE}, the x-axis is plotted on log2 scale.
}
}
}
\value{
Produces a base R plot. Invisibly returns \code{NULL}.
}
\description{
Plots sensitivity-versus-awareness (SvA) or sensitivity-versus-metacognition
(SvM) curves from a fitted GRTapas model, including bootstrap confidence bands
for conditional sensitivity.
}
\details{
If \code{sub = "group"}, the group model SvA curve is plotted. If \code{sub}
is a participant index, the function delegates to \code{sva_curve_ind} and
plots that participant's curve using the selected \code{curve} mode.
}
