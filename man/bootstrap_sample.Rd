\name{bootstrap_sample}
\alias{bootstrap_sample}
\title{Obtain Bootstrap Samples for a Fitted GRT-wIND Model}
\usage{
bootstrap_sample(x, ...)
}
\arguments{
\item{x}{An object returned by \code{\link{grt_wind_fit}} or
\code{\link{grt_wind_fit_parallel}}.}

\item{\dots}{Additional arguments passed to methods. For
\code{bootstrap_sample.grt_wind_fit}, these include:
\itemize{
\item \code{N_rows}: number of trials per row. It can be a single number, a
length-4 vector, or a list with one entry per participant.
\item \code{N}: number of bootstrap samples (default \code{1000}).
\item \code{alpha}: error probability used for two-sided confidence intervals
(default \code{0.05}).
\item \code{n_cores}: number of cores to use; \code{0} means auto-detect.
}
}
}
\value{
An object of class \code{"model_samples"} containing:
\itemize{
\item \code{samples}: bootstrap parameter samples.
\item \code{cis}: confidence intervals computed from bootstrap samples.
}
}
\description{
Produces bootstrap samples for parameter estimates of a fitted GRT-wIND model
and returns bootstrap confidence intervals.
}
\details{
This procedure uses parametric bootstrap sampling from a fitted GRT-wIND model.
Each bootstrap sample is obtained by simulating data from the fitted model,
refitting the model, and storing estimated parameters.
}
\examples{
# Create list with confusion matrices.
cmats <- list(matrix(c(161,41,66,32,24,147,64,65,37,41,179,43,14,62,22,202), nrow = 4, byrow = TRUE))
cmats[[2]] <- matrix(c(126,82,67,25,8,188,54,50,34,75,172,19,7,103,14,176), nrow = 4, byrow = TRUE)
cmats[[3]] <- matrix(c(117,64,89,30,11,186,69,34,21,81,176,22,10,98,30,162), nrow = 4, byrow = TRUE)
cmats[[4]] <- matrix(c(168,57,47,28,15,203,33,49,58,54,156,32,9,96,9,186), nrow = 4, byrow = TRUE)
cmats[[5]] <- matrix(c(169,53,53,25,34,168,69,29,38,48,180,34,19,44,60,177), nrow = 4, byrow = TRUE)

# Fit model and bootstrap.
fitted_model <- grt_wind_fit(cmats)
N_rows <- lapply(cmats, rowSums)
bootstrap_samples <- bootstrap_sample(fitted_model, N_rows = N_rows, N = 10, n_cores = 1)

# Inspect confidence intervals.
bootstrap_samples$cis
}
