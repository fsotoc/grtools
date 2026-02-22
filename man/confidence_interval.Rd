\name{confidence_interval}
\alias{confidence_interval}
\title{Calculate Confidence Intervals from Parameter Samples}
\usage{
confidence_interval(input_matrix, mle_pars = NULL, col_selection = NULL,
  alpha = 0.05)
}
\arguments{
\item{input_matrix}{An \eqn{N \times K} matrix of parameter samples with
parameters in columns and samples in rows.}

\item{mle_pars}{Optional vector of parameter estimates for selected columns.
If omitted, estimates are computed from \code{input_matrix}.}

\item{col_selection}{Optional vector of selected columns (indexes or names). If
omitted, all columns are used.}

\item{alpha}{Error probability for two-sided intervals. Default is
\code{0.05}.}
}
\value{
A matrix with three rows for each selected parameter: estimate, lower bound,
and upper bound.
}
\description{
Computes percentile-based confidence intervals for parameter samples.
}
\details{
The percentile method is used to calculate confidence intervals from a set of
samples.
}
\references{
Carpenter, J., \& Bithell, J. (2000). Bootstrap confidence intervals: when,
which, what? A practical guide for medical statisticians. Statistics in
medicine, 19(9), 1141-1164.
}
\examples{
data <- matrix(c(161,41,66,32,24,147,64,65,37,41,179,43,14,62,22,202,
26,82,67,25,8,188,54,50,34,75,172,19,7,103,14,176,
117,64,89,30,11,186,69,34,21,81,176,22,10,98,30,162,
168,57,47,28,15,203,33,49,58,54,156,32,9,96,9,186,
169,53,53,25,34,168,69,29,38,48,180,34,19,44,60,177), nrow = 5, byrow = TRUE)

ci_data <- confidence_interval(data, alpha = 0.05)
ci_data
}
