\name{data_sample}
\alias{data_sample}
\title{Sample Confusion Matrices from a Fitted Model}
\usage{
data_sample(x, ...)
}
\arguments{
\item{x}{A fitted model object. For \code{grt_wind_fit} objects, use
\code{data_sample(x, N_row = ...)}.}

\item{\dots}{Additional arguments passed to methods. For
\code{data_sample.grt_wind_fit}, \code{N_row} controls trial counts per row.}
}
\value{
A list of confusion matrices sampled from the fitted model.
}
\description{
Generic function to sample data from fitted models.
}
