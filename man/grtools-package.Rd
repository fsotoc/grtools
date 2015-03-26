\name{grtools-package}
\alias{grtools-package}
\alias{grtools}
\docType{package}
\title{
General recognition theory tools for the analysis of perceptual independence
}
\description{
Statistical tools from General Recognition Theory for psychophysical data analyses aimed at determining independent processing of perceptual dimensions
}
\details{
\tabular{ll}{
Package: \tab grtools\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2015-03-26\cr
License: \tab GPL (>= 2)\cr
}
\code{grtools} provides functions for the following analyses using general recognition theory:

1. Model-based analyses of separability and independence with GRT-wIND for the 2x2 identification experiment (Soto et al., 2015). See \code{\link{grt_wind_fit}} and \code{\link{grt_wind_fit_parallel}}

2. Model-based analyses of separability and independence with traditional GRT models for the 2x2 identification experiment (Ashby & Soto, 2015). See \code{\link{grt_hm_fit}}

3. Summary statistics analysis (i.e. Kadlec's MDSDA; see Kadlec & Townsend, 1992) for the 2x2 identification experiment. See \code{\link{sumstats_micro}} and \code{\link{sumstats_macro}}

4. Summary statistic analysis for the 2x2 Garner filtering task (Ashby & Maddox, 1994). See \code{\link{sumstats_garner}}
}

\author{
Fabian Soto, Emily Zheng

Maintainer: Fabian Soto <fabian.soto@psych.ucsb.edu>
}
\references{
Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal detection theory. In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A. Eidels (Eds.), \emph{Oxford handbook of computational and mathematical psychology} (pp. 13-34). Oxford University Press: New York, NY.

Ashby, F. G., & Maddox, W. T. (1994). A response time theory of separability and integrality in speeded classification. \emph{Journal of Mathematical Psychology, 38}(4), 423-466.

Kadlec, H., & Townsend, J. T. (1992). Signal detection analyses of multidimensional interactions. In F. G. Ashby (Ed.), \emph{Multidimensional models of perception and cognition} (pp. 181–231). Hillsdale, NJ: Erlbaum.

Soto, F. A., Musgrave, R., Vucovich, L., & Ashby, F. G. (2015). General recognition theory with individual differences: A new method for examining perceptual and decisional interactions with an application to face perception. \emph{Psychonomic Bulletin & Review, 22}(1), 88-111.
  
}

\keyword{ package, general recognition theory, grt, mdsdt, signal detection theory, 
perceptual separability, perceptual independence, separability }
\seealso{
For applications of General Recognition Theory to perceptual categorization experiments, see \code{\href{http://cran.r-project.org/web/packages/grt/index.html}{grt}
}
}