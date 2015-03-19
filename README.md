# grtools
**R package for the analysis of perceptual independence using general recognition theory.**

**grtools** provides functions for the following analyses:

1. Model-based analyses of separability and independence with GRT-wIND (Soto et al., 2015) for the 2x2 identification experiment.
2. Model-based analyses of separability and independence with traditional GRT models for the 2x2 identification experiment (Ashby & Soto, 2015).
3. Summary statistics analysis (i.e. Kadlec's MSDA; see Kadlec & Townsend, 1992) for the 2x2 identification experiment.
4. Summary statistics analysis for the 2x2 Garner filtering task (Ashby & Maddox, 1994).

The easiest way to instal **grtools** and its dependencies is using devtools. In the R console, type:

```R
install.packages("devtools")
devtools::install_github("fsotoc/grtools", dependencies="Imports")
```

After installation. Type the following in the R console:

```R
library(grtools)
?grtools
```

This will open a document that includes links to help documentation for each of the main analyses included in **grtools** (including examples).

**grtools** requires [Rcpp](http://cran.r-project.org/web/packages/Rcpp/index.html) to work, which in turns requires a development environment with a suitable compiler. If you are getting compiling errors, see the [Rcpp FAQ](http://cran.r-project.org/web/packages/Rcpp/vignettes/Rcpp-FAQ.pdf), particularly points 1.2 and 1.3.

**Note that this package is still under development. This is a pre-release version that has not been extensively tested. We welcome your comments, feature requests, bug reports, etc.**

References
----------
Ashby, F. G., & Maddox, W. T. (1994). A response time theory of separability and integrality in speeded classification. *Journal of Mathematical Psychology, 38*(4), 423-466.

Ashby, F. G., & Soto, F. A. (2015). Multidimensional signal detection theory. In J. R. Busemeyer, J. T. Townsend, Z. J. Wang, & A. Eidels (Eds.), *Oxford handbook of computational and mathematical psychology* (pp. 13-34). Oxford University Press: New York, NY.

Kadlec, H., & Townsend, J. T. (1992). Signal detection analyses of multidimensional interactions. In F. G. Ashby
(Ed.), *Multidimensional models of perception and cognition* (pp. 181–231). Hillsdale, NJ: Erlbaum.

Soto, F. A., Musgrave, R., Vucovich, L., & Ashby, F. G. (2015). General recognition theory with individual differences: A new method for examining perceptual and decisional interactions with an application to face perception. *Psychonomic Bulletin & Review, 22*(1), 88-111.

