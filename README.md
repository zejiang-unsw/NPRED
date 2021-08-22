# NPRED
 
Predictor Identifier: Nonparametric PREDiction (NPRED)
Partial informational correlation (PIC) is used to identify the meaningful predictors to the response from a large set of potential predictors.

The initial version of NPRED is at [Hydrology@UNSW](https://www.hydrology.unsw.edu.au/download/software/npred). This is a new version of NPRED without calling Fortran codes.

Applications of this package can be found in: 
* Jiang, Z., Sharma, A., & Johnson, F. (2020). <doi:10.1029/2019WR026962> 
* Jiang, Z., Rashid, M. M., Johnson, F., & Sharma, A. (2020). <doi:10.1016/j.envsoft.2020.104907>
* Jiang, Z., Sharma, A., & Johnson, F. (2021). <doi:10.1016/J.JHYDROL.2021.126816>

## Installation
You can install the package via devtools from [GitHub](https://github.com/) with:

``` r
devtools::install_github("zejiang-unsw/NPRED")
```

or via CRAN with: 

``` r
install.packages("NPRED")
```

## Citations
Sharma, A., Mehrotra, R. (2014). An information theoretic alternative to model a natural system using observational information alone. Water Resources Research, 50(1): 650-660.

Galelli S., Humphrey G.B., Maier H.R., Castelletti A., Dandy G.C. and Gibbs M.S. (2014). An evaluation framework for input variable selection algorithms for environmental data-driven models, Environmental Modelling and Software, 62, 33-51.

Sharma, A., Mehrotra, R., Li, J., & Jha, S. (2016). A programming tool for nonparametric system prediction using Partial Informational Correlation and Partial Weights. Environmental Modelling & Software, 83, 271-275. 

Mehrotra, R., & Sharma, A. (2006). Conditional resampling of hydrologic time series using multiple predictor variables: A K-nearest neighbour approach. Advances in Water Resources, 29(7), 987-999. 
