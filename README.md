# randomForestCI

&#x1F534; This package is deprecated. Please use one of the following packages instead: &#x1F534;

- <a href="https://github.com/swager/grf">grf</a>, which has built-in support for resampling-based confidence intervals, or
- <a href="https://github.com/imbs-hl/ranger">ranger</a>, which has an actively maintained version of the infinitesimal jackknife for random forests.

&#x1F534; Both packages are available from CRAN. &#x1F534;

Confidence intervals for random forests using the infinitesimal jackknife, as developed by Efron (2014) and Wager et al. (2014).

To install this package in R, run the following commands:

```R
install.packages("devtools")
library(devtools) 
install_github("swager/randomForestCI")
```

Example usage:

```R
library(randomForestCI)

# Make some data...
n = 250
p = 100
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
  
#  Run the method
rf = randomForest(X, Y, keep.inbag = TRUE)
ij = randomForestInfJack(rf, X, calibrate = TRUE)

plot(ij)
```

#### References
Efron, Bradley. <b>Estimation and accuracy after model selection.</b> <i>Journal of the American Statistical Association</i>, 109(507), 2014. [<a href="http://statweb.stanford.edu/~ckirby/brad/papers/2013ModelSelection.pdf">link</a>]

Wager, Stefan, Trevor Hastie, and Bradley Efron. <b>Confidence intervals for random forests: The jackknife and the infinitesimal jackknife.</b> <i>The Journal of Machine Learning Research</i>, 15(1), 2014. [<a href="http://jmlr.csail.mit.edu/papers/volume15/wager14a/wager14a.pdf">link</a>]
