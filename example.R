library(devtools)

install_github("swager/randomForest")
install_github("swager/randomForestCI")

library(randomForest)
library(randomForestCI)

n = 1000
p = 100
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)

rf = randomForest(X, Y, keep.inbag = TRUE)
ij = randomForestInfJack(rf, X, calibrate = TRUE)

plot(ij)