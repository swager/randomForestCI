library(devtools)

install_github("swager/randomForest")
install_github("swager/randomForestCI")

library(randomForest)
library(randomForestCI)

X = matrix(rnorm(200 * 100), 200, 100)
Y = rnorm(200)

rf = randomForest(X, Y, keep.inbag = TRUE)
ij = randomForestInfJack(rf, X)

plot(ij)