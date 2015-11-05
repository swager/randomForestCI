library(devtools)

install_github("swager/randomForestCI")

library(randomForest)
library(randomForestCI)

#
# Simple simulation
#

n = 250
p = 100
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)

rf = randomForest(X, Y, keep.inbag = TRUE)
ij = randomForestInfJack(rf, X, calibrate = TRUE)

plot(ij)

#
# Spambase example
#

spam.data = read.csv("http://archive.ics.uci.edu/ml/machine-learning-databases/spambase/spambase.data",
                     header = FALSE)

# Separate into predictors and response
X = spam.data[, 1:57]
Y = spam.data[, 58]

# Note -- some variance estimates with y.hat ~~ 0.5 may be close to 0;
# this is a Monte Carlo noise effect
rf.spam = randomForest(X, factor(Y), keep.inbag = TRUE, ntree = 500)
ij.spam = randomForestInfJack(rf.spam, X, calibrate = TRUE)
plot(ij.spam)

# With 2000 trees, training takes longer, but the variance estimates are much better
rf.spam2 = randomForest(X, factor(Y), keep.inbag = TRUE, ntree = 2000)
ij.spam2 = randomForestInfJack(rf.spam2, X, calibrate = TRUE)
plot(ij.spam2)
