set.seed(1987)

library(randomForest)

n = 250
p = 100
X = matrix(rnorm(n * p), n, p)
Y = rnorm(n)
rf = randomForest(X, Y, keep.inbag = TRUE)
pred = predict(rf, newdata = X, predict.all = TRUE)

rf.small = randomForest(X[1:20, ], Y[1:20], keep.inbag = TRUE)

test_that("sanity", {
  ij = randomForestInfJack(rf, X, calibrate = FALSE)
  expect_that(all(is.finite(ij$var.hat)), is_true())
  
  ij.calibrated = randomForestInfJack(rf, X, calibrate = TRUE)
  expect_that(all(ij.calibrated$var.hat > 0), is_true())

  expect_warning(randomForestInfJack(rf, X[1:20, ], calibrate = TRUE))

  ij.standalone = infJack(pred$individual, rf$inbag, calibrate = FALSE)
  expect_equal(ij.standalone, ij)

  ij.standalone.calibrated = infJack(pred$individual, rf$inbag, calibrate = TRUE)
  expect_equal(ij.standalone.calibrated, ij.calibrated, tolerance = .1)
})
