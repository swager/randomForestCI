#' The infinitesimal jackknife for random forests
#'
#' @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
#' @param newdata A set of test points at which to evaluate standard errors

randomForestInfJack = function(rf, newdata, calibrate = TRUE, used.trees = NULL) {
	
	if (is.null(rf$inbag)) {
		stop("Random forest must be trained with keep.inbag = TRUE")
	}
	
	if (length(levels(factor(colSums(rf$inbag)))) > 1) {
		stop("The keep.inbag field must store the number of times each observation was used")
	}
	
	if (is.null(used.trees)) {
		used.trees = 1:rf$ntree
	}
	
	B = length(used.trees)
	n = length(rf$y)
	s = sum(rf$inbag) / rf$ntree
	
	predictions = predict(rf, newdata, predict.all = TRUE)
	pred = predictions$individual[, used.trees]
	y.hat = rowMeans(pred)
	pred.centered = pred - rowMeans(pred)
	
	N = Matrix::Matrix(rf$inbag[, used.trees], sparse = TRUE)
	N.avg = Matrix::rowMeans(N)
	
	C = N %*% t(pred.centered) - Matrix::Matrix(N.avg, nrow(N), 1) %*% Matrix::Matrix(rowSums(pred.centered), 1, nrow(pred.centered))	
	raw.IJ = Matrix::colSums(C^2) / B^2
	
	N.var = mean(Matrix::rowMeans(N^2) - Matrix::rowMeans(N)^2)
	boot.var = rowSums(pred.centered^2) / B
	bias.correction = n * N.var * boot.var / B
	vars = raw.IJ - bias.correction

	results = data.frame(y.hat=y.hat, var.hat=vars)
	
	if (nrow(results) <= 100) {
		calibrate = FALSE
		warning("No calibration with n <= 100")
	}
	
	if (calibrate) {
		calibration.ratio = 2
		n.sample = ceiling(B / calibration.ratio)
		results.ss = randomForestInfJack(rf, newdata, calibrate = FALSE, used.trees = sample(used.trees, n.sample))
		sigma2.ss = var(results.ss$var.hat - results$var.hat)
		delta = n.sample / B
		sigma2 = (delta^2 + (1 - delta)^2) / (2 * (1 - delta)^2) * sigma2.ss
		
		vars.calibrated = calibrateEB(vars, sigma2)
		results$var.hat = vars.calibrated
	}
	
	return (results)
}
