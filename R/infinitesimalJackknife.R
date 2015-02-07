#' The infinitesimal jackknife for random forests
#'
#' @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
#' @param newdata A set of test points at which to evaluate standard errors

randomForestInfJack = function(rf, newdata) {
	
	if (is.null(rf$inbag)) {
		stop("Random forest must be trained with keep.inbag = TRUE")
	}
	
	if (length(levels(factor(colSums(rf$inbag)))) > 1) {
		stop("The keep.inbag field must store the number of times each observation was used")
	}
	
	B = rf$ntree
	n = length(rf$y)
	s = sum(rf$inbag) / B
	
	predictions = predict(rf, newdata, predict.all = TRUE)
	y.hat = predictions$aggregate
	pred = predictions$individual
	pred.centered = pred - rowMeans(pred)
	
	N = Matrix(rf$inbag, sparse = TRUE)
	N.avg = rowMeans(N)
	
	C = N %*% t(pred.centered) - Matrix(N.avg, nrow(N), 1) %*% Matrix(rowSums(pred.centered), 1, nrow(pred.centered))	
	raw.IJ = colSums(C^2) / B^2
	
	N.var = mean(rowMeans(N^2) - rowMeans(N)^2)
	boot.var = rowSums(pred.centered^2) / B
	bias.correction = n * N.var * boot.var / B
	vars = raw.IJ - bias.correction;

	data.frame(y.hat=y.hat, var.hat=vars);
}