#' The infinitesimal jackknife for random forests with subsampling
#'
#' @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
#' @param newdata A set of test points at which to evaluate standard errors

randomForestInfJackSubsample = function(rf, newdata) {
	
	if (is.null(rf$inbag)) {
		stop("Random forest must be trained with keep.inbag = TRUE")
	}
	
	B = rf$ntree
	n = length(rf$y)
	s = sum(rf$inbag) / B
	
	obs.seen = sum(rf$inbag) + sum(rf$oob.times)
	if (obs.seen != B * n | s %% 1 != 0) {
		stop("Random forest must be trained with subsampling (i.e., replace = FALSE)")
	}
	
	predictions = predict(rf, newdata, predict.all = TRUE)
	agg_pred = predictions$aggregate
	pred = predictions$individual
	pred.centered = pred - rowMeans(pred)
	
	N = Matrix(rf$inbag, sparse = TRUE)
	
	raw.IJ = colSums((N %*% t(pred.centered) - s/n * Matrix(1, nrow(N), 1) %*% Matrix(colSums(pred.centered), 1, nrow(pred.centered)))^2) / B^2
	boot.var = rowSums(pred.centered^2) / B
	bias.correction = s * (n - s) / n * boot.var / B
	vars = raw.IJ - bias.correction;

	data.frame(y.hat=agg_pred, var.hat=vars);
}
