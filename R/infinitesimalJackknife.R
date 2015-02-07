#' The infinitesimal jackknife for random forests with subsampling
#'
#' @param rf A random forest trained with replace = TRUE and keep.inbag = TRUE
#' @param newdata A set of test points at which to evaluate standard errors

randomForestInfJackSubsample = function(rf, newdata) {
	predictions = predict(rf, newdata, predict.all = TRUE);
	agg_pred = predictions$aggregate;
	pred = predictions$individual;
	N = rf$inbag;
	Nc = N - rowMeans(N);
	
	cent_pred = pred - rowMeans(pred);
	raw.vars = colSums((Nc %*% t(cent_pred))^2) / rf$ntree/(rf$ntree - 1);
	bias.corr = colSums(Nc^2 %*% t(cent_pred)^2) / rf$ntree/(rf$ntree - 1);
	vars = raw.vars - bias.corr;
	
	var.mc = sapply(1:nrow(pred), function(rr)var(pred[rr,])) / rf$ntree
	data.frame(predictions=agg_pred, variance=vars, var.mc = var.mc);
}