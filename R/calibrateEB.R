#' Fit an empirical Bayes prior in the hierarchical model
#'     mu ~ G, X ~ N(mu, sigma^2)
#'
#' @param X a vector of observations
#' @param sigma noise estimate
#' @param p tuning parameter -- number of parameters used to fit G
#' @param nbin tuning parameter -- number of bins used for discrete approximation
#' @param unif.fraction tuning parameter -- fraction of G modeled as "slab"
#'
#' @return posterior density estimate g
#'
#' @section References:
#' For more details about "g-estimation", see: B Efron. Two modeling strategies for
#' empirical Bayes estimation. Stat. Sci., 29(2): 285–301, 2014.

gfit = function(X, sigma, p = 2, nbin = 1000, unif.fraction = 0.1) {
	
	xvals = seq(min(min(X) - 2 * sd(X), 0), max(max(X) + 2 * sd(X), sd(X)), length.out = nbin)
	binw = xvals[2] - xvals[1]
	
	zero.idx = max(which(xvals <= 0))
	noise.kernel = dnorm(xvals / sigma) * binw / sigma
	
	if (zero.idx > 1) {
		noise.rotate = noise.kernel[c(zero.idx:length(xvals), 1:(zero.idx - 1))]
	} else {
		noise.rotate = noise.kernel
	}
	
	XX = sapply(1:p, function(j) xvals^j * as.numeric(xvals >= 0))
	neg.loglik = function(eta) {
		g.eta.raw = exp(XX %*% eta) * as.numeric(xvals >= 0)
		if ((sum(g.eta.raw) == Inf) | (sum(g.eta.raw) <= 100 * .Machine$double.eps)) {
			return (1000 * (length(X) + sum(eta^2)))
		}
		g.eta.main = g.eta.raw / sum(g.eta.raw)
		g.eta = (1 - unif.fraction) * g.eta.main +
			unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
		f.eta = convolve(g.eta, noise.rotate)
		sum(approx(xvals, -log(pmax(f.eta, 0.0000001)), X)$y)
	}
	
	eta.hat = nlm(neg.loglik, rep(-1, p))$estimate
	g.eta.raw = exp(XX %*% eta.hat) * as.numeric(xvals >= 0)
	g.eta.main = g.eta.raw / sum(g.eta.raw)
	g.eta = (1 - unif.fraction) * g.eta.main +
		unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
	
	return(data.frame(x=xvals, g=g.eta))
}

#' Bayes posterior estimation with Gaussian noise
#'
#' @param x0 an obsevation
#' @param g.est a prior density, as returned by gfit
#' @param sigma noise estimate
#'
#' @return posterior estimate E[mu | x0]

gbayes = function(x0, g.est, sigma) {
	Kx = dnorm((g.est$x - x0) / sigma)
	post = Kx * g.est$g
	post = post / sum(post)
	sum(post * g.est$x)
}

#' Empirical Bayes calibration of noisy variance estimates
#'
#' @param vars list of variance estimates
#' @param sigma2 estimate of the Monte Carlo noise in vars
#'
#' @return calibrated variance estimates

calibrateEB = function(vars, sigma2) {
	
	if(sigma2 <= 0 | min(vars) == max(vars)) {
		return(pmax(vars, 0))
	}
	
	sigma = sqrt(sigma2)
	eb.prior = gfit(vars, sigma)
	
	if (length(vars >= 200)) {
		# If there are many test points, use interpolation to speed up computations
		calib.x = quantile(vars, q = seq(0, 1, by = 0.02))
		calib.y = sapply(calib.x, function(xx) gbayes(xx, eb.prior, sigma))
		calib.all = approx(x=calib.x, y=calib.y, xout=vars)$y
	} else {
		calib.all = sapply(vars, function(xx) gbayes(xx, eb.prior, sigma))
	}
	return(calib.all)
}
	