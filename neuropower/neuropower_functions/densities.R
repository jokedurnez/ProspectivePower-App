#################################
# NEUROPOWER: DENSITY FUNCTIONS #
#################################

# Family of functions that compute the pdf and cdf of
# null and alternative distributions in the mixture model.
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
# inputs:
# 	- x: a vector of peak values
#   - par: a vector with two values: c(mu,sigma) of the
#          estimated alternative distribution
#   - pi0: the proportion of null peaks referring to the
#          proportion of non-active peaks
#		- u: excursion threshold

# 1) Function: mix.sum.log
# ------------------------
# This function computes the negative of the sum of the logs of the
# likelihoods of the mixture distribution.
# When this function is minimised over the variables 'par',
# the MLE of the alternative distribution is found.

mix.sum.log <- function(x,par,pi0,u){
	f.x <- mix.pdf(x,par,pi0,u)
	return(-sum(log(f.x)))
}

# 2) Function: mix.pdf
# ------------------------
# This function computes the probability density function
# of the mixture of null and alternative distributions with weights pi0.

mix.pdf <- function(x,par,pi0,u){
	f.x <- pi0*null.pdf(x,u) + (1-pi0)*alt.pdf(x,par,u)
	return(f.x)
}

# 3) Function: null.pdf
# ------------------------
# This function computes the probability density function
# of the null distribution under Random Field Theory.

null.pdf <- function(x,u){
	f.x <- u*exp(-u*(x-u))
	return(f.x)
}

# 3) Function: alt.pdf
# ------------------------
# This function computes the probability density function
# of the alternative distribution of peaks: a truncated normal.

alt.pdf <- function(x,par,u){
  mu <- par[1]; sigma <- par[2]
  f.x.num <- 1/sigma * dnorm((x-mu)/sigma)
	f.x.den <- 1-pnorm((u-mu)/sigma)
  f.x <- f.x.num/f.x.den
  return(f.x)
}

# 3) Function: alt.cdf
# ------------------------
# This function computes the cumulative density function
# of the alternative distribution to get values for power.

alt.cdf <- function(x,par,u){
	mu <- par[1]; sigma <- par[2]
  ksi <- (x-mu)/sigma
  alpha <- (u-mu)/sigma
  F.x <- (pnorm(ksi)-pnorm(alpha))/(1-pnorm(alpha))
  return(F.x)
}
