#################################
# NEUROPOWER: MODELFIT FUNCTION #
#################################

# Function: modelestimator
# ---------------------------
# This function estimates the mixture of the null and alternative distribution
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
# 	- x: statistical values of peaks
#   - p: pvalues of peaks
#   - u: the excursion threshold
#   - df: number of degrees of freedom (n-1 or n-2 depending on design)
#   - n: number of subjects
#   - plot: Boolean for execution of the plot

modelestimator <- function(x,p,u,df,n,plot){
	options(warn=-1)
	# pi0 estimate using BioNet
	try(BUM <- bumOptim(p,starts=1))
	if(BUM=="BUM model could not be fitted to data"){
		cat("\n BUM model could not be fitted to the data")
		stop
	}
	pi0e <- BUM$lambda + (1-BUM$lambda)*BUM$a
	if(round(pi0e,digits=6)==1){
		cat("\n Estimated prevalence of activation is 0.  Power cannot be estimated.")
		stop
	}
	# effect estimate using MLE (L-BFGS-B algoritm)
	mixture <- optim(par=c(5,0.5),method="L-BFGS-B",lower=c(2.5,0.1),fn=mix.sum.log,x=x,pi0=pi0e,u=u)
	estimates <- list(pi0e,mixture$par,BUM$a)
	names(estimates) <- c("pi0","Ha","a")
	names(estimates$Ha) <- c("Delta","Sigma")
	if(plot==TRUE){model.plot(x,p,pi0e,estimates$a,estimates$Ha,u,n)}
	return(estimates)
}
