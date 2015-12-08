######################################
# NEUROPOWER: PLOT MODELFIT FUNCTION #
######################################

# Function: modelfitplot
# ----------------------
# This function creates a plot with a histogram of the observed
# distribution and overlays of the distributions using the MLE estimator
# of the null and alternative distribution.
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
# 	- x: statistical values of peaks
#		- p: p-values of peaks
#   - pi0: the estimated proportion of null peaks
#   - mu: the estimated centrality parameter of the alternative distribution
#   - sigma: the estimated variance parameter of the alternative distribution
#   - u: the excursion threshold
#   - n: number of subjects

model.plot <- function(x,p,pi0,a,par,u,n){
  # set parameters for visualisation
	par(mfrow=c(1,2),mar=c(4,5,5,1))
	col <- c("#A6CEE3","#1F78B4","#41AB5D","#006D2C")

	# left plot: histogram of uncorrected p-values for peaks with fitted distributions

  # histogram
  breakshist <- seq(from=0,to=1,by=0.1)
	histogram <- hist(p,freq=FALSE,col=col[1],border=col[1],breaks=breakshist,xlab="peak p-values",main=paste("Distribution of",length(p),"peak p-values"))
  # overlay lines
	xn <- seq(from=0,to=1,length=1000)
	lines(xn,rep(pi0,1000),col=col[3],lwd=3)
	lines(xn,(dbeta(xn,a,1)*(1-pi0)+pi0),col=col[2],lwd=3)
  # title and legend
	mtext(bquote(pi[1] ~ "=" ~ .(round(1-pi0,digits=2))))
	legend(0.3,max(histogram$density),c("Estimated null","Estimated total"),fill=col[c(3,2)],border=col[c(3,2)],box.lwd=0,box.col="white")

	# right plot: histogram of peak heights with fitted distributions

  # histogram
	breakshist <- seq(from=u,to=15,by=0.40)
	histo <- hist(x,freq=FALSE,col=col[1],border=col[1],breaks=breakshist,xlab="peak heights",main="Distribution of peak heights",xlim=c(u,8))
  # lines
  x <- seq(from=u,to=15,length=10000)
	y.tot <- mix.pdf(x,par,pi0,u)
	y.nul <- null.pdf(x,u)
	y.alt <- alt.pdf(x,par,u)
	lines(x,y.alt*(1-pi0),col=col[3],lwd=3)
	lines(x,y.nul*pi0,col=col[4],lwd=3)
	lines(x,y.tot,col=col[2],lwd=3)
  # title and legend
  mtext(bquote(pi[1] ~ "=" ~ .(round(1-pi0,digits=2)) ~ "-" ~ delta ~ "=" ~ .(round(par[1]/sqrt(n),digits=2))))
	legend(3.7,max(histo$density),c("Estimated alternative","Estimated null","Estimated total"),fill=col[c(3,4,2)],border=col[c(3,4,2)],box.lwd=0,box.col="white")
}
