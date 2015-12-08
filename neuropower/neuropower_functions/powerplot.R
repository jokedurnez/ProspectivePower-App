##########################################
# NEUROPOWER: PLOT POWER CURVES FUNCTION #
##########################################

# Function: power.plot
# ----------------------
# This function creates a plot with power curves for
# different MTP corrections
#
# written by Joke Durnez (joke.durnez@gmail.com)
#

# function that generates the plot with power curves
power.plot <- function(newsubjects,powers,newsamplesize,newpower,power,subjects,MCP,alpha,samplesize){
	options(warn=-1)

	# aesthetic considerations
	cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")
	layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.7,0.3),heights=c(0.7,0.3))
	par(mar=c(4,4,2,1),oma=c(0,0,0,0))

	# power curves
	plot(newsubjects,powers$power.UN,type="l",col=cols[1],ylim=c(-0.1,1),xlim=c(min(newsubjects),max(newsubjects)),xlab="Subjects",ylab="Average power",lwd=2,main="Power curves")
	lines(newsubjects,powers$power.BF,col=cols[2],lwd=2)
	lines(newsubjects,powers$power.RFT,col=cols[3],lwd=2)
	lines(newsubjects,powers$power.BH,col=cols[4],lwd=2)

	# power in *this* study
	T <- c(powers$power.UN[1],powers$power.BF[1],powers$power.RFT[1],powers$power.BH[1])

	# dots at current study
	for(i in 1:4){points(subjects,T[i],pch=16,col=cols[i])}
	# lines at new study
	for(i in 1:4){lines(c(newsamplesize[i],newsamplesize[i]),c(0,newpower[i]),col=cols[i],lty=2)}
	for(i in 1:4){lines(c(min(newsubjects),newsamplesize[i]),c(newpower[i],newpower[i]),col=cols[i],lty=2)}
	# dots at required sample size for new study
	for(i in 1:4){points(newsamplesize[i],newpower[i],pch=16,col=cols[i])}
	# numeric value of new sample size
	for(i in 1:4){text(newsamplesize[i],-0.05,newsamplesize[i],col=cols[i])}

	# legend
	par(mar=c(0,0,0,0))
	plot(1:3, 1:3, col="white",axes=FALSE,xlab="",ylab="")
	legend(1,2.8,c("UN","BF","RFT","BH"),col=cols[1:4],lty=rep(1,4),lwd=2,box.lwd=0,bty="n",border="white",cex=1.5)

	# text with required sample size
	par(mar=c(0,0,0,0))
	plot(1:3, 1:3, col="white",axes=FALSE,xlab="",ylab="")
	samplesize <- ifelse(is.na(samplesize),">50",samplesize)
	text(1,2.5,pos=4,labels=paste("To obtain a power level of ",power,"with ",MCP,"\n control at level",alpha,", the minimal sample size is",samplesize,"."),cex=1.5)
}
