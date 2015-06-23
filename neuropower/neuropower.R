#DEPENDENCIES: BioNet, AnalyzeFMRI, QVALUE

# cluster function: extracts peaks with a 26-point clustering algoritm
cluster <- function(zmap,u){
	# enlarge the map first for easier calculation at borders
	options(warn=-1)
	dimZ <- dim(zmap)
	Zstar <- array(0,dim=dimZ+c(2,2,2))  
	Zstar[2:(dimZ[1]+1),2:(dimZ[2]+1),2:(dimZ[3]+1)] <- zmap
	Zstar[is.na(Zstar)] <- -100
	# empty container of peaks and coordinates
	peaks <- xco <- yco <- zco <- c()  
	for(m in 2:(dim(Zstar)[1]-1)){
		incProgress(1/m,detail=paste("Percentage finished:",100*round(m/(dim(Zstar)[1]-1),digits=2),"%")) # shows progression bar
		for(n in 2:(dim(Zstar)[2]-1)){
			for(o in 2:(dim(Zstar)[3]-1)){
				if(Zstar[m,n,o]<u){next}
				# creates a container with all values around voxel, if none are higher, it's a peak!
				surroundings <- c(Zstar[m-1,n+1,o],
				Zstar[m,n+1,o],
				Zstar[m+1,n+1,o],
				Zstar[m-1,n,o],
				Zstar[m+1,n,o],
				Zstar[m-1,n-1,o],
				Zstar[m,n-1,o],
				Zstar[m+1,n-1,o],

				Zstar[m-1,n+1,o-1],
				Zstar[m,n+1,o-1],
				Zstar[m+1,n+1,o-1],
				Zstar[m-1,n,o-1],
				Zstar[m+1,n,o-1],
				Zstar[m-1,n-1,o-1],
				Zstar[m,n-1,o-1],
				Zstar[m,n,o-1],
				Zstar[m+1,n-1,o-1],

				Zstar[m-1,n+1,o+1],
				Zstar[m,n+1,o+1],
				Zstar[m+1,n+1,o+1],
				Zstar[m-1,n,o+1],
				Zstar[m+1,n,o+1],
				Zstar[m-1,n-1,o+1],
				Zstar[m,n,o+1],
				Zstar[m,n-1,o+1],
				Zstar[m+1,n-1,o+1])

				if(Zstar[m,n,o]>max(surroundings)){
					peaks <- c(peaks,Zstar[m,n,o])
					xco <- c(xco,m-1)
					yco <- c(yco,n-1)
					zco <- c(zco,o-1)          
				}
			}
		}
	}
	results <- data.frame(peaks,xco,yco,zco)        
	return(results)
}

# the (negative) sum of the log of the likelihood: 
# when minimising this function over "par" we find the MLE
sumlogtruncdensBIS <- function(x,par,pi0,a){
	mu.2 <- par[1]; sigma.2 <- par[2]
	x2 <- x-a
	f.x.0 <- a*exp(-a*(x2))

	f.x.num.1 <- 1/sigma.2 * dnorm((x-mu.2)/sigma.2)
	F.x.den.1 <- 1-pnorm((a-mu.2)/sigma.2)

	f.x <- pi0*f.x.0 + (1-pi0)*(f.x.num.1/F.x.den.1)
	return(-sum(log(f.x)))
}

# density function as above, for plots 
truncdensBIS <- function(x,par,pi0,a){
	mu.2 <- par[1]; sigma.2 <- par[2]
	x2 <- x-a
	f.x.0 <- a*exp(-a*(x2))

	f.x.num.1 <- 1/sigma.2 * dnorm((x-mu.2)/sigma.2)
	F.x.den.1 <- 1-pnorm((a-mu.2)/sigma.2)

	f.x <- pi0*f.x.0 + (1-pi0)*f.x.num.1/F.x.den.1
	return(f.x)
}

# density of null peaks, for plots
truncdensNUL <- function(x,a){
	f.x.0 <- a*exp(-a*(x-a))
	return(f.x.0)
}

# density of alternative peaks (truncated normal), for plots
CDFtrunc <- function(x,mu,sigma,a){
	ksi <- (x-mu)/sigma
	alpha <- (a-mu)/sigma
	F.x <- (pnorm(ksi)-pnorm(alpha))/(1-pnorm(alpha))
	return(F.x)
}

# function returns plot with estimated distribution for \hat(pi_0) and effect 
NPplot <- function(estimates,peaklist){
	options(warn=-1)
	pi0e <- estimates[[1]]
	mixture <- estimates[[2]]
	u <- estimates$u
	par(mfrow=c(1,2),mar=c(4,5,5,1))
	cols <- c("#A6CEE3","#1F78B4","#41AB5D","#006D2C")
	breakshist <- seq(from=0,to=1,by=0.1)
	# left plot: histogram of uncorrected p-values for peaks with fitted distributions
	histo <- hist(peaklist$pvalue,freq=FALSE,col=cols[1],border=cols[1],breaks=breakshist,xlab="peak p-values",main=paste("Distribution of",length(peaklist$pvalue),"peak p-values"))
	xn <- seq(from=0,to=1,length=1000)
	lines(xn,rep(pi0e,1000),col=cols[3],lwd=3)
	lines(xn,(dbeta(xn,estimates$a,1)*(1-pi0e)+pi0e),col=cols[2],lwd=3)
	mtext(bquote(pi[1] ~ "=" ~ .(round(1-pi0e,digits=2))))
	legend(0.3,max(histo$density),c("Estimated null","Estimated total"),fill=cols[c(3,2)],border=cols[c(3,2)],box.lwd=0,box.col="white")
	
	# right plot: histogram of peak heights with fitted distributions
	x <- seq(from=2,to=15,length=10000)
	yt <- truncdensBIS(x,mixture,pi0=pi0e,u)
	yN <- truncdensBIS(x,mixture,pi0=1,u)
	yA <- truncdensBIS(x,mixture,pi0=0,u)
	breakshist <- seq(from=2,to=15,by=0.40)
	histo <- hist(peaklist$peaks,freq=FALSE,col=cols[1],border=cols[1],breaks=breakshist,xlab="peak heights",main="Distribution of peak heights",xlim=c(u,8))
	mtext(bquote(pi[1] ~ "=" ~ .(round(1-pi0e,digits=2)) ~ "-" ~ delta ~ "=" ~ .(round(mixture[1]/sqrt(estimates$subjects),digits=2))))
	lines(x,yA*(1-pi0e),col=cols[3],lwd=3)
	lines(x,yN*pi0e,col=cols[4],lwd=3)
	lines(x,yt,col=cols[2],lwd=3)
	legend(3.7,max(histo$density),c("Estimated alternative","Estimated null","Estimated total"),fill=cols[c(3,4,2)],border=cols[c(3,4,2)],box.lwd=0,box.col="white")
}

# function estimates pi0 and effect size  
NPestimate <- function(peaklist,STAT,u,df,subs,plot){
	options(warn=-1)
	# pi0 estimate using BioNet
	try(BUM <- bumOptim(peaklist$pvalue,starts=1))
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
	mixture <- optim(par=c(5,0.5),method="L-BFGS-B",lower=c(2.5,0.1),fn=sumlogtruncdensBIS,x=peaklist$peaks,pi0=pi0e,a=u)
	estimates <- list(pi0e,mixture$par)
	names(estimates) <- c("pi1","Ha")
	names(estimates$Ha) <- c("Delta","Sigma")
	estimates$u <- u
	estimates$a <- BUM$a
	estimates$lambda <- BUM$lambda
	estimates$peaks <- peaklist
	estimates$subjects <- subs
	if(plot==TRUE){NPplot(estimates,peaklist)}
	return(estimates)
}

# function estimates post-hoc power  
NPposthoc <- function(estimates,subjects,MCP,u,alpha,resels){
	options(warn=-1)
	# for a range of values: what is the probability of X>u
	thresh <- seq(from=u,to=15,length=100)
	cdfN <- exp(-u*(thresh-u))
	cdfN_RFT <- resels*exp(-thresh^2/2)*thresh^2
	# FDR cutoff estimation
	ps <- estimates$peaks$pvalue
	pvalms <- sort(ps)
	orderpvalms <- rank(ps)
	FDRqval <- (orderpvalms/length(ps))*alpha
	pr <- ifelse(pvalms[orderpvalms] < FDRqval,1,0)
	FDRc <- ifelse(sum(pr)==0,0,max(ps[pr==1]))
	# cutoffs for all MCP's
	cutoff.BH <- ifelse(FDRc==0,NA,thresh[min(which(cdfN<FDRc))])
	Q <- qvalue(ps,fdr.level=alpha)
	cutoff.Q <- ifelse(!is.list(Q),NA,ifelse(sum(Q$significant)==0,NA,min(estimates$peaks$peaks[Q$significant==TRUE])))
	cutoff.UN <- thresh[min(which(cdfN<alpha))]
	cutoff.FWE <- thresh[min(which(cdfN<(alpha/length(estimates$peaks))))]
	cutoff.RFT <- thresh[min(which(cdfN_RFT<alpha))]

	# power computed with CDF of alternative distribution
	power.BH <- 1-CDFtrunc(cutoff.BH,estimates$Ha[1],estimates$Ha[2],2)
	power.Q <- 1-CDFtrunc(cutoff.Q,estimates$Ha[1],estimates$Ha[2],2)
	power.UN <- 1-CDFtrunc(cutoff.UN,estimates$Ha[1],estimates$Ha[2],2)
	power.FWE <- 1-CDFtrunc(cutoff.FWE,estimates$Ha[1],estimates$Ha[2],2)
	power.RFT <- 1-CDFtrunc(cutoff.RFT,estimates$Ha[1],estimates$Ha[2],2)

	# posthoc power for requested MCP
	if(MCP == "BH"){
		power <- power.BH
	} else if (MCP == "Q"){
		power <- power.Q
	} else if (MCP == "FWE"){
		power <- power.FWE
	} else if (MCP == "UN"){
		power <- power.UN
	} else if (MCP == "RFT"){
	power <- power.RFT
	} else {cat("Unknown MCP"); stop}

	power <- ifelse(is.na(power),0,power)
	powertext <- paste("\n The power of this study with",subjects,"subjects with",MCP,"control at level",alpha,"equals",round(power,digits=3),": \n the average chance of detecting an activated peak is",100*round(power,digits=3),"%.")
	return(powertext)

}

# this function computes power and minimal sample size for future study

NPsamplesize <- function(estimates,subjects,maxsubjects,MCP,u,alpha,power,resels,plot){
	options(warn=-1)

	# compute cutoffs in original analysis

	# compute standardised effect
	eff.sta <- estimates$Ha[1]/sqrt(subjects)
	# compute FDR threshold
	thresh <- seq(from=u,to=15,length=100)
	cdfN <- exp(-u*(thresh-u))
	cdfN_RFT <- resels*exp(-thresh^2/2)*thresh^2
	ps <- estimates$peaks$pvalue
	pvalms <- sort(ps)
	orderpvalms <- rank(ps)
	FDRqval <- (orderpvalms/length(ps))*alpha
	pr <- ifelse(pvalms[orderpvalms] < FDRqval,1,0)
	FDRc <- ifelse(sum(pr)==0,0,max(FDRqval[pr==1]))
	cutoff.BH <- ifelse(FDRc==0,NA,thresh[min(which(cdfN<FDRc))])
	# compute Qvalue threshold
	Q <- qvalue(ps,fdr.level=alpha)
	cutoff.Q <- ifelse(!is.list(Q),NA,ifelse(sum(Q$significant)==0,NA,min(estimates$peaks$peaks[Q$significant==TRUE])))
	# compute threshold for uncorrected, FWE and RFT
	cutoff.UN <- thresh[min(which(cdfN<alpha))]
	cutoff.FWE <- thresh[min(which(cdfN<(alpha/length(estimates$peaks))))]
	cutoff.RFT <- thresh[min(which(cdfN_RFT<alpha))]

	# compute new power for a range of number of subjects
	power.BH <- power.Q <- power.UN <- power.FWE <- power.RFT <- c()
	newsubjects <- seq(from=subjects,to=maxsubjects)
	for(p in 1:length(newsubjects)){
		proj.eff <- eff.sta*sqrt(newsubjects[p])
		newpar <- c(proj.eff,estimates$Ha[2])      
		power.BH[p] <- 1-CDFtrunc(cutoff.BH,proj.eff,estimates$Ha[2],u)
		power.Q[p] <- 1-CDFtrunc(cutoff.Q,proj.eff,estimates$Ha[2],u)
		power.UN[p] <- 1-CDFtrunc(cutoff.UN,proj.eff,estimates$Ha[2],u)
		power.FWE[p] <- 1-CDFtrunc(cutoff.FWE,proj.eff,estimates$Ha[2],u)
		power.RFT[p] <- 1-CDFtrunc(cutoff.RFT,proj.eff,estimates$Ha[2],u)
	}

	# compute minimal samplesize for required power
	newsamplesize <- c(newsubjects[min(which(power.BH>power))],newsubjects[min(which(power.Q>power))],newsubjects[min(which(power.UN>power))],newsubjects[min(which(power.FWE>power))],newsubjects[min(which(power.RFT>power))])

	# combine powers for plot
	powers <- data.frame(power.BH,power.Q,power.UN,power.FWE,power.RFT)
	newpower <- c(power.BH[min(which(power.BH>power))],power.Q[min(which(power.Q>power))],power.UN[min(which(power.UN>power))],power.FWE[min(which(power.FWE>power))],power.RFT[min(which(power.RFT>power))]) # EXACT power larger than required power

	if(MCP == "BH"){
		samplesize <- newsamplesize[1]
	} else if (MCP == "Q"){
		samplesize <- newsamplesize[2]
	} else if (MCP == "UN"){
		samplesize <- newsamplesize[3]
	} else if (MCP == "FWE"){
		samplesize <- newsamplesize[4]
	} else if (MCP == "RFT"){
		samplesize <- newsamplesize[4]  
	} else {cat("Unknown MCP"); stop}

	if(plot==TRUE){plotje <- NPssplot(newsubjects,powers,newsamplesize,newpower,power,subjects,MCP,alpha,samplesize)}

	res <- list()
	res$plotje <- plotje
	return(res)
} 
  
# function that generates the plot with power curves
NPssplot <- function(newsubjects,powers,newsamplesize,newpower,power,subjects,MCP,alpha,samplesize){
	options(warn=-1)
	cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00")

	layout(matrix(c(1,2,3,3),2,2,byrow=TRUE),widths=c(0.7,0.3),heights=c(0.7,0.3))
	par(mar=c(4,4,2,1),oma=c(0,0,0,0))

	# power curves
	plot(newsubjects,powers$power.BH,type="l",col=cols[1],ylim=c(-0.1,1),xlim=c(min(newsubjects),max(newsubjects)),xlab="Subjects",ylab="Average power",lwd=2,main="Power curves")
	lines(newsubjects,powers$power.Q,col=cols[2],lwd=2) 
	lines(newsubjects,powers$power.UN,col=cols[3],lwd=2)
	lines(newsubjects,powers$power.FWE,col=cols[4],lwd=2)
	lines(newsubjects,powers$power.RFT,col=cols[5],lwd=2)

	# power in *this* study
	T <- c(powers$power.BH[1],powers$power.Q[1],powers$power.UN[1],powers$power.FWE[1],powers$power.RFT[1])

	for(i in 1:5){lines(c(subjects,subjects),c(0,T[i]),col=cols[i])} # volle lijnen van huidige studie
	for(i in 1:5){points(subjects,T[i],pch=16,col=cols[i])} # dikke bollen op huidige studie
	for(i in 1:5){lines(c(newsamplesize[i],newsamplesize[i]),c(0,newpower[i]),col=cols[i],lty=2)} # volle vertikale lijnen van nieuwe studie
	for(i in 1:5){lines(c(min(newsubjects),newsamplesize[i]),c(newpower[i],newpower[i]),col=cols[i],lty=2)} # volle horizontale lijnen van nieuwe studie
	for(i in 1:5){points(newsamplesize[i],newpower[i],pch=16,col=cols[i])} # dikke bollen op huidige studie
	for(i in 1:5){text(newsamplesize[i],-0.05,newsamplesize[i],col=cols[i])}


	par(mar=c(0,0,0,0))
	plot(1:3, 1:3, col="white",axes=FALSE,xlab="",ylab="")
	legend(1,2.8,c("FDR (BH)","aFDR (Q)","UN","FWER","RFT"),col=cols[1:5],lty=rep(1,5),lwd=2,box.lwd=0,bty="n",border="white",cex=1.5)

	par(mar=c(0,0,0,0))
	plot(1:3, 1:3, col="white",axes=FALSE,xlab="",ylab="")
	text(1,2.5,pos=4,labels=paste("To obtain a power level of ",power,"with ",MCP,"\n control at level",alpha,", the minimal sample size is",samplesize,"."),cex=1.5)
}





   
  
  