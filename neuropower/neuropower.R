#DEPENDENCIES: BioNet, AnalyzeFMRI, QVALUE


cluster <- function(zmap,u){
  options(warn=-1)
  cat("Searching for peaks")
#  zmap <- ifelse(length(dim(zmap))==4,zmap[,,,1],zmap)
  dimZ <- dim(zmap)
  Zstar <- array(0,dim=dimZ+c(2,2,2))  
  Zstar[2:(dimZ[1]+1),2:(dimZ[2]+1),2:(dimZ[3]+1)] <- zmap
  Zstar[is.na(Zstar)] <- -100
  peaks <- xco <- yco <- zco <- c()  
  locmax <- array(NA,dim=dimZ+c(2,2,2))
  for(m in 2:(dim(Zstar)[1]-1)){
      incProgress(1/m,detail=paste("Percentage finished:",100*round(m/(dim(Zstar)[1]-1),digits=2),"%"))
     for(n in 2:(dim(Zstar)[2]-1)){
      for(o in 2:(dim(Zstar)[3]-1)){
        if(Zstar[m,n,o]<u){next}
        surroundings <- c(Zstar[m-1,n+1,o],
                          Zstar[m,n+1,o],
                          Zstar[m+1,n+1,o],
                          Zstar[m-1,n,o],
                          Zstar[m+1,n,o],
                          Zstar[m-1,n-1,o],
                          Zstar[m,n-1,o],
                          Zstar[m+1,n-1,o],
                          
                          Zstar[m,n+1,o-1],
                          Zstar[m-1,n,o-1],
                          Zstar[m+1,n,o-1],
                          Zstar[m,n-1,o-1],
                          Zstar[m+1,n+1,o-1],
                          Zstar[m+1,n-1,o-1],
                          Zstar[m-1,n+1,o-1],
                          Zstar[m-1,n-1,o-1],
                          
                          Zstar[m,n+1,o+1],
                          Zstar[m-1,n,o+1],
                          Zstar[m+1,n,o+1],
                          Zstar[m,n-1,o+1],
                          Zstar[m+1,n+1,o+1],
                          Zstar[m+1,n-1,o+1],
                          Zstar[m-1,n+1,o+1],
                          Zstar[m-1,n-1,o+1])
        
        locmax[m,n,o] <- ifelse(Zstar[m,n,o]>max(surroundings),Zstar[m,n,o],0)        
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
  
sumlogtruncdensBIS <- function(x,par,pi0,a){
  mu.2 <- par[1]; sigma.2 <- par[2]
  x2 <- x-a
  f.x.0 <- a*exp(-a*(x2))
  
  f.x.num.1 <- 1/sigma.2 * dnorm((x-mu.2)/sigma.2)
  F.x.den.1 <- 1-pnorm((a-mu.2)/sigma.2)
  
  f.x <- pi0*f.x.0 + (1-pi0)*(f.x.num.1/F.x.den.1)
  return(-sum(log(f.x)))
}

truncdensBIS <- function(x,par,pi0,a){
  mu.2 <- par[1]; sigma.2 <- par[2]
  x2 <- x-a
  f.x.0 <- a*exp(-a*(x2))
  
  f.x.num.1 <- 1/sigma.2 * dnorm((x-mu.2)/sigma.2)
  F.x.den.1 <- 1-pnorm((a-mu.2)/sigma.2)
  
  f.x <- pi0*f.x.0 + (1-pi0)*f.x.num.1/F.x.den.1
  return(f.x)
}

truncdensNUL <- function(x,a){
  #PDF = u*exp(-u*(x-u))
  #CDF = 1-exp(-u*(x-u))
  f.x.0 <- a*exp(-a*(x-a))
  return(f.x.0)
}


CDFtrunc <- function(x,mu,sigma,a){
  ksi <- (x-mu)/sigma
  alpha <- (a-mu)/sigma
  F.x <- (pnorm(ksi)-pnorm(alpha))/(1-pnorm(alpha))
  return(F.x)
}


NPplot <- function(estimates,peaklist){
    options(warn=-1)
    pi0e <- estimates[[1]]
    mixture <- estimates[[2]]
    u <- estimates$u
    par(mfrow=c(1,2),mar=c(4,5,5,1))
    cols <- c("#A6CEE3","#1F78B4","#41AB5D","#006D2C")
    breakshist <- seq(from=0,to=1,by=0.1)
    histo <- hist(peaklist$pvalue,freq=FALSE,col=cols[1],border=cols[1],breaks=breakshist,xlab="peak p-values",main=paste("Distribution of",length(peaklist$pvalue),"peak p-values"))
    xn <- seq(from=0,to=1,length=1000)
    lines(xn,rep(pi0e,1000),col=cols[3],lwd=3)
    lines(xn,(dbeta(xn,estimates$a,1)*(1-pi0e)+pi0e),col=cols[2],lwd=3)
    mtext(bquote(pi[1] ~ "=" ~ .(round(1-pi0e,digits=2))))
    legend(0.3,max(histo$density),c("Estimated null","Estimated total"),fill=cols[c(3,2)],border=cols[c(3,2)],box.lwd=0,box.col="white")

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
    
NPestimate <- function(peaklist,STAT,u,df,subs,plot){
  options(warn=-1)
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
  	mixture <- optim(par=c(5,0.5),method="L-BFGS-B",lower=c(2.5,0.1),fn=sumlogtruncdensBIS,x=peaklist$peaks,pi0=pi0e,a=u)
  estimates <- list(1-pi0e,mixture$par)
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
  
NPposthoc <- function(estimates,subjects,MCP,u,alpha){
  options(warn=-1)
  thresh <- seq(from=u,to=15,length=100)
  cdfN <- exp(-u*(thresh-u))
  ps <- estimates$peaks$pvalue
  pvalms <- sort(ps)
  orderpvalms <- rank(ps)
  FDRqval <- (orderpvalms/length(ps))*alpha
  pr <- ifelse(pvalms[orderpvalms] < FDRqval,1,0)
  FDRc <- 1
  cutoff.BH <- ifelse(FDRc==0,NA,thresh[min(which(cdfN<FDRc))])
  Q <- qvalue(ps,fdr.level=alpha)
  cutoff.Q <- ifelse(!is.list(Q),NA,ifelse(sum(Q$significant)==0,NA,min(estimates$peaks$peaks[Q$significant==TRUE])))
  cutoff.UN <- thresh[min(which(cdfN<alpha))]
  cutoff.FWE <- thresh[min(which(cdfN<(alpha/length(estimates$peaks))))]

  power.BH <- 1-CDFtrunc(cutoff.BH,estimates$Ha[1],estimates$Ha[2],2)
  power.Q <- 1-CDFtrunc(cutoff.Q,estimates$Ha[1],estimates$Ha[2],2)
  power.UN <- 1-CDFtrunc(cutoff.UN,estimates$Ha[1],estimates$Ha[2],2)
  power.FWE <- 1-CDFtrunc(cutoff.FWE,estimates$Ha[1],estimates$Ha[2],2)

  if(MCP == "BH"){
	 	power <- power.BH
  } else if (MCP == "Q"){
	 	power <- power.Q
  } else if (MCP == "FWE"){
	 	power <- power.FWE
  } else if (MCP == "uncorrected"){
	 	power <- power.UN
  } else {cat("Unknown MCP"); stop}
	
  powertext <- paste("The power of this study with",subjects,"subjects with",MCP,"control at level",alpha,"equals",round(power,digits=3),": \n the average chance of detecting an activated peak is",100*round(power,digits=3),"%.")
  powers <- c(power.BH,power.Q,power.UN,power.FWE)
  names(powers) <- "power"
  return(powertext)
 
}

NPsamplesize <- function(estimates,subjects,maxsubjects,MCP,u,alpha,power,plot){
  options(warn=-1)
  
# compute cutoffs in original analysis

  # compute standardised effect
  eff.sta <- estimates$Ha[1]/sqrt(subjects)
  # compute FDR threshold
  thresh <- seq(from=u,to=15,length=100)
  cdfN <- exp(-u*(thresh-u))
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
  # compute threshold for uncorrected and FWE
  cutoff.UN <- thresh[min(which(cdfN<alpha))]
  cutoff.FWE <- thresh[min(which(cdfN<(alpha/length(estimates$peaks))))]
  
# compute new power for a range of number of subjects
  power.BH <- power.Q <- power.UN <- power.FWE <- c()
  newsubjects <- seq(from=subjects,to=maxsubjects)
  for(p in 1:length(newsubjects)){
	  proj.eff <- eff.sta*sqrt(newsubjects[p])
	  newpar <- c(proj.eff,estimates$Ha[2])      
	  power.BH[p] <- 1-CDFtrunc(cutoff.BH,proj.eff,estimates$Ha[2],u)
	  power.Q[p] <- 1-CDFtrunc(cutoff.Q,proj.eff,estimates$Ha[2],u)
	  power.UN[p] <- 1-CDFtrunc(cutoff.UN,proj.eff,estimates$Ha[2],u)
  	  power.FWE[p] <- 1-CDFtrunc(cutoff.FWE,proj.eff,estimates$Ha[2],u)
	
  }

# compute minimal samplesize for required power
  newsamplesize <- c(newsubjects[min(which(power.BH>power))],newsubjects[min(which(power.Q>power))],newsubjects[min(which(power.UN>power))],newsubjects[min(which(power.FWE>power))])
  
# error message (not implemented) when power is impossible to obtain
  if(MCP == "BH"){if(sum(power.BH>power)==0){samplesize <- NA} else {samplesize <- newsamplesize[1]}
  } else if (MCP == "Q"){if(sum(power.Q>power)==0){samplesize <- NA} else {samplesize <- newsamplesize[2]}
  } else if (MCP == "FWE"){if(sum(power.FWE>power)==0){samplesize <- NA} else {samplesize <- newsamplesize[3]}
  } else if (MCP == "UN"){if(sum(power.UN>power)==0){samplesize <- NA} else {samplesize <- newsamplesize[4]}
  } else {stop}
	
  powers <- data.frame(power.BH,power.Q,power.UN,power.FWE)
  newpower <- c(power.BH[min(which(power.BH>power))],power.Q[min(which(power.Q>power))],power.UN[min(which(power.UN>power))],power.FWE[min(which(power.FWE>power))]) # EXACT power larger than required power
  
  tekstje <- paste("\n The minimal samplesize of this study with with",MCP,"control at level",alpha," for a power of",power,":",samplesize,"\n\n")
  
  if(plot==TRUE){plotje <- NPssplot(newsubjects,powers,newsamplesize,newpower,power,subjects)}
  
  res <- list()
  res$plotje <- plotje
  res$tekstje <- tekstje

  return(res)
} 
  
 
NPssplot <- function(newsubjects,powers,newsamplesize,newpower,power,subjects){
  options(warn=-1)
  cols <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3")

 layout(matrix(c(1,2),1,2,byrow=TRUE),widths=c(0.7,0.3))
par(mar=c(4,4,2,1),oma=c(0,0,0,0))

# power curves
plot(newsubjects,powers$power.BH,type="l",col=cols[1],ylim=c(-0.1,1),xlim=c(min(newsubjects),max(newsubjects)),xlab="Subjects",ylab="Average power",lwd=2,main="Power with varying sample size")
 lines(newsubjects,powers$power.Q,col=cols[2],lwd=2) 
 lines(newsubjects,powers$power.UN,col=cols[3],lwd=2)
 lines(newsubjects,powers$power.FWE,col=cols[4],lwd=2)

# power in *this* study
T <- c(powers$power.BH[1],powers$power.Q[1],powers$power.UN[1],powers$power.FWE[1])

for(i in 1:4){lines(c(subjects,subjects),c(0,T[i]),col=cols[i])} # volle lijnen van huidige studie
for(i in 1:4){points(subjects,T[i],pch=16,col=cols[i])} # dikke bollen op huidige studie

for(i in 1:4){lines(c(newsamplesize[i],newsamplesize[i]),c(0,newpower[i]),col=cols[i],lty=2)} # volle vertikale lijnen van nieuwe studie
for(i in 1:4){lines(c(min(newsubjects),newsamplesize[i]),c(newpower[i],newpower[i]),col=cols[i],lty=2)} # volle horizontale lijnen van nieuwe studie
for(i in 1:4){points(newsamplesize[i],newpower[i],pch=16,col=cols[i])} # dikke bollen op huidige studie
for(i in 1:4){text(newsamplesize[i],-0.05,newsamplesize[i],col=cols[i])}


 par(mar=c(0,0,0,0))
plot(1:3, 1:3, col="white",axes=FALSE,xlab="",ylab="")
legend(1,2.8,c("FDR (BH)","aFDR (Q)","UN","FWER"),col=cols[1:4],lty=rep(1,4),lwd=2,box.lwd=0,bty="n",cex=1,border="white")
}





   
  
  