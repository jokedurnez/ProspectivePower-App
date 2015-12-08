#########################################
# NEUROPOWER: PREDICTIVE POWER FUNCTION #
#########################################

# Function: power.predictive
# ---------------------------
# This function computes the predicted power.
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
#   - p: pvalues of peaks
#   - resels: resel count
#   - u: the excursion threshold
#   - alpha: level of false positive control
#   - n: samplesize
#   - nmax: what is the maximum sample size
#   - MCP: multiple comparison control:
#       o "UN": uncorrected
#       o "BF": bonferroni
#       o "RFT": random field Theory
#       o "BH": FDR control (Benjamini-Hochberg)
#   - estimates: list with values:
#       o estimates$pi0: estimates pi0
#       o estimates$Ha$Delta: estimated mu
#       o estimates$Ha$Sigma: estimated sigma
#   - plot: Boolean to show power curves

power.predictive <- function(p,resels,u,alpha,n,nmax,MCP,estimates,power,plot){
  # compute cutoffs
  cutoffs <- cutoff(p,resels,u,alpha)

  # compute standardised effect size
  cohenD <- estimates$Ha[1]/sqrt(n)

  # compute new power for a range of number of subjects
	power.BH <- power.UN <- power.BF <- power.RFT <- c()
	newsubjects <- seq(from=n,to=nmax)
	for(p in 1:length(newsubjects)){
		proj.eff <- cohenD*sqrt(newsubjects[p])
		power.UN[p] <- 1-alt.cdf(cutoffs$UN,c(proj.eff,estimates$Ha[2]),u)
		power.BF[p] <- 1-alt.cdf(cutoffs$BF,c(proj.eff,estimates$Ha[2]),u)
		power.RFT[p] <- 1-alt.cdf(cutoffs$RFT,c(proj.eff,estimates$Ha[2]),u)
    power.BH[p] <- 1-alt.cdf(cutoffs$BH,c(proj.eff,estimates$Ha[2]),u)
	}

  powers <- data.frame(power.UN,power.BF,power.RFT,power.BH)
  powers[powers>1] <- 1

  # compute minimal samplesize for required power
	newsamplesize <- c(newsubjects[min(which(powers$power.UN>power))],
                    newsubjects[min(which(powers$power.BF>power))],
                    newsubjects[min(which(powers$power.RFT>power))],
                    newsubjects[min(which(powers$power.BH>power))])
	newpower <- c(power.UN[min(which(powers$power.UN>power))],
                power.BF[min(which(powers$power.BF>power))],
                power.RFT[min(which(powers$power.RFT>power))],
                power.BH[min(which(powers$power.BH>power))]) # EXACT power larger than required power

  if(MCP == "UN"){samplesize <- newsamplesize[1]
	} else if (MCP == "BF"){samplesize <- newsamplesize[2]
	} else if (MCP == "RFT"){samplesize <- newsamplesize[3]
	} else if (MCP == "BH"){samplesize <- newsamplesize[4]
	} else {cat("Unknown MCP"); stop}

  if(plot==TRUE){powercurve <- power.plot(newsubjects,powers,newsamplesize,newpower,power,n,MCP,alpha,samplesize)}

  return(powercurve)
}
