######################################
# NEUROPOWER: POSTHOC POWER FUNCTION #
######################################

# Function: power.posthoc
# ---------------------------
# This function computes the posthoc power as described in Durnez et al. (2014)
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
#   - p: pvalues of peaks
#   - resels: resel count
#   - u: the excursion threshold
#   - alpha: level of false positive control
#   - n: samplesize
#   - MCP: multiple comparison control:
#       o "UN": uncorrected
#       o "BF": bonferroni
#       o "RFT": random field Theory
#       o "BH": FDR control (Benjamini-Hochberg)
#   - estimates: list with values:
#       o estimates$pi0: estimates pi0
#       o estimates$Ha$Delta: estimated mu
#       o estimates$Ha$Sigma: estimated sigma

power.posthoc <- function(p,resels,u,alpha,n,MCP,estimates){
	options(warn=-1)
	# compute cutoff
	cutoffs <- cutoff(p,resels,u,alpha)

	# power computed with CDF of alternative distribution
	power.UN <- 1-alt.cdf(cutoffs$UN,estimates$Ha,u)
	power.BF <- 1-alt.cdf(cutoffs$BF,estimates$Ha,u)
	power.RFT <- 1-alt.cdf(cutoffs$RFT,estimates$Ha,u)
	power.BH <- 1-alt.cdf(cutoffs$BH,estimates$Ha,u)

	# posthoc power for requested MCP
	if(MCP == "UN"){power <- power.UN
	} else if (MCP == "BF"){power <- power.BF
	} else if (MCP == "RFT"){power <- power.RFT
	} else if (MCP == "BH"){power <- power.BH
	} else {cat("Unknown MCP"); stop}

	power <- ifelse(ifelse(is.na(power),0,power)>1,1,power)
	powertext <- paste("\n The power of this study with",n,"subjects with",MCP,"control at level",alpha,"equals",round(power,digits=3),": \n the average chance of detecting an activated peak is",100*round(power,digits=3),"%.")
	return(powertext)

}
