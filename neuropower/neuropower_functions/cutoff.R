################################
# NEUROPOWER: CUT-OFF FUNCTION #
################################

# Function: cutoff
# ---------------------------
# This function computes the cutoff for different Multiple Comparison Procedures
#
# written by Joke Durnez (joke.durnez@gmail.com)
#
#  inputs:
#   - p: pvalues of peaks
#   - u: the excursion rangeold
#   - resels: resel count
#   - alpha: level of false positive control

cutoff <- function(p,resels,u,alpha){
  range <- seq(from=u,to=15,length=100)
  # UN and BF
  cdfN <- exp(-u*(range-u))
  UN <- range[min(which(cdfN<alpha))]
  BF <- range[min(which(cdfN<(alpha/length(p))))]
  # RFT
  cdfN_RFT <- resels*exp(-range^2/2)*range^2
  RFT <- range[min(which(cdfN_RFT<alpha))]
  # FDR cutoff estimation
  p.order <- rank(p)
  FDRqval <- (p.order/length(p))*alpha
  p.sig <- ifelse(p < FDRqval,1,0)
  FDRc <- ifelse(sum(p.sig)==0,0,max(p[p.sig==1]))
  BH <- ifelse(FDRc==0,NA,range[min(which(cdfN<FDRc))])
  cutoffs <- list(UN,BF,RFT,BH)
  names(cutoffs) <- c("UN","BF","RFT","BH")
  return(cutoffs)
}
