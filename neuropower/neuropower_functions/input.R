##############################
# NEUROPOWER: FUNCTION INPUT #
##############################

# Function: neuropower.input
# ---------------------------
# This function handles the input and converts it to peak table and other values.
#
# written by Joke Durnez (joke.durnez@gmail.com)
#

neuropower.input <- function(input){
  # read in data with oro.nifti
  file <- paste(input$MapName$datapath[1],sep="")
  folder.new <- paste(file,".nii",sep="")
  true <- file.rename(file,folder.new)
  data <- readNIfTI(folder.new,reorient=FALSE)@.Data
  mask <- ifelse(is.na(data) | data==0,0,1)

  # number of subjects and degrees of freedom
  n <- as.numeric(input$Subjects) # number of subjects
  df <- ifelse(input$OneTwoSample=="One-Sample",n-1,n-2)

  # compute excursion threshold
  u <- ifelse(input$torp==2,as.numeric(input$PeakThres),-qt(as.numeric(input$PeakThres),df))

  # compute clusters and peaks with progression bar
  withProgress(message = 'Calculation in progress',
         detail = paste("Percentage finished:",0), value = 0,expr={
    peaks <- cluster(data,u)
  })

  # compute or estimate smoothness and calculate resels
  if(input$Smoothchoice == 1){
    sigma <- SmoothEst(data,mask=mask,voxdim=c(1,1,1))
    FWHM_voxel <- diag(sqrt(sigma*(8*log(2))))
  } else {
    FWHM.stripped <- sub("\\]","",sub("\\[","",input$FWHMmm))
    FWHM_mm <- as.numeric(unlist(strsplit(FWHM.stripped,",")))
    vsize.stripped <- sub("\\]","",sub("\\[","",input$vsize))
    FWHM_voxel <- FWHM_mm*as.numeric(unlist(strsplit(vsize.stripped,",")))
  }
  resels <- sum(mask)/prod(FWHM_voxel)


  # adjust peak list: compute P, transform to z-values, change p-values = 0 or 1
  peaks$pvalue <- exp(-u*(peaks$peaks-u))
  if(input$TorZ == "T"){peaks$pvalue <- -qnorm(pt(-peaks$pvalue,df))}
  peaks$pvalue[peaks$pvalue==0] <- 10^(-6)
  peaks$pvalue[peaks$pvalue==1] <- 1-10^(-6)

  # create output
  out <- list(peaks,n,u,df,resels)
  names(out) <- c("peaks","n","u","df","resels")
  return(out)
}
