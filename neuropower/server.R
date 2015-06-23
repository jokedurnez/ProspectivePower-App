# read in libraries

library(shiny)
library(BioNet)
library(oro.nifti)
library(qvalue)
library(AnalyzeFMRI)
source("neuropower.R")

# start the reactive output

shinyServer(function(input,output){
  
  ##########################
  # INPUT and calculations #
  ##########################
  
  # read in data and extract peaks
  
	data <- reactive({
		if(input$Compute1==0) {return()}
    	# read in data with oro.nifti
		folder <- input$MapName$datapath
		folder.new <- paste(folder,".nii",sep="")
		true <- file.rename(folder,folder.new)
		data <- readNIfTI(folder.new,reorient=FALSE)@.Data
		mask <- ifelse(is.na(data) | data==0,0,1)	

    	# define other specified parameters
		sub <- as.numeric(input$Subjects) # number of subjects
		df <- ifelse(input$OneTwoSample=="One-Sample",sub-1,sub-2) # degrees of freedom (depending on 1- or 2-sample test)
		if(input$torp==2){
		  u <- as.numeric(input$PeakThres)
		} else {
		  u <- -qt(as.numeric(input$PeakThres),df)
		}
		if(input$Smoothchoice == 1){
		  sigma <- SmoothEst(data,mask=mask,voxdim=c(1,1,1))
		  FWHM_voxel <- diag(sqrt(sigma*(8*log(2))))
		} else {
		  help <- sub("\\]","",sub("\\[","",input$FWHMmm))
		  FWHM_mm <- as.numeric(unlist(strsplit(help,",")))
		  help2 <- sub("\\]","",sub("\\[","",input$vsize))
		  FWHM_voxel <- FWHM_mm*as.numeric(unlist(strsplit(help2,",")))
		}
		resels <- sum(mask)/prod(FWHM_voxel)

		# compute clusters and peaks with progression bar
		withProgress(message = 'Calculation in progress',
				   detail = paste("Percentage finished:",0), value = 0,expr={
			peaks <- cluster(data,u)})
		}  
		 
		# adjust peak list
		peaks <- peaks[peaks$peaks>u,]
		peaks$pvalue <- exp(-u*(peaks$peaks-u)) 
		if(input$TorZ == "T"){peaks$pvalue <- -qnorm(pt(-peaks$pvalue,df))
		peaks$pvalue[peaks$pvalue==0] <- 10^(-6)
		peaks$pvalue[peaks$pvalue==1] <- 1-10^(-6)

		out <- list()
		out$peaks <- peaks
		out$TorZ <- input$TorZ
		out$subs <- sub
		out$u <- u
		out$df <- df
		out$resels <- resels

		out
	})

	# estimate model
	
	np.est <- reactive({
		if(input$Compute2==0){return()}
		npest <- NPestimate(data()$peaks,data()$TorZ,data()$u,data()$df,data()$subs,plot=FALSE)
		plot <- NPplot(npest,data()$peaks)
		out <- list()
		out$npest <- npest
		out$plot <- plot

		out
	})

	# estimate post-hoc power
	
	np.posthoc <- reactive({
		if(input$Compute3==0){return()}
		npposthoc <- NPposthoc(np.est()$npest,data()$subs,input$MCP,data()$u,as.numeric(input$alpha),data()$resels)

		npposthoc
	})

	# estimate samplesize
	
	np.samplesize <- reactive({
		if(input$Compute4==0){return()}
		npss <- NPsamplesize(np.est()$npest,data()$subs,50,input$MCP,data()$u,as.numeric(input$alpha),as.numeric(input$power),data()$resels,plot=TRUE)

		npss
	})

	
  ##########
  # OUTPUT #
  ##########
  
  	# table with peaks
	output$nTable <- renderTable({
		input$Compute1
		m <- data()$peaks
		m
	})

	# fit plot
	output$NPplot <- renderPlot({
		input$Compute2
		m <- np.est()$plot
		m
	})

	# post-hoc power
	output$Posthoc <- renderText({
		input$Compute3
		m <- np.posthoc()
		m
	})

	# sample size calculations
	output$power <- renderPlot({
		input$Compute5
		n <- np.samplesize()$plotje
		n
	})
})

 