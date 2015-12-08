# read in libraries

library(shiny)
library(BioNet)
library(oro.nifti)
library(qvalue)
library(AnalyzeFMRI)
source("neuropower_functions/cluster.R")
source("neuropower_functions/cutoff.R")
source("neuropower_functions/densities.R")
source("neuropower_functions/input.R")
source("neuropower_functions/modelestimator.R")
source("neuropower_functions/modelfitplot.R")
source("neuropower_functions/powerplot.R")
source("neuropower_functions/powerposthoc.R")
source("neuropower_functions/powerpredictive.R")

# start the reactive output

shinyServer(function(input,output){

  ############
  # 1. INPUT #
  ############

	Data <- reactive({
		if(input$Compute1==0) {return()}
    data <- neuropower.input(input)
    return(data)
		})

  # OUTPUT TAB 1
	output$PeakTable <- renderTable({
		input$Compute1
		m <- Data()$peaks
		m
	})

  #####################
  # 2. ESTIMATE MODEL #
  #####################

  ModelEst <- reactive({
		if(input$Compute2==0){return()}
		est <- modelestimator(Data()$peaks$peaks,
                Data()$peaks$pvalue,
                Data()$u,
                Data()$df,
                Data()$n,
                plot=FALSE
                )
		plot <- model.plot(Data()$peaks$peaks,
                Data()$peaks$pvalue,
                est$pi0,
								est$a,
                est$Ha,
                Data()$u,
                Data()$n
                )
		out <- list(est,plot)
    names(out) <- c("est","plot")
		return(out)
	})

  # OUTPUT TAB 2
  output$ModelPlot <- renderPlot({
		input$Compute2
		m <- ModelEst()$plot
		m
	})

  #####################
  # 3. POST HOC POWER #
  #####################

	PowerPostHoc <- reactive({
		if(input$Compute3==0){return()}
		powerposthoc <- power.posthoc(
      Data()$peaks$pvalue,
      Data()$resels,
      Data()$u,
      as.numeric(input$alpha),
      Data()$n,
      input$MCP,
      ModelEst()$est
      )

		return(powerposthoc)
	})

  #output TAB 3

  output$Posthoc <- renderText({
		input$Compute3
		m <- PowerPostHoc()
		m
	})

  ########################
  # 3. PROSPECTIVE POWER #
  ########################

	ProspectivePower <- reactive({
		if(input$Compute4==0){return()}
		prospow <- power.predictive(
      Data()$peaks$pvalue,
      Data()$resels,
      Data()$u,
      as.numeric(input$alpha),
      Data()$n,
      50,
      input$MCP,
      ModelEst()$est,
      input$power,
      plot=TRUE
    )
		return(prospow)
	})

  #output tab 4: power curves

	output$Power <- renderPlot({
		input$Compute5
		n <- ProspectivePower()$powercurve
		n
	})
})
