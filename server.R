library(shiny)
library(BioNet)
library(oro.nifti)
library(qvalue)
source("neuropower.R")


shinyServer(
  function(input,output){
    data <- reactive({
      if(input$Compute1==0) {return()}
      folder <- input$MapName$datapath
      folder.new <- paste(folder,".nii",sep="")
      true <- file.rename(folder,folder.new)
      data <- readNIfTI(folder.new,reorient=FALSE)@.Data
      withProgress(message = 'Calculation in progress',
                   detail = paste("Percentage finished:",0), value = 0,expr={
      peaks <- cluster(data)})
      sub <- as.numeric(input$Subjects)
      df <- ifelse(input$OneTwoSample=="One-Sample",sub-1,sub-2)
      u <- as.numeric(input$PeakThres)
      peaks <- peaks[peaks$peaks>u,]
      peaks$pvalue <- exp(-u*(peaks$peaks-u)) 
      if(input$TorZ == "T"){peaks$pvalue <- -qnorm(pt(-peaks$pvalue,df))}   
      peaks$pvalue[peaks$pvalue==0] <- 10^(-6)
      peaks$pvalue[peaks$pvalue==1] <- 1-10^(-6)
      
      out <- list()
      out$peaks <- peaks
      out$TorZ <- input$TorZ
      out$subs <- sub
      out$u <- u
      out$df <- df
      
     out
      })
    
    np.est <- reactive({
      if(input$Compute2==0){return()}
      npest <- NPestimate(data()$peaks,data()$TorZ,data()$u,data()$df,plot=FALSE)
      plot <- NPplot(npest,data()$peaks)
      out <- list()
      out$npest <- npest
      out$plot <- plot
      
      out
    })
    
    np.posthoc <- reactive({
      if(input$Compute3==0){return()}
      npposthoc <- NPposthoc(np.est()$npest,data()$subs,input$MCP,data()$u,as.numeric(input$alpha))
      
      npposthoc
    })
    
    np.samplesize <- reactive({
      if(input$Compute4==0){return()}
      npss <- NPsamplesize(np.est()$npest,data()$subs,100,input$MCP,data()$u,as.numeric(input$alpha),as.numeric(input$power),plot=TRUE)
      
      npss
      
    })
    
     output$nTable <- renderTable({
       input$Compute1
       m <- data()$peaks
       m
     })
    
    output$NPplot <- renderPlot({
      input$Compute2
      m <- np.est()$plot
      m
    })
    
    output$Posthoc <- renderText({
      input$Compute3
      m <- np.posthoc()
      m
    })
    
     output$power <- renderPlot({
      input$Compute4
      m <- np.samplesize()$plotje
      m
    })
})

 