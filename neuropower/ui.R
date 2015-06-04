library(shiny)
library(shinyIncubator)

shinyUI(fluidPage(
  titlePanel(h1("NeuroPower")),
  sidebarLayout(
  sidebarPanel(
    tabsetPanel(
      tabPanel("1. Load data, extract peaks and compute peak p-values",
                 br(),
                 fileInput("MapName","Select your statistical parametric map for a certain contrast (T or Z) in nifti format (.nii).",multiple=FALSE),
                 selectInput("TorZ","Are the values Z- or T-values?",choices=c("Z","T")),
                 textInput("PeakThres",label="Type here a peakforming threshold",value="2"),
               textInput("Subjects",label="How many subjects?",value="10"),
               selectInput("OneTwoSample",label="Is the study a one- or two-sample test?",choices=c("One-sample","Two-sample")),
               
               
                  actionButton("Compute1",label="Compute")
    ),
 tabPanel("2. Compute power",
          h3("Fit mixture model to peaks."),
           actionButton("Compute2",label="Compute"),
          h3("Compute post-hoc power"),
          selectInput("MCP",label="Which multiple testing procedure has been used?",choices=list("Benjamini-Hochberg"="BH","Bonferroni"="FWE","Uncorrected"="UN","Storey (Q-value)"="Q")),
          textInput("alpha",label="At which significance level was this control?",value="0.05"),
          actionButton("Compute3",label="Compute"),   
          h3("Compute minimal sample size for a certain power."),
          sliderInput("power", "What level of power do you want to obtain?", min = 0, max = 1, value = 0.8, step= 0.05),
          actionButton("Compute4",label="Compute")
 )  
 
    )),
    
    mainPanel(
        tabsetPanel(
        tabPanel("Peaks",tableOutput("nTable")),
        tabPanel("Fit",plotOutput("NPplot")),
        tabPanel("Post-hoc",textOutput("Posthoc")),
        tabPanel("Power",plotOutput("power"))
        )
    )
)))
  