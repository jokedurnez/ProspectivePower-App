library(shiny)
library(shinyIncubator)

shinyUI(fluidPage(
  titlePanel("NeuroPower"),
  sidebarLayout(
  sidebarPanel(
    tabsetPanel(
      tabPanel("1. Load data",
               br(),
          fileInput("MapName","Select your statistical parametric map for a certain contrast (T or Z) in nifti format (.nii).",multiple=FALSE),
          selectInput("TorZ","Are the values Z- or T-values?",choices=c("Z","T")),
          textInput("PeakThres",label="Type here a peakforming threshold",value="2"),
          textInput("Subjects",label="How many subjects?",value="10"),
          selectInput("OneTwoSample",label="Is the study a one- or two-sample test?",choices=c("One-sample","Two-sample")),
          actionButton("Compute1",label="Extract peaks")
    ),  

    tabPanel("2. Estimate model",
          h4("Fit mixture model to peaks."),
            actionButton("Compute2",label="Estimate")
 ),
 tabPanel("3. Power",
          br(),
          h4("Post-hoc power"),
          selectInput("MCP",label="Which multiple testing procedure has been used?",choices=list("Benjamini-Hochberg"="BH","Bonferroni (peaks)"="FWE","Uncorrected"="UN","Storey (Q-value)"="Q","Random Field Theory"="RFT")),
          textInput("alpha",label="At which significance level was this control?",value="0.05"),
          actionButton("Compute3",label="Compute"),   
          h4("Sample size calculation."),
          sliderInput("power", "What level of power do you want to obtain?", min = 0, max = 1, value = 0.8, step= 0.05),
          actionButton("Compute4",label="Compute"),
          actionButton("Compute4",label="Display Plot")
 )  
 
    )),
    
 
 mainPanel(
        tabsetPanel(
        tabPanel("Welcome",
        h2("Introduction"),
        p("This power calculation tool was first introduced in this",
          a("OHBM 2015 poster",href = "http://users.ugent.be/~jdurnez/ProsPeakPow_1406_OHBM.pdf"),"."),
        p("This app applies functions to compute post-hoc power in a fMRI pilot study  as described in",a("Durnez, Moerkerke & Nichols, 2014",href="http://www.ncbi.nlm.nih.gov/pubmed/23927901"),"and to compute the minimal sample size needed to obtain a preferred power rate in a subsequent fMRI study based on that pilot study."),
        h2("Use"),
        p("In a first step, you can upload your statistical map (z-map or t-map) to extract peak information.  Note that the data are NOT stored on any server once the peak information is extracted.  This step usually takes less than 1 minute.  The result of this step can be seen in the tab",strong("Peaks"),", where a table will be displayed with all local maxima above the preferred excursion threshold.  The displayed x-,y- and z-coordinates are in voxel space."),
        p("Subsequently, the mixture model of inactive and active peaks is fit to the data.  This model fitting includes an estimation of both the prevalence of active peaks, the mean (delta) and the standard deviation (sigma) of the distribution of active peaks.  The fit of the model can be inpected in the tab",strong("Model fit"),".  The left hand panel will show a histogram of the peak p-values with their estimated distribution, and the right hand panel will show a histogram of the peak heights (t or z) with their estimated null, alternative and combined distributions."),
        p("Finally, the power can be calculated.  When entering the preferred thresholding procedure and level of significance, the post-hoc power will be estimated for the current data (result desplayed in tab",strong("Post-hoc"),",as well as the minimal sample size for a certain level of power in a new dataset.")),
        tabPanel("Peaks",br(),tableOutput("nTable")),
        tabPanel("Model fit",plotOutput("NPplot")),
        tabPanel("Post-hoc",br(),textOutput("Posthoc")),
        tabPanel("Power",plotOutput("power"))
        )
    )
)))
  