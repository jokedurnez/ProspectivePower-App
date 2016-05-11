library(shiny)
library(shinyIncubator)

shinyUI(fluidRow(
	titlePanel("NeuroPower"),
  fluidRow(
    column(4,
            wellPanel(
              "Hi,",
              br(),
              p("During the past months, we've been working really hard for a better and faster way to perform power analyses for fMRI. The new version can be found on ",
                a("neuropowertools.org",href = "http://www.neuropowertools.org"),"."),
              p("This shiny app is no longer supported or updated.  We very much would like to welcome you to the new application, with lots of new features and better visualisations.")
              )
            )
    )
	)
)
