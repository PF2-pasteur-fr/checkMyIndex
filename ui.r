library(shiny)

shinyUI(fluidPage(theme = "bootstrap.min.css",
  
  titlePanel("Search for a set of compatible indexes for your sequencing experiment"),
  
  sidebarLayout(
    
    # parameters
    sidebarPanel(
      fileInput("inputFile", label="Select your tab-delimited file containing the index ids and sequences", accept="text"),
      conditionalPanel(condition="output.indexUploaded", {uiOutput("nbSamples")}),
      conditionalPanel(condition="output.indexUploaded", {uiOutput("multiplexingRate")}),
      checkboxInput("uniqueCombinations", "Use each combination only once", value=TRUE),
      conditionalPanel(condition="input.uniqueCombinations", {checkboxInput("uniqueIndexes", "Use each index only once", value=FALSE)}),
      selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:5), selected=10),
      actionButton("go", label="Search for a solution"),
      br(),br(),
      p("Contact: "), a("hugo.varet@pasteur.fr")
    ),
    
    # output
    mainPanel(
      tabsetPanel(
        # 1st panel: input
        tabPanel("Input indexes", dataTableOutput("inputIndex")),
        
        # 2nd panel: results
        tabPanel("Proposed flowcell design",
                 p(textOutput("textNbCompatibleIndexes")),
                 p(textOutput("textDescribingSolution")),
                 dataTableOutput("solution"),
                 br(""),
                 downloadButton("downloadData", "Download"))
      )
    )
  )
))
