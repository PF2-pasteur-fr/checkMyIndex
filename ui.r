library(shiny)

shinyUI(fluidPage(theme = "bootstrap.min.css",
                  
                  titlePanel("Search for a set of compatible indexes for your sequencing experiment"),
                  
                  sidebarLayout(
                    
                    # parameters
                    sidebarPanel(
                      fileInput("inputFile", label="Select your tab-delimited file containing the index ids and sequences", accept="text"),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("nbSamples")}),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("multiplexingRate")}),
                      selectInput("unicityConstraint", label="Constraint on the indexes", choices=c("None","Use each combination only once","Use each index only once"), selected="None"),
                      selectInput("advancedOptions", label="Optional parameters", choices=c("Hide","Show"), selected="Hide"),
                      conditionalPanel(condition="input.advancedOptions=='Show' & output.indexUploaded", {uiOutput("minRedGreen")}),
                      conditionalPanel(condition="input.advancedOptions=='Show'", {selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:4), selected=10)}),
                      actionButton("go", label="Search for a solution")
                    ),
                    
                    # output
                    mainPanel(
                      tabsetPanel(
                        # 1st panel: input
                        tabPanel("Input indexes", dataTableOutput("inputIndex")),
                        
                        # 2nd panel: results
                        tabPanel("Proposed flowcell design",
                                 p(textOutput("textNbCombinations")),
                                 p(textOutput("textDescribingSolution")),
                                 dataTableOutput("solution"),
                                 p(textOutput("textDescribingMinRedGreen")),
                                 br(""),
                                 downloadButton("downloadData", "Download")),
                        
                        # 3rd panel: help
                        tabPanel("Help",
                                 
                                 h3("Input indexes file"),
                                 p("The user must provide the list of its available indexes as a two-column tab-delimited text file (without header). 
                                   Index ids are in the first column and the corresponding sequences in the second. An example of such a file is available ", 
                                   a("here", href="https://github.com/PF2-pasteur-fr/checkMyIndex/blob/master/inputIndexesExample.txt")," to test the application."),
                                 
                                 h3("Mandatory parameters"),
                                 p(strong("Total number of samples"), "in your experiment (can be greater than the number of available indexes)."),
                                 p(strong("Multiplexing rate"), "i.e. number of samples per lane (only divisors of the total number of samples are proposed)."),
                                 p(strong("Constraint on the indexes"), "to avoid having two samples or two lanes with the same index(es)."),
                                 h3("Optional parameters"),
                                 p(strong("Minimal number of red/green lights"), "required at each position is equal to 1 by default to have compatible indexes but can be increased."),
                                 p(strong("Maximum number of trials"), "can be increased if a solution is difficult to find with the parameters chosen."),
                                 
                                 h3("About"),
                                 p("This application has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet. Feel free to send an e-mail to", 
                                   a("hugo.varet@pasteur.fr"), "for any suggestion or bug report."),
                                 p("Source code and instructions to run it locally are available on", a("GitHub", href="https://github.com/PF2-pasteur-fr/checkMyIndex"), "."))
                        
                      )
                    )
                  )
))
