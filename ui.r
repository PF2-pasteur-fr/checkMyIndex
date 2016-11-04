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
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("minRedGreen")}),
                      checkboxInput("completeLane", "Directly look for a solution with the desired multiplexing rate", value=FALSE),
                      checkboxInput("selectCompIndexes", "Select compatible indexes before looking for a solution", value=FALSE),
                      selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:4), selected=10),
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
                                 p(textOutput("textNbCompCombinations")),
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
                                 
                                 h3("How the algorithm works"),
                                 p("There can be many combinations of indexes to check according to the number of input indexes and the multiplexing rate. Thus, testing for 
                                    the compatibility of all the combinations may be long or even impossible. The trick is to find a partial solution with the desired number 
                                    of lanes but with a lower multiplexing rate and then to complete each lane with some of the remaining indexes to reach the desired multiplexing 
                                    rate (if N indexes are compatibles then the same N indexes + 1 are compatible). Briefly, a lower multiplexing rate generates a lower number of 
                                    combinations to test and thus make the research of a partial solution very fast. Adding some indexes to each lane to reach the speficied
                                    multiplexing rate is fast too."),
                                 p("Unfortunately, the research of a final solution might become impossible as the astuteness reduces the number of combinations of indexes.
                                    Note that one can look for a solution using directly the desired multiplexing rate (see parameters), the only risk being to increase 
                                    the computational time."),
                                 
                                 h3("Parameters"),
                                 p(strong("Total number of samples"), "in your experiment (can be greater than the number of available indexes)."),
                                 p(strong("Multiplexing rate"), "i.e. number of samples per lane (only divisors of the total number of samples are proposed)."),
                                 p(strong("Constraint on the indexes"), "to avoid having two samples or two lanes with the same index(es)."),
                                 p(strong("Minimal number of red and green lights"), "required at each position is equal to 1 by default to have compatible indexes but can be increased."),
                                 p(strong("Directly look for a solution with the desired multiplexing rate"), "instead of looking for a solution with a few samples per lane and add 
                                           some of the remaining indexes to reach the desired multiplexing rate."),
                                 p(strong("Select compatible indexes"), "before looking for a solution can take some time but then speed up the algorithm. 
                                           It also gives the possibility to get the number of compatible combinations of indexes if the algorithm directly 
                                           look for a solution with the desired multiplexing rate."),
                                 p(strong("Maximum number of trials"), "can be increased if a solution is difficult to find with the parameters chosen."),

                                 h3("About"),
                                 p("This application has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet. Feel free to send an e-mail to", 
                                   a("hugo.varet@pasteur.fr"), "for any suggestion or bug report."),
                                 p("Source code and instructions to run it locally are available on", a("GitHub", href="https://github.com/PF2-pasteur-fr/checkMyIndex"), "."))
                        
                      )
                    )
                  )
))
