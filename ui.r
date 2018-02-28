library(shiny)

shinyUI(fluidPage(theme = "bootstrap.min.css",

                  titlePanel(title=div(img(src="dna.png", width=50), strong("Search for a set of compatible indexes for your sequencing experiment")), windowTitle="checkMyIndex"),
                  
                  sidebarLayout(
                    
                    # parameters
                    sidebarPanel(
                      fileInput("inputFile", label="Select your tab-delimited file containing the index ids and sequences", accept="text"),
                      selectizeInput("chemistry", label="Illumina chemistry", 
                                     choices=c("Four-channels (HiSeq & MiSeq)" = 4, "Two-channels (NovaSeq, NextSeq & MiniSeq)" = 2)),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("nbSamples")}),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("multiplexingRate")}),
                      selectInput("unicityConstraint", label="Constraint on the indexes", 
                                  choices=c("None" = "none", "Use each combination only once" = "lane", "Use each index only once" = "index")),
                      checkboxInput("completeLane", "Directly look for a solution with the desired multiplexing rate", value=FALSE),
                      checkboxInput("selectCompIndexes", "Select compatible indexes before looking for a solution", value=FALSE),
                      selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:4)),
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
                                 br(""),
                                 downloadButton("downloadData", "Download")),
                        
                        # 3rd panel: plot results
                        tabPanel("Visualization of the design",
                                 uiOutput("heatmapindex2")),
                        
                        # 4th panel: help
                        tabPanel("Help",
                                 
                                 h3("Input indexes file"),
                                 p("The user must provide the list of its available indexes as a two-column tab-delimited text file (without header). 
                                   Index ids are in the first column and the corresponding sequences in the second. An example of such a file is available ", 
                                   a("here", href="https://github.com/PF2-pasteur-fr/checkMyIndex/blob/master/inputIndexesExample.txt")," to test the application."),
                                 
                                 h3("How the algorithm works"),
                                 p("There can be many combinations of indexes to check according to the number of input indexes and the multiplexing rate. Thus, testing for 
                                    the compatibility of all the combinations may be long or even impossible. The trick is to find a partial solution with the desired number 
                                    of pools/lanes but with fewer samples than asked and then to complete each pool/lane with some of the remaining indexes to reach the desired 
                                    multiplexing rate. Indeed, adding indexes to a combination of compatible indexes will give a compatible combination for sure. Briefly, a lower 
                                    number of samples per pool/lane generates a lower number of combinations to test and thus makes the research of a partial solution very fast. 
                                    Adding some indexes to complete each pool/lane is fast too and gives the final solution."),
                                 p("Unfortunately, the research of a final solution might become impossible as the astuteness reduces the number of combinations of indexes.
                                    In such a case, one can look for a solution using directly the desired multiplexing rate (see parameters), the only risk is to increase 
                                    the computational time."),
                                 
                                 h3("Parameters"),
                                 p(strong("Illumina chemistry"), "can be either four-channels (HiSeq & MiSeq) or two-channels (NovaSeq, NextSeq & MiniSeq). With the four-channel
                                   chemistry A/C are red and G/T are green and indexes are compatible if there are at least one red light and one green light at each position.
                                   With the two-channel chemistry G has no color, A is orange, C is red and T is green and indexes are compatible if there is at least one
                                   color at each position. Note that indexes starting with GG are not compatible with the two-channel chemistry.
                                   Please refer to the Illumina documentation for more details."),
                                 p(strong("Total number of samples"), "in your experiment (can be greater than the number of available indexes)."),
                                 p(strong("Multiplexing rate"), "i.e. number of samples per pool/lane (only divisors of the total number of samples are proposed)."),
                                 p(strong("Constraint on the indexes"), "to avoid having two samples or two pools/lanes with the same index(es)."),
                                 p(strong("Directly look for a solution with the desired multiplexing rate"), "instead of looking for a partial solution with a few samples per pool/lane 
                                           and then add some of the remaining indexes to reach the desired multiplexing rate."),
                                 p(strong("Select compatible indexes"), "before looking for a (partial) solution can take some time but then speed up the algorithm. 
                                           If one looks for a solution directly with the desired multiplexing rate, turning on this options allows to get the number of 
                                           compatible combinations of indexes."),
                                 p(strong("Maximum number of trials"), "can be increased if a solution is difficult to find with the parameters chosen."),

                                 h3("About"),
                                 p("This application has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet. Feel free to send an e-mail to", 
                                   a("hugo.varet@pasteur.fr"), "for any suggestion or bug report."),
                                 p("Source code and instructions to run it locally are available on", a("GitHub", href="https://github.com/PF2-pasteur-fr/checkMyIndex"), ". 
                                   Please note that checkMyIndex is provided without any guarantees as to its accuracy."),
                                 div(img(src="logo_c3bi_citech.jpg", width=300), style="text-align: center;"))
                        
                      )
                    )
                  )
))
