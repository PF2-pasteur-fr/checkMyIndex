library(shiny)
library(shinyjs)

shinyUI(fluidPage(theme = "bootstrap.min.css", shinyjs::useShinyjs(),

                  titlePanel(title=div(img(src="dna.png", width=50), strong("Search for a set of compatible indexes for your sequencing experiment")), windowTitle="checkMyIndex"),
                  
                  sidebarLayout(
                    
                    # parameters
                    sidebarPanel(div(id="allParameters",
                      
                      # input indexes
                      conditionalPanel(condition="!output.testdataProvided", {
                        fileInput("inputFile", label="Select your tab-delimited file containing the index 1 (i7) ids and sequences", accept="text")
                      }),
                      conditionalPanel(condition="!output.testdataProvided", {
                        fileInput("inputFile2", label="Optional index 2 (i5) file for dual-indexing", accept="text")
                      }),
                      conditionalPanel(condition="!output.inputFileProvided & !output.inputFile2Provided", {
                        selectizeInput("testdata", label="Load test indexes",
                                       choices=c("None"="none",
                                                 "24 indexes 1 (i7)"="simple",
                                                 "6 indexes 1 (i7) and 4 indexes 2 (i5)"="dual"))
                      }),
                      
                      # parameters
                      selectizeInput("chemistry", label="Illumina chemistry", 
                                     choices=c("Four-channels (HiSeq & MiSeq)" = 4, 
                                               "Two-channels (NovaSeq, NextSeq & MiniSeq)" = 2,
                                               "One-channel (iSeq 100)" = 1)),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("nbSamples")}),
                      conditionalPanel(condition="output.indexUploaded", {uiOutput("multiplexingRate")}),
                      conditionalPanel(condition="!output.indexUploaded2", {selectInput("unicityConstraint", label="Constraint on the indexes (single-indexing only)", 
                                                                                      choices=c("None" = "none", "Use each combination only once" = "lane", "Use each index only once" = "index"),
                                                                                      selected="none")}),
                      conditionalPanel(condition="output.indexUploaded & !output.indexUploaded2", {checkboxInput("completeLane", "Directly look for a solution with the desired multiplexing rate", value=FALSE)}),
                      conditionalPanel(condition="output.indexUploaded & !output.indexUploaded2", {checkboxInput("selectCompIndexes", "Select compatible indexes before looking for a solution", value=FALSE)}),
                      selectInput("nbMaxTrials", label="Maximum number of trials to find a solution", choices=10^(1:4)),
                      
                      # go and reset buttons
                      actionButton("go", label="Search for a solution"),
                      actionButton("reset", "Reset parameters")
                      
                    )),
                    
                    # output
                    mainPanel(
                      
                      tabsetPanel(id="mainPanel",
                        # 1st panel: input
                        tabPanel("Input indexes",
                                 value="inputIndexes",
                                 p(textOutput("textIndex")),
                                 dataTableOutput("inputIndex"),
                                 p(textOutput("textIndex2")),
                                 dataTableOutput("inputIndex2")),
                        
                        # 2nd panel: results
                        tabPanel("Proposed flowcell design",
                                 value="proposedSolution",
                                 shinyjs::hidden(div(
                                   id="proposedSolution",
                                   p(textOutput("textDescribingSolution")),
                                   dataTableOutput("solution"),
                                   br(""),
                                   uiOutput("downloadButton")
                                 ))
                                 ),
                        
                        # 3rd panel: plot results
                        tabPanel("Visualization of the design",
                                 shinyjs::hidden(div(
                                   id="visualization",
                                   p(textOutput("textDescribingHeatmap")),
                                   uiOutput("heatmapindex2")
                                 ))
                                 ),
                        
                        # 4th panel: help
                        tabPanel("Help",
                                 
                                 h3("Input indexes file(s)"),
                                 p("The user must provide its available indexes as two-column tab-delimited text file(s) without header: index ids are in the first column
                                   and corresponding sequences in the second. An example of such a file is available into the GitHub repository and also ", 
                                   a("here", href="inputIndexesExample.txt", target="blank", download="inputIndexesExample.txt")," to test the application.
                                   Note that for dual-indexing sequencing experiments the first file corresponds to indexes 1 (i7) and the second file to indexes 2 (i5).
                                   Example files for dual-indexing are available ",
                                   a("here", href="index24-i7.txt", target="blank", download="index24-i7.txt"), " (index 1) and ", 
                                   a("here", href="index24-i5.txt", target="blank", download="index24-i5.txt"), " (index 2)."),
                                 
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
                                 p(strong("Illumina chemistry"), "can be either four-channels (HiSeq & MiSeq), two-channels (NovaSeq, NextSeq & MiniSeq) or one-channel (iSeq 100).
                                   With the four-channel chemistry, a red laser detects A/C bases and a green laser detects G/T bases and the indexes are compatible if there is at 
                                   least one red light and one green light at each position. With the two-channel chemistry, G bases have no color, A bases are orange, C bases are
                                   red and T bases are green and indexes are compatible if there is at least one color at each position. Note that indexes starting with GG are not
                                   compatible with the two-channel chemistry. With the one-channel chemistry, compatibility cannot be defined with colors and indexes are compatible
                                   if there is at least one A or C or T base at each position. Please refer to the Illumina documentation for more detailed information on the different
                                   chemistries."),
                                 p(strong("Total number of samples"), "in your experiment (can be greater than the number of available indexes)."),
                                 p(strong("Multiplexing rate"), "i.e. number of samples per pool/lane (only divisors of the total number of samples are proposed)."),
                                 p(strong("Constraint on the indexes"), "(only for single-indexing) to avoid having two samples or two pools/lanes with the same index(es)."),
                                 p(strong("Directly look for a solution with the desired multiplexing rate"), "(only for single-indexing) instead of looking for a partial solution 
                                           with a few samples per pool/lane and then add some of the remaining indexes to reach the desired multiplexing rate."),
                                 p(strong("Select compatible indexes"), "(only for single-indexing) before looking for a (partial) solution can take some time but then speed up the algorithm."),
                                 p(strong("Maximum number of trials"), "can be increased if a solution is difficult to find with the parameters chosen."),

                                 h3("About"),
                                 p("This application has been developed at the Biomics pole of the Institut Pasteur by Hugo Varet and an ", a("Application Note", href="https://doi.org/10.1093/bioinformatics/bty706"), 
                                   "describing it has been published in 2018 in Bioinformatics. Feel free to send an e-mail to", a("hugo.varet@pasteur.fr"), "for any suggestion or bug report."),
                                 p("Source code and instructions to run it locally are available the", a("GitHub", href="https://github.com/PF2-pasteur-fr/checkMyIndex"), "repository. 
                                   Please note that checkMyIndex is provided without any guarantees as to its accuracy."),
                                 
                                 h3("Version"),
                                 p(paste0("This website executes checkMyIndex version ", checkMyIndexVersion, ".")),
                                 div(img(src="logo_c3bi_citech.jpg", width=300), style="text-align: center;"))
                        
                      )
                    )
                  )
))
