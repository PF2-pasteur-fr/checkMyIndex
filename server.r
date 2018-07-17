library(shiny)

options(shiny.sanitize.errors = FALSE,   # to display informative error messages
        shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output, session) {
  
  # list of input indexes 1
  inputIndex <- reactive({
    if (input$testdata %in% c("simple", "dual")){
      file <- ifelse(input$testdata=="simple", "www/inputIndexesExample.txt", "www/index24-i7.txt")
    } else{
      if (!is.null(input$inputFile)) file <- input$inputFile$datapath else return(NULL)
    }
    index <- tryCatch({readIndexesFile(file)}, 
                      error = function(e) stop("An error occured when loading index 1 file, please check its structure."))
    index <- addColors(index, input$chemistry)
    index$score <- scores(index$sequence)
    return(index)
  })
  output$indexUploaded <- reactive({!is.null(inputIndex())})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  output$inputIndex <- renderDataTable({inputIndex()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  textIndex <- reactive({
    index <- inputIndex()
    if (is.null(index)){
      "No index file loaded yet, use the left panel to select an input file."
    } else{
      paste0("The table below shows the ", nrow(index), " indexes 1 (i7) uploaded with the colors corresponding
             to the chosen Illumina chemistry and the minimum number of mismatches with the other indexes.
             Note that the smallest number of mismatches between two indexes of this list is ", min(index$score), ".")
    }
  })
  output$textIndex <- renderText({tryCatch({textIndex()}, error = function(e) NULL)})
  
  # list of input indexes 2
  inputIndex2 <- reactive({
    if (input$testdata == "simple") return(NULL)
    if (input$testdata == "dual"){
      file2 <- "www/index24-i5.txt"
    } else{
      if (!is.null(input$inputFile2)){
        if (is.null(inputIndex())) stop("Please load indexes 1 (i7) first.")
        file2 <- input$inputFile2$datapath
      } else{
        return(NULL)
      }
    }
    index2 <- tryCatch({readIndexesFile(file2)}, 
                       error = function(e) stop("An error occured when loading index 2 file, please check its structure."))
    index2 <- addColors(index2, input$chemistry)
    index2$score <- scores(index2$sequence)
    return(index2)
  })
  output$indexUploaded2 <- reactive({!is.null(inputIndex2())})
  outputOptions(output, "indexUploaded2", suspendWhenHidden=FALSE)
  output$inputIndex2 <- renderDataTable({inputIndex2()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  textIndex2 <- reactive({
    index2 <- inputIndex2()
    if (is.null(index2)){
      ""
    } else{
      paste0("The table below shows the ", nrow(index2), " indexes 2 (i5) uploaded with the colors corresponding
             to the chosen Illumina chemistry and the minimum number of mismatches with the other indexes.
             Note that the smallest number of mismatches between two indexes of this list is ", min(index2$score), ".")
      
    }
  })
  output$textIndex2 <- renderText({tryCatch({textIndex2()}, error = function(e) NULL)})
  
  # propose both the possible nb samples and multiplexing rates according to the input list of indexes
  output$nbSamples <- renderUI({
    index <- tryCatch({inputIndex()}, error = function(e) NULL)
    if (is.null(index)){
      return("")
    } else{
      nr <- ifelse(input$chemistry == "2", nrow(index[substr(index$sequence, 1, 2) != "GG",]), nrow(index))
      index2 <- tryCatch({inputIndex2()}, error = function(e) NULL)
      if (is.null(index2)){
        nr2 <- 1
      } else{
        nr2 <- ifelse(input$chemistry == "2", nrow(index2[substr(index2$sequence, 1, 2) != "GG",]), nrow(index2))
      }
      numericInput("nbSamples", label="Total number of samples in the experiment", value=nr*nr2, min=2, step=1)
    }
  })
  output$multiplexingRate <- renderUI({
    index <- tryCatch({inputIndex()}, error = function(e) NULL)
    if (is.null(index) | is.null(input$nbSamples)){
      return("")
    } else{
      nbSamples <- as.numeric(input$nbSamples)
      if (is.na(nbSamples) || nbSamples %% 1 != 0) stop("Number of samples must be an integer.")
      if (nbSamples <= 1) stop("Number of samples must be greater than 1.")
      mr <- 1:nbSamples
      choices <- mr[sapply(mr, function(x) nbSamples %% x == 0)]
      nr <- ifelse(input$chemistry == "2", nrow(index[substr(index$sequence, 1, 2) != "GG",]), nrow(index))
      index2 <- tryCatch({inputIndex2()}, error = function(e) NULL)
      if (is.null(index2)){
        nr2 <- 1
      } else{
        nr2 <- ifelse(input$chemistry == "2", nrow(index2[substr(index2$sequence, 1, 2) != "GG",]), nrow(index2))
      }
      choices <- choices[choices <= nr*nr2]
      selectInput("multiplexingRate", label="Multiplexing rate (i.e. # of samples per pool)", 
                  choices=choices, selected=ifelse(input$chemistry == "2", choices[length(choices)], choices[2]))
    }
  })
  
  # generate list(s) of indexes
  generateList <- reactive({
    return(generateListOfIndexesCombinations(index = inputIndex(),
                                             nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                             completeLane = input$completeLane,
                                             selectCompIndexes = input$selectCompIndexes,
                                             chemistry = input$chemistry))
  })
  generateList2 <- reactive({
    index2 <- inputIndex2()
    if (is.null(index2)){
      return(NULL)
    } else{
      return(generateListOfIndexesCombinations(index = index2,
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry))
    }
  })
  
  # automatically go to the proposed solution when pressing "search for a solution"
  observeEvent(input$go, {updateTabsetPanel(session, "mainPanel", selected = "proposedSolution")})
  
  # text describing the solution
  textDescribingSolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      "Please load indexes before pressing the \"Search for a solution\" button."
    } else{
      paste("Below is a solution for", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate),
            "pool(s) of", input$multiplexingRate, "samples using the parameters specified. The table contains
            one row per sample to be sequenced and several columns: pool/lane labels, index ids, index sequences,
            the corresponding colors according to the chosen Illumina chemistry and a score equal to the minimum
            number of mismatches with the other indexes of the pool/lane.")
    }
  })
  output$textDescribingSolution <- renderText({tryCatch({textDescribingSolution()}, error = function(e) NULL)})
  
  # display the solution
  displaySolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      return(NULL)
    } else{
      withProgress({
        return(findSolution(indexesList = generateList(),
                            index = inputIndex(),
                            indexesList2 = generateList2(),
                            index2 = inputIndex2(),
                            nbSamples = as.numeric(input$nbSamples),
                            multiplexingRate = as.numeric(input$multiplexingRate),
                            unicityConstraint = ifelse(is.null(inputIndex2()), input$unicityConstraint, "none"),
                            nbMaxTrials = as.numeric(input$nbMaxTrials),
                            completeLane = input$completeLane,
                            selectCompIndexes = input$selectCompIndexes,
                            chemistry = input$chemistry))
      }, message="R is looking for a solution...", max=0)
      
    }
  })
  output$solution <- renderDataTable({displaySolution()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  
  # get the number samples of the solution to avoid resizing the heatmap when modifying the input number of samples
  getNbSamples <- eventReactive(input$go, {nrow(displaySolution())})
  # text describing the solution
  textDescribingHeatmap <- eventReactive(input$go, {
    if (is.null(displaySolution())){
      ""
    } else{
      paste0("The plot below allows to vizualize the proposed solution. Samples (in rows) are grouped by pool/lane
             and each nucleotide of each index is displayed with a color according to the chosen Illumina chemistry. 
             One can thus quickly check whether each color is used at each position. Note that sample ids (from 1 to ",
             getNbSamples(), ") are printed on the left while index ids are printed on the right.")
    }
  })
  output$textDescribingHeatmap <- renderText({textDescribingHeatmap()})
  # plot of the solution
  output$heatmapindex <- renderPlot({heatmapindex(displaySolution())}, res=90)
  output$heatmapindex2 <- renderUI({
    if (!is.null(tryCatch({displaySolution()}, error = function(e) NULL))){
      plotOutput("heatmapindex", width=900, height=220+20*getNbSamples())
    }
  })
  
  # download the solution
  output$downloadData <- downloadHandler(
    filename = "chosenIndexes.txt",
    content = function(file) write.table(displaySolution(), file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  output$downloadButton <- renderUI({
    if (!is.null(tryCatch({displaySolution()}, error = function(e) NULL))) downloadButton("downloadData", "Download")
  })
  
})
