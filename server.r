library(shiny)

options(shiny.sanitize.errors = FALSE,   # to display informative error messages
        shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output) {
  
  # list of input indexes
  inputIndex <- reactive({
    if(is.null(input$inputFile)){
      return(NULL)
    } else{
      index <- readIndexesFile(file=input$inputFile$datapath)
      index <- addColors(index, input$chemistry)
      return(index)
    }
  })
  output$indexUploaded <- reactive({return(!is.null(inputIndex()))})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  output$inputIndex <- renderDataTable({
    index <- inputIndex()
    startGG <- sapply(index$sequence, substr, 1, 2) == "GG"
    if (any(startGG)){
      index$comment <- ""
      index[startGG, "comment"] <- "Can't be used with the two-channel chemistry"
    }
    return(index)
  }, options=list(paging=FALSE, searching=FALSE))
  
  # propose both the possible nb samples and multiplexing rates according to the input list of indexes
  output$nbSamples <- renderUI({
    if (!is.null(inputIndex())){
      index <- inputIndex()
      nr <- ifelse(input$chemistry == "4", nrow(index), nrow(index[substr(index$sequence, 1, 2) != "GG",]))
      numericInput("nbSamples", label="Total number of samples in the experiment", value=nr, min=2, step=1)
    } else{
      ""
    }
  })
  output$multiplexingRate <- renderUI({
    if (!is.null(input$nbSamples)){
      nbSamples <- as.numeric(input$nbSamples)
      if (is.na(nbSamples) || nbSamples %% 1 != 0) stop("Number of samples must be an integer")
      if (nbSamples <= 1) stop("Number of samples must be greater than 1.")
      mr <- 1:nbSamples
      choices <- mr[sapply(mr, function(x) nbSamples %% x == 0)]
      choices <- choices[choices <= nrow(inputIndex())]
      selectInput("multiplexingRate", label="Multiplexing rate (i.e. # of samples per pool)", 
                  choices=choices, selected=ifelse(input$chemistry == "2", choices[length(choices)], choices[2]))
    } else{
      ""
    }
  })
  
  # number of candidate combinations of indexes
  textNbCombinations <- reactive({
    if (is.null(inputIndex()) | is.null(input$multiplexingRate)){
      ""
    } else{
      paste("In the input list of", nrow(inputIndex()), "indexes: there are", choose(n=nrow(inputIndex()), k=as.numeric(input$multiplexingRate)), 
            "possible combinations of", as.numeric(input$multiplexingRate), "indexes (not necessarily compatible).")
    }
  })
  output$textNbCombinations <- renderText({textNbCombinations()})
  
  # generate list of indexes
  generateList <- reactive({
    generateListOfIndexesCombinations(index = inputIndex(),
                                      nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                      completeLane = input$completeLane,
                                      selectCompIndexes = input$selectCompIndexes,
                                      chemistry = input$chemistry)
  })
  
  # number of compatible combinations of indexes
  textNbCompCombinations <- eventReactive(input$go, {
    if (is.null(inputIndex()) | is.null(input$multiplexingRate)){
      ""
    } else{
      if (nrow(generateList()[[1]]) == as.numeric(input$multiplexingRate) & input$selectCompIndexes){
        paste("Among them", length(generateList()), "contain compatible indexes.")
      } else{
        paste("Note: the number of combinations containing compatible indexes cannot be calculated with the parameters used.")
      }
    }
  })
  output$textNbCompCombinations <- renderText({textNbCompCombinations()})
  
  # text describing the solution
  textDescribingSolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      ""
    } else{
      paste("Below is a solution for", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate), 
            "pool(s) of", input$multiplexingRate, "samples using the parameters specified:")
    }
  })
  output$textDescribingSolution <- renderText({textDescribingSolution()})
  
  displaySolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      ""
    } else{
      return(findSolution(indexesList = generateList(),
                          index = inputIndex(),
                          nbSamples = as.numeric(input$nbSamples),
                          multiplexingRate = as.numeric(input$multiplexingRate),
                          unicityConstraint = input$unicityConstraint,
                          nbMaxTrials = as.numeric(input$nbMaxTrials),
                          completeLane = input$completeLane,
                          selectCompIndexes = input$selectCompIndexes,
                          chemistry = input$chemistry))
    }
  })
  output$solution <- renderDataTable({displaySolution()}, options=list(paging=FALSE, searching=FALSE))
  
  getNbSamples <- eventReactive(input$go, {nrow(displaySolution())})
  output$heatmapindex <- renderPlot({heatmapindex(displaySolution())}, res=90)
  output$heatmapindex2 <- renderUI({
    plotOutput("heatmapindex", width=900, height=220+20*getNbSamples())
  })
  
  output$downloadData <- downloadHandler(
    filename = "chosenIndexes.txt",
    content = function(file) write.table(displaySolution(), file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  
})
