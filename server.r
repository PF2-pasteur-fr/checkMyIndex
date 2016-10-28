library(shiny)

options(shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output) {
  
  inputIndex <- reactive({
    if(is.null(input$inputFile)){
      return(NULL)
    } else{
      index <- readIndexesFile(input$inputFile$datapath)
      checkInputIndexes(index)
      return(index)
    }
  })
  output$indexUploaded <- reactive({return(!is.null(inputIndex()))})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  
  # propose both the possible nb samples and multiplexing rates according to the input list of indexes
  output$nbSamples <- renderUI({if (!is.null(inputIndex())){
    nr <- nrow(inputIndex())
    numericInput("nbSamples", label="Total number of samples in the experiment", value=nr, min=2, step=1)
  } else{
    ""
  }})
  output$multiplexingRate <- renderUI({if (!is.null(input$nbSamples)){
    nbSamples <- as.numeric(input$nbSamples)
    if (is.na(nbSamples) || nbSamples %% 1 != 0) stop("Number of samples must be an integer")
    if (nbSamples <= 1) stop("Number of samples must be greater than 1.")
    mr <- 1:nbSamples
    choices <- mr[sapply(mr, function(x) nbSamples %% x == 0)]
    selectInput("multiplexingRate", label="Multiplexing rate (i.e. number of samples per lane)", choices=choices, selected=choices[2])
  } else{
    ""
  }})
  # propose some minimal number of red/green per position according to the multiplexing rate
  output$minRedGreen <- renderUI({if (!is.null(input$multiplexingRate)){
    maxRedGreen <- max(trunc(as.numeric(input$multiplexingRate)/2), 1)
    selectInput("minRedGreen", label="Minimal number of red/green lights per position", choices=1:maxRedGreen, selected=1)
  } else{
    ""
  }})
  
  # list of input indexes
  output$inputIndex <- renderDataTable({inputIndex()},options=list(paging=FALSE,searching=FALSE))
  
  textNbCombinations <- reactive({
    if (is.null(input$multiplexingRate) | is.null(input$minRedGreen) | is.null(inputIndex())){
      ""
    } else{
      paste("In the input list of", nrow(inputIndex()), "indexes: there are", choose(n=nrow(inputIndex()), k=as.numeric(input$multiplexingRate)), 
            "possible combinations of", as.numeric(input$multiplexingRate), "indexes (not necessarily compatible).")
    }
  })
  output$textNbCombinations <- renderText({textNbCombinations()})
  
  # generate list of indexes
  generateList <- reactive({
    generateListOfIndexesCombinations(inputIndex(), as.numeric(input$multiplexingRate), as.numeric(input$minRedGreen), input$selectCompIndexes)
  })
  
  textNbCompatibleIndexes <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(input$minRedGreen) | is.null(inputIndex()) | !input$selectCompIndexes){
      ""
    } else{
      paste("Among them", length(generateList()), "contain compatible indexes (i.e. at least", input$minRedGreen, "red/green light(s) per position).")
    }
  })
  output$textNbCompatibleIndexes <- renderText({textNbCompatibleIndexes()})

  textDescribingSolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      ""
    } else{
      paste0("Below is a solution for ", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate), 
             " lanes of ", input$multiplexingRate, " samples using the parameters specified:")
    }
  })
  output$textDescribingSolution <- renderText({textDescribingSolution()})
  
  displaySolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(input$minRedGreen) | is.null(inputIndex())){
      ""
    } else{
      unicityConstraint <- ifelse(input$unicityConstraint=="None", "none",
                                  ifelse(input$unicityConstraint=="Use each combination only once", 
                                         "lane", "index"))
      return(findSolution(generateList(), inputIndex(), as.numeric(input$nbSamples), as.numeric(input$multiplexingRate), 
                          unicityConstraint, as.numeric(input$minRedGreen), as.numeric(input$nbMaxTrials), input$selectCompIndexes))
    }
  })
  output$solution <- renderDataTable({displaySolution()},options=list(paging=FALSE,searching=FALSE))
  
  output$downloadData <- downloadHandler(
    filename = "chosenIndexes.txt",
    content = function(file) write.table(displaySolution(), file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  
})
