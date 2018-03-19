library(shiny)

options(shiny.sanitize.errors = FALSE,   # to display informative error messages
        shiny.maxRequestSize = 5*1024^2) # limit size for the input file to upload (5Mo here)

shinyServer(function(input, output) {
  
  # list of input indexes 1
  inputIndex <- reactive({
    if (is.null(input$inputFile)){
      return(NULL)
    } else{
      index <- readIndexesFile(file=input$inputFile$datapath)
      index <- addColors(index, input$chemistry)
      return(index)
    }
  })
  output$indexUploaded <- reactive({return(!is.null(inputIndex()))})
  outputOptions(output, "indexUploaded", suspendWhenHidden=FALSE)
  output$inputIndex <- renderDataTable({inputIndex()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  textIndex <- reactive({
    if (is.null(inputIndex())) "" else paste("The table below shows the", nrow(inputIndex()), "i7 indexes uploaded:")
  })
  output$textIndex <- renderText({textIndex()})
  
  # list of input indexes 2
  inputIndex2 <- reactive({
    if (is.null(input$inputFile2)){
      return(NULL)
    } else{
      if (is.null(inputIndex())) stop("Please load i7 indexes first.")
      index2 <- readIndexesFile(file=input$inputFile2$datapath)
      index2 <- addColors(index2, input$chemistry)
      return(index2)
    }
  })
  output$indexUploaded2 <- reactive({return(!is.null(inputIndex2()))})
  outputOptions(output, "indexUploaded2", suspendWhenHidden=FALSE)
  output$inputIndex2 <- renderDataTable({inputIndex2()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  textIndex2 <- reactive({
    if (is.null(inputIndex2())) "" else paste("The table below shows the", nrow(inputIndex2()), "i5 indexes uploaded:")
  })
  output$textIndex2 <- renderText({textIndex2()})
  
  # propose both the possible nb samples and multiplexing rates according to the input list of indexes
  output$nbSamples <- renderUI({
    if (is.null(inputIndex())){
      return("")
    } else{
      index <- inputIndex()
      nr <- ifelse(input$chemistry == "4", nrow(index), nrow(index[substr(index$sequence, 1, 2) != "GG",]))
      if (!is.null(inputIndex2())){
        index2 <- inputIndex2()
        nr2 <- ifelse(input$chemistry == "4", nrow(index2), nrow(index2[substr(index2$sequence, 1, 2) != "GG",]))
      } else{
        nr2 <- 1
      }
      numericInput("nbSamples", label="Total number of samples in the experiment", value=nr*nr2, min=2, step=1)
    }
  })
  output$multiplexingRate <- renderUI({
    if (is.null(input$nbSamples)){
      return("")
    } else{
      nbSamples <- as.numeric(input$nbSamples)
      if (is.na(nbSamples) || nbSamples %% 1 != 0) stop("Number of samples must be an integer")
      if (nbSamples <= 1) stop("Number of samples must be greater than 1.")
      mr <- 1:nbSamples
      choices <- mr[sapply(mr, function(x) nbSamples %% x == 0)]
      index <- inputIndex()
      nr <- ifelse(input$chemistry == "4", nrow(index), nrow(index[substr(index$sequence, 1, 2) != "GG",]))
      if (!is.null(inputIndex2())){
        index2 <- inputIndex2()
        nr2 <- ifelse(input$chemistry == "4", nrow(index2), nrow(index2[substr(index2$sequence, 1, 2) != "GG",]))
      } else{
        nr2 <- 1
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
    if (is.null(inputIndex2())){
      return(NULL)
    } else{
      return(generateListOfIndexesCombinations(index = inputIndex2(),
                                               nbSamplesPerLane = as.numeric(input$multiplexingRate),
                                               completeLane = input$completeLane,
                                               selectCompIndexes = input$selectCompIndexes,
                                               chemistry = input$chemistry))
    }
  })
  
  # text describing the solution
  textDescribingSolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      return("")
    } else{
      return(paste("Below is a solution for", as.numeric(input$nbSamples)/as.numeric(input$multiplexingRate),
                   "pool(s) of", input$multiplexingRate, "samples using the parameters specified:"))
    }
  })
  output$textDescribingSolution <- renderText({textDescribingSolution()})
  
  # print the solution
  displaySolution <- eventReactive(input$go, {
    if (is.null(input$multiplexingRate) | is.null(inputIndex())){
      return("")
    } else{
      return(findSolution(indexesList = generateList(),
                          index = inputIndex(),
                          indexesList2 = generateList2(),
                          index2 = inputIndex2(),
                          nbSamples = as.numeric(input$nbSamples),
                          multiplexingRate = as.numeric(input$multiplexingRate),
                          unicityConstraint = input$unicityConstraint,
                          nbMaxTrials = as.numeric(input$nbMaxTrials),
                          completeLane = input$completeLane,
                          selectCompIndexes = input$selectCompIndexes,
                          chemistry = input$chemistry))
    }
  })
  output$solution <- renderDataTable({displaySolution()}, options=list(paging=FALSE, searching=FALSE, info=FALSE))
  
  # get the number samples of the solution to avoid resizing the heatmap when modifying the input number of samples
  getNbSamples <- eventReactive(input$go, {nrow(displaySolution())})
  # plot of the solution
  output$heatmapindex <- renderPlot({heatmapindex(displaySolution())}, res=90)
  output$heatmapindex2 <- renderUI({plotOutput("heatmapindex", width=900, height=220+20*getNbSamples())})
  
  # download the solution
  output$downloadData <- downloadHandler(
    filename = "chosenIndexes.txt",
    content = function(file) write.table(displaySolution(), file, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  )
  
})
