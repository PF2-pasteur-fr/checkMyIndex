# this file contains the functions used in both the "checkMyIndex" script and shiny application

checkInputIndexes <- function(index){
  if (any(duplicated(index$id))) stop("Input index ids are not unique, please check your input file.")
  if (!(all(unlist(strsplit(index$seq,"")) %in% c("A","C","T","G")))) stop("Input indexes contain other characters than A, T, C and G, please check your input file.")
  if (length(table(sapply(index$seq, nchar)))>1) stop("Input indexes do not have all the same length, please check your input file.")
}

isIndexesCombinationCompatible <- function(index){
  # return TRUE if the input indexes are compatible (i.e. can be used within the same lane)
  index$seqColor <- gsub("A|C", "R", index$seq)             # A and C are red
  index$seqColor <- gsub("G|T", "G", index$seqColor)        # G and T are green
  matColors <- matrix(unlist(sapply(index$seqColor, strsplit, "")), byrow=TRUE, nrow=nrow(index))  
  percentageRed <- apply(matColors, 2, function(x) mean(x=="R"))
  !any(percentageRed == 0 | percentageRed == 1) | nrow(index)==1
}

searchCompatibleIndexes <- function(index, nbSamplesPerLane){
  # generate the list of all the possible combinations
  possibleCombinations <- combn(x=index$id, m=nbSamplesPerLane, simplify=FALSE)
  indexesCombinations <- lapply(possibleCombinations, function(x) index[index$id %in% x,])
  # select only those compatible
  indexesCombinations[sapply(indexesCombinations, isIndexesCombinationCompatible)]
}

searchOneSolution <- function(compatibleIndexes, index, nbSamples, multiplexingRate, unicityConstraint="none"){
  # goal: look for a solution (i.e. a combination of combination of indexes) such that
  #  - each index is used only once if required (unicityConstraint parameter)
  #  - each combination of indexes is used only once if required (unicityConstraint parameter)
  # this function can return NULL if no solution is found (need to re-run in that case)
  nbLanes <- nbSamples/multiplexingRate
  compatibleCombinations <- vector(mode="list", length=nbLanes)
  k <- 1
  while (k <= nbLanes){
    if (length(compatibleIndexes)==0){
      return(NULL)
    } else{
      i <- sample(1:length(compatibleIndexes), 1, FALSE)
      compatibleCombinations[[k]] <- compatibleIndexes[[i]]
      # remove either the combination used or all the combinations for which an index has already been selected
      if (unicityConstraint=="index") compatibleIndexes <- compatibleIndexes[!sapply(compatibleIndexes, function(tab) any(tab$id %in% compatibleCombinations[[k]]$id))]
      if (unicityConstraint=="lane") compatibleIndexes <- compatibleIndexes[-i]
      k <- k+1
    }
  }
  data.frame(sample=1:nbSamples, lane=rep(1:nbLanes, each=multiplexingRate), do.call("rbind", compatibleCombinations))
}

findSolution <- function(compatibleIndexes, index, nbSamples, multiplexingRate, unicityConstraint="none", nbMaxTrials=1000){
  # this function run searchOneSolution() searchOneSolution times until finding a solution based on the parameters defined
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indexes: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  if (unicityConstraint!="none" && length(compatibleIndexes)<nbSamples/multiplexingRate){
    stop("Number of combinations of compatible indexes lower than the number of lanes desired. You can remove the index or lane unicity constraint.")
  }
  nbTrials <- 1
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(compatibleIndexes, index, nbSamples, multiplexingRate, unicityConstraint)
    if (!is.null(solution)) return(solution) else nbTrials <- nbTrials + 1
  }
  stop(paste("No solution found after", nbMaxTrials, "trials, you can increase this number in the parameters."))
}
