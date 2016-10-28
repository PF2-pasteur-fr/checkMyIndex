# this file contains the functions used in both the "checkMyIndex" script and shiny application
library(parallel)

readIndexesFile <- function(file) read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","sequence"))

checkInputIndexes <- function(index){
  if (any(duplicated(index$id))) stop("Input index ids are not unique, please check your input file.")
  if (any(duplicated(index$sequence))) stop("Input indexes are not unique, please check your input file.")
  if (!(all(unlist(strsplit(index$sequence,"")) %in% c("A","C","T","G")))) stop("Input indexes contain other characters than A, T, C and G, please check your input file.")
  if (length(table(sapply(index$sequence, nchar)))>1) stop("Input indexes do not have all the same length, please check your input file.")
}

areIndexesCompatible <- function(index, minRedGreen){
  # return TRUE if the input indexes are compatible (i.e. can be used within the same lane)
  if (nrow(index)==1) return(TRUE)
  matColors <- matrix(unlist(sapply(index$color, strsplit, "")), byrow=TRUE, nrow=nrow(index))  
  sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
  sumGreen <- nrow(matColors) - sumRed
  return(all(sumRed >= minRedGreen & sumGreen >= minRedGreen))
}

generateListOfIndexesCombinations <- function(index, multiplexingRate, minRedGreen, selectCompIndexes){
  index$color <- gsub("G|T", "G", gsub("A|C", "R", index$sequence)) # A and C are red and G and T are green
  # generate the list of all the possible combinations
  Clust=makeCluster(max(c(detectCores(logical=FALSE)-1, 1)))
  possibleCombinations <- combn(x=index$id, m=multiplexingRate, simplify=FALSE)
  indexesCombinations <- parLapply(cl=Clust, X=possibleCombinations, fun=function(x, index) index[index$id %in% x,], index=index)
  # select only combinations of compatible indexes before searching for a solution
  if (selectCompIndexes) indexesCombinations <- indexesCombinations[parSapply(cl=Clust, X=indexesCombinations, FUN=areIndexesCompatible, minRedGreen=minRedGreen)]
  stopCluster(Clust)
  return(indexesCombinations)
}

searchOneSolution <- function(indexesList, nbSamples, multiplexingRate, unicityConstraint, minRedGreen, selectCompIndexes){
  # goal: look for a solution (i.e. a combination of combination of indexes) such that
  #  - each index is used only once if required (unicityConstraint = index)
  #  - each combination of indexes is used only once if required (unicityConstraint = lane)
  # this function can return NULL if no solution is found (need to re-run in that case)
  # indexesList is either the list of all the compatible indexes (selectCompIndexes) or all the possible combinations
  nbLanes <- nbSamples/multiplexingRate
  compatibleCombinations <- vector(mode="list", length=nbLanes)
  k <- 1
  while (k <= nbLanes){
    if (length(indexesList)==0){
      # no available combination of indexes anymore, try again!
      return(NULL)
    } else{
      i <- sample(1:length(indexesList), 1, FALSE)
      # if selectCompIndexes, the indexes a compatible for sure, no need to test the compatibility
      if (selectCompIndexes || areIndexesCompatible(indexesList[[i]], minRedGreen)){
        compatibleCombinations[[k]] <- indexesList[[i]]
        # remove either the combination used or all the combinations for which an index has already been selected
        if (unicityConstraint=="index") indexesList <- indexesList[!sapply(indexesList, function(tab) any(tab$id %in% indexesList[[i]]$id))]
        if (unicityConstraint=="lane") indexesList <- indexesList[-i]
        k <- k+1
      } else{
        indexesList <- indexesList[-i]
      }
    }
  }
  return(data.frame(sample=1:nbSamples, lane=rep(1:nbLanes, each=multiplexingRate), do.call("rbind", compatibleCombinations)))
}

findSolution <- function(indexesList, index, nbSamples, multiplexingRate, unicityConstraint, minRedGreen, nbMaxTrials, selectCompIndexes){
  # this function run searchOneSolution() nbMaxTrials times until finding a solution based on the parameters defined
  if (!unicityConstraint %in% c("none","index","lane")) stop("unicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indexes: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  nbLanes <- nbSamples/multiplexingRate
  if (selectCompIndexes & unicityConstraint!="none" & length(indexesList)<nbLanes){
    stop("There are only ", length(indexesList), " combinations of compatible indexes to fill ", nbLanes, " lanes. You can remove the index or lane unicity constraint.")
  }

  nbTrials <- 1
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(indexesList, nbSamples, multiplexingRate, unicityConstraint, minRedGreen, selectCompIndexes)
    if (!is.null(solution)){
      checkProposedSolution(solution, unicityConstraint, minRedGreen)
      return(solution)
    } else{
      nbTrials <- nbTrials + 1
    }
  }
  stop(paste("No solution found after", nbMaxTrials, "trials, you can increase this number in the parameters."))
}

checkProposedSolution <- function(solution, unicityConstraint, minRedGreen){
  # indexes not unique
  err1 <- unicityConstraint=="index" & any(duplicated(solution$id))
  # several lanes with the same combination of indexes
  err2 <- unicityConstraint=="lane" & any(duplicated(split(solution$id, solution$lane)))
  # different number of samples on the lanes
  err3 <- length(table(table(solution$lane)))>1
  # indexes not compatible on at least one lane
  err4 <- any(!sapply(split(solution, solution$lane), areIndexesCompatible, minRedGreen))
  err <- c(err1, err2, err3, err4)
  if (any(err)){
    stop("The solution proposed by the algorithm does not satisfy the constraints (",
         paste(paste0("err", 1:length(err))[err], collapse=", "), 
         "). Thanks to report this error to Hugo Varet <hugo.varet@pasteur.fr>")
  }
}
