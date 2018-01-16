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
  matColors <- do.call("rbind", strsplit(index$color, ""))
  sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
  sumGreen <- nrow(matColors) - sumRed
  return(all(sumRed >= minRedGreen & sumGreen >= minRedGreen))
}

generateListOfIndexesCombinations <- function(index, nbSamplesPerLane, minRedGreen, completeLane, selectCompIndexes){
  # optimize nbSamplesPerLane to generate less combinations of indexes
  while (!completeLane & choose(nrow(index), nbSamplesPerLane)>1e5 & nbSamplesPerLane>2 & nbSamplesPerLane-1>=2*minRedGreen){
    nbSamplesPerLane <- nbSamplesPerLane - 1
  }
  if (choose(nrow(index), nbSamplesPerLane)>1e9){
    stop("Too many candidate combinations of indexes to easily find a solution, please use different parameters.")
  }
  index$color <- gsub("G|T", "G", gsub("A|C", "R", index$sequence)) # A and C are red and G and T are green
  # generate the list of all the possible combinations
  Clust=makeCluster(max(c(detectCores(logical=FALSE)-1, 1)))
  possibleCombinations <- combn(x=index$id, m=nbSamplesPerLane, simplify=FALSE)
  indexesCombinations <- parLapply(cl=Clust, X=possibleCombinations, fun=function(x, index) index[index$id %in% x,], index=index)
  # select only combinations of compatible indexes before searching for a solution
  if (selectCompIndexes) indexesCombinations <- indexesCombinations[parSapply(cl=Clust, X=indexesCombinations, FUN=areIndexesCompatible, minRedGreen=minRedGreen)]
  stopCluster(Clust)
  return(indexesCombinations)
}

searchOneSolution <- function(indexesList, index, nbLanes, multiplexingRate, unicityConstraint, minRedGreen){
  # goal: look for a solution (i.e. a combination of combination of indexes) such that
  #  - each index is used only once if required (unicityConstraint = index)
  #  - each combination of indexes is used only once if required (unicityConstraint = lane)
  # two steps:
  #  1) fill the lanes with as many samples per lane as in indexesList (may be lower than multiplexingRate for optimization purposes)
  #  2) if this number of lower than the desired multiplexing rate, complete the solution returned in 1) adding indexes
  # this function can return NULL if no solution is found (need to re-run in that case)
  inputNbSamplesPerLane <- nrow(indexesList[[1]])
  compatibleCombinations <- vector(mode="list", length=nbLanes)
  k <- 1
  while (k <= nbLanes){
    if (length(indexesList)==0){
      # no available combination of indexes anymore, try again!
      return(NULL)
    } else{
      i <- sample(1:length(indexesList), 1, FALSE)
      if (areIndexesCompatible(indexesList[[i]], minRedGreen)){
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
  solution <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane), lane=rep(1:nbLanes, each=inputNbSamplesPerLane), do.call("rbind", compatibleCombinations))
  if (multiplexingRate > inputNbSamplesPerLane){
    # only a partial solution has been found, need to complete it
    solution <- completeSolution(solution, index, multiplexingRate, unicityConstraint)
  }
  return(solution)
}

completeSolution <- function(partialSolution, index, multiplexingRate, unicityConstraint){
  nbSamplesToAdd <- multiplexingRate - nrow(partialSolution)/max(partialSolution$lane) # to each lane
  for (l in unique(partialSolution$lane)){
    if (unicityConstraint=="index"){
      # remove all the indexes already used
      index2 <- index[!(index$id %in% partialSolution$id),]
    } else{
      # remove the indexes already used in the current lane
      index2 <- index[!(index$id %in% partialSolution$id[partialSolution$lane==l]),]
    }
    if (nbSamplesToAdd > nrow(index2)) return(NULL) # not enough remaining indexes to complete the solution
    indexesToAdd <- index2[sample(1:nrow(index2), nbSamplesToAdd, FALSE),]
    indexesToAdd$lane <- l
    partialSolution <- rbind.data.frame(partialSolution[,c("lane","id","sequence")], indexesToAdd)
  }
  finalSolution <- data.frame(sample=1:nrow(partialSolution), partialSolution[order(partialSolution$lane, partialSolution$id),])
  finalSolution$color <- gsub("G|T", "G", gsub("A|C", "R", finalSolution$sequence)) # A and C are red and G and T are green
  return(finalSolution)
}

findSolution <- function(indexesList, index, nbSamples, multiplexingRate, unicityConstraint, minRedGreen, nbMaxTrials, completeLane, selectCompIndexes){
  # this function run searchOneSolution() nbMaxTrials times until finding a solution based on the parameters defined
  if (!unicityConstraint %in% c("none","index","lane")) stop("unicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indexes: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  nbLanes <- nbSamples/multiplexingRate
  if (unicityConstraint!="none" & completeLane & selectCompIndexes & length(indexesList) < nbLanes){
    stop("There are only ", length(indexesList), " combinations of compatible indexes to fill ", nbLanes, " lanes.")
  }
  nbTrials <- 1
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(indexesList, index, nbLanes, multiplexingRate, unicityConstraint, minRedGreen)
    if (!is.null(solution)){
      checkProposedSolution(solution, unicityConstraint, minRedGreen)
      return(solution)
    } else{
      nbTrials <- nbTrials + 1
    }
  }
  stop(paste("No solution found after", nbMaxTrials, "trials using these parameters, you can modify them or increase the number of trials."))
}

checkProposedSolution <- function(solution, unicityConstraint, minRedGreen){
  # indexes not unique
  if (unicityConstraint=="index" & any(duplicated(solution$id))) stop("The solution proposed uses some indexes several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  # indexes not unique within a lane
  if (any(sapply(split(solution$id, solution$lane), function(x) any(duplicated(x))))) stop("The solution proposed uses some indexes several times within a lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  # several lanes with the same combination of indexes
  if (unicityConstraint=="lane" & any(duplicated(split(solution$id, solution$lane)))) stop("The solution proposed uses some combinations of indexes several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  # different number of samples on the lanes
  if (length(table(table(solution$lane)))>1) stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  # indexes not compatible on at least one lane
  if (any(!sapply(split(solution, solution$lane), areIndexesCompatible, minRedGreen))) stop("The solution proposed uses incompatible indexes. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
}

calculateFinalMinRedGreen <- function(solution){
  tmpfun <- function(solution.lane){
    matColors <- matrix(unlist(sapply(solution.lane$color, strsplit, "")), byrow=TRUE, nrow=nrow(solution.lane))
    sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
    sumGreen <- nrow(matColors) - sumRed
    return(min(c(sumRed, sumGreen)))
  }
  min(sapply(split(solution, solution$lane), tmpfun))
}

heatmapindex <- function(solution){
  splitsol <- split(solution, solution$lane)
  # build a matrix containing the index bases
  seqmat <- lapply(splitsol, function(sol) matrix(unlist(strsplit(sol$sequence, "")), nrow=nrow(sol), byrow=TRUE))
  seqmat <- do.call("rbind", lapply(seqmat, function(x) rbind(x, rep(NA, ncol(x)))))
  seqmat <- seqmat[-nrow(seqmat),]
  # build a matrix containing the colors
  seqcol <- lapply(splitsol, function(sol) matrix(unlist(strsplit(sol$color, "")), nrow=nrow(sol), byrow=TRUE))
  seqcol <- do.call("rbind", lapply(seqcol, function(x) rbind(x, rep(NA, ncol(x)))))
  seqcol <- seqcol[-nrow(seqcol),]
  # plot
  par(mar=c(2, 6, 1, 6) + 0.1, xpd=TRUE, xaxs="i", yaxs="i")
  plot.new()
  plot.window(xlim=c(-1.5, ncol(seqmat)+1.5), ylim=c(-2, nrow(seqmat)))
  for (i in 1:nrow(seqmat)){
    for (j in 1:ncol(seqmat)){
      rect(xleft=j-1, ybottom=nrow(seqmat)-(i-1), xright=j, ytop=nrow(seqmat)-i, 
           col=ifelse(seqcol[i,j]=="R", "orangered", "darkseagreen4"),
           border=ifelse(seqcol[i,j]=="R", "orangered", "darkseagreen4"))
      text(x=j-0.5, y=nrow(seqmat)-(i-0.5), labels=seqmat[i,j])
    }
  }
  # print positions
  text(x=1:ncol(seqmat)-0.5, y=-0.75, labels=1:ncol(seqmat))
  text(x=mean(1:ncol(seqmat))-0.5, y=-2, labels="Position")
  # extract and print sample ids
  tmp <- unlist(lapply(splitsol, function(sol) c(sol$sample, NA)))
  text(x=-0.25, y=nrow(seqmat):1 - 0.5, labels=tmp[-length(tmp)], pos=2)
  text(x=-1.5, y=nrow(seqmat)/2, labels="Sample", srt=90)
  # extract and print index ids
  tmp <- unlist(lapply(splitsol, function(sol) c(sol$id, NA)))
  text(x=ncol(seqmat)+0.25, y=nrow(seqmat):1 - 0.5, labels=tmp[-length(tmp)], pos=4)
  text(x=ncol(seqmat)+1.5, y=nrow(seqmat)/2, labels="Index", srt=90)
}
