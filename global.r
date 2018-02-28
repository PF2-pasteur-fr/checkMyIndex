# this file contains the functions used in both the "checkMyIndex" script and shiny application
library(parallel)

readIndexesFile <- function(file){
  index <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","sequence"))
  checkInputIndexes(index)
  return(index)
}

addColors <- function(index, chemistry){
  # convert bases into colors
  if (chemistry == "2"){
    # G has no color, A is orange, C is red and T is green
    index$color <- gsub("T", "G", gsub("C", "R", gsub("A", "O", gsub("G", "-", index$sequence))))
  } else{
    # A/C are red and G/T are green
    index$color <- gsub("G|T", "G", gsub("A|C", "R", index$sequence))
  }
  return(index)
}

checkInputIndexes <- function(index){
  if (any(duplicated(index$id))) stop("Input index ids are not unique, please check your input file.")
  if (any(duplicated(index$sequence))) stop("Input indexes are not unique, please check your input file.")
  if (!(all(unlist(strsplit(index$sequence,"")) %in% c("A","C","T","G")))) stop("Input indexes contain other characters than A, T, C and G, please check your input file.")
  if (length(table(sapply(index$sequence, nchar))) > 1) stop("Input indexes do not have all the same length, please check your input file.")
}

areIndexesCompatible <- function(index, chemistry){
  # return TRUE if the input indexes are compatible (i.e. can be used within the same pool/lane)
  if (nrow(index)==1) return(TRUE)
  matColors <- do.call("rbind", strsplit(index$color, ""))
  if (chemistry == "2"){
    sumNoColor <- apply(matColors, 2, function(x) sum(x=="-"))
    return(all(sumNoColor < nrow(index)))
  } else{
    sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
    sumGreen <- nrow(matColors) - sumRed
    return(all(sumRed >= 1 & sumGreen >= 1))
  }
}

generateListOfIndexesCombinations <- function(index, nbSamplesPerLane, completeLane, selectCompIndexes, chemistry){
  # remove indexes starting with GG if two-channel chemistry
  if (chemistry == "2") index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]
  # optimize nbSamplesPerLane to generate a reduced number of combinations of indexes
  while (!completeLane & choose(nrow(index), nbSamplesPerLane) > 1e5 & nbSamplesPerLane > 2) nbSamplesPerLane <- nbSamplesPerLane - 1
  # stop if too many combinations to test
  if (choose(nrow(index), nbSamplesPerLane) > 1e9) stop("Too many candidate combinations of indexes to easily find a solution, please use different parameters.")
  # generate the list of all the possible combinations
  Clust <- makeCluster(max(c(detectCores(logical=FALSE)-1, 1)))
  possibleCombinations <- combn(x=index$id, m=nbSamplesPerLane, simplify=FALSE)
  indexesCombinations <- parLapply(cl=Clust, X=possibleCombinations, fun=function(x, index) index[index$id %in% x,], index=index)
  # select only combinations of compatible indexes before searching for a solution
  if (selectCompIndexes) indexesCombinations <- indexesCombinations[parSapply(cl=Clust, X=indexesCombinations, FUN=areIndexesCompatible, chemistry=chemistry)]
  stopCluster(Clust)
  return(indexesCombinations)
}

searchOneSolution <- function(indexesList, index, nbLanes, multiplexingRate, unicityConstraint, chemistry){
  # goal: look for a solution (i.e. a combination of combination of indexes) such that
  #  - each index is used only once if required (unicityConstraint = index)
  #  - each combination of indexes is used only once if required (unicityConstraint = lane)
  # two steps:
  #  1) fill the lanes with as many samples per lane as in indexesList (may be lower than multiplexingRate for optimization purposes)
  #  2) if this number is lower than the desired multiplexing rate, complete the solution returned at step 1 adding indexes
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
      if (areIndexesCompatible(indexesList[[i]], chemistry)){
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
  solution <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane), 
                         pool=rep(1:nbLanes, each=inputNbSamplesPerLane), 
                         do.call("rbind", compatibleCombinations))
  if (multiplexingRate > inputNbSamplesPerLane){
    # only a partial solution has been found, need to complete it
    if (chemistry == "2") index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]
    solution <- completeSolution(partialSolution = solution,
                                 index = index,
                                 multiplexingRate = multiplexingRate,
                                 unicityConstraint = unicityConstraint)
  }
  return(solution)
}

completeSolution <- function(partialSolution, index, multiplexingRate, unicityConstraint){
  nbSamplesToAdd <- multiplexingRate - nrow(partialSolution)/max(partialSolution$pool) # to each lane
  for (l in unique(partialSolution$pool)){
    if (unicityConstraint == "index"){
      # remove all the indexes already used
      index2 <- index[!(index$id %in% partialSolution$id),]
    } else{
      # remove the indexes already used in the current lane
      index2 <- index[!(index$id %in% partialSolution$id[partialSolution$pool==l]),]
    }
    if (nbSamplesToAdd > nrow(index2)) return(NULL) # not enough remaining indexes to complete the solution
    indexesToAdd <- index2[sample(1:nrow(index2), nbSamplesToAdd, FALSE),]
    indexesToAdd$pool <- l
    print(partialSolution)
    print(indexesToAdd)
    partialSolution <- rbind.data.frame(partialSolution[, c("pool","id","sequence","color")],
                                        indexesToAdd[, c("pool","id","sequence","color")])
  }
  finalSolution <- data.frame(sample=1:nrow(partialSolution), partialSolution[order(partialSolution$pool, partialSolution$id),])
  return(finalSolution)
}

findSolution <- function(indexesList, index, nbSamples, multiplexingRate, unicityConstraint, 
                         nbMaxTrials, completeLane, selectCompIndexes, chemistry){
  # this function run searchOneSolution() nbMaxTrials times until finding a solution based on the parameters defined
  if (!unicityConstraint %in% c("none", "index", "lane")) stop("unicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indexes: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  nbLanes <- nbSamples/multiplexingRate
  if (unicityConstraint!="none" & completeLane & selectCompIndexes & length(indexesList) < nbLanes){
    stop("There are only ", length(indexesList), " combinations of compatible indexes to fill ", nbLanes, " lanes.")
  }
  nbTrials <- 1
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(indexesList = indexesList,
                                  index = index,
                                  nbLanes = nbLanes,
                                  multiplexingRate = multiplexingRate,
                                  unicityConstraint = unicityConstraint,
                                  chemistry = chemistry)
    if (!is.null(solution)){
      checkProposedSolution(solution = solution, unicityConstraint = unicityConstraint, chemistry = chemistry)
      return(solution)
    } else{
      nbTrials <- nbTrials + 1
    }
  }
  stop(paste("No solution found after", nbMaxTrials, "trials using these parameters, you can modify them or increase the number of trials."))
}

checkProposedSolution <- function(solution, unicityConstraint, chemistry){
  # indexes not unique
  if (unicityConstraint=="index" & any(duplicated(solution$id))){
    stop("The solution proposed uses some indexes several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indexes not unique within a pool/lane
  if (any(sapply(split(solution$id, solution$pool), function(x) any(duplicated(x))))){
    stop("The solution proposed uses some indexes several times within a lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # several pools/lanes with the same combination of indexes
  if (unicityConstraint=="lane" & any(duplicated(split(solution$id, solution$pool)))){
    stop("The solution proposed uses some combinations of indexes several times. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # different number of samples on the pools/lanes
  if (length(table(table(solution$pool)))>1){
    stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indexes not compatible on at least one pool/lane
  if (any(!sapply(split(solution, solution$pool), areIndexesCompatible, chemistry))){
    stop("The solution proposed uses incompatible indexes. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # two-channel chemistry and indexes starting with GG
  if (chemistry == "2" && any(sapply(solution$sequence, substr, 1, 2) == "GG")){
    stop("Indexes starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
}

heatmapindex <- function(solution){
  splitsol <- split(solution, solution$pool)
  # build a matrix containing the index bases and NAs to separate the pools/lanes
  seqmat <- lapply(splitsol, function(splitsol.pool) do.call("rbind", strsplit(splitsol.pool$sequence, "")))
  seqmat <- do.call("rbind", lapply(seqmat, function(seqmat.pool) rbind(seqmat.pool, rep(NA, ncol(seqmat.pool)))))
  seqmat <- seqmat[-nrow(seqmat),]
  # build a matrix containing the colors and NAs to separate the pools/lanes
  seqcol <- lapply(splitsol, function(splitsol.pool) do.call("rbind", strsplit(splitsol.pool$color, "")))
  seqcol <- do.call("rbind", lapply(seqcol, function(seqcol.pool) rbind(seqcol.pool, rep(NA, ncol(seqcol.pool)))))
  seqcol <- seqcol[-nrow(seqcol),]
  # plot
  par(mar=c(2, 6, 1, 6) + 0.1, xpd=TRUE, xaxs="i", yaxs="i")
  plot.new()
  plot.window(xlim=c(-1.5, ncol(seqmat)+1.5), ylim=c(-2, nrow(seqmat)))
  for (i in 1:nrow(seqmat)){
    for (j in 1:ncol(seqmat)){
      rect(xleft=j-1, ybottom=nrow(seqmat)-(i-1), xright=j, ytop=nrow(seqmat)-i,
           col=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white"),
           border=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white"))
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
