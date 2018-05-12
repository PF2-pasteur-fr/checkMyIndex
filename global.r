# this file contains the functions used in both the "checkMyIndex" script and shiny application
library(parallel)

readIndexesFile <- function(file){
  index <- read.table(file, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","sequence"))
  checkInputIndexes(index)
  return(index)
}

addColors <- function(index, chemistry){
  # convert bases into colors
  if (chemistry == "4"){
    # A/C are red and G/T are green
    index$color <- gsub("G|T", "G", gsub("A|C", "R", index$sequence))
  }
  if (chemistry == "2"){
    # G has no color, A is orange, C is red and T is green
    index$color <- gsub("T", "G", gsub("C", "R", gsub("A", "O", gsub("G", "-", index$sequence))))
  } 
  if (chemistry == "1"){
    # G has no color, A is orange, C is red and T is green
    index$color <- gsub("A", "G", gsub("C", "R", gsub("T", "O", gsub("G", "-", index$sequence))))
  } 
  return(index)
}

checkInputIndexes <- function(index){
  if (any(duplicated(index$id))){
    stop("Input index ids are not unique, please check your input file.")
  }
  if (any(duplicated(index$sequence))){
    stop("Input indexes are not unique, please check your input file.")
  }
  if (!(all(unlist(strsplit(index$sequence,"")) %in% c("A","C","T","G")))){
    stop("Input indexes contain other characters than A, T, C and G, please check your input file.")
  }
  if (length(table(sapply(index$sequence, nchar))) > 1){
    stop("Input indexes do not have all the same length, please check your input file.")
  }
}

areIndexesCompatible <- function(index, chemistry){
  # return TRUE if the input indexes are compatible (i.e. can be used within the same pool/lane)
  if (nrow(index)==1) return(TRUE)
  matColors <- do.call("rbind", strsplit(index$color, ""))
  if (chemistry == "4"){
    sumRed <- apply(matColors, 2, function(x) sum(x=="R"))
    sumGreen <- nrow(matColors) - sumRed
    return(all(sumRed >= 1 & sumGreen >= 1))
  }
  if (chemistry %in% c("1","2")){
    sumNoColor <- apply(matColors, 2, function(x) sum(x=="-"))
    return(all(sumNoColor < nrow(index)))
  }
}

generateListOfIndexesCombinations <- function(index, nbSamplesPerLane, completeLane, selectCompIndexes, chemistry){
  # remove indexes starting with GG if two-channel chemistry
  if (chemistry == "2") index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]
  nbSamplesPerLane <- min(nrow(index), nbSamplesPerLane)
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

searchOneSolution <- function(indexesList, index, indexesList2=NULL, index2=NULL,
                              nbLanes, multiplexingRate, unicityConstraint, chemistry){
  # goal: look for a solution (i.e. a combination of combination of indexes) such that
  #  - each index is used only once if required (unicityConstraint = index)
  #  - each combination of indexes is used only once if required (unicityConstraint = lane)
  # two steps:
  #  1) fill the lanes with as many samples per lane as in indexesList (may be lower than multiplexingRate for optimization purposes)
  #  2) if this number is lower than the desired multiplexing rate, complete the solution returned at step 1 adding indexes
  # this function can return NULL if no solution is found (need to re-run in that case)
  inputNbSamplesPerLane <- nrow(indexesList[[1]])
  compatibleCombinations <- vector(mode="list", length=nbLanes)
  # single-indexing
  if (is.null(index2) | is.null(indexesList2)){
    k <- 1
    while (k <= nbLanes){
      if (length(indexesList) == 0){
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
  # dual-indexing
  } else {
    # reminder: unicityConstraint has been set to "none" to simplify the algorithm
    # and also because it is not so important with dual-indexing
    inputNbSamplesPerLane2 <- nrow(indexesList2[[1]])
    compatibleCombinations2 <- vector(mode="list", length=nbLanes)
    k <- 1
    while (k <= nbLanes){
      if (length(indexesList) == 0 | length(indexesList2) == 0){
        # no available combination of indexes anymore, try again!
        return(NULL)
      } else{
        i <- sample(1:length(indexesList), 1, FALSE)
        if (areIndexesCompatible(indexesList[[i]], chemistry)){
          compatibleCombinations[[k]] <- indexesList[[i]]
          for (i2 in sample(1:length(indexesList2), 1, FALSE)){
            if (areIndexesCompatible(indexesList2[[i2]], chemistry)){
              compatibleCombinations2[[k]] <- indexesList2[[i2]]
              k <- k+1
              break
            } else{
              indexesList2 <- indexesList2[-i2]
            }
          }
        } else{
          indexesList <- indexesList[-i]
        }
      }
    }
    solution <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane), 
                           pool=rep(1:nbLanes, each=inputNbSamplesPerLane), 
                           do.call("rbind", compatibleCombinations))
    solution2 <- data.frame(sample=1:(nbLanes*inputNbSamplesPerLane2), 
                            pool=rep(1:nbLanes, each=inputNbSamplesPerLane2), 
                            do.call("rbind", compatibleCombinations2))
    # remove indexes starting with GG before completing the solutions if four-channel chemistry
    if (chemistry == "2"){
      index <- index[!sapply(index$sequence, substr, 1, 2) == "GG",]
      index2 <- index2[!sapply(index2$sequence, substr, 1, 2) == "GG",]
    }
    # only a partial solution has been found, need to complete it
    if (min(multiplexingRate, nrow(index)) > inputNbSamplesPerLane){
      solution <- completeSolution(partialSolution = solution,
                                   index = index,
                                   multiplexingRate = min(multiplexingRate, nrow(index)),
                                   unicityConstraint = unicityConstraint)
    }
    if (min(multiplexingRate, nrow(index2)) > inputNbSamplesPerLane2){
      solution2 <- completeSolution(partialSolution = solution2,
                                    index = index2,
                                    multiplexingRate = min(multiplexingRate, nrow(index2)),
                                    unicityConstraint = unicityConstraint)
    }
    names(solution) <- paste0(names(solution), "1")
    names(solution2) <- paste0(names(solution2), "2")
    # generate all the possible couple of indexes
    solution.merged <- merge(solution, solution2, by.x="pool1", by.y="pool2")
    solution <- list()
    for (pool in unique(solution.merged$pool1)){
      # look for the most diverse couple of indexes can beslightly difficult
      # we thus encapsulate it into a while() loop
      solution.pool.OK <- FALSE
      while (!solution.pool.OK){
        solution.merged.pool <- solution.merged[which(solution.merged$pool1 == pool),]
        counts.id1 <- numeric(length(unique(solution.merged.pool$id1)))
        names(counts.id1) <- unique(solution.merged.pool$id1)
        counts.id2 <- numeric(length(unique(solution.merged.pool$id2)))
        names(counts.id2) <- unique(solution.merged.pool$id2)
        solution.tmp.pool <- NULL
        for (i in 1:multiplexingRate){
          for (k in sample(1:nrow(solution.merged.pool), nrow(solution.merged.pool), FALSE)){
            counts.id1.k <- counts.id1 + I(names(counts.id1) == solution.merged.pool[k, "id1"])
            counts.id2.k <- counts.id2 + I(names(counts.id2) == solution.merged.pool[k, "id2"])
            maxdiff <- function(x) max(x) - min(x)
            # make sure no index is used n+2 times and another n times
            # i.e. indexes must be used a homogeneous number of times
            if (maxdiff(counts.id1.k) <= 1 & maxdiff(counts.id2.k) <= 1){
              counts.id1 <- counts.id1.k
              counts.id2 <- counts.id2.k
              solution.tmp.pool <- rbind(solution.tmp.pool, solution.merged.pool[k,])
              solution.merged.pool <- solution.merged.pool[-k,]
              break
            }
          }
        }
        if (nrow(solution.tmp.pool) == multiplexingRate) solution.pool.OK <- TRUE
      }
      solution[[pool]] <- solution.tmp.pool
    }
    solution <- do.call("rbind", solution)
    solution <- solution[order(solution$pool1, solution$id1, solution$id2),]
    solution <- data.frame(sample=1:(nbLanes*multiplexingRate),
                           pool=solution$pool1,
                           id1=solution$id1,
                           sequence1=solution$sequence1,
                           color1=solution$color1,
                           id2=solution$id2,
                           sequence2=solution$sequence2,
                           color2=solution$color2,
                           stringsAsFactors = FALSE)
    return(solution)
  }
}

completeSolution <- function(partialSolution, index, multiplexingRate, unicityConstraint){
  nbSamplesToAdd <- multiplexingRate - nrow(partialSolution)/max(partialSolution$pool) # to each lane
  for (l in unique(partialSolution$pool)){
    if (unicityConstraint == "index"){
      # remove all the indexes already used
      index.remaining <- index[!(index$id %in% partialSolution$id),]
    } else{
      # remove the indexes already used in the current lane
      index.remaining <- index[!(index$id %in% partialSolution$id[partialSolution$pool==l]),]
    }
    if (nbSamplesToAdd > nrow(index.remaining)) return(NULL) # not enough remaining indexes to complete the solution
    indexesToAdd <- index.remaining[sample(1:nrow(index.remaining), nbSamplesToAdd, FALSE),]
    indexesToAdd$pool <- l
    partialSolution <- rbind.data.frame(partialSolution[, c("pool","id","sequence","color")],
                                        indexesToAdd[, c("pool","id","sequence","color")])
  }
  finalSolution <- data.frame(sample=1:nrow(partialSolution), partialSolution[order(partialSolution$pool, partialSolution$id),])
  return(finalSolution)
}

findSolution <- function(indexesList, index, indexesList2=NULL, index2=NULL,
                         nbSamples, multiplexingRate, unicityConstraint, nbMaxTrials, 
                         completeLane, selectCompIndexes, chemistry){
  # this function run searchOneSolution() nbMaxTrials times until finding a solution based on the parameters defined
  if (!unicityConstraint %in% c("none", "index", "lane")) stop("unicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
  if (unicityConstraint=="index" & nbSamples > nrow(index)) stop("More samples than available indexes: cannot use each index only once.")
  if (nbSamples %% multiplexingRate != 0) stop("Number of samples must be a multiple of the multiplexing rate.")
  if (I(!is.null(index2) | !is.null(indexesList2)) & unicityConstraint != "none") stop("unicityConstraint parameter must be equal to 'none' with dual-indexing.")
  nbLanes <- nbSamples/multiplexingRate
  if (unicityConstraint!="none" & completeLane & selectCompIndexes & length(indexesList) < nbLanes){
    stop("There are only ", length(indexesList), " combinations of compatible indexes to fill ", nbLanes, " lanes.")
  }
  nbTrials <- 1
  while (nbTrials <= nbMaxTrials){
    solution <- searchOneSolution(indexesList = indexesList,
                                  index = index,
                                  indexesList2 = indexesList2,
                                  index2 = index2,
                                  nbLanes = nbLanes,
                                  multiplexingRate = multiplexingRate,
                                  unicityConstraint = unicityConstraint,
                                  chemistry = chemistry)
    if (!is.null(solution)){
      # check solution only for single-indexing as dual-indexing has no constraint
      if (is.null(indexesList2) & is.null(index2)){
        checkProposedSolution(solution = solution, unicityConstraint = unicityConstraint, chemistry = chemistry)
      } else{
        checkProposedSolution2(solution = solution, chemistry = chemistry)
      }
      return(solution)
    } else{
      nbTrials <- nbTrials + 1
    }
  }
  stop(paste("No solution found after", nbMaxTrials, "trials using these parameters, you can modify them or increase the number of trials."))
}

# single-indexing solution checking
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
  if (length(table(table(solution$pool))) > 1){
    stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indexes not compatible in at least one pool/lane
  if (any(!sapply(split(solution, solution$pool), areIndexesCompatible, chemistry))){
    stop("The solution proposed uses incompatible indexes. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # two-channel chemistry and indexes starting with GG
  if (chemistry == "2" && any(sapply(solution$sequence, substr, 1, 2) == "GG")){
    stop("Indexes starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # one-channel chemistry and all indexes starting with GG
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indexes of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
}

# dual-indexing solution checking
checkProposedSolution2 <- function(solution, chemistry){
  # indexes not unique within a pool/lane
  if (any(sapply(split(paste(solution$id1, solution$id2, sep="-"), solution$pool), function(x) any(duplicated(x))))){
    stop("The solution proposed uses some dual-index combinations several times within a lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # different number of samples on the pools/lanes
  if (length(table(table(solution$pool))) > 1){
    stop("The solution proposed uses different numbers of samples per lane. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # indexes not compatible in at least one pool/lane
  sol1 <- data.frame(id=solution$id1, sequence=solution$sequence1, color=solution$color1, stringsAsFactors = FALSE)
  if (any(!sapply(split(sol1, solution$pool), areIndexesCompatible, chemistry))){
    stop("The solution proposed uses incompatible indexes. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  sol2 <- data.frame(id=solution$id2, sequence=solution$sequence2, color=solution$color2, stringsAsFactors = FALSE)
  if (any(!sapply(split(sol2, solution$pool), areIndexesCompatible, chemistry))){
    stop("The solution proposed uses incompatible indexes. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # two-channel chemistry and indexes starting with GG
  if (chemistry == "2" && any(sapply(solution$sequence1, substr, 1, 2) == "GG")){
    stop("Indexes starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  if (chemistry == "2" && any(sapply(solution$sequence2, substr, 1, 2) == "GG")){
    stop("Indexes starting with GG are not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  # one-channel chemistry and all indexes starting with GG
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence1, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indexes of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
  if (chemistry == "1" && any(sapply(lapply(split(solution$sequence2, solution$pool), substr, 1, 2), function(x) all(x=="GG")))){
    stop("Having all the indexes of a pool starting with GG is not compatible with the chosen chemistry. Thanks to report this error to Hugo Varet (hugo.varet@pasteur.fr)")
  }
}

heatmapindex <- function(solution){
  if ("id1" %in% names(solution)){
    # dual indexing, need to format the solution data.frame as for single-indexing
    solution <- data.frame(sample=solution$sample, pool=solution$pool,
                           id=paste(solution$id1, solution$id2, sep="-"),
                           sequence=paste(solution$sequence1, solution$sequence2),
                           color=paste(solution$color1, solution$color2),
                           stringsAsFactors = FALSE)
  }
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
           col=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white", " "="white"),
           border=switch(seqcol[i,j], "R"="orangered", "G"="darkseagreen4", "-"="white", "O"="orange", "NA"="white", " "="white"))
      text(x=j-0.5, y=nrow(seqmat)-(i-0.5), labels=seqmat[i,j])
    }
  }
  # extract and print sample ids
  tmp <- unlist(lapply(splitsol, function(sol) c(sol$sample, NA)))
  text(x=-0.25, y=nrow(seqmat):1 - 0.5, labels=tmp[-length(tmp)], pos=2)
  text(x=-1.5, y=nrow(seqmat)/2, labels="Sample", srt=90)
  # extract and print index ids
  tmp <- unlist(lapply(splitsol, function(sol) c(sol$id, NA)))
  text(x=ncol(seqmat)+0.25, y=nrow(seqmat):1 - 0.5, labels=tmp[-length(tmp)], pos=4)
}
