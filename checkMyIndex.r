#! /usr/bin/Rscript

# Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 8 -m 2
# Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 6 -m 2 -u index
# Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 12 -m 3 -u lane
# Rscript checkMyIndex.r --inputFile7=inputIndexesExample.txt -C 4 -n 24 -m 4 -u index
# Rscript checkMyIndex.r --inputFile7=index24-i7.txt --inputFile5=index24-i5.txt -C 2 -n 24 -m 12

library(optparse)

option_list <- list( 
  make_option("--inputFile7",
              type="character",
              dest="inputFile",
              help="two columns tab-delimited text file without header containing the available i7 index ids and sequences"),
  
  make_option("--inputFile5",
              type="character",
              default=NULL,
              dest="inputFile2",
              help="optional two columns tab-delimited text file without header containing the available i5 index ids
              and sequences (for dual-indexing sequencing experiments"),
  
  make_option(c("-n","--nbSamples"),
              type="integer",
              dest="nbSamples",
              help="total number of samples in the experiment (can be greater than the number of available indexes)"),
  
  make_option(c("-C","--chemistry"),
              type="integer",
              dest="chemistry",
              help="Illumina chemistry: either 1 (iSeq 100), 2 (NovaSeq, NextSeq & MiniSeq) or 4 (HiSeq & MiSeq) channel chemistry.
              With the four-channel chemistry A/C are red and G/T are green and indexes are compatible if there are at
              least one red light and one green light at each position. With the two-channel chemistry G has no color, 
              A is orange, C is red and T is green and indexes are compatible if there is at least one color at each
              position. Note that indexes starting with GG are not compatible with the two-channel chemistry. With the
              one-channel chemistry compatibility cannot be defined with colors and indexes are compatible is there is
              at least one A or C or T at each position.Please refer to the Illumina documentation for more details."),
  
  make_option(c("-m","--multiplexingRate"),
              type="integer",
              dest="multiplexingRate",
              help="multiplexing rate, i.e. number of samples per sequencing pool (must be a divisor of the total number of samples)"),
  
  make_option(c("-u","--unicityConstraint"),
              type="character",
              default="none",
              dest="unicityConstraint",
              help="(only for single-indexing) set to 'lane' to use each combination of compatible indexes only once or
              'index' to use each index only once. Note that it is automatically set to 'none' in case of dual-indexing [default: %default]"),
  
  make_option(c("-o","--outputFile"),
              type="character",
              default=NULL,
              dest="outputFile",
              help="name of the output file to save the proposed sequencing design [default: %default (no output file created)]"),

  make_option(c("-c","--complete"),
              type="logical",
              default=FALSE,
              action="store_true",
              dest="completeLane",
              help="(only for single-indexing) directly look for a solution with the desired multiplexing rate instead of
              looking for a solution with a few samples per pool/lane and add some of the remaining indexes to reach the
              desired multiplexing rate. Note that it is automatically set to FALSE in case of dual-indexing [default: %default]"),
  
  make_option(c("-s","--selectCompIndexes"),
              type="logical",
              default=FALSE,
              action="store_true",
              dest="selectCompIndexes",
              help="(only for single-indexing) select compatible indexes before looking for a solution (can take some time
              but then speed up the algorithm). Note that it is automatically set to FALSE in case of dual-indexing [default: %default]"),
  
  make_option(c("-b","--nbMaxTrials"),
              default=10,
              dest="nbMaxTrials", 
              help="maximum number of iterations to find a solution [default: %default]")
)

description <- c("Search for a set of compatible indexes for your sequencing experiment.

There can be many combinations of indexes to check according to the number of input indexes and the multiplexing rate.
Thus, testing for the compatibility of all the combinations may be long or even impossible. The trick is to find a partial
solution with the desired number of pools/lanes but with a lower multiplexing rate and then to complete each pool/lane with
some of the remaining indexes to reach the desired multiplexing rate. Indeed, adding indexes to a combination of compatible
indexes will give a compatible combination for sure. Briefly, a lower multiplexing rate generates a lower number of combinations
to test and thus make the research of a partial solution very fast. Adding some indexes to complete each pool/lane is fast too
and give the final solution.

Unfortunately, the research of a final solution might become impossible as the astuteness reduces the number of combinations
of indexes. In such a case, one can look for a solution using directly the desired multiplexing rate (see parameters),
the only risk is to increase the computational time.")

parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description=description,
                       epilogue="For comments, suggestions, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

inputFile <- opt$inputFile
inputFile2 <- opt$inputFile2
multiplexingRate <- as.numeric(opt$multiplexingRate)
nbSamples <- as.numeric(opt$nbSamples)
chemistry <- opt$chemistry
unicityConstraint <- opt$unicityConstraint
outputFile <- opt$outputFile
completeLane <- as.logical(opt$completeLane)
selectCompIndexes <- as.logical(opt$selectCompIndexes)
nbMaxTrials <- as.numeric(opt$nbMaxTrials)

source("global.r")

nbLanes <- nbSamples/multiplexingRate
index <- readIndexesFile(inputFile)
index <- addColors(index, chemistry)
nr <- nrow(index)
# dual-indexing
if (!is.null(inputFile2)){
  index2 <- readIndexesFile(inputFile2)
  index2 <- addColors(index2, chemistry)
  nr2 <- nrow(index2)
  unicityConstraint <- "none"
  completeLane <- FALSE
  selectCompIndexes <- FALSE
} else{
  nr2 <- 1
  index2 <- NULL
}

# some basic checkings
if (length(chemistry) != 1 || !I(chemistry %in% c("1", "2", "4"))) stop("\nChemistry must be equal to 1 (iSeq 100), 2 (NovaSeq, NextSeq & MiniSeq) or 4 (HiSeq & MiSeq).")
if (nbSamples %% 1 != 0 || nbSamples <= 1) stop("\nNumber of samples must be an integer greater than 1.")
if (nbSamples %% multiplexingRate != 0) stop("\nNumber of samples must be a multiple of the multiplexing rate.")
if (multiplexingRate > nr*nr2) stop("\nMultiplexing rate can't be higher than the number of input indexes.")

cat("--------------- Parameters ---------------\n")
cat("Input file i7: ", inputFile, " (", nrow(index), " indexes)\n", sep="")
cat("Input file i5: ", inputFile2, 
    ifelse(!is.null(index2), paste0("(", nrow(index2), " indexes)"), ""), "\n")
cat("Multiplexing rate:", multiplexingRate, "\n")
cat("Number of samples:", nbSamples, "\n")
cat("Chemistry:", switch(chemistry,
                         "1" = "one-channel (iSeq 100)",
                         "2" = "two-channels (NovaSeq, NextSeq & MiniSeq)", 
                         "4" = "four-channels (HiSeq & MiSeq)"), "\n")
cat("Number of pools/lanes:", nbLanes, "\n")
cat("Constraint:", ifelse(unicityConstraint=="none","none",
                          ifelse(unicityConstraint=="lane", "use each combination of compatible indexes only once", 
                                 "use each index only once")), "\n")
cat("Directly look for complete pools/lanes:", completeLane, "\n")
cat("Select compatible indexes before looking for a solution:", selectCompIndexes, "\n")
cat("Maximum number of iterations to find a solution:", nbMaxTrials, "\n")
cat("Output file:", outputFile, "\n")
cat("------------------------------------------\n")

# generate the list of indexes
indexesList <- generateListOfIndexesCombinations(index = index,
                                                 nbSamplesPerLane = multiplexingRate,
                                                 completeLane = completeLane,
                                                 selectCompIndexes = selectCompIndexes,
                                                 chemistry = chemistry)

if (!is.null(index2)){
  indexesList2 <- generateListOfIndexesCombinations(index = index2,
                                                   nbSamplesPerLane = multiplexingRate,
                                                   completeLane = completeLane,
                                                   selectCompIndexes = selectCompIndexes,
                                                   chemistry = chemistry)
} else{
  indexesList2 <- NULL
}

cat("Let's try to find a solution for", nbLanes, "pool(s) of", multiplexingRate, "samples using the specified parameters:\n\n")
print(solution <- findSolution(indexesList = indexesList,
                               index = index,
                               indexesList2 = indexesList2,
                               index2 = index2,
                               nbSamples = nbSamples,
                               multiplexingRate = multiplexingRate,
                               unicityConstraint = unicityConstraint,
                               nbMaxTrials = nbMaxTrials,
                               completeLane = completeLane,
                               selectCompIndexes = selectCompIndexes,
                               chemistry = chemistry), row.names=FALSE)

if (!is.null(outputFile)){
  write.table(solution, outputFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  cat(paste0("\nThe proposed sequencing design has been exported into ", outputFile))
}

cat("\nRun the program again to obtain another solution!\n")
