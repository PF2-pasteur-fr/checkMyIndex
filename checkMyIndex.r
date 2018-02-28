#! /usr/bin/Rscript

# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -C 4 -n 8 -m 2
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -C 4 -n 6 -m 2 -u index
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -C 4 -n 12 -m 3 -u lane
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -C 4 -n 24 -m 4 -u index

library(optparse)

option_list <- list( 
  make_option(c("-i","--inputFile"),
              type="character",
              dest="inputFile",
              help="two columns tab-delimited text file without header containing the available index ids and sequences"),
  
  make_option(c("-n","--nbSamples"),
              type="integer",
              dest="nbSamples",
              help="total number of samples in the experiment (can be greater than the number of available indexes)"),
  
  make_option(c("-C","--chemistry"),
              type="integer",
              dest="chemistry",
              help="Illumina chemistry: either 2 (NovaSeq, NextSeq & MiniSeq) or 4 (HiSeq & MiSeq) channel chemistry.
              With the four-channel chemistry A/C are red and G/T are green and indexes are compatible if there are at
              least one red light and one green light at each position. With the two-channel chemistry G has no color, 
              A is orange, C is red and T is green and indexes are compatible if there is at least one color at each
              position. Note that indexes starting with GG are not compatible with the two-channel chemistry. Please
              refer to the Illumina documentation for more details."),
  
  make_option(c("-m","--multiplexingRate"),
              type="integer",
              dest="multiplexingRate",
              help="multiplexing rate, i.e. number of samples per sequencing pool (must be a divisor of the total number of samples)"),
  
  make_option(c("-u","--unicityConstraint"),
              type="character",
              default="none",
              dest="unicityConstraint",
              help="set to 'lane' to use each combination of compatible indexes only once or 'index' to use each index only once [default: %default]"),
  
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
              help="directly look for a solution with the desired multiplexing rate instead of looking for a solution with 
              a few samples per pool/lane and add some of the remaining indexes to reach the desired multiplexing rate [default: %default]"),
  
  make_option(c("-s","--selectCompIndexes"),
              type="logical",
              default=FALSE,
              action="store_true",
              dest="selectCompIndexes",
              help="select compatible indexes before looking for a solution (can take some time but then speed up the algorithm) [default: %default]"),
  
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

# some basic checkings
if (length(chemistry) != 1 || !I(chemistry %in% c("2", "4"))) stop("\nChemistry must be equal to 2 (NovaSeq, NextSeq & MiniSeq) or 4 (HiSeq & MiSeq).")
if (nbSamples %% 1 != 0 || nbSamples <= 1) stop("\nNumber of samples must be an integer greater than 1.")
if (nbSamples %% multiplexingRate != 0) stop("\nNumber of samples must be a multiple of the multiplexing rate.")
if (multiplexingRate > nrow(index)) stop("\nMultiplexing rate can't be higher than the number of input indexes.")

cat("--------------- Parameters ---------------\n")
cat("Input file:", inputFile,"\n")
cat("Multiplexing rate:", multiplexingRate,"\n")
cat("Number of samples:", nbSamples,"\n")
cat("Chemistry:", ifelse(chemistry==2, "two-channels (NovaSeq, NextSeq & MiniSeq)", "four-channels (HiSeq & MiSeq)"), "\n")
cat("Number of pools/lanes:", nbLanes,"\n")
cat("Constraint:", ifelse(unicityConstraint=="none","none",
                          ifelse(unicityConstraint=="lane", "use each combination of compatible indexes only once", 
                                 "use each index only once")),"\n")
cat("Directly look for complete pools/lanes:", completeLane, "\n")
cat("Select compatible indexes before looking for a solution:", selectCompIndexes, "\n")
cat("Maximum number of iterations to find a solution:", nbMaxTrials,"\n")
cat("Output file:", outputFile,"\n")
cat("------------------------------------------\n")

# how many possible combinations (combinatory logic using C^n_k)
cat("\nIn the input list of", nrow(index), "indexes:", choose(n=nrow(index), k=multiplexingRate), "possible combinations of", multiplexingRate, "indexes (not necessarily compatible)\n")

# generate the list of indexes
indexesList <- generateListOfIndexesCombinations(index = index,
                                                 nbSamplesPerLane = multiplexingRate,
                                                 completeLane = completeLane,
                                                 selectCompIndexes = selectCompIndexes,
                                                 chemistry = chemistry)

# print the number of compatible combinations
if (nrow(indexesList[[1]]) == multiplexingRate & selectCompIndexes){
  cat("Among them", length(indexesList), "contain compatible indexes.\n")
} else{
  cat("Note: the number of combinations containing compatible indexes cannot be calculated with the parameters used.\n")
}

cat("Let's try to find a solution for", nbLanes, "pool(s) of", multiplexingRate, "samples using the specified parameters:\n\n")
print(solution <- findSolution(indexesList = indexesList,
                               index = index,
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
