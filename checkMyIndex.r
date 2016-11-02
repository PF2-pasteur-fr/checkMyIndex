#! /usr/bin/Rscript

# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -n 8 -m 2 -u
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -n 6 -m 2 -u index
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -n 12 -m 3 -u lane
# Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt -n 24 -m 4 -u index -b 20

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
  
  make_option(c("-m","--multiplexingRate"),
              type="integer",
              dest="multiplexingRate",
              help="multiplexing rate, i.e. number of samples per lane (must be a divisor of the total number of samples)"),
  
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
  
  make_option(c("-r","--minRedGreen"),
              type="integer",
              default=1,
              dest="minRedGreen",
              help="minimum number of red/green lights required at each position [default: %default]"),
  
  make_option(c("-s","--selectCompIndexes"),
              default=FALSE,
              action="store_true",
              dest="selectCompIndexes", 
              help="use this option to select the compatible indexes before looking for a solution (can take some time but be informative and then speed up the research of a solution)"),
  
  make_option(c("-b","--nbMaxTrials"),
              default=10000,
              dest="nbMaxTrials", 
              help="maximum number of iterations to find a solution [default: %default]")
)

parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description="Search for a set of compatible indexes for your sequencing experiment.",
                       epilogue="For comments, suggestions, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

inputFile <- opt$inputFile
multiplexingRate <- as.numeric(opt$multiplexingRate)
nbSamples <- as.numeric(opt$nbSamples)
unicityConstraint <- opt$unicityConstraint
minRedGreen <- as.numeric(opt$minRedGreen)
outputFile <- opt$outputFile
nbMaxTrials <- as.numeric(opt$nbMaxTrials)
selectCompIndexes <- as.logical(opt$selectCompIndexes)

source("global.r")

nbLanes <- nbSamples/multiplexingRate
index <- readIndexesFile(inputFile)

# some basic checkings
checkInputIndexes(index)
if (nbSamples %% 1 != 0 || nbSamples <= 1) stop("\nNumber of samples must be an integer greater than 1.")
if (nbSamples %% multiplexingRate != 0) stop("\nNumber of samples must be a multiple of the multiplexing rate.")
if (minRedGreen > multiplexingRate/2) stop("\nMinimal number of red/green lights per position can't be higher than the multiplexing rate divided by 2.")

cat("--------------- Parameters ---------------\n")
cat("Input file:", inputFile,"\n")
cat("Multiplexing rate:", multiplexingRate,"\n")
cat("Number of samples:", nbSamples,"\n")
cat("Number of lanes:", nbLanes,"\n")
cat("Constraint:", ifelse(unicityConstraint=="none","none",
                          ifelse(unicityConstraint=="lane", "use each combination of compatible indexes only once", 
                                 "use each index only once")),"\n")
cat("Minimal number of red/green lights per position:", minRedGreen,"\n")
cat("Output file:", outputFile,"\n")
cat("Maximum number of iterations to find a solution:", nbMaxTrials,"\n")
cat("Prior selection of the compatible indexes:", selectCompIndexes,"\n")
cat("------------------------------------------\n")

# how many possible combinations (combinatory logic using C^n_k)
cat("\nIn the input list of", nrow(index), "indexes:", choose(n=nrow(index), k=multiplexingRate), "possible combinations of", multiplexingRate, "indexes (not necessarily compatible)\n")

# generate the list of indexes
indexesList <- generateListOfIndexesCombinations(index, multiplexingRate, minRedGreen, selectCompIndexes)
if (selectCompIndexes & nrow(indexesList[[1]])==multiplexingRate) cat("Among them", length(indexesList), "contain compatible indexes (i.e. at least", minRedGreen, "red/green light(s) per position)\n")

cat("Let's try to find a solution for", nbLanes, "lanes of", multiplexingRate, "samples using the specified parameters\n\n")
print(solution <- findSolution(indexesList, index, nbSamples, multiplexingRate, unicityConstraint, minRedGreen, nbMaxTrials, selectCompIndexes), row.names=FALSE)

if (!is.null(outputFile)){
  write.table(solution, outputFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  cat(paste0("\nThe proposed sequencing design has been exported into ", outputFile))
}

cat("\nRun the program again to obtain another solution!\n")
