#! /usr/bin/Rscript

# Rscript checkMyIndex.r --inputFile=inputIndexExamples/TruSeqRNA_v2_24index.txt --multiplexingRate=2 --nbSamples=6 --uniqueIndexes --uniqueCombinations
# Rscript checkMyIndex.r --inputFile=inputIndexExamples/TruSeqRNA_v2_24index.txt --multiplexingRate=2 --nbSamples=8
# Rscript checkMyIndex.r --inputFile=inputIndexExamples/TruSeqRNA_v2_24index.txt --multiplexingRate=3 --nbSamples=12 --uniqueIndexes --uniqueCombinations
# Rscript checkMyIndex.r --inputFile=inputIndexExamples/TruSeqRNA_v2_24index.txt --multiplexingRate=4 --nbSamples=24 --uniqueIndexes --uniqueCombinations --nbMaxTrials=5

library(optparse)

option_list <- list( 
  make_option("--inputFile",
              type="character",
              dest="inputFile",
              help="two columns tab-delimited text file without header containing the index ids and sequences"),
  
  make_option("--multiplexingRate",
              type="integer",
              dest="multiplexingRate",
              help="multiplexing rate, i.e. number of samples per lane"),
  
  make_option("--nbSamples",
              type="integer",
              dest="nbSamples",
              help="total number of samples in the experiment"),

  make_option("--uniqueIndexes",
              type="logical",
              action="store_true",
              default=FALSE,
              dest="uniqueIndexes", 
              help="use each index only once"),
  
  make_option("--uniqueCombinations",
              type="logical",
              action="store_true",
              default=FALSE,
              dest="uniqueCombinations",
              help="use each combination only once"),
  
  make_option("--nbMaxTrials",
              default=1000,
              dest="nbMaxTrials", 
              help="maximum number of iterations to find a solution [default: %default]")
)

# now parse the command line to check which option is given and get associated values
parser <- OptionParser(usage="usage: %prog [options]",
                       option_list=option_list, 
                       description="Search for a set of compatible indexes for your sequencing experiment.",
                       epilogue="For comments, bug reports etc... please contact Hugo Varet <hugo.varet@pasteur.fr>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=TRUE), positional_arguments=0)$options

inputFile <- opt$inputFile
multiplexingRate <- as.numeric(opt$multiplexingRate)
nbSamples <- as.numeric(opt$nbSamples)
uniqueIndexes <- as.logical(opt$uniqueIndexes)
uniqueCombinations <- as.logical(opt$uniqueCombinations)
nbMaxTrials <- as.numeric(opt$nbMaxTrials)

source("global.r")

nbLanes <- nbSamples/multiplexingRate
index <- read.table(inputFile, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","seq"))

# some basic checkings
checkInputIndexes(index)
if (nbSamples <= 1) stop("\nNumber of samples must be greater than 1.")
if (nbSamples > nrow(index)) stop("\nMore samples than available indexes.")
if (nbSamples %% multiplexingRate != 0) stop("\nNumber of samples must be a multiple of the multiplexing rate.")
if (uniqueIndexes & !uniqueCombinations){
  uniqueCombinations <- TRUE
  warning("uniqueCombinations has been set to TRUE as uniqueIndexes is TRUE")
}

cat("--------------- Parameters ---------------\n")
cat("Input file:", inputFile,"\n")
cat("List of input indexes:\n")
print(index, row.names=FALSE)
cat("Multiplexing rate:", multiplexingRate,"\n")
cat("Number of samples:", nbSamples,"\n")
cat("Number of lanes:", nbLanes,"\n")
cat("Use only unique indexes:", uniqueIndexes,"\n")
cat("Use unique combinations:", uniqueCombinations,"\n")
cat("------------------------------------------\n")

# how many possible combinations (combinatory logic using C^n_k)
cat("\nIn the input list of", nrow(index), "indexes:", choose(n=nrow(index), k=multiplexingRate), "possible combinations of", multiplexingRate, "indexes (not necessarily compatible)\n")

compatibleIndexes <- searchCompatibleIndexes(index, nbSamplesPerLane=multiplexingRate)
cat("In the input list of", nrow(index), "indexes:", length(compatibleIndexes), "combinations of", multiplexingRate, "compatibles indexes\n\n")

cat(paste0("Let's try to find a solution for ", nbLanes, " lanes of ", multiplexingRate, " samples using:\n - ",
           nbSamples, ifelse(!uniqueIndexes," non",""), " unique indexes\n",
           " - each combination ", ifelse(uniqueCombinations, "only once", "potentially several times"), "\n\n"))
print(findSolution(compatibleIndexes, index, nbSamples, multiplexingRate, uniqueIndexes, uniqueCombinations, nbMaxTrials), row.names=FALSE)

cat("\nRun the program again to obtain another solution!\n")
