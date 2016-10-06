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
outputFile <- opt$outputFile
nbMaxTrials <- as.numeric(opt$nbMaxTrials)

source("global.r")

nbLanes <- nbSamples/multiplexingRate
index <- read.table(inputFile, header=FALSE, sep="\t", stringsAsFactors=FALSE, col.names=c("id","seq"))

# some basic checkings
checkInputIndexes(index)
if (nbSamples <= 1) stop("\nNumber of samples must be greater than 1.")
if (!unicityConstraint %in% c("none","index","lane")) stop("\nunicityConstraint parameter must be equal to 'none', 'lane' or 'index'.")
if (nbSamples %% multiplexingRate != 0) stop("\nNumber of samples must be a multiple of the multiplexing rate.")

cat("--------------- Parameters ---------------\n")
cat("Input file:", inputFile,"\n")
cat("List of input indexes:\n")
print(index, row.names=FALSE)
cat("Multiplexing rate:", multiplexingRate,"\n")
cat("Number of samples:", nbSamples,"\n")
cat("Number of lanes:", nbLanes,"\n")
cat("Constraint:", ifelse(unicityConstraint=="none","none",
                          ifelse(unicityConstraint=="lane", "use each combination of compatible indexes only once", "use each index only once")),"\n")
cat("Output file:", outputFile,"\n")
cat("Maximum number of iterations to find a solution:", nbMaxTrials,"\n")
cat("------------------------------------------\n")

# how many possible combinations (combinatory logic using C^n_k)
cat("\nIn the input list of", nrow(index), "indexes:", choose(n=nrow(index), k=multiplexingRate), "possible combinations of", multiplexingRate, "indexes (not necessarily compatible)\n")

compatibleIndexes <- searchCompatibleIndexes(index, nbSamplesPerLane=multiplexingRate)
cat("In the input list of", nrow(index), "indexes:", length(compatibleIndexes), "combinations of", multiplexingRate, "compatibles indexes\n\n")

cat(paste0("Let's try to find a solution for ", nbLanes, " lanes of ", multiplexingRate, " samples using:\n - ",
           nbSamples, ifelse(unicityConstraint!="index"," non",""), " unique indexes\n",
           " - each combination ", ifelse(unicityConstraint!="none", "only once", "potentially several times"), "\n\n"))
print(solution <- findSolution(compatibleIndexes, index, nbSamples, multiplexingRate, unicityConstraint, nbMaxTrials), row.names=FALSE)

if (!is.null(outputFile)){
  write.table(solution, outputFile, col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
  cat(paste0("\nThe proposed sequencing design has been exported into ", outputFile))
}

cat("\nRun the program again to obtain another solution!\n")
