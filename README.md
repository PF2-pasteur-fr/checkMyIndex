# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment according to:

* the number of samples
* the desired multiplexing rate (i.e. number of samples per lane)
* the constraint on the indexes (none, use each one or each combination only once)

## Input indexes file

The list of the available indexes must be stored in a text file containing two tab-separated columns (without header): index ids are in the first column and the corresponding sequences in the second. `inputIndexesExample.txt` is an example of such a file and can be used to test *checkMyIndex*.

## Shiny application

### Public website

Click [here](http://shiny01.hosting.pasteur.fr/checkMyIndex/) to use the shiny interface of *checkMyIndex*.

### Locally

One can use the application locally running the two following lines in R:

`library(shiny)`

`runGitHub("PF2-pasteur-fr/checkMyIndex", launch.browser=TRUE)`

**Note**: the *shiny* R package must be available, run `install.packages(shiny)` in R to install it.

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory.

**Note**: the *optparse* R package must be available to interpret the input parameters, run `install.packages(optparse)` in R to install it. 

Here are 3 examples using the input example file `inputIndexesExample.txt`:

* List of 9 indexes for 9 samples distributed on 3 lanes:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --nbSamples=9 --multiplexingRate=3`

* List of 12 indexes for 12 samples distributed on 4 lanes using each lane combination only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --nbSamples=12 --multiplexingRate=3 --unicityConstraint=lane`

* List of 12 indexes for 12 samples distributed on 4 lanes using each index only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --nbSamples=12 --multiplexingRate=3 --unicityConstraint=index`

The help page of the script can be displayed running the following command: 

`Rscript checkMyIndex.r --help`

## About checkMyIndex

This tool has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet (<hugo.varet@pasteur.fr>).
