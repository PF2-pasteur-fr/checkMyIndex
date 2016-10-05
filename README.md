# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment according to the number of samples and the desired multiplexing rate (i.e. number of samples per lane).

## Input indexes file

The list of the available indexes must be stored in a text file containing two tab-separated columns (without header): index ids are in the first column and the corresponding sequences in the second. `inputIndexesExample.txt` is an example of such a file and can be used to test *checkMyIndex*.

## Shiny application

**Public website**

The last version of *checkMyIndex* is available on <http://shiny01.hosting.pasteur.fr/checkMyIndex/>.

**Locally**

The *shiny* R package must be available to run the *checkMyIndex* interface locally, run `install.packages(shiny)` in R to install it.

Within R, run the following commands to launch the application:

`library(shiny)

runGitHub("PF2-pasteur-fr/checkMyIndex", launch.browser=TRUE)`

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory.

**Note**: the *optparse* R package must be available to interpret the input parameters, run `install.packages(optparse)` in R to install it. 

**Examples**:

* List of 9 indexes for 9 samples distributed on 3 lanes:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9`

* List of 12 indexes for 12 samples distributed on 4 lanes using each lane combination only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9 --uniqueCombination`

* List of 12 indexes for 12 samples distributed on 4 lanes using each index only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9 --uniqueIndexes`

**To get some help**:

Run the following command line to print the help page of the script: 

`Rscript checkMyIndex.r --help`

## About checkMyIndex

This tool has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet (<hugo.varet@pasteur.fr>).
