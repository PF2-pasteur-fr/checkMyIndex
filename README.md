# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment

## Usage

### Shiny application

**Public website**

The last version of *checkMyIndex* is available on <http://shiny01.hosting.pasteur.fr/checkMyIndex/>.

**Locally**

The `shiny` R package must be available, run `install.packages(shiny)` in R to install it.

Within R, run the following commands to launch the shiny application:

`library(shiny)`
`runGitHub("PF2-pasteur-fr/checkMyIndex", launch.browser=TRUE)`

### Command line

*checkMyIndex* be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources the `global.r`, both files must be placed in the same directory.

**Example**:

List of 9 indexes for 9 samples distributed on 3 lanes:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9`

List of 12 indexes for 12 samples distributed on 4 lanes using each lane combination only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9 --uniqueCombination`

List of 12 indexes for 12 samples distributed on 4 lanes using each index only once:

`Rscript checkMyIndex.r --inputFile=inputIndexesExample.txt --multiplexingRate=3 --nbSamples=9 --uniqueIndexes`

**To get some help**:

`Rscript checkMyIndex.r --help`

**Note**: the `optparse` R package must be available to interpret the input parameters, run `install.packages(optparse)` in R to install it. 

## About checkMyIndex

This tool has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet (<hugo.varet@pasteur.fr>).
