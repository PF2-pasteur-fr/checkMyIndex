# checkMyIndex

Search for a set of compatible indexes for your sequencing experiment according to:

* the number of samples
* the desired multiplexing rate (i.e. number of samples per lane)
* the constraint on the indexes (none, use each one or each combination only once)
* the minimum number of red/green lights required at each position

## Input indexes file

The list of the available indexes must be stored in a text file containing two tab-separated columns (without header): index ids are in the first column and the corresponding sequences in the second. `inputIndexesExample.txt` is an example of such a file and can be used to test *checkMyIndex*.

## Shiny application

### Public website

Click [here](http://shiny01.hosting.pasteur.fr/checkMyIndex/) to use the shiny interface of *checkMyIndex*.

### Locally

One can use the application locally running the two following lines in R:

`library(shiny)`

`runGitHub("PF2-pasteur-fr/checkMyIndex", launch.browser=TRUE)`

## Rscript command line

*checkMyIndex* can be executed calling `checkMyIndex.r` with `Rscript`. As `checkMyIndex.r` sources `global.r`, both files must be placed in the same directory. Here are 3 examples using the input example file `inputIndexesExample.txt`:

* List of 9 indexes for 9 samples distributed on 3 lanes:

`Rscript checkMyIndex.r -i inputIndexesExample.txt -n 9 -m 3`

* List of 12 indexes for 12 samples distributed on 4 lanes using each lane combination only once:

`Rscript checkMyIndex.r -i inputIndexesExample.txt -n 12 -m 3 -u lane`

* List of 12 indexes for 12 samples distributed on 4 lanes using each index only once:

`Rscript checkMyIndex.r -i inputIndexesExample.txt -n 12 -m 3 -u index`

The help page of the script can be displayed running the following command: 

`Rscript checkMyIndex.r --help`

## Requirements

Here is the list of the R packages needed to run *checkMyIndex*:

* *shiny* to run the shiny application locally

* *optparse* to interpret the input parameters when using the Rscript command line

* *parallel* to speed up the calculations

One can install each of these packages running `install.packages(packageName)` in R.

## About checkMyIndex

This tool has been developed at the Transcriptome & Epigenome Platform of the Biomics pole by Hugo Varet (<hugo.varet@pasteur.fr>).
