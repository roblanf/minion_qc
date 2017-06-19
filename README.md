# Quality control for MinION sequencing data

## What?

An R script to do some basic QC on data from Oxford Nanopore's MinION sequencer, using the `sequencing_summary.txt` file from Albacore as input.


## Quick start

```
Rscript minion_QC.R sequencing_summary.txt output_directory
```

*minion_QC.R*: path to this script
*sequencing_summary.txt*: path to a sequencing_summary.txt file from Albacore
*output_directory*: path to an output directory. Files will be overwritten.

## Dependencies
A recent version of R, and install the following:

```
install.packages("ggplot2")
install.packages("viridis")
install.packages("reshape2")
install.packages("plyr")
```

