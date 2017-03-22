# Quality control for MinION sequencing data

## What?

We're just getting started with MinION data, so we thought we'd share our QC scripts

## Why? 

MinION data is great.

## Getting started

1. Install the following software: 

	* *fastqc*: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
	* *poretools*: http://poretools.readthedocs.io/en/latest/content/installation.html

2. Run the scripts like this:

```
sh minion_qc.sh
```

The scripts will write a bunch of stuff to the output directories you specify in them.

## fast5_qc.sh

This just uses poretools to look at the raw fast5 files that come from metrichor

## fastq_qc.sh

So far all this does is create a fastq file with poretools, then use fastqc to look at it... more to come.