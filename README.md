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

This just uses poretools to look at the raw fast5 files that come from metrichor. Please note that if you are running this remotely on a server, then you will need to edit a couple of things in poretools first. Specifically:

1. Go to ```/poretools/poretools/hist.py```, edit this line:

```
#matplotlib.use('Agg') # Must be called before any other matplotlib calls
```
to remove the comment, so it looks like this:
```
matplotlib.use('Agg') # Must be called before any other matplotlib calls
```

2. Do the same for ```/poretools/poretools/yield_plot.py```

3. Add the following two lines to the top of ```/poretools/poretools/qual_v_pos.py```, right after the ```import pandas``` line.

```
import matplotlib
matplotlib.use('Agg') # Must be called before any other matplotlib calls
```                                                                            

4. Ditto 3 for ```/poretools/poretools/occupancy.py``` (make sure it goes right after the pandas line)

5. cd to /poretools and run ```python setup.py install```

This sets things up so that you can write PDFs on your server while logged in remotely.

## fastq_qc.sh

So far all this does is create a fastq file with poretools, then use fastqc to look at it... more to come.