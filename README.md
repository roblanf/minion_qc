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

## Output details
This script is designed as a fairly quick-and-dirty way to examine the quality and quantity of your MinION output. More details of the rationale behind the plots is here: robertlanfear.com/blog/files/nanopore_performance.html, below I provide a brief overview and example of the files produced (which are all the `/example_output` folder of this repository).

For convenience, a lot of the plots and summary statistics are repeated on both all of the data, and on just those reads with a mean Q score >10. We found this to be a convenient cutoff to get a rough idea of how much *good* data is in a run, rather than just how much data there is overall.

### summary.txt

This txt file has some simple summary statistics for the data. Entries should be self exlpanatory for the most part. The `reads` and `bases` parts are the end are just the number of reads and number of bases from reads greater than the cutoffs in the table. E.g. in the example below, for the data with reads >Q10, we have 205 reads >200Kb that sum to a total of 22624337 bases. I like this stat, because I can quickly know that I have ~4% coverage of my genome (which is ~500Kb) from these very long reads. So if I get ~25 of these flowcells I can expect ~1x coverage of ultra-long reads.

```
Summary stats from input file sequencing_summary.txt 


all.reads.summary 
 	  	 
total.bases 	 3843373607 
 	  	 
N50.length 	 34354 
 	  	 
mean.length 	 16150.05 
 	  	 
median.length 	 7176 
 	  	 
max.length 	 826249 
 	  	 
mean.q 	 9.502917 
 	  	 
median.q 	 10.131 
 	 >20kb 	 >50kb 	 >100kb 	 >200kb 	 >500kb 	 >1m 	 
reads 	 83614 14651 253 10 1 0 
 	 >20kb 	 >50kb 	 >100kb 	 >200kb 	 >500kb 	 >1m 	 
bases 	 3113987988 933314126 30331350 3045299 826249 0 



q10.reads.summary 
 	  	 
total.bases 	 3276103734 
 	  	 
N50.length 	 36010 
 	  	 
mean.length 	 27284.04 
 	  	 
median.length 	 24669 
 	  	 
max.length 	 210971 
 	  	 
mean.q 	 12.71787 
 	  	 
median.q 	 12.858 
 	 >20kb 	 >50kb 	 >100kb 	 >200kb 	 >500kb 	 >1m 	 
reads 	 75226 13456 205 1 0 0 
 	 >20kb 	 >50kb 	 >100kb 	 >200kb 	 >500kb 	 >1m 	 
bases 	 2811649966 854005042 22624337 210971 0 0 
```

### length_histogram.png
Read length, on a log10 scale, on the X axis, and counts on the Y axis.
![length_histogram](example_output/length_histogram.png)

### q_histogram.png
Mean Q score for a read on the X axis, and counts on the Y axis. 
![q_histogram](example_output/q_histogram.png)

### epb_histogram.png
Events per base (i.e. numbe of events for each read, divided by the number of bases called for that read) on the X axis, and counts on the Y axis. We have found this measure to be potentially useful in finding dodgy reads, see robertlanfear.com/blog/files/nanopore_performance.html for more.
![epb_histogram](example_output/epb_histogram.png)

### length_vs_q.png
Read length (log10 scale) on the X axis, mean Q score on the Y axis. Points are coloured by the events per base, scaled such that all events per base >10 are recorded as a 10. The latter makes it easier to distinguish 'good' and 'bad' reads. 
![length_vs_q](example_output/length_vs_q.png)

### yield_summary.png
Minimum read length on the X axis, and the yield of bases with reads at least that long on the Y axis. This is just like the 'reads' table in the `summary.txt` output, but done across all read lengths up to 100KB. I cut off at 100KB because you (probably) don't have most of your data at those lenghts. Good on you if you do though.
![yield_summary](example_output/yield_summary.png)

### channel_summary.png
Histograms of total bases, total reads, mean read length, and median read length that show the variance across the 512 available channels. Repeated for all data and reads with Q>10.
![channel_summary](example_output/channel_summary.png)

### flowcell_channels_epb.png
This one's busy, but hopefully useful. The 512 channels are laid out as on the R9.5 flowcell. Then each sub-panel of the plot simply plots out the time of the run in hours on the X axis, and the events per base (log scale, cut off at 10 events per base) on the Y axis. This gives a little insight into exactly what was going on in each of your channels over the course of the run. Blow it up big!
![flowcell_channels_epb](example_output/flowcell_channels_epb.png)
