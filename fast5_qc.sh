# A script by Rob Lanfear and Miriam Schalamun to do basic QC on MinIon data

# set up directories for input and output
inputf="/disks/dacelo/data/test_data_minion"
outputf="/disks/dacelo/data/QC/test_data_minion"

# calculate basic stats and raw data and save them
poretools stats $inputf > $outputf/stats.txt
poretools nucdist $inputf > $outputf/nucdist.txt
poretools qualdist $inputf > $outputf/qualdist.txt

# plot histograms and save the pdfs
poretools hist --saveas $outputf/readlength_hist.pdf --theme-bw $inputf
poretools qualpos --saveas $outputf/quality_boxplot.pdf $inputf

# yield plots
poretools yield_plot --plot-type basepairs --saveas $outputf/yield_basepairs.pdf --theme-bw $inputf
poretools yield_plot --plot-type reads --saveas $outputf/yield_reads.pdf --theme-bw $inputf

# occupancy plot
poretools occupancy --saveas $outputf/occupancy.pdf $inputf


