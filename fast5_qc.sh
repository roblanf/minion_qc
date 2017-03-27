# A script by Rob Lanfear and Miriam Schalamun to do basic QC on MinIon data

# set up directories for input and output
inputf="/disks/dacelo/data/raw_data/tree_EG1/1003_Miriam_data_downloaded/pass/"
outputf="/disks/dacelo/data/QC/1003_Miriam_data_downloaded/"

# calculate basic stats and raw data and save them
poretools stats $inputf > $outputf/stats.txt
poretools nucdist $inputf > $outputf/nucdist.txt
poretools qualdist $inputf > $outputf/qualdist.txt

# plot histograms and save the pdfs
poretools hist --saveas $outputf/readlength_hist.pdf --theme-bw $inputf

# yield plots
poretools yield_plot --plot-type basepairs --saveas $outputf/yield_basepairs.pdf --theme-bw $inputf
poretools yield_plot --plot-type reads --saveas $outputf/yield_reads.pdf --theme-bw $inputf

# occupancy plot
poretools occupancy --saveas $outputf/occupancy.pdf $inputf


