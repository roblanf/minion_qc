# A script by Rob Lanfear and Miriam Schalamun to do basic QC on MinIon data


# set up directories for input and output
input_dir="/Users/roblanfear/Desktop/minion_test"
output_dir="/Users/roblanfear/Desktop/minion_qc"
mkdir $output_dir

# calculate basic stats and raw data and save them
poretools stats $input_dir > $output_dir/stats.txt
poretools nucdist $input_dir > $output_dir/nucdist.txt
poretools qualdist $input_dir > $output_dir/qualdist.txt

# plot histograms and save the pdfs
poretools hist --saveas $output_dir/readlength_hist.pdf --theme-bw $input_dir
poretools qualpos --saveas $output_dir/quality_boxplot.pdf $input_dir

# yield plots
poretools yield_plot --plot-type basepairs --saveas $output_dir/yield_basepairs.pdf --theme-bw $input_dir
poretools yield_plot --plot-type reads --saveas $output_dir/yield_reads.pdf --theme-bw $input_dir

# occupancy plot
poretools occupancy --saveas $output_dir/occupancy.pdf $input_dir