# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/disks/dacelo/data/test_data_minion/"
outputbase="/disks/dacelo/data/QC/test_data_minion/"
ref="/disks/dacelo/data/raw_data/active_refs/Emel.fa.gz" # reference file as a fasta
gff="/disks/dacelo/data/raw_data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=10 # number of threads to use


# set up directories for output
outputrawqc=$outputbase"rawqc/"
mkdir $outputbase
mkdir $outputrawqc


# make the fastq file from the fast5 files using poretools
fname=$(basename "$input_dir")
poretools fastq $inputf | gzip > $outputbase'/'$fname'.fastq.gz'

# TODO fix up the hard-coded filename
# run fastqc on all the raw and trimmed data files
fastqc $outputbase'minion_test.fastq.gz' -o $outputrawqc -t $threads

