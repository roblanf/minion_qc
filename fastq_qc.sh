# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/Users/roblanfear/Desktop/minion_test/"
outputbase="/Users/roblanfear/Desktop/minion_qc/"
ref="/disks/dacelo/data/raw_data/active_refs/Emel.fa.gz" # reference file as a fasta
gff="/disks/dacelo/data/raw_data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=10 # number of threads to use


outputrawqc=$outputbase"rawqc/"
echo $outputrawqc

# set up dirs
mkdir $outputbase
mkdir $outputbase"ngm"
mkdir $outputrawqc


# make the fastq file from the fast5 files using poretools
fname=$(basename "$input_dir")
poretools fastq $inputf | gzip > $outputbase'/'$fname'.fastq.gz'


# run fastqc on all the raw and trimmed data files
echo "Running fastqc"
date
find $inputf -name '*.fastq.gz' | xargs fastqc  -o $outputrawqc -t $threads
echo "Done running fastqc"
date

