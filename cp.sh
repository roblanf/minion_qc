# shell script to try and assemble a chloroplast genome from raw minion reads


# A few things to set before you go
inputf="/disks/dacelo/data/raw_data/tree_EG1/1003_Miriam_data_downloaded/pass/"
outputf="/disks/dacelo/data/mincan/"
ref="/disks/dacelo/data/active_refs/Egra.fa.gz" # reference file as a fasta
threads=10 # number of threads to use

# make the fastq file from the fast5 files using poretools
# uncoment both the nanopolish and poretools lines to compare them
# just leave one of them there to use that one - they do the same
echo "Creating fastq file with poretools"
fastqf=$outputf'reads.fastq.gz'
time poretools fastq $inputf | gzip > $fastqf

# TODO: trim adaptors

# map with BWA mem, sort
# the -M flag marks shorter split hits as secondary, so we can find them later
outBWAMEM=$outputf"BWAMEM/"
mkdir $outBWAMEM
cd $outBWAMEM
bwa index $ref
echo "Mapping with BWAMEM"
date
time bwa mem -x ont2d -M -t $threads $ref $fastq_file > out.sam
echo "Done Mapping with BWAMEM"
date
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
rm out.sam

# get just the things that map to the chloroplast
samtools view -b out.bam NC_014570.1 > cp.bam
bamToFastq -i out.bam -fq cp.fq


# now assemble those reads with CANU
# we know the genome size is ~163672
# since we have ~200x coverage we could also try corOutCoverage parameter
# E.g. to use all the reads we could set corOutCoverage=1000
# and to use the longest set of reads to get ~100x coverage we would use corOutCoverage=100
# the manual suggests that when increasing the coverage we can 
# since we have R9 1D reads, we should set errorRate=0.025 (according to the manual)
# but if we increase coverage >60x, the manual suggests decreasing this to errorRate=0.013 which should use just the better reads.
canu -p cp -d cp genomeSize=164k -nanopore-raw cp.fq maxThreads=10

# trim and circularise with circlator as per: https://github.com/marbl/canu/issues/87

