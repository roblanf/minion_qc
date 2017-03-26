# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/disks/dacelo/data/raw_data/tree_EG1/1003_Miriam_data_downloaded/pass/"
outputbase="/disks/dacelo/data/QC/1003_Miriam_data_downloaded/"
ref="/disks/dacelo/data/active_refs/Egra.fa.gz" # reference file as a fasta
gff="/disks/dacelo/data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=20 # number of threads to use

# set up directories for output
outputrawqc=$outputbase"rawqc/"
mkdir $outputbase
mkdir $outputrawqc

# make the fastq file from the fast5 files using poretools
# uncoment both the nanopolish and poretools lines to compare them
# just leave one of them there to use that one - they do the same
echo "Creating fastq file with poretools"
fastq_file=$outputbase'reads.fastq.gz'
#time poretools fastq $inputf | gzip > $fastq_file
time nanopolish extract --fastq $inputf | gzip > $fastq_file

# run fastqc on the raw fastq data
echo "Running fastqc"
fastqc $fastq_file -o $outputrawqc -t $threads

# map with BWA mem, pipe, sort
# TODO: make sure BWA prints proper secondary mapping flags 
outBWAMEM=$outputbase"BWAMEM/"
mkdir $outBWAMEM
cd $outBWAMEM
bwa index $ref
echo "Mapping with BWAMEM"
date
time bwa mem -x ont2d -t $threads $ref $fastq_file > out.sam
echo "Done Mapping with BWAMEM"
date
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_gff/" -gff $gff -nt $threads -c

# graphmap...
outgm=$outputbase"gm/"
mkdir $outgm
cd $outgm
echo "Mapping with graphmap"
date
time graphmap align -t $threads -r $ref -d $fastq_file -o out.sam --extcigar
echo "Done Mapping with graphmap"
date
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outgm"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outgm"qualimap_gff/" -gff $gff -nt $threads -c


# ngmlr
outngmlr=$outputbase"ngmlr/"
mkdir $outngmlr
cd $outngmlr
echo "Mapping with ngmlr"
date
time ngmlr -t $threads -r $ref -q $fastq_file -o out.sam -x ont
echo "Done Mapping with ngmlr"
date
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outngmlr"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outngmlr"qualimap_gff/" -gff $gff -nt $threads -c

# TODO - consider turning off secondary read mapping in BWA and NGMLR to make them more comparable.

