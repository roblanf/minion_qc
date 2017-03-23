# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/disks/dacelo/data/test_data_minion/"
outputbase="/disks/dacelo/data/QC/test_data_minion/"
ref="/disks/dacelo/data/active_refs/Emel.fa.gz" # reference file as a fasta
gff="/disks/dacelo/data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=10 # number of threads to use

# set up directories for output
outputrawqc=$outputbase"rawqc/"
mkdir $outputbase
mkdir $outputrawqc

# make the fastq file from the fast5 files using poretools
fastq_file=$outputbase'reads.fastq.gz'
poretools fastq $inputf | gzip > $fastq_file

# run fastqc on all the raw and trimmed data files
fastqc $fastq_file -o $outputrawqc -t $threads

# map with BWA mem, pipe, sort
outBWAMEM=$outputbase"BWAMEM/"
mkdir $outBWAMEM
cd $outBWAMEM
bwa index $ref
bwa mem -x ont2d -t $threads $ref $fastq_file > out.sam
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outBWAMEM"qualimap_gff/" -gff $gff -nt $threads -c

# graphmap...
outgm=$outputbase"gm/"
mkdir $outgm
cd $outgm
graphmap align -t $threads -r $ref -d $fastq_file -o out.sam --extcigar
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outgm"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outgm"qualimap_gff/" -gff $gff -nt $threads -c


# ngmlr
outngmlr=$outputbase"ngmlr/"
mkdir $outngmlr
cd $outngmlr
ngmlr -t $threads -r $ref -q $fastq_file -o out.sam -x ont
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outngmlr"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outngmlr"qualimap_gff/" -gff $gff -nt $threads -c

# nanoBLASTer doesn't seem to work.
outnano=$outputbase"nano/"
mkdir $outnano
cd $outnano
nanoblaster -C10 -r $ref -i $fastq_file -o out
samtools view -bS -@ $threads out.sam > out.bam
samtools sort -@ $threads out.bam -o out.bam
samtools index out.bam
qualimap bamqc -bam out.bam -outdir $outnnano"qualimap_all/" -nt $threads -c
qualimap bamqc -bam out.bam -outdir $outnnano"qualimap_gff/" -gff $gff -nt $threads -c

# let's see what we get 
cd outputbase
multiqc . --force