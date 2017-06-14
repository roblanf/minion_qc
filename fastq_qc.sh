# Basic quality control for mapping MinION data to a reference

# Rob Lanfear and Miriam Schalamun

# A few things to set before you go
inputf="/data/nanopore/RB7_A2/20170609_0424_Epauciflora_A2/"
outputbase="/data/nanopore/RB7_A2/"
ref="/data/active_refs/Epau.fa.gz" # reference file as a fasta
gff="/data/active_refs/Egrandis_genes_chr1_to_chr11.gff3"
threads=55 # number of threads to use
mem_size='50G' # memory size for Qualimap
flowcellID="FLO-MIN107"
kidID="SQK-LSK108"

mkdir $outputbase

# basecall with albacore
read_fast5_basecaller.py -i $inputf -t $threads -s $outputbase -f FLO-MIN107 -k SQK-RAD002 -r -o fastq

# cat together the fastq's to one big fastq
fastq_file=$outputbase"reads.fastq"
cat $outputbase/workspace/*.fastq > $fastq_file

# map with ngmlr
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
rm out.sam
qualimap bamqc -bam out.bam -outdir $outputbase"qualimap_all/" -nt $threads -c --java-mem-size=$mem_size
qualimap bamqc -bam out.bam -outdir $outputbase"qualimap_gff/" -gff $gff -nt $threads -c --java-mem-size=$mem_size

# run nanoplot
echo "Running nanoplot"
outnano=$outputbase"nanoplot/"
mkdir $outnano
NanoPlot --fastq_rich $fastq_file --outdir $outnano --threads $threads --loglength
NanoPlot --bam out.bam --outdir $outnano --threads $threads --loglength --prefix bam

# stats on reads > various length (thanks to @gringer here: https://bioinformatics.stackexchange.com/questions/678/get-the-mapping-statistics-of-a-single-read-from-a-bam-file/696#696)
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>1000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_1k.txt
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>2000)){print $0}}' |  samtools stats   | grep '^SN' | cut -f 2- > stats_2k.txt
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>10000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_10k.txt
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>20000)){print $0}}' | samtools stats   | grep '^SN' | cut -f 2- > stats_20k.txt
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>100000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_100k.txt
samtools view -h ngmlr/out.bam |     awk -F'\t' '{if((/^@/) || (length($10)>200000)){print $0}}' | samtools stats  | grep '^SN' | cut -f 2- > stats_200k.txt 

