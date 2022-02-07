#!/bin/bash
set -e

source /home/jingyaq/Minn/data/ifnar_epigenome/scripts/processing/activate_process_RNA.sh

# cutadapt version 2.9
# salmon version 1.2.0
# picard.jar version 2.23.3
# samtools version 1.3.1
# deeptools version 3.5.0

dirpath="/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/RNA"

sample=$1

echo $sample

##############################
### Trim adapter sequences ###
##############################

mkdir -p /home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/trimmed/logs

echo "Trim adapters"

r1_adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
r2_adapter="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

fastq1=/home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/${sample}.R1.fastq.gz
fastq2=/home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/${sample}.R2.fastq.gz
fastq1_trimmed=/home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/trimmed/${sample}_trimmed_R1.fastq.gz
fastq2_trimmed=/home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/trimmed/${sample}_trimmed_R2.fastq.gz

~/.local/bin/cutadapt -m 30 -e 0.2 -q 10 -O 3 --cores=16 -a ${r1_adapter} -A ${r2_adapter} -o ${fastq1_trimmed} -p ${fastq2_trimmed} ${fastq1} ${fastq2} > /home/jingyaq/Minn/data/B16_R499_Stat1KO_RNA/data/fastq/trimmed/logs/${sample}.log

#############################
### Salmon quantification ###
#############################

echo "Gene quantification with salmon"

mkdir -p ${dirpath}/salmon

source activate salmon

index_dir="/home/jingyaq/Minn/resources/salmon/GRCm38_vM24_decoyaware"

salmon_out=${dirpath}/salmon/${sample}

salmon quant -i ${index_dir} -l A --validateMappings --rangeFactorizationBins 4 --numBootstraps 3 --seqBias --gcBias -p 16 -1 ${fastq1_trimmed} -2 ${fastq2_trimmed} -o ${salmon_out}

####################
### STAR mapping ###
####################

echo "Align reads with STAR"

mkdir -p ${dirpath}/star/logs

genome_dir="/home/jingyaq/Minn/resources/STAR/GRCm38_vM24"
star_out=${dirpath}/star/${sample}_

STAR --runThreadN 32 --genomeDir ${genome_dir} --readFilesCommand zcat --readFilesIn ${fastq1_trimmed} ${fastq2_trimmed} --outFileNamePrefix ${star_out} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMtype BAM SortedByCoordinate

bam=${dirpath}/star/${sample}_Aligned.sortedByCoord.out.bam
samtools index ${bam}
mv ${dirpath}/star/${sample}_Log* ${dirpath}/star/logs

### Gather mapping summary statistics ###

# Sequencing depth, alignment rate, number of mappable fragments, duplication rate, unique library size, fragment size distribution #

########################
### Filter bam files ###
########################

# Retain only properly paired, unique reads
# Remove:
# 1. unmapped reads (-F 1804)
# 2. reads that are not properly paired (-f 2)
# 3. MAPQ (-q 5)

echo "Filter bam"

mkdir -p ${dirpath}/star/filtered

bam_filtered=${dirpath}/star/filtered/${sample}_filtered.bam

samtools view -h -F 1804 -f 2 -q 5 -b ${bam} > ${bam_filtered}
samtools index ${bam_filtered}

# #############################
# ### Remove PCR duplicates ###
# #############################

# echo "Remove duplicates"

# bam_deduped=${dirpath}/star/filtered/${sample}_filtered_deduped.bam
# deduped_metrics=${dirpath}/star/logs/${sample}.deduped.log

# java -jar /home/jingyaq/picard.jar MarkDuplicates -INPUT ${bam_filtered} -OUTPUT ${bam_deduped} -METRICS_FILE ${deduped_metrics} -REMOVE_DUPLICATES FALSE

# samtools index ${bam_deduped}
# rm ${bam_filtered}*

#####################
### bam to bigwig ###
#####################

echo "Generate bigwig files"

mkdir -p ${dirpath}/bw

bw=${dirpath}/bw/${sample}_filtered_deduped.bw

# CHECK "chrM" notation is correct!
~/.local/bin/bamCoverage -b ${bam_deduped} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2407883318 --extendReads --ignoreForNormalization chrM -o ${bw}


