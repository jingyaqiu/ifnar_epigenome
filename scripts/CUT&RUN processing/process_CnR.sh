#!/bin/bash
set -e

source /home/jingyaq/Minn/data/R499_CnR/scripts/activate_process_CnR.sh

# cutadapt version 2.9
# bowtie2 version 2.3.4.1
# picard.jar version 2.23.3
# samtools version 1.11
# bedtools version 2.29.2
# R 3.5.1

DIRPATH=$(pwd)

NAME=$1

echo ${NAME}

###################################
### Adapter trimming at 3' ends ###
###################################

mkdir -p ${DIRPATH}/fastq/trimmed/logs

echo "Trim adapters"

# NEBNext adapters are the same as TruSeq kits
# [page 30] https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/experiment-design/illumina-adapter-sequences-1000000002694-14.pdf

adapter1="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
adapter2="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

fastq1=${DIRPATH}/fastq/${NAME}_R1_001.fastq.gz
fastq2=${DIRPATH}/fastq/${NAME}_R2_001.fastq.gz
fastq1_trimmed=${DIRPATH}/fastq/trimmed/${NAME}_trimmed_R1.fastq.gz
fastq2_trimmed=${DIRPATH}/fastq/trimmed/${NAME}_trimmed_R2.fastq.gz
trim_log=${DIRPATH}/fastq/trimmed/logs/${NAME}.log

~/.local/bin/cutadapt -m 5 -e 0.1 -q 10 --cores=16 -a ${adapter1} -A ${adapter2} -o ${fastq1_trimmed} -p ${fastq2_trimmed} ${fastq1} ${fastq2} > ${trim_log}

### Adapter trimming statistics ###

log_summary=${DIRPATH}/fastq/trimmed/${NAME}.summary

echo ${NAME} >> ${log_summary}
grep ${trim_log} -e "Total read pairs processed:" >> ${log_summary}
grep ${trim_log} -e "Read 1 with adapter:" >> ${log_summary}
grep ${trim_log} -e "Read 2 with adapter:" >> ${log_summary}
grep ${trim_log} -e "Pairs written (passing filters):" >> ${log_summary}
cat ${log_summary} >> ${DIRPATH}/fastq/trimmed/logs/combined_log.summary
rm ${log_summary}

################################
### Align reads with Bowtie2 ###
################################

echo "Align reads"

mkdir -p ${DIRPATH}/bowtie2/logs

ref="/home/jingyaq/Minn/resources/bowtie2/mm10/mm10"

sam=${DIRPATH}/bowtie2/${NAME}.sam
bam=${DIRPATH}/bowtie2/${NAME}.sorted.bam
bowtie2_log=${DIRPATH}/bowtie2/logs/${NAME}.bowtie2.log

# For >25 bp reads
bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 16 \
	-x ${ref} \
	-1 ${fastq1_trimmed} -2 ${fastq2_trimmed} \
	-S ${sam} 2> ${bowtie2_log}

# For 25 bp reads
# bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p 16 -x ${ref} -1 ${dirpath}/fastq/trimmed/${name}_trimmed_R1.fastq.gz -2 ${dirpath}/fastq/trimmed/${name}_trimmed_R2.fastq.gz -S ${dirpath}/bowtie2/${name}.sam 2> ${dirpath}/bowtie2/logs/${name}.bowtie2.log

# Align untrimmed reads
# bowtie2 --local --very-sensitive-local --no-unal --no-mixed --no-discordant --phred33 -I 10 -X 700 -x ${ref} -1 ${dirpath}/fastq/${name}_R1_001.fastq.gz -2 ${dirpath}/fastq/${name}_R2_001.fastq.gz -S ${dirpath}/bowtie2/${name}.sam 2> ${dirpath}/bowtie2/logs/${name}.bowtie2.log

/home/jingyaq/samtools-1.11/samtools view -u ${sam} | /home/jingyaq/samtools-1.11/samtools sort - > ${bam}
/home/jingyaq/samtools-1.11/samtools index ${bam}

# Remove intermediate files

rm ${sam}

########################
### Filter bam files ###
########################

# Retain only properly paired, unique reads
# Remove:
# 1. unmapped reads (-F 1804)
# 2. reads that are not properly paired (-f 2)
# 3. Low MAPQ (-q 5) # ???

echo "Filter bam"

bam_filtered=${DIRPATH}/bowtie2/${NAME}.sorted.filtered.bam

# /home/jingyaq/samtools-1.11/samtools view -h -F 1804 -f 2 -q 5 -b ${bam} > ${bam_filtered}
/home/jingyaq/samtools-1.11/samtools view -h -F 1804 -f 2 -b ${bam} > ${bam_filtered}
/home/jingyaq/samtools-1.11/samtools index ${bam_filtered}

###########################
### MARK PCR duplicates ###
###########################

## For CUT&TAG, not recommended to remove duplicates because Tn5 bias often results in cutting at the exact same sites - easy to confuse true fragments as duplicates. Removing duplicates only recommended when there is very small amounts of material, or when PCR duplication is suspected. 
## https://yezhengstat.github.io/CUTTag_tutorial/#iv_alignment_results_filtering_and_file_format_conversion
## CUT&RUN protocol: "Use Picard 'MarkDuplicates' to mark presumed PCR duplicates for removal from low-cell number data"
## https://www.nature.com/articles/nprot.2018.015#Sec26

echo "Mark duplicates"

bam_dupMarked=${DIRPATH}/bowtie2/${NAME}.sorted.filtered.dupMarked.bam
bam_dupMarked_metrics=${DIRPATH}/bowtie2/logs/${NAME}.dupMarked.log

java -jar /home/jingyaq/picard.jar MarkDuplicates -INPUT ${bam_filtered} -OUTPUT ${bam_dupMarked} -METRICS_FILE ${bam_dupMarked_metrics} -REMOVE_DUPLICATES FALSE
/home/jingyaq/samtools-1.11/samtools index $bam_dupMarked

############################
### Filter out blacklist ###
############################

echo "Filter blacklist reads"

BLACKLIST=/home/jingyaq/Minn/references/GRCm38/mm10.blacklist.bed
# BLACKLIST=/home/jingyaq/Minn/references/GRCm38/mm10.blacklist_modified.bed

bam_blfilt=${DIRPATH}/bowtie2/${NAME}.sorted.filtered.dupMarked.blfilt.bam

/home/jingyaq/bedtools2/bin/intersectBed -a ${bam_dupMarked} -b ${BLACKLIST} -v > ${bam_blfilt}
/home/jingyaq/samtools-1.11/samtools index ${bam_blfilt}

rm ${bam_filtered}* ${bam_dupMarked}*

#####################
### bam to bigwig ###
#####################

echo "Generate bigwig files"

mkdir -p ${DIRPATH}/bw

bw_blfilt=${DIRPATH}/bw/${NAME}.sorted.filtered.dupMarked.blfilt.bw

# Multi-mapped reads excluded so effective genome size is lower
# https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
~/.local/bin/bamCoverage -b ${bam_blfilt} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2308125349 --extendReads --ignoreForNormalization chrM -o ${bw_blfilt}

##################################
### Fragment size distribution ###
##################################

echo "Calculate fragment lengths"

# 9th column of bam file is fragment length
samtools view ${bam_blfilt} | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${DIRPATH}/bowtie2/logs/${NAME}_fragment_lengths.txt

#######################
### bam to bedgraph ###
#######################

echo "Generate bedgraph files for SEACR peak calling"

# Prepare input bedgraph files for SEACR
# https://github.com/FredHutch/SEACR

mkdir -p ${DIRPATH}/bowtie2/bedgraph

bed=${DIRPATH}/bowtie2/bedgraph/${NAME}.bed
bed_clean=${DIRPATH}/bowtie2/bedgraph/${NAME}.clean.bed
bed_fragments=${DIRPATH}/bowtie2/bedgraph/${NAME}.fragments.bed
bedgraph=${DIRPATH}/bowtie2/bedgraph/${NAME}.fragments.bedgraph

# Convert bam to bed file
/home/jingyaq/bedtools2/bin/bamToBed -i ${bam_blfilt} -bedpe > ${bed}

# Keep read pairs on the same chromosome and fragment length <1000 bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bed} > ${bed_clean}

# Only extract fragment-related columns
cut -f 1,2,6 ${bed_clean} | sort -k1,1 -k2,2n -k3,3n > ${bed_fragments}

/home/jingyaq/bedtools2/bin/genomeCoverageBed -bg -i ${bed_fragments} -g "/home/jingyaq/Minn/STAR_genome_indices/STAR_2.7.1a/GRCm38/chrNameLength.txt" > ${bedgraph}

rm ${bed} ${bed_fragments}
