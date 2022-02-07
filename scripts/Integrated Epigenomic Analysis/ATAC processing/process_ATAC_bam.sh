#!/bin/bash
set -e

source /home/jingyaq/Minn/data/ifnar_epigenome/scripts/processing/activate_process_ATAC_bam.sh

# samtools version 1.3.1
# picard.jar version 2.23.3
# deeptools version 3.1.3
# bedtools version 2.29.2
# R 3.5.1

dirpath="/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC"

sample=$1

echo $sample

###################################################################
### Further processing of bam file produced by ENCODE pipeline ####
###################################################################

mkdir -p ${dirpath}/bam/filtered/logs
mkdir -p ${dirpath}/bw/filtered

bam=${dirpath}/bam/FGC1587_${sample}.R1.trim_pooled.PE2SE.nodup.bam
bam_noChrM=${dirpath}/bam/filtered/${sample}.R1.trim_pooled.PE2SE.nodup_noChrM.bam
bam_tn5=${dirpath}/bam/filtered/${sample}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift.bam
bam_tn5_sorted=${dirpath}/bam/filtered/${sample}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam
bw=${dirpath}/bw/filtered/${sample}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bw

###################################################################

### Already done by ENCODE pipeline ###

# # Retain only properly paired reads with MAPQ > 30
# samtools view -h -F 1804 -f 2 -q 30 -b ${bam} > ${bam_filt}
# samtools index ${bam_filt}

# # Remove duplicate reads
# java -jar /home/jingyaq/picard.jar MarkDuplicates I={bam_filt} O=${bam_deduped} M=${deduped_metrics} REMOVE_DUPLICATES=true CREATE_INDEX=true

###################################################################

echo "Remove mitochondrial reads"
samtools idxstats ${bam} | cut -f 1 | grep -v "chrM" | xargs samtools view -b ${bam} | samtools sort -@ 12 -O bam -o ${bam_noChrM}
samtools index ${bam_noChrM}

echo "tn5 shift"
alignmentSieve -b ${bam_noChrM} --ATACshift -p 16 -o ${bam_tn5}
samtools sort -o ${bam_tn5_sorted} -@ 12 ${bam_tn5}
samtools index ${bam_tn5_sorted}

# Remove intermediate bam files
rm ${bam_noChrM}* ${bam_tn5}*

echo "Generate bigwigs"
~/.local/bin/bamCoverage -b ${bam_tn5_sorted} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${bw}

