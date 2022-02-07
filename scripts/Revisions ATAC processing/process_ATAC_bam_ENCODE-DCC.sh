#!/bin/bash
set -e

source /home/jingyaq/Minn/data/ifnar_epigenome/scripts/processing/activate_process_ATAC_bam.sh

# samtools version 1.11 (old 1.3.1)
# picard.jar version 2.23.3
# deeptools version 3.1.5 (old 3.1.3)
# bedtools version 2.29.2
# R 3.5.1

DIRPATH=$(pwd)

SAMPLE=$1

echo $SAMPLE

###################################################################
### Further processing of bam file produced by ENCODE pipeline ####
###################################################################

# bam=${DIRPATH}/bam/${SAMPLE}_*_R1_001.trim.srt.nodup.no_chrM_MT.bam
bam=${DIRPATH}/bam/${SAMPLE}_R1_001.trim.srt.nodup.no_chrM_MT.bam
bam_tn5=${DIRPATH}/bam/${SAMPLE}.trim.srt.nodup.no_chrM_MT_tn5shift.bam
bam_tn5_sorted=${DIRPATH}/bam/${SAMPLE}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam
bw=${DIRPATH}/bw/${SAMPLE}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bw

###################################################################

### Already done by ENCODE pipeline ###

# mkdir -p ${DIRPATH}/bam/filtered/logs
# mkdir -p ${DIRPATH}/bw/filtered

# # Retain only properly paired reads with MAPQ > 30
# samtools view -h -F 1804 -f 2 -q 30 -b ${bam} > ${bam_filt}
# samtools index ${bam_filt}

# # Remove duplicate reads
# java -jar /home/jingyaq/picard.jar MarkDuplicates I={bam_filt} O=${bam_deduped} M=${deduped_metrics} REMOVE_DUPLICATES=true CREATE_INDEX=true

# echo "Remove mitochondrial reads"
# samtools idxstats ${bam} | cut -f 1 | grep -v "chrM" | xargs samtools view -b ${bam} | samtools sort -@ 12 -O bam -o ${bam_noChrM}
# samtools index ${bam_noChrM}

###################################################################

/home/jingyaq/samtools-1.11/samtools index ${bam}

echo "tn5 shift"
~/.local/bin/alignmentSieve -b ${bam} --ATACshift -p 16 -o ${bam_tn5}
/home/jingyaq/samtools-1.11/samtools sort -o ${bam_tn5_sorted} -@ 12 ${bam_tn5}
/home/jingyaq/samtools-1.11/samtools index ${bam_tn5_sorted}

# Remove intermediate bam files
rm ${bam_tn5}*

echo "Generate bigwigs"
~/.local/bin/bamCoverage \
	-b ${bam_tn5_sorted} \
	-of bigwig \
	--binSize 1 \
	-p 16 \
	--normalizeUsing RPGC \
	--effectiveGenomeSize 2652783500 \
	--extendReads \
	--centerReads \
	-o ${bw}


