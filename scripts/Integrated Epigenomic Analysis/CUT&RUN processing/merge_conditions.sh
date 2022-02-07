#!/bin/bash
set -e

# samtools version 1.11
# deeptools version 3.1.3

DIRPATH=$(pwd)

COND=$1
HISTONE=$2
NAME=${COND}_${HISTONE}

echo $NAME

mkdir -p ${DIRPATH}/bowtie2/merged
mkdir -p ${DIRPATH}/bw/merged

BAM=(${DIRPATH}/bowtie2/${COND}*${HISTONE}.sorted.filtered.dupMarked.blfilt.bam)

# Remove poor quality samples
# remove_samples=(R499_SKO_H3K27Ac_rep1 R499_SKO_H3K4me3_rep2 R499_SKO_H3K4me1_rep2)
# for sample in ${remove_samples[@]}
# do
# 	file=${DIRPATH}/bowtie2/${sample}.sorted.filtered.dupMarked.blfilt.bam
# 	BAM=("${BAM[@]/$file}")
# done
# echo ${BAM[@]}

#######################
### Merge bam files ###
#######################

BAM_MERGED=${DIRPATH}/bowtie2/merged/${NAME}.sorted.filtered.dupMarked.blfilt_merged.bam
BW_MERGED=${DIRPATH}/bw/merged/${NAME}.sorted.filtered.dupMarked.blfilt_merged.bw

# Merge bam files
/home/jingyaq/samtools-1.11/samtools merge -@ 16 ${BAM_MERGED} ${BAM[@]}
/home/jingyaq/samtools-1.11/samtools index ${BAM_MERGED}

# Generate merged bigwigs
~/.local/bin/bamCoverage -b ${BAM_MERGED} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${BW_MERGED}

###############################
### Convert bam to bedgraph ###
###############################

mkdir -p bowtie2/merged/bedgraph

bed=${DIRPATH}/bowtie2/merged/bedgraph/${NAME}.bed
bed_clean=${DIRPATH}/bowtie2/merged/bedgraph/${NAME}.clean.bed
bed_fragments=${DIRPATH}/bowtie2/merged/bedgraph/${NAME}.fragments.bed
bedgraph=${DIRPATH}/bowtie2/merged/bedgraph/${NAME}.fragments.bedgraph

# Fragment size distribution
samtools view ${BAM_MERGED} | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${DIRPATH}/bowtie2/logs/${NAME}_fragment_lengths.txt

# Convert bam to bed file
/home/jingyaq/bedtools2/bin/bamToBed -i ${BAM_MERGED} -bedpe > ${bed}

# Keep read pairs on the same chromosome and fragment length <1000 bp
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bed} > ${bed_clean}

# Only extract fragment-related columns
cut -f 1,2,6 ${bed_clean} | sort -k1,1 -k2,2n -k3,3n > ${bed_fragments}

/home/jingyaq/bedtools2/bin/genomeCoverageBed -bg -i ${bed_fragments} -g "/home/jingyaq/Minn/STAR_genome_indices/STAR_2.7.1a/GRCm38/chrNameLength.txt" > ${bedgraph}

rm ${bed} ${bed_fragments}

#########################################
### Peak calling on merged conditions ###
#########################################

source /home/jingyaq/Minn/data/R499_CnR/scripts/activate_peak_calling_CnR.sh

mkdir -p macs2/merged/logs seacr/merged

CONTROL_BAM=${DIRPATH}/bowtie2/merged/${COND}_IgG.sorted.filtered.dupMarked.blfilt_merged.bam
TREATMENT_BAM=${DIRPATH}/bowtie2/merged/${NAME}.sorted.filtered.dupMarked.blfilt_merged.bam
MACS2_OUTDIR=${DIRPATH}/macs2/merged
MACS2_LOG=${DIRPATH}/macs2/merged/logs/${NAME}.log

echo "Control sample:" ${CONTROL_BAM}
echo "Treatment sample:" ${TREATMENT_BAM}

### macs2 ###

# Stricter cutoff (q=0.1 or p=1e-5)
macs2 callpeak -t ${TREATMENT_BAM} \
	-c ${CONTROL_BAM} \
	-f BAMPE \
	-g 1.87e9 \
	-n ${NAME}_1e-5 \
	-p 1e-5 \
	--broad \
	--broad-cutoff 1e-5 \
	--outdir ${MACS2_OUTDIR} \
	--keep-dup auto \
	--SPMR 2> ${MACS2_LOG}_1e-5

### 08/25/21: H3K27Ac and H3K4me3 should use NARROW (not broad) peak calling ###

macs2 callpeak -t ${TREATMENT_BAM} \
    -c ${CONTROL_BAM} \
    -f BAMPE \
    -g 1.87e9 \
    -n ${NAME}_1e-5 \
    -q 0.05 \
    --outdir ${MACS2_OUTDIR} \
    --keep-dup auto > ${MACS2_LOG}1e-5_narrowPeaks
