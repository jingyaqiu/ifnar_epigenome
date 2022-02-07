#!/bin/bash
set -e

# macs2 version 2.2.7.1

source /home/jingyaq/Minn/data/R499_CnR/scripts/activate_peak_calling_CnR.sh

DIRPATH=$(pwd)

NAME=$1
CONTROL_MACS2=$2
CONTROL_SEACR=$3

# COND=$1
# HISTONE=$2
# REP=$3
# NAME=${COND}_${HISTONE}_${REP}

###############################
### Peak calling with MACS2 ###
###############################

echo "Peak calling with MACS2"

mkdir -p ${DIRPATH}/macs2/logs

# For MACS2 peak calling, parameters used were macs2 callpeak – t input_file –p 1e-5 –f BEDPE/BED(Paired End vs. Single End sequencing data) –keep-dup all –n out_name. 
# Found peaks using macs2 2.1.1.20160309 with parameters: callpeak -t fragments_not_- mask.bed -f BED -g hs –keep-dup all -p 1e-5 -n not_mask –SPMR.

CONTROL_BAM=${CONTROL_MACS2}
TREATMENT_BAM=${DIRPATH}/bowtie2/${NAME}.sorted.filtered.dupMarked.blfilt.bam
MACS2_OUTDIR=${DIRPATH}/macs2
MACS2_LOG=${DIRPATH}/macs2/logs/${NAME}.log

echo "Control sample:" ${CONTROL_BAM}
echo "Treatment sample:" ${TREATMENT_BAM}

# Relax significance threshold to run IDR later
# "To run IDR the recommendation is to run MACS2 less stringently. The IDR algorithm requires sampling of both signal and noise distributions to separate the peaks into two groups, so having a more liberal threshold allows us to bring in some noise."
# https://hbctraining.github.io/Intro-to-ChIPseq/lessons/07_handling-replicates-idr.html

# Merged IgG file??
macs2 callpeak -t ${TREATMENT_BAM} -c ${CONTROL_BAM} -p 1e-3 -f BAMPE --broad --broad-cutoff 1e-3 --keep-dup auto -n ${NAME}_1e-3 --outdir ${MACS2_OUTDIR} -g 1.87e9 --SPMR 2> ${MACS2_LOG}_1e-3

macs2 callpeak -t ${TREATMENT_BAM} -c ${CONTROL_BAM} -p 1e-4 -f BAMPE --broad --broad-cutoff 1e-4 --keep-dup auto -n ${NAME}_1e-4 --outdir ${MACS2_OUTDIR} -g 1.87e9 --SPMR 2> ${MACS2_LOG}_1e-4

# Stricter cutoff (q=0.1 or p=1e-5)
macs2 callpeak -t ${TREATMENT_BAM} -c ${CONTROL_BAM} -p 1e-5 -f BAMPE --broad --broad-cutoff 1e-5 --keep-dup auto -n ${NAME}_1e-5 --outdir ${MACS2_OUTDIR} -g 1.87e9 --SPMR 2> ${MACS2_LOG}_1e-5

# Sort peaks by -log10(pval)
sort -k8,8nr ${DIRPATH}/macs2/${NAME}_1e-3_peaks.broadPeak > ${DIRPATH}/macs2/${NAME}_1e-3_peaks_sorted.broadPeak

###############################
### Peak calling with SEACR ###
###############################

echo "Peak calling with SEACR"

mkdir -p ${DIRPATH}/seacr
cd ${DIRPATH}/seacr

CONTROL_BG=${CONTROL_SEACR}
TREATMENT_BG=${DIRPATH}/bowtie2/bedgraph/${NAME}.fragments.bedgraph

echo "Control sample:" $CONTROL_BG
echo "Treatment sample:" $TREATMENT_BG

bash /home/jingyaq/SEACR_1.3.sh ${TREATMENT_BG} ${CONTROL_BG} norm stringent ${DIRPATH}/seacr/${NAME}_seacr.peaks
bash /home/jingyaq/SEACR_1.3.sh ${TREATMENT_BG} 0.01 norm stringent ${DIRPATH}/seacr/${NAME}_seacr_top0.01.peaks
