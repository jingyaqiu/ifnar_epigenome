## Revisions ATAC-seq

*Last updated 11/10/20*

#### Dataset

Our revisions experiments consisted of paired-end ATAC-seq generated from sorted mouse melanoma cells from B16, R499, B16 STAT1 KO, R499 STAT1 KO, R499 IRF3 KO, and R499 STAT1/IRF3 DKO (24 libraries total, n=4 biological replicates each).

#### Demultiplex

bcl2fastq2 version 2.20.0

```
module load bcl2fastq2/v2.20.0.422

cd /home/jingyaq/Minn/data/revisions_NatCancer

mkdir -p data/201019 data/201020 data/201021

# 201019
bsub -o output_logs/bcl2fastq_atac_revisions_201019.out -e output_logs/bcl2fastq_atac_revisions_201019.err -M 100000 -R "span[hosts=1] rusage [mem=100000]" \
    bcl2fastq --runfolder-dir /home/jingyaq/Minn/data/revisions_NatCancer/201019_NB551627_0192_AH2CKWBGXH \
    -p 12 \
    --output-dir /home/jingyaq/Minn/data/revisions_NatCancer/data/201019/fastq \
    --sample-sheet /home/jingyaq/Minn/data/revisions_NatCancer/data/sample_sheets/201019_ATAC_sample_sheet.csv \
    --no-lane-splitting
    
# 201020
bsub -o output_logs/bcl2fastq_atac_revisions_201020.out -e output_logs/bcl2fastq_atac_revisions_201020.err -M 100000 -R "span[hosts=1] rusage [mem=100000]" \
    bcl2fastq --runfolder-dir /home/jingyaq/Minn/data/revisions_NatCancer/201020_NB551627_0193_AH2CTKBGXH \
    -p 12 \
    --output-dir /home/jingyaq/Minn/data/revisions_NatCancer/data/201020/fastq \
    --sample-sheet /home/jingyaq/Minn/data/revisions_NatCancer/data/sample_sheets/201020_ATAC_sample_sheet.csv \
    --no-lane-splitting
    
# 201021
bsub -o output_logs/bcl2fastq_atac_revisions_201021.out -e output_logs/bcl2fastq_atac_revisions_201021.err -M 100000 -R "span[hosts=1] rusage [mem=100000]" \
    bcl2fastq --runfolder-dir /home/jingyaq/Minn/data/revisions_NatCancer/201021_NB551627_0194_AH25YYBGXH \
    -p 12 \
    --output-dir /home/jingyaq/Minn/data/revisions_NatCancer/data/201021/fastq \
    --sample-sheet /home/jingyaq/Minn/data/revisions_NatCancer/data/sample_sheets/201021_ATAC_sample_sheet.csv \
    --no-lane-splitting
```
#### FastQC

```
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/201019/fastq
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/201020/fastq
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/201021/fastq

mkdir fastqc

for file in *fastq.gz;  do /home/jingyaq/FastQC/fastqc $file -o fastqc; done

# Aggregate fastqc files

cd /home/jingyaq/Minn/data/revisions_NatCancer/data
mkdir processed/qc/fastqc

cp 201019/fastq/fastqc/* processed/qc/fastqc
cp 201020/fastq/fastqc/* processed/qc/fastqc
cp 201021/fastq/fastqc/* processed/qc/fastqc

cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/qc
export PATH="/home/jingyaq/anaconda3-new/bin:$PATH"
multiqc fastqc
mv multiqc* fastqc

```
### NEW ENCODE-DCC ATAC-SEQ PIPELINE

ATAC experiments for revisions were processed with the updated [ENCODE-DCC ATAC-seq pipeline](https://github.com/ENCODE-DCC/atac-seq-pipeline). 

[ENCODE standards](https://www.encodeproject.org/atac-seq/) <br>
[Github](https://github.com/ENCODE-DCC/atac-seq-pipeline) <br>

#### Run

*B16 SKO 2, B16 SKO 3, and R499 DKO 2 were removed because libraries were empty (no fastq reads).*

```
export PATH="/home/jingyaq/miniconda3/bin/:$PATH"
source activate /home/jingyaq/miniconda3/envs/encode-atac-seq-pipeline

conditions=(B16_WT B16_SKO R499_WT R499_SKO R499_IRF3KO R499_DKO)

for CONDITION in ${conditions[@]}
do
	echo ${CONDITION}

	WDL=/home/jingyaq/atac-seq-pipeline/atac.wdl
	INPUT_JSON=/home/jingyaq/Minn/data/revisions_NatCancer/data/input_json/${CONDITION}_ATAC.json
	
	OUTPUT_DIR=/home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}
	mkdir ${OUTPUT_DIR} # make a separate directory for each workflow
	
	cd ${OUTPUT_DIR}
	bsub -n 16 -M 150000 -R "span[hosts=1] rusage [mem=150000]" \
		-o /home/jingyaq/Minn/data/revisions_NatCancer/output_logs/encode_${CONDITION}.out \
		-e /home/jingyaq/Minn/data/revisions_NatCancer/output_logs/encode_${CONDITION}.err \
		caper run ${WDL} -i ${INPUT_JSON}
done

rm -r ~/atac-seq-pipeline/caper-local-dir/*
```
#### Gather results

```
export PATH="/home/jingyaq/miniconda3/bin/:$PATH"
source activate /home/jingyaq/miniconda3/envs/encode-atac-seq-pipeline
pip install qc2tsv

conditions=(B16 B16_SKO R499_WT R499_SKO R499_IRF3KO R499_DKO)

for CONDITION in ${conditions[@]}
do
	echo ${CONDITION}
	
	METADAT_JSON=/home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/atac/*/metadata.json

	cd /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}
	croo ${METADAT_JSON}
done

# QC metrics
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/
qc2tsv B16/qc/qc.json \
	B16_SKO/qc/qc.json \
	R499_WT/qc/qc.json \
	R499_SKO/qc/qc.json \
	R499_IRF3KO/qc/qc.json \
	R499_DKO/qc/qc.json > qc.metrics.tsv
```
**Organize output**

Move filtered bam, bigwigs, and called peaks (raw individual, IDR) output files to processed/ directory

```
cd /home/jingyaq/Minn/data/revisions_NatCancer/data
mkdir -p processed/bam processed/peaks/idr processed/mat processed/bw processed/qc

conditions=(B16 B16_SKO R499_WT R499_SKO R499_IRF3KO R499_DKO)

for CONDITION in ${conditions[@]}
do
	echo ${CONDITION}
	cd /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/align
	
	reps=($(ls -d rep*))
	
	for REP in ${reps[@]}
	do
		echo ${REP}
		cp /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/align/${REP}/*.trim.srt.nodup.no_chrM_MT.bam /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam
		cp /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/peak/${REP}/*.trim.srt.nodup.no_chrM_MT.tn5.pval0.01.300K.narrowPeak.gz /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/peaks
		cp /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/peak/idr_reproducibility/idr.optimal_peak.narrowPeak.gz /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/peaks/idr/${CONDITION}_idr.optimal_peak.narrowPeak.gz
		cp /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/peak/idr_reproducibility/idr.conservative_peak.narrowPeak.gz /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/peaks/idr/${CONDITION}_idr.conservative_peak.narrowPeak.gz
		cp -r /home/jingyaq/Minn/data/revisions_NatCancer/data/${CONDITION}/qc /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/qc/${CONDITION}
	done
done

# Combine IDR peaks
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/peaks/idr

gunzip *
cat *.optimal_peak.narrowPeak > combined_idr.optimal_peak.narrowPeak

```
### tn5 shift and generate bigwigs

**process\_ATAC\_bam\_ENCODE-DCC.sh**

bam files are already filtered for properly paired, deduped, and non-chrM reads.

1. tn5 shift: In order to perform DNA footprinting analysis, we reconstruct Tn5 transposition events by shifting the 5' position of the positive strand +4 bp and the negative strand -5 bp.

2. Generate bigwigs: For visualization, we use bamCoverage from the deepTools python package to convert our final, processed bam files to normalized coverage files.

```
cd /home/jingyaq/Minn/data/revisions_NatCancer

samples=(B16_WT-1 B16_WT-2 B16_WT-3 B16_WT-4 B16_SKO-1 B16_SKO-4 R499_WT-1 R499_WT-2 R499_WT-3 R499_WT-4 R499_SKO-1 R499_SKO-2 R499_SKO-3 R499_SKO-4 R499_IRF3KO-1 R499_IRF3KO-2 R499_IRF3KO-3 R499_IRF3KO-4 R499_DKO-1 R499_DKO-3 R499_DKO-4)

for id in ${samples[@]}
do
	echo ${id}
	bsub -M 90000 -R "span[hosts=1] rusage [mem=90000]" -o output_logs/process_ATAC_bam/${id}.out -e output_logs/process_ATAC_bam/${id}.err \
	bash /home/jingyaq/Minn/data/ifnar_epigenome/scripts/processing/process_ATAC_bam_ENCODE-DCC.sh ${id}
done
```
### Merge bam files by condition

*R499 SKO 3 removed because poor mapping rate (<10%). B16 WT 3 removed because lacking di-nucleosome peak, NFR too high, likely over-transposed?.*

```
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed

mkdir -p bam/merged bw/merged

conditions=(B16_WT B16_SKO R499_WT R499_SKO R499_IRF3KO R499_DKO)

for CONDITION in ${conditions[@]}
do
	echo $CONDITION
	bam=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/${CONDITION}*.trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam
	bam_merged=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/merged/${CONDITION}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted_merged.bam
	bw_merged=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bw/merged/${CONDITION}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted_merged.bw
	
	# Merge bam files
	echo ${bam[@]}
	samtools merge -@ 12 ${bam_merged} ${bam[@]}
	samtools index ${bam_merged}
	
	# Generate merged bigwigs
	 ~/.local/bin/bamCoverage -b ${bam_merged} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${bw_merged}
done
```
**Count number of reads**

```
# Individual samples
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam
for file in *.trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done

# Merged samples
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/merged
for file in *.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done
```

### Downsample

Downsample individual replicates to 65 million reads per sample.

```
cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed

mkdir -p ${dirpath}/bam/downsample

samples=(B16_WT-1 B16_WT-2 B16_WT-3 B16_WT-4 B16_SKO-1 B16_SKO-4 R499_WT-1 R499_WT-2 R499_WT-3 R499_WT-4 R499_SKO-1 R499_SKO-2 R499_SKO-4 R499_IRF3KO-1 R499_IRF3KO-2 R499_IRF3KO-3 R499_IRF3KO-4 R499_DKO-1 R499_DKO-3 R499_DKO-4)

for SAMPLE in ${samples[@]}
do
	bam=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/${SAMPLE}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam
	bam_downsample=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/downsample/${SAMPLE}_downsample_65M.bam
	
	echo ${bam}
	
	frac=$( samtools idxstats ${bam} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=65000000/total; if (frac > 1) {print 1} else {print frac}}' )
	
	samtools view -bs ${frac} ${bam} > ${bam_downsample}
	samtools index ${bam_downsample}
done
```
```
# Check
cd bam/downsample
for file in *_downsample_65M.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done
```
**Downsample merged files to 170 million reads**

```
# Check
cd bam/merged
for file in *.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done

cd /home/jingyaq/Minn/data/revisions_NatCancer/data/processed
mkdir -p bam/merged/downsample bw/merged/downsample

conditions=(B16_WT B16_SKO R499_WT R499_SKO R499_IRF3KO R499_DKO)

for CONDITION in ${conditions[@]}
do
	echo $CONDITION
	
	bam=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/merged/${CONDITION}.trim.srt.nodup.no_chrM_MT_tn5shift_sorted_merged.bam
	bam_downsample=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bam/merged/downsample/${CONDITION}_downsample_170M_merged.bam
	
	echo ${bam}
	
	frac=$( samtools idxstats ${bam} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=170000000/total; if (frac > 1) {print 1} else {print frac}}' )
	
	/home/jingyaq/samtools-1.11/samtools view -bs ${frac} ${bam} > ${bam_downsample}
	/home/jingyaq/samtools-1.11/samtools index ${bam_downsample}
	
	bw_downsample=/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/bw/merged/downsample/${CONDITION}_downsample_170M_merged.bw
	
	# Generate merged bigwigs
	bsub -M 100000 -R "span[hosts=1] rusage [mem=100000]" ~/.local/bin/bamCoverage -b ${bam_downsample} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${bw_downsample}
done

```
### Generate consensus ATAC count matrix

**generate\_consensus\_ATAC\_matrix.R**

1. Generate consensus peakset: Use readNarrowpeaks() function from chromVAR package to read in IDR peaks (from running MACS2 in ENCODE pipeline), resize to a fixed width around centers (n=750 bp), and remove peaks that overlap peak with stronger signal.

2. Count Tn5 insertions: Use countOverlaps() from GRanges package to count insertion counts (corrected for Tn5 offset) at each consensus peak.

3. Normalize count matrix: We use the regularized-log transform rlog() function from DESeq2 to normalize the insertion count matrix.

4. Differential accessibility analysis: Use Diffbind package to identify differentially accessible peaks for contrasts of interest. *The assumption that there are no global differences in chromatin accessibility between B16 and Res499 is likely not appropriate. Use full library size normalization (not reads in peaks) to capture global accessibility differences.*

5. Identify newly opened (fold change > 1, pval < 0.05) and constitutive peaks.

```
Rscript generate_consensus_ATAC_matrix_ENCODE-DCC.R
### ```