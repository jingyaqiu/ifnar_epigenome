## ATAC-seq processing

*Last updated 09/12/20*

We processed ATAC-seq data by running the ENCODE ATAC-seq pipeline from Anshul Kundaje's lab. The [version](https://github.com/kundajelab/atac_dnase_pipelines) of the pipeline we ran in January 2018 has since been deprecated and updated for ENCODE-DCC.

### Dataset

Our original dataset consists of paired-end ATAC-seq reads generated from day 15 sorted mouse melanoma cells from B16, R499, B16 STAT1 KO, and R499 STAT1 KO cell lines (20 libraries total, n=5 biological replicates for each condition). Each sample was resequenced twice and reads from the same sample were pooled to improve sequencing depth.

### Install ENCODE ATAC-seq/DNase-seq pipeline

Requires BigDataScript 0.99999e and Java 8 (jdk >= 1.8)

```
export PATH=$PATH:$HOME/.bds
export PATH=/home/jingyaq/jdk1.8.0_144/bin:$PATH
```
Install pipeline and mm10 genome data

```
git clone https://github.com/kundajelab/atac_dnase_pipelines --recursive
bash install_dependencies.sh
bash install_genome_data.sh mm10 atac_dnase_pipelines/mm10
```

### Run pipeline

This pipeline accepts raw fastq files as input, trims adapters with cutadapt, aligns with Bowtie2, calls peaks with macs2, filters out peaks in blacklist regions, identifies IDR peaks, and reports library and sequencing QC metrics.

We run this pipeline aligning to mm10 and the --enable_idr flag.

```
conditions=(B16_cas B16_SKO R499_cas R499_SKO)
for cond in ${conditions[@]}
do
	bsub -M 120000 -R "span[hosts=1] rusage [mem=120000]" -o output_logs/ENCODE_${cond}.out -e output_logs/ENCODE_${cond}.err \
	bash scripts/processing/ATAC_ENCODE_job_submits/${cond}_ATAC.sh
done
```
### Filter chrM reads and tn5 shift

**process\_ATAC\_bam.sh**

1. Remove chrM: The ENCODE pipeline removes unmapped reads and improperly paired reads (-F 1804 -f 2), as well as duplicate reads with Picard's MarkDuplicates. We further process the output bam files by filtering out non-unique alignments (-q 30) and mitochondrial reads. 

2. tn5 shift: In order to perform DNA footprinting analysis, we reconstruct Tn5 transposition events by shifting the 5' position of the positive strand +4 bp and the negative strand -5 bp.

3. Generate bigwigs: For visualization, we use bamCoverage from the deepTools python package to convert our final, processed bam files to normalized coverage files.

```
cd /home/jingyaq/Minn/data/ifnar_epigenome

samples=(B16_cas_1 B16_cas_2 B16_cas_3 B16_cas_4 B16_cas_5 B16_SKO_1 B16_SKO_2 B16_SKO_3 B16_SKO_4 B16_SKO_5 R499_cas_1 R499_cas_2 R499_cas_3 R499_cas_4 R499_cas_5 R499_SKO_1 R499_SKO_2 R499_SKO_3 R499_SKO_4 R499_SKO_5)

for id in ${samples[@]}
do
	echo ${id}
	bsub -M 120000 -R "span[hosts=1] rusage [mem=120000]" -o output_logs/process_ATAC_bam_${id}.out -e output_logs/process_ATAC_bam_${id}.err \
	bash scripts/processing/process_ATAC_bam.sh ${id}
done
```
### Merge bam files by condition

```
dirpath="/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC"

mkdir -p ${dirpath}/bam/filtered/merged
mkdir -p ${dirpath}/bw/filtered/merged

conditions=(B16_cas B16_SKO R499_cas R499_SKO)
for cond in ${conditions[@]}
do
	echo $cond
	bam=(${dirpath}/bam/filtered/${cond}*.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam)
	bam_merged=${dirpath}/bam/filtered/merged/${cond}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted_merged.bam
	bw_merged=${dirpath}/bw/filtered/merged/${cond}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted_merged.bw
	
	# Merge bam files
	samtools merge ${bam_merged} ${bam[@]}
	samtools index ${bam_merged}
	
	# Generate merged bigwigs
	 ~/.local/bin/bamCoverage -b ${bam_merged} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${bw_merged}
done
```
**Count number of reads**

```
# Individual samples
cd /home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/bam/filtered
for file in *.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done

# Merged samples
cd /home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/bam/filtered/merged
for file in *.bam;  do echo $file; samtools idxstats $file | awk '{s+=$3+$4} END {print s}'; done
```

### Downsample

Downsample individual replicates to 10 million reads per sample.

```
DIRPATH=/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC
cd ${DIRPATH}

mkdir -p ${DIRPATH}/bam/filtered/downsample ${DIRPATH}/bw/filtered/downsample

samples=(B16_cas_1 B16_cas_2 B16_cas_3 B16_cas_4 B16_cas_5 B16_SKO_1 B16_SKO_2 B16_SKO_3 B16_SKO_4 B16_SKO_5 R499_cas_1 R499_cas_2 R499_cas_3 R499_cas_4 R499_cas_5 R499_SKO_1 R499_SKO_2 R499_SKO_3 R499_SKO_4 R499_SKO_5)

for SAMPLE in ${samples[@]}
do
	bam=/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/bam/filtered/${SAMPLE}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam
	bam_downsample=/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/bam/filtered/downsample/${SAMPLE}_downsample_10M.bam
	bw_downsample=/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/bw/filtered/downsample/${SAMPLE}_downsample_10M.bw
	echo ${bam}
	
	frac=$( samtools idxstats ${bam} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=10000000/total; if (frac > 1) {print 1} else {print frac}}' )
	
	samtools view -bs ${frac} ${bam} > ${bam_downsample}
	samtools index ${bam_downsample}
	
	# Generate merged bigwigs
	bsub -M 90000 -R "span[hosts=1] rusage [mem=90000]" ~/.local/bin/bamCoverage -b ${bam_downsample} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2652783500 --extendReads --centerReads -o ${bw_downsample}
done
```
Downsample merged files to 53 million reads

```
DIRPATH=/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC
cd ${DIRPATH}

mkdir -p ${DIRPATH}/bam/filtered/merged/downsample bw/filtered/merged/downsample

conditions=(B16_cas B16_SKO R499_cas R499_SKO)
for CONDITION in ${conditions[@]}
do
	echo ${CONDITION}
	bam=${DIRPATH}/bam/filtered/merged/${CONDITION}.R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted_merged.bam
	bam_downsample=${DIRPATH}/bam/filtered/merged/downsample/${CONDITION}_downsample_53M.bam
	
	frac=$( samtools idxstats ${bam} | cut -f3 | awk 'BEGIN {total=0} {total += $1} END {frac=53000000/total; if (frac > 1) {print 1} else {print frac}}' )
	
	samtools view -bs ${frac} ${bam} > ${bam_downsample}
	samtools index ${bam_downsample}
	
	bw_downsample=${DIRPATH}/bw/filtered/merged/downsample/${CONDITION}_downsample_53M_merged.bw
	
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
Rscript generate_consensus_ATAC_matrix.R

```