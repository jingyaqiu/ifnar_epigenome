## RNA-seq processing

*Last updated 09/06/20*

### Dataset

Our RNA-seq libraries are paired to original ATAC-seq libraries (same samples). Paired-end, total RNA-seq reads were generated from day 15 sorted mouse melanoma cells from B16, R499, B16 STAT1 KO, and R499 STAT1 KO cell lines (16 libraries total).

### Processing

**process\_RNA.sh**

1. Trim adapters
2. Quantification with salmon (version 1.2.0)
3. For visualization, we align reads using STAR and convert properly and uniquely aligned reads (-F 1804 -q 5) to normalized bigwigs using bamCoverage from the deepTools python package.

```
cd /home/jingyaq/Minn/data/ifnar_epigenome

samples=(B16_cas_1 B16_cas_4 B16_cas_5 B16_SKO_1 B16_SKO_2 B16_SKO_3 B16_SKO_4 B16_SKO_5 R499_cas_1 R499_cas_2 R499_cas_4 R499_SKO_1 R499_SKO_2 R499_SKO_3 R499_SKO_4 R499_SKO_5)

for id in ${samples[@]}
do
	echo ${id}
	bsub -M 120000 -R "span[hosts=1] rusage [mem=120000]" -o output_logs/process_RNA_${id}.out -e output_logs/process_RNA_${id}.err \
	bash scripts/processing/process_RNA.sh ${id}
done
```
### Merge bam files by condition

```
cd /home/jingyaq/Minn/data/ifnar_epigenome/

dirpath="/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/RNA"

mkdir -p ${dirpath}/star/filtered/merged
mkdir -p ${dirpath}/bw/filtered/merged

conditions=(B16_cas B16_SKO R499_cas R499_SKO)
for cond in ${conditions[@]}
do
	echo $cond
	bam=(${dirpath}/star/filtered/${cond}*_filtered_deduped.bam)
	bam_merged=${dirpath}/star/filtered/merged/${cond}_filtered_deduped_merged.bam
	bw_merged=${dirpath}/bw/filtered/merged/${cond}_filtered_deduped_merged.bw
	
	# Merge bam files
	samtools merge ${bam_merged} ${bam[@]}
	samtools index ${bam_merged}
	
	# Generate merged bigwigs
	 ~/.local/bin/bamCoverage -b ${bam_merged} -of bigwig --binSize 1 -p 16 --normalizeUsing RPGC --effectiveGenomeSize 2407883318 --extendReads --centerReads -o ${bw_merged}
done
```

### Create normalized RNA count matrix

generate\_RNA\_count\_matrix.R

1. Use tximport package to load in salmon-quantified transcript counts, reduce to gene level.
2. Normalize count matrix using rlog() function from DESeq2.

```
Rscript generate_RNA_count_matrix.R
```

