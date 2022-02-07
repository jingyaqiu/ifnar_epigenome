## Processing pipeline for CUT&RUN data

*Last updated 02/10/21*

### Dataset

This experiment consists of genome-wide profiling of histone modifications (H3K4me3, H3K27Ac, and H3K4me1) and a control IgG antibody by CUT&RUN. The data were generated from sorted mouse melanoma cells from B16 WT, R499 WT, B16 STAT1 KO, and R499 STAT1 KO conditions. Each condition for each assayed histone modification consisted of 2 biological replicates for a total of 24 libraries. Some libraries were sequenced twice to increase sequencing depth.

### Demultiplex

bcl2fastq2 version 2.20.0

```
module load bcl2fastq2/v2.20.0.422

cd /home/jingyaq/Minn/data/R499_CnR

# Run FINAL1 #
mkdir -p /home/jingyaq/Minn/data/R499_CnR/data/run_FINAL1
bsub -o output_logs/bcl2fastq_run_final1.out -e output_logs/bcl2fastq_run_final1.err -M 100000 -R "span[hosts=1] rusage [mem=100000]" \
    bcl2fastq --runfolder-dir /home/jingyaq/Minn/data/R499_CnR/201121_NB551627_0206_AHYG35BGXC \
    -p 16 \
    --output-dir /home/jingyaq/Minn/data/R499_CnR/data/run_FINAL1/fastq \
    --sample-sheet /home/jingyaq/Minn/data/R499_CnR/data/sample_sheets/CnR_112120_FINAL1.csv \
    --no-lane-splitting

# Run FINAL2 #
mkdir -p /home/jingyaq/Minn/data/R499_CnR/data/run_FINAL2
bsub -o output_logs/bcl2fastq_run_final2.out -e output_logs/bcl2fastq_run_final2.err -M 100000 -R "span[hosts=1] rusage [mem=100000]" \
    bcl2fastq --runfolder-dir /home/jingyaq/Minn/data/R499_CnR/201122_NB551627_0207_AH7TLHBGXG \
    -p 16 \
    --output-dir /home/jingyaq/Minn/data/R499_CnR/data/run_FINAL2/fastq \
    --sample-sheet /home/jingyaq/Minn/data/R499_CnR/data/sample_sheets/CnR_112120_FINAL2.csv \
    --no-lane-splitting
```
Rename fastq files

```
rm(list=ls())

sampleID <- paste0(rep(c("B16_WT_H3K4me1", "B16_WT_H3K4me3", "B16_WT_H3K27Ac", "B16_WT_IgG", "B16_SKO_H3K4me1", "B16_SKO_H3K4me3", "B16_SKO_IgG", "B16_SKO_H3K27Ac", "R499_WT_H3K4me1", "R499_WT_H3K4me3", "R499_WT_H3K27Ac", "R499_WT_IgG", "R499_SKO_H3K4me1", "R499_SKO_H3K4me3", "R499_SKO_H3K27Ac", "R499_SKO_IgG"), each=2), "_rep", rep(1:2, 12))

# run_FINAL1
files1 <- list.files("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL1/fastq", pattern=".gz", full.names=T)
files1 <- files1[-grep(files1, pattern="Undetermined")]
files1.renamed <- gsub("S\\d+_R","R",files1)

# Check
samp1 <- sapply(strsplit(files1, split="/|_"), function(x) paste(x[12:15], collapse="_"))
samp1.renamed <- sapply(strsplit(files1.renamed, split="/|_"), function(x) paste(x[12:15], collapse="_"))
all(samp1 %in% sampleID)
all(samp1 == samp1.renamed)

for(i in 1:length(files1.renamed)) {
    print(files1[i])
    print(files1.renamed[i])
    system(paste("mv", files1[i], files1.renamed[i], sep=" "))
}
                        
# run_FINAL2
files2 <- list.files("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL2/fastq", pattern=".gz", full.names=T)
files2 <- files2[-grep(files2, pattern="Undetermined")]
files2.renamed <- gsub("S\\d+_R","R",files2)

# Check
samp2 <- sapply(strsplit(files2, split="/|_"), function(x) paste(x[12:15], collapse="_"))
samp2.renamed <- sapply(strsplit(files2.renamed, split="/|_"), function(x) paste(x[12:15], collapse="_"))
all(samp2 %in% sampleID)
all(samp2 == samp2.renamed)

for(i in 1:length(files2.renamed)) {
    print(samp2[i])
    print(files2[i])
    print(files2.renamed[i])
    system(paste("mv", files2[i], files2.renamed[i], sep=" "))
}
```
#### FastQC

```
DIRPATH="/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL1"
DIRPATH="/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL2"

cd ${DIRPATH}/fastq
mkdir fastqc

for FILE in *fastq.gz
do 
    bsub /home/jingyaq/FastQC/fastqc $FILE -o fastqc
done

# multiqc
export PATH="/home/jingyaq/anaconda3-new/bin:$PATH"
multiqc fastqc
mv multiqc* fastqc
```
#### Combine reads from technical replicates

Exclude poor quality samples:

R499\_SKO\_H3K27Ac\_rep1 (run2) <br>
R499\_SKO\_H3K4me3\_rep2 (run1/run2) <br>
R499\_SKO\_H3K4me1\_rep2 (run2)

```
cd /home/jingyaq/Minn/data/R499_CnR/data
mkdir -p /home/jingyaq/Minn/data/R499_CnR/data/run_FINAL/fastq

###########################################

rm(list=ls())

sampleID <- paste0(rep(c("B16_WT_H3K4me1", "B16_WT_H3K4me3", "B16_WT_H3K27Ac", "B16_WT_IgG", "B16_SKO_H3K4me1", "B16_SKO_H3K4me3", "B16_SKO_IgG", "B16_SKO_H3K27Ac", "R499_WT_H3K4me1", "R499_WT_H3K4me3", "R499_WT_H3K27Ac", "R499_WT_IgG", "R499_SKO_H3K4me1", "R499_SKO_H3K4me3", "R499_SKO_H3K27Ac", "R499_SKO_IgG"), each=2), "_rep", rep(1:2, 12))

files1 <- list.files("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL1/fastq", pattern=".gz", full.names=T)
files2 <- list.files("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL2/fastq", pattern=".gz", full.names=T)

for(i in 1:length(sampleID)) {
    files.R1 <- c(grep(files1, pattern=paste0(sampleID[i], "_R1"), value=T), grep(files2, pattern=paste0(sampleID[i], "_R1"), value=T))
    files.R2 <- c(grep(files1, pattern=paste0(sampleID[i], "_R2"), value=T), grep(files2, pattern=paste0(sampleID[i], "_R2"), value=T))
    print(sampleID[i])
    print(files.R1)
    print(files.R2)
    
    file.out.R1 <- paste0("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL/fastq/", sampleID[i], "_R1_001.fastq.gz")
    file.out.R2 <- paste0("/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL/fastq/", sampleID[i], "_R2_001.fastq.gz")
    
    system(paste("cat", paste(files.R1, collapse=" "), ">", file.out.R1, sep=" "))
    system(paste("cat", paste(files.R2, collapse=" "), ">", file.out.R2, sep=" "))
}
```
### Processing

**process_CnR.sh**

- Adapter trimming with cutadapt
- Align reads with bowtie2
- Filter reads: remove unmapped reads and improperly paired reads (-F 1804 -f 2), low-quality alignments (-q 5). (duplicate reads kept)
- Generate bigwigs: For visualization, we use bamCoverage from the deepTools python package to convert our final, processed bam files to normalized coverage files.
- Generate bedgraphs and binned genome coverage files for SEACR peak calling and assessing replicate reproducibility

*Blacklist reads filtered*

```
### SUBMIT JOBS ###

DIRPATH="/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL"
cd ${DIRPATH}

# Run final (n=32)
ids=(B16_WT_IgG_rep1 B16_WT_IgG_rep2 \
     B16_WT_H3K4me1_rep1 B16_WT_H3K4me1_rep2 \
     B16_WT_H3K4me3_rep1 B16_WT_H3K4me3_rep2 \
     B16_WT_H3K27Ac_rep1 B16_WT_H3K27Ac_rep2 \
     B16_SKO_IgG_rep1 B16_SKO_IgG_rep2 \
     B16_SKO_H3K4me1_rep1 B16_SKO_H3K4me1_rep2 \
     B16_SKO_H3K4me3_rep1 B16_SKO_H3K4me3_rep2 \
     B16_SKO_H3K27Ac_rep1 B16_SKO_H3K27Ac_rep2 \
     R499_WT_IgG_rep1 R499_WT_IgG_rep2 \
     R499_WT_H3K4me1_rep1 R499_WT_H3K4me1_rep2 \
     R499_WT_H3K4me3_rep1 R499_WT_H3K4me3_rep2 \
     R499_WT_H3K27Ac_rep1 R499_WT_H3K27Ac_rep2 \
     R499_SKO_IgG_rep1 R499_SKO_IgG_rep2 \
     R499_SKO_H3K4me1_rep1 R499_SKO_H3K4me1_rep2 \
     R499_SKO_H3K4me3_rep1 R499_SKO_H3K4me3_rep2 \
     R499_SKO_H3K27Ac_rep1 R499_SKO_H3K27Ac_rep2)

for ID in ${ids[@]}
do
    echo $ID
    bsub -n 16 -M 80000 -R "span[hosts=1] rusage [mem=80000]" -o ${DIRPATH}/output_logs/${ID}_run_FINAL.out -e ${DIRPATH}/output_logs/${ID}_run_FINAL.err bash /home/jingyaq/Minn/data/R499_CnR/scripts/process_CnR.sh $ID
done
```
**peak\_calling\_CnR.sh**

*Didn't use these peaks for downstream analysis because low depth in some libraries led to bad peak calls*

- macs2 broad peak calling (default threshold p=1e-05, also try relaxing to p=1e-04, p=1e-03) 
- SEACR peak calling (specific to CUT&RUN data)

```
### SUBMIT JOBS ###

DIRPATH="/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL"
cd ${DIRPATH}
mkdir -p output_logs/peak_calling
rm output_logs/peak_calling/*
rm -r macs2

### Use correct IgG control ###

conditions=(B16_WT B16_SKO R499_WT R499_SKO)
histones=(H3K27Ac H3K4me1 H3K4me3 IgG)
reps=(rep1 rep2)

for COND in ${conditions[@]}
do
    for HIST in ${histones[@]}
    do
        for REP in ${reps[@]}
        do
            NAME=${COND}_${HIST}_${REP}
            echo ${NAME}
        
            bsub -n 12 -M 50000 -R "span[hosts=1] rusage [mem=50000]" -o output_logs/peak_calling/${NAME}.out -e output_logs/peak_calling/${NAME}.err bash /home/jingyaq/Minn/data/R499_CnR/scripts/peak_calling_CnR.sh ${COND} ${HIST} ${REP}
        done
    done
done
```
**merge_conditions.sh**

- Merge final processed bam reads for replicates
- Generate merged bigwigs
- Peak calling on merged bam files - **use for generating normalized count matrix, downstream analysis**

macs2 parameters:
- broad peak calling (narrow for H3K27ac, H3K4me3?)
- -p 1e-05, --broad-cutoff = 1e-05
- --keep-dup auto
- --SPMR 2

```
DIRPATH="/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL"
cd ${DIRPATH}

names=(B16_WT_IgG B16_WT_H3K4me1 B16_WT_H3K4me3 B16_WT_H3K27Ac \
       B16_SKO_IgG B16_SKO_H3K4me1 B16_SKO_H3K4me3 B16_SKO_H3K27Ac \
       R499_WT_IgG R499_WT_H3K4me1 R499_WT_H3K4me3 R499_WT_H3K27Ac \
       R499_SKO_IgG R499_SKO_H3K4me1 R499_SKO_H3K4me3 R499_SKO_H3K27Ac)

conditions=(B16_WT B16_SKO R499_WT R499_SKO)
histones=(H3K27Ac H3K4me1 H3K4me3 IgG)

for COND in ${conditions[@]}
do
    for HIST in ${histones[@]}
    do
        echo ${COND}_${HIST}

        bsub -o output_logs/merge_conditions/${COND}_${HIST}.out -e output_logs/merge_conditions/${COND}_${HIST}.err -n 16 -M 120000 -R "span[hosts=1] rusage [mem=120000]" bash /home/jingyaq/Minn/data/R499_CnR/scripts/merge_conditions.sh ${COND} ${HIST}
    done
done
```
**generate\_consensus\_count\_matrix.R**

1. Create consensus peaksets for each histone mark
2. Count Tn5 insertions
3. Normalize count matrix
4. Differential region analysis

Consensus peaksets are combined broad peak calls from B16 and Res499 (no STAT1 KO peaks included). Removed peaks that overlap peak with a stronger signal.

```
~/Dropbox/Minn/R499_CnR/scripts/generate_consensus_count_matrix.R
```
