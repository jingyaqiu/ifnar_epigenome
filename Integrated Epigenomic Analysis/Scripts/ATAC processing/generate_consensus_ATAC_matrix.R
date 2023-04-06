# Center and resize IDR peaks to fixed-width, non-overlapping peaks
# Count number of insertions in consensus IDR peakset

# Use tn5 shifted bam files

# Run on pmacs with R/3.5.1

rm(list=ls())

library(rtracklayer) # version 1.40.6
library(chromVAR) # version 1.4.1
library(GenomicRanges) # version 1.34.0
library(Rsamtools) # version 1.32.3
library(DESeq2) # version 1.22.2
library(DiffBind) # version 2.10.0
library(Biostrings) # version 2.50.2
library(BSgenome.Mmusculus.UCSC.mm10) # version 1.4.0

dir.path <- "/home/jingyaq/Minn/data/ifnar_epigenome/ENCODE_pipeline_data/ATAC/"

conditions <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
sampleID <- paste0(rep(conditions, each=5), "_", 1:5)

peaks.use <- "IDR" # IDR, pooled, overlap
bam_type <- "all" # all, downsample

##################################
### Generate consensus peakset ###
##################################

# Read in narrowPeaks file (from macs2), resize to fixed width around summit, remove peaks that overlap peak with stronger signal
if(peaks.use == "IDR") {
	gr.consensus <- readNarrowpeaks(paste0(dir.path, "peaks/idr/combined_ppr.IDR0.1.filt.narrowPeak.gz"), width=750, non_overlapping=TRUE)
} else if(peaks.use == "pooled") {
	gr.consensus <- readNarrowpeaks(paste0(dir.path, "peaks/pooled_rep/combined_R1.trim_pooled.PE2SE.nodup.tn5_pooled.pf.pval0.01.300K.filt.narrowPeak.gz"), width=750, non_overlapping=TRUE)
}

# Write out consensus peakset
write.table(data.frame(gr.consensus)[,1:3], file=paste0(dir.path, "mat/consensus_peaks_", peaks.use, "_fixed750.bed"), sep="\t", quote=F, row.names=F, col.names=F)

# Write out in fasta format
fa.consensus <- getSeq(Mmusculus, gr.consensus)
names(fa.consensus) <- paste0(as.character(gr.consensus@seqnames), ":", start(gr.consensus), "-", end(gr.consensus))
writeXStringSet(fa.consensus, file=paste0(dir.path, "mat/consensus_peaks_", peaks.use, "_fixed750.fasta"))

############################
### Count Tn5 insertions ###
############################

# Using countOverlaps, following Corces et al. 2018 methods

### Convert bam files to GRanges ###

# Load in bam files for each sample
if(bam_type == "all") {
	bamReads <- paste0(dir.path, "bam/filtered/", sampleID, ".R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam")
} else if(bam_type == "downsample") {
	bamReads <- paste0(dir.path, "bam/filtered/downsample/", sampleID, "_downsample_10M.bam")
}

bamFiles_list <- lapply(bamReads, function(x) BamFile(x))
bam_list <- lapply(bamFiles_list, function(x) scanBam(x))

gr.bam_list <- lapply(bam_list, function(x) GRanges(seqnames=Rle(x[[1]]$rname), ranges=IRanges(x[[1]]$pos, end=x[[1]]$pos+x[[1]]$qwidth-1), strand=Rle(x[[1]]$strand)))
names(gr.bam_list) <- sampleID
saveRDS(gr.bam_list, file=paste0(dir.path, "mat/gr.bam_list_", bam_type, ".rds"))

### Count raw insertion counts at each peak (not single-base resolution) ###

insertions_in_peaks_list <- lapply(gr.bam_list, function(x) countOverlaps(gr.consensus, x))
mat.counts_consensus <- do.call(cbind, insertions_in_peaks_list)
colnames(mat.counts_consensus) <- sampleID
rownames(mat.counts_consensus) <- paste(as.character(gr.consensus@seqnames), start(gr.consensus), end(gr.consensus), sep="_")
write.table(mat.counts_consensus, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_", peaks.use, "_", bam_type, ".txt"), quote=F, col.names=T, row.names=T, sep="\t")

##############################
### Normalize count matrix ###
##############################

mat.counts_consensus <- read.table(paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_", peaks.use, "_", bam_type, ".txt"), header=T, sep="\t", stringsAsFactors=F)

# Set up data for DESeq2
coldata <- data.frame(sampleID=sampleID, condition=rep(conditions, each=5))
rownames(coldata) <- colnames(mat.counts_consensus)

# DESeq2
dds <- DESeqDataSetFromMatrix(countData=mat.counts_consensus, colData=coldata, design= ~ condition)
dds <- DESeq(dds)

# Count data transformation and normalization
vst <- vst(dds, blind=FALSE)
mat.atac_vst <- assay(vst)

rld <- rlog(dds, blind=FALSE)
mat.atac_rld <- assay(rld)

# Write out
saveRDS(dds, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_", peaks.use, "_dds_", bam_type, ".rds"))
write.table(mat.atac_vst, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_", peaks.use, "_vst_", bam_type, ".txt"), col.names=T, row.names=T, quote=F, sep="\t")
write.table(mat.atac_rld, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_", peaks.use, "_rlog_", bam_type, ".txt"), col.names=T, row.names=T, quote=F, sep="\t")

###########################################
### Differential accessibility analysis ###
###########################################

gr.consensus <- import.bed(paste0(dir.path, "mat/consensus_peaks_", peaks.use, "_fixed750.bed"))

### Using DiffBind size factors based on full library size ###

sampleID <- paste0(rep(conditions, each=5), "_", 1:5)
Condition <- rep(conditions, each=5)
Replicate <- rep(1:5, 4)

# Sample sheet for individual replicates
if(bam_type == "all") {
	bamReads <- paste0(dir.path, "bam/filtered/", sampleID, ".R1.trim_pooled.PE2SE.nodup_noChrM_tn5shift_sorted.bam")
} else if(bam_type == "downsample") {
	bamReads <- paste0(dir.path, "bam/filtered/downsample/", sampleID, "_downsample_10M.bam")
}

Peaks <- list.files(paste0(dir.path, "peaks"), pattern=".narrowPeak.gz$", full.names=TRUE)
PeakCaller <- rep("narrow", length(sampleID))
ss <- data.frame(sampleID, Condition, Replicate, bamReads, Peaks, PeakCaller)

remove_samples <- c("B16_SKO_3")
ss <- ss[!ss$sampleID %in% remove_samples, ]

# Load in individual sample sheet
db <- dba(sampleSheet=ss)

# Generate normalized counts for each peak in consensus set
db.counts <- dba.count(db, peaks=gr.consensus, score=DBA_SCORE_TMM_READS_FULL, bParallel=TRUE)

# Differential peak analysis
db.DE <- dba.contrast(db.counts, categories=DBA_CONDITION)
db.DE <- dba.analyze(db.DE, bFullLibrarySize=TRUE, method=c(DBA_DESEQ2), bSubControl=TRUE, bParallel=TRUE)
saveRDS(db.DE, file=paste0(dir.path, "diffbind/consensus_peaks_", peaks.use, "_Diffbind_DE_", bam_type, ".rds"))

db.DE <- dba.contrast(db.counts, categories=DBA_CONDITION)
db.DE <- dba.analyze(db.DE, bFullLibrarySize=TRUE, method=c(DBA_EDGER), bSubControl=FALSE, bParallel=TRUE)
db.DE$config$AnalysisMethod <- factor("edgeR")
saveRDS(db.DE, file=paste0(dir.path, "diffbind/consensus_peaks_", peaks.use, "_Diffbind_DE_", bam_type, ".rds"))

###################################################
### Write out differentially accessible regions ###
###################################################

source("~/Minn/data/ifnar_epigenome/scripts/Integrated Epigenomic Analysis/Integrated Epigenomic Analysis Functions.R")

### Original ATAC ###

db.DE <- readRDS("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_original_IDR_Diffbind_DE_all.rds")
contrasts.idx <- c(1,6,2,5)
contrasts.name <- c(
	"B16_cas_vs_B16_SKO", "R499_cas_vs_R499_SKO", 
	"B16_cas_vs_R499_cas", "B16_SKO_vs_R499_SKO")
lapply(1:length(contrasts.idx), function(x) {
	gr <- dba.report(db.DE, contrast = contrasts.idx[x], th = 1)
	saveRDS(gr, file=paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_original/", contrasts.name[x], "_DiffBind.rds"))
})

### Revisions ATAC ###

db.DE <- readRDS("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_revisions_IDR_Diffbind_DE_all.rds")
db.DE$config$AnalysisMethod <- factor("edgeR")
contrasts.idx <- c(2,7,1,10,11,12)
contrasts.name <- c(
	"B16_WT_vs_R499_WT", "B16_SKO_vs_R499_SKO", 
	"B16_WT_vs_B16_SKO", 
	"R499_WT_vs_R499_SKO", "R499_WT_vs_R499_IRF3KO", "R499_WT_vs_R499_DKO")
lapply(1:length(contrasts.idx), function(x) {
	gr <- dba.report(db.DE, contrast = contrasts.idx[x], th = 1)
	saveRDS(gr, file=paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_revisions/", contrasts.name[x], "_DiffBind.rds"))
})

### H3K4me1 ###

db.DE <- readRDS(paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K4me1_pooled_WT_Diffbind_DE.rds"))
contrasts.idx <- c(1,2,3,4)
contrasts.name <- c(
	"R499_WT_vs_B16_WT", "B16_WT_vs_B16_SKO", "R499_WT_vs_R499_SKO", "R499_SKO_vs_B16_SKO")
lapply(1:length(contrasts.idx), function(x) {
	db.sub <- db.DE[[contrasts.idx[x]]]
	gr <- dba.report(db.sub, contrast = 1, th = 1)
	saveRDS(gr, file=paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K4me1/", contrasts.name[x], "_DiffBind.rds"))
})

### H3K27ac ###

db.DE <- readRDS(paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K27Ac_pooled_WT_Diffbind_DE.rds"))
contrasts.idx <- c(1,2,3,4)
contrasts.name <- c(
	"R499_WT_vs_B16_WT", "B16_WT_vs_B16_SKO", "R499_WT_vs_R499_SKO", "R499_SKO_vs_B16_SKO")
lapply(1:length(contrasts.idx), function(x) {
	db.sub <- db.DE[[contrasts.idx[x]]]
	gr <- dba.report(db.sub, contrast = 1, th = 1)
	saveRDS(gr, file=paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K27Ac/", contrasts.name[x], "_DiffBind.rds"))
})

### H3K4me3 ###

db.DE <- readRDS(paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K4me3_pooled_WT_Diffbind_DE.rds"))
contrasts.idx <- c(1,2,3,4)
contrasts.name <- c(
	"R499_WT_vs_B16_WT", "B16_WT_vs_B16_SKO", "R499_WT_vs_R499_SKO", "R499_SKO_vs_B16_SKO")
lapply(1:length(contrasts.idx), function(x) {
	db.sub <- db.DE[[contrasts.idx[x]]]
	gr <- dba.report(db.sub, contrast = 1, th = 1)
	saveRDS(gr, file=paste0("~/Minn/data/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/H3K4me3/", contrasts.name[x], "_DiffBind.rds"))
})
