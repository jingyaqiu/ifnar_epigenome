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

# Revisions ATAC
dir.path <- "/home/jingyaq/Minn/data/revisions_NatCancer/data/processed/"
conditions <- c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO")
sampleID <- c(paste0("B16_WT-", c(1,2,4)), paste0("B16_SKO-", c(1,4)), paste0("R499_WT-", 1:4), paste0("R499_SKO-", c(1,2,4)), paste0("R499_IRF3KO-", 1:4), paste0("R499_DKO-", c(1,3,4)))

peaks.use <- "IDR" # IDR, pooled, overlap
bam_type <- "downsample" # all, downsample

##################################
### Generate consensus peakset ###
##################################

# Read in combined IDR peak file (from macs2), resize to fixed width around summit, remove peaks that overlap peak with stronger signal
gr.consensus <- readNarrowpeaks(paste0(dir.path, "peaks/idr/combined_idr.optimal_peak.narrowPeak"), width=750, non_overlapping=TRUE)

# Write out consensus peakset
write.table(data.frame(gr.consensus)[,1:3], file=paste0(dir.path, "mat/consensus_peaks_IDR_fixed750.bed"), sep="\t", quote=F, row.names=F, col.names=F)

# Write out in fasta format
fa.consensus <- getSeq(Mmusculus, gr.consensus)
names(fa.consensus) <- paste0(as.character(gr.consensus@seqnames), ":", start(gr.consensus), "-", end(gr.consensus))
writeXStringSet(fa.consensus, file=paste0(dir.path, "mat/consensus_peaks_IDR_fixed750.fasta"))

# gr.consensus <- import.bed(paste0(dir.path, "mat/consensus_peaks_IDR_fixed750.bed"))

############################
### Count Tn5 insertions ###
############################

# Using countOverlaps, following Corces et al. 2018 methods

### Convert bam files to GRanges ###

# Load in bam files for each sample
if(bam_type == "all") {
	bamReads <- paste0(dir.path, "bam/", sampleID, ".trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam")
} else if(bam_type == "downsample") {
	bamReads <- paste0(dir.path, "bam/downsample/", sampleID, "_downsample_51M.bam")
}

bamFiles_list <- lapply(bamReads, function(x) BamFile(x, asMates=T))

gr.bam_list <- lapply(1:length(bamFiles_list), function(x) {
	print(sampleID[x])
	bam <- scanBam(bamFiles_list[[x]])
	gr <- GRanges(seqnames=Rle(bam[[1]]$rname), ranges=IRanges(bam[[1]]$pos, end=bam[[1]]$pos+bam[[1]]$qwidth-1), strand=Rle(bam[[1]]$strand))
	gr <- keepStandardChromosomes(gr, pruning.mode="coarse") # Remove non-standard chromosomes
	saveRDS(gr, file=paste0(dir.path, "mat/gr_", sampleID[x], "_", bam_type, ".rds"))
	return(gr)
})
names(gr.bam_list) <- sampleID
saveRDS(gr.bam_list, file=paste0(dir.path, "mat/gr.bam_list_", bam_type, ".rds"))

### Count raw insertion counts at each peak (not single-base resolution) ###

insertions_in_peaks_list <- lapply(gr.bam_list, function(x) countOverlaps(gr.consensus, x))
mat.counts_consensus <- do.call(cbind, insertions_in_peaks_list)
colnames(mat.counts_consensus) <- sampleID
rownames(mat.counts_consensus) <- paste(as.character(gr.consensus@seqnames), start(gr.consensus), end(gr.consensus), sep="_")
write.table(mat.counts_consensus, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_IDR_", bam_type, ".txt"), quote=F, col.names=T, row.names=T, sep="\t")

##############################
### Normalize count matrix ###
##############################

mat.counts_consensus <- read.table(paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_IDR_", bam_type, ".txt"), header=T, sep="\t", stringsAsFactors=F)

# # Remove outlier sample
# remove_sample <- "B16_WT.3"
# mat.counts_consensus <- mat.counts_consensus[ ,!colnames(mat.counts_consensus) %in% remove_sample]

# Set up data for DESeq2
coldata <- data.frame(sampleID=colnames(mat.counts_consensus), condition=sapply(strsplit(colnames(mat.counts_consensus), split="_"), function(x) paste(x[1:2], collapse="_")))
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
saveRDS(dds, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_IDR_dds_", bam_type, ".rds"))
write.table(mat.atac_vst, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_IDR_vst_", bam_type, ".txt"), col.names=T, row.names=T, quote=F, sep="\t")
write.table(mat.atac_rld, file=paste0(dir.path, "mat/consensus_mat_tn5_insertion_counts_IDR_rlog_", bam_type, ".txt"), col.names=T, row.names=T, quote=F, sep="\t")

###########################################
### Differential accessibility analysis ###
###########################################

gr.consensus <- import.bed(paste0(dir.path, "mat/consensus_peaks_IDR_fixed750.bed"))

### Using DiffBind size factors based on full library size ###

Condition <- sapply(strsplit(sampleID, split="_"), function(x) paste(x[1:2], collapse="_"))
Replicate <- sapply(strsplit(sampleID, split="_"), function(x) x[[3]])

# Sample sheet for individual replicates
if(bam_type == "all") {
	bamReads <- paste0(dir.path, "bam/", sampleID, ".trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam")
} else if(bam_type == "downsample") {
	bamReads <- paste0(dir.path, "bam/downsample/", sampleID, "_downsample_51M.bam")
}

Peaks <- sapply(1:length(sampleID), function(x) list.files(paste0(dir.path, "peaks/"), pattern=sampleID[x], full.names=T))
PeakCaller <- rep("narrow", length(sampleID))
ss <- data.frame(sampleID, Condition, Replicate, bamReads, Peaks, PeakCaller)

# Load in individual sample sheet
db <- dba(sampleSheet=ss)

# Generate normalized counts for each peak in consensus set
db.counts <- dba.count(db, peaks=gr.consensus, score=DBA_SCORE_TMM_READS_FULL, bParallel=TRUE)
saveRDS(db.counts, file=paste0(dir.path, "diffbind/consensus_peaks_IDR_Diffbind_counts_", bam_type, ".rds"))

# db.counts <- readRDS(paste0(dir.path, "diffbind/consensus_peaks_IDR_Diffbind_counts_", bam_type, ".rds"))

# Differential peak analysis
db.DE <- dba.contrast(db.counts, categories=DBA_CONDITION, minMembers=2)
db.DE <- dba.analyze(db.DE, bFullLibrarySize=TRUE, method=c(DBA_EDGER), bSubControl=TRUE, bTagwise=TRUE, bParallel=TRUE)
saveRDS(db.DE, file=paste0(dir.path, "diffbind/consensus_peaks_IDR_Diffbind_DE_", bam_type, ".rds"))

#####################
### Identify DARs ###
#####################

db.DE <- readRDS(paste0(dir.path, "diffbind/consensus_peaks_IDR_Diffbind_DE_", bam_type, ".rds"))
db.DE$config$AnalysisMethod <- factor("edgeR")

contrasts <- c(2,1,10,7,11,12,13,14,15)
contrast.names <- c("B16_v_R499", "B16_v_B16_SKO", "R499_v_R499_SKO", "B16_SKO_v_R499_SKO", "R499_v_R499_IRF3KO", "R499_v_R499_DKO", "R499_SKO_v_R499_IRF3KO", "R499_SKO_v_R499_DKO", "R499_IRF3KO_v_R499_DKO")

# Write out DARs

make_HOMER_bed <- function(gr, file=NULL) {
	bed <- data.frame(gr)[,1:3]
	bed$ID <- paste(bed[,1], bed[,2], bed[,3], sep="_")
	bed$col5 <- NA
	bed$strand <- as.character(strand(gr))
	if(!is.null(file)) write.table(bed, file=file, sep="\t", row.names=F, col.names=F, quote=F)
	return(bed)
}

write_fasta <- function(gr, file=NULL) {
	fa <- getSeq(Mmusculus, gr)
	names(fa) <- paste("chromosome:GRCm38:", as.character(gr@seqnames), ":", start(gr), ":", end(gr), sep="")
	writeXStringSet(fa, filepath=file, append=FALSE)
}

sig <- 0.05
# max_peaks <- 8000
DARs_list <- lapply(1:length(contrasts), function(x) {
	gr <- dba.report(db.DE, contrast=contrasts[x], bUsePval=F, th=1)

	gr.sig1 <- gr[gr$FDR < sig & gr$Fold > 0]
	gr.sig2 <- gr[gr$FDR < sig & gr$Fold < 0]
	print(paste0(length(gr.sig1), ", ", length(gr.sig2)))

	# if(length(gr.sig1) > max_peaks) gr.sig1 <- gr.sig1[1:max_peaks]
	# if(length(gr.sig2) > max_peaks) gr.sig2 <- gr.sig2[1:max_peaks]

	make_HOMER_bed(gr=gr.sig1, file=paste0(dir.path, "diffbind/DAR/DAR_", contrast.names[x], "_cond1.bed"))
	make_HOMER_bed(gr=gr.sig2, file=paste0(dir.path, "diffbind/DAR/DAR_", contrast.names[x], "_cond2.bed"))
	write_fasta(gr=gr.sig1, file=paste0(dir.path, "diffbind/DAR/DAR_", contrast.names[x], "_cond1.fasta"))
	write_fasta(gr=gr.sig2, file=paste0(dir.path, "diffbind/DAR/DAR_", contrast.names[x], "_cond2.fasta"))
	
	return(list(gr.sig1, gr.sig2))
})
names(DARs_list) <- contrast.names
sapply(DARs_list, function(x) c(length(x[[1]]), length(x[[2]])))

####################################################
### Identify constitutive and newly opened peaks ###
####################################################

sig <- 0.1
fc <- 1

# Newly opened (FDR < sig, FC > fc)
newlyOpened_peaks <- lapply(1:length(contrasts), function(x) {
	gr <- dba.report(db.DE, contrast=contrasts[x], bUsePval=F, th=1)
	if(contrast.names[x] %in% c("B16_v_R499", "B16_SKO_v_R499_SKO", "R499_SKO_v_R499_IRF3KO")) {
		gr.filt <- gr[gr$FDR < sig & gr$Fold < -fc]
	} else if(contrast.names[x] %in% c("B16_v_B16_SKO", "R499_v_R499_SKO", "R499_v_R499_IRF3KO", "R499_v_R499_DKO", "R499_SKO_v_R499_DKO", "R499_IRF3KO_v_R499_DKO")) {
		gr.filt <- gr[gr$FDR < sig & gr$Fold > fc]
	}
	return(gr.filt)
})
sapply(newlyOpened_peaks, function(x) length(x))

# Constitutive (FDR > sig, top 5000)
constitutive_peaks <- lapply(1:length(contrasts), function(x) {
	gr <- dba.report(db.DE, contrast=contrasts[x], bUsePval=F, th=1)
	gr.filt <- gr[(length(gr)-4999):length(gr)]
	return(gr.filt)
})
sapply(constitutive_peaks, function(x) length(x))

# Write out
lapply(1:length(contrasts), function(x) {
	make_HOMER_bed(newlyOpened_peaks[[x]], file=paste0(dir.path, "diffbind/DAR/newlyOpened_peaks_", contrast.names[x], ".bed"))
	make_HOMER_bed(constitutive_peaks[[x]], file=paste0(dir.path, "diffbind/DAR/constitutive_peaks_", contrast.names[x], ".bed"))
	write_fasta(gr=newlyOpened_peaks[[x]], file=paste0(dir.path, "diffbind/DAR/newlyOpened_peaks_", contrast.names[x], ".fasta"))
	write_fasta(gr=constitutive_peaks[[x]], file=paste0(dir.path, "diffbind/DAR/constitutive_peaks_", contrast.names[x], ".fasta"))
})

#############
### GREAT ###
#############

library(rGREAT)

great_list <- lapply(1:length(contrasts), function(x) {
	gr1 <- DARs_list[[x]][[1]]
	gr2 <- DARs_list[[x]][[2]]

	great1 <- submitGreatJob(gr1, species="mm10")
	tb1 <- getEnrichmentTables(great1, ontology=c("GO Biological Process"))[[1]]

	great2 <- submitGreatJob(gr2, species="mm10")
	tb2 <- getEnrichmentTables(great2, ontology=c("GO Biological Process"))[[1]]
	return(list(tb1[tb1$Binom_Adjp_BH < 0.05, ], tb2[tb2$Binom_Adjp_BH < 0.05, ]))
})
names(great_list) <- contrast.names
saveRDS(great_list, file=paste0(dir.path, "diffbind/DAR/great/DARs_GREAT_GO_BP.rds"))

# GREAT on newly opened peaks
great_newlyOpened <- lapply(1:length(contrasts), function(x) {
	great <- submitGreatJob(newlyOpened_peaks[[x]], species="mm10")
	tb <- getEnrichmentTables(great, ontology=c("GO Biological Process"))[[1]]
	return(tb[tb$Binom_Adjp_BH < 0.05, ])
})
names(great_newlyOpened) <- contrast.names
saveRDS(great_newlyOpened, file=paste0(dir.path, "diffbind/DAR/great/newlyOpened_GREAT_GO_BP.rds"))

# ##############################################

# # R499 > B16 - all peaks
# db.DE <- readRDS(paste0(dir.path, "diffbind/consensus_peaks_", peaks.use, "_Diffbind_DE.rds"))
# gr <- dba.report(db.DE, contrast=2, bUsePval=F, th=1)
# gr.sig <- gr[gr$FDR < 0.05 & gr$Fold < 0]
# make_HOMER_bed(gr.sig, file=paste0(dir.path, "diffbind/DAR/DAR_B16_v_R499_group2_FULL.bed"))
# fa.sig <- getSeq(Mmusculus, gr.sig)
# names(fa.sig) <- paste("chromosome:GRCm38:", as.character(gr.sig@seqnames), ":", start(gr.sig), ":", end(gr.sig), sep="")
# writeXStringSet(fa.sig, filepath=paste0(dir.path, "diffbind/DAR/DAR_B16_v_R499_group2_FULL.fasta"), append=FALSE)

# ########################
# ### csaw DA analysis ###
# ########################

# library(csaw) # version 1.16.1
# library(edgeR) # version 3.24.3

# # Load in bam files
# if(bam_type == "all") {
# 	bam.files <- paste0(dir.path, "bam/", sampleID, ".trim.srt.nodup.no_chrM_MT_tn5shift_sorted.bam")
# } else if(bam_type == "downsample") {
# 	bam.files <- paste0(dir.path, "bam/downsample/", sampleID, "_downsample_65M.bam")
# }

# param <- readParam(minq=20, pe="both", dedup=TRUE, BPPARAM=MulticoreParam(workers=16))
# data <- windowCounts(bam.files, ext=110, width=10, param=param)

# # Filter out uninteresting regions
# keep <- aveLogCPM(asDGEList(data)) >= -1
# data <- data[keep,]

# # Calculate normalization factors
# binned <- windowCounts(bam.files, bin=TRUE, width=10000, param=param)
# data <- normFactors(binned, se.out=data)

# # Identify DA windows

# sampleID <- c(paste0("B16_WT-", 1:4), paste0("B16_SKO-", c(1,4)), paste0("R499_WT-", 1:4), paste0("R499_SKO-", c(1,2,4)), paste0("R499_IRF3KO-", 1:4), paste0("R499_DKO-", c(1,3,4)))
# condition <- sapply(strsplit(sampleID, split="-"), function(x) x[[1]])

# design <- model.matrix(~factor(condition, levels=c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO")))
# # colnames(design) <- c("intercept", "condition")

# y <- asDGEList(data)
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design, robust=TRUE)
# results <- glmQLFTest(fit)

# # Correct for multiple testing
# merged <- mergeResults(data, results$table, tol=1000L)

####################

s1 <- readNarrowpeaks(paste0("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/peaks/idr/R499_WT_idr.optimal_peak.narrowPeak"))
s2 <- readNarrowpeaks(paste0("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/peaks/idr/R499_cas_ppr.IDR0.1.filt.narrowPeak.gz"))
s3 <- readNarrowpeaks(paste0("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/peaks/idr/B16_cas_ppr.IDR0.1.filt.narrowPeak.gz"))
s4 <- readNarrowpeaks(paste0("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/peaks/idr/R499_IRF3KO_idr.optimal_peak.narrowPeak"))

db.DE <- readRDS("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/diffbind/consensus_peaks_IDR_Diffbind_DE_all.rds")
db.DE$config$AnalysisMethod <- factor("edgeR")

dar1 <- import.bed("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/diffbind/DAR/DAR_B16_v_R499_cond2.bed")
# dar1 <- dba.report(db.DE, contrast=4, th=0.05)
dar2 <- import.bed("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/diffbind/DAR/DAR_B16_v_R499_group2_FULL.bed")

ol <- findOverlaps(dar1, dar2)
length(unique(from(ol))) / length(dar1)
length(unique(to(ol))) / length(dar2)


gr <- dba.report(db.DE, contrast=4, bUsePval=F, th=1)
gr.sig <- gr[gr$FDR < 0.05 & gr$Fold < 0]
dar1 <- gr.sig

dar1 <- import.bed("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/diffbind/DAR/DAR_B16_v_B16_SKO_cond1.bed")
dar2 <- import.bed("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/diffbind/DAR/DAR_B16_v_B16_SKO_group1.bed")
ol <- findOverlaps(dar1, dar2)
length(unique(from(ol))) / length(dar1)
length(unique(to(ol))) / length(dar2)

#########################
### Load in ATAC data ###
#########################

# Convert peak IDs (chr_start_end) to GRanges
peakids2GRanges <- function(peakids, delim="_") {
	df <- data.frame(chr=sapply(strsplit(peakids, split=delim), function(x) x[1]), start=as.numeric(sapply(strsplit(peakids, split=delim), function(x) x[2])), end=as.numeric(sapply(strsplit(peakids, split=delim), function(x) x[3])))
	gr <- makeGRangesFromDataFrame(df, starts.in.df.are.0based=FALSE)
	return(gr)
}

# ATAC1 data
adat1 <- read.table("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/mat/consensus_mat_tn5_insertion_counts_IDR_vst_downsample.txt", sep="\t", stringsAsFactors=F, header=T)
# gr.atac_IDR <- import.bed(paste0(dir.path, "mat/consensus_peaks_IDR_fixed750.bed"))
gr.atac1 <- peakids2GRanges(rownames(adat1), delim="_")
gr.atac1$idx <- 1:length(gr.atac1)

# ATAC2
adat2 <- read.table("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/mat/consensus_mat_tn5_insertion_counts_IDR_rlog.txt", sep="\t", stringsAsFactors=F, header=T)
gr.atac2 <- peakids2GRanges(peakids=rownames(adat2), delim="_")

# DARs
dar1 <- import.bed("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/diffbind/DAR/DAR_B16_v_R499_cond2.bed")
dar1 <- import.bed("/Users/jingyaq/Dropbox/Minn/revisions_NatCancer/data/processed/diffbind/DAR/DAR_R499_v_R499_IRF3KO_cond2.bed")
dar1 <- dba.report(db.DE, contrast=4, bUsePval=F, th=1)
dar1 <- dar1[dar1$FDR < 0.05 & dar1$Fold < 0]
dar1 <- dar1[dar1$FDR < 0.05 & dar1$Fold > 0]
dar2 <- import.bed("/Users/jingyaq/Dropbox/Minn/ifnar_epigenome/data/ATAC/diffbind/DAR/DAR_B16_v_R499_group2_FULL.bed")

# Test revisions DARs
idx1 <- unique(from(findOverlaps(gr.atac1, dar1)))
idx2 <- unique(from(findOverlaps(gr.atac2, dar1)))

mat <- t(scale(t(adat1[idx1, ])))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
pheatmap(mat, show_rownames = F, angle_col = 45, cellwidth = 20, cluster_cols = F, gaps_col = c(3,5,9,12,16), color=rev(brewer.pal(9, "RdBu")))

mat <- t(scale(t(adat2[idx2, ])))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
pheatmap(mat, show_rownames = F, angle_col = 45, cellwidth = 20, cluster_cols = F, gaps_col = c(5,10,15), color=rev(brewer.pal(9, "RdBu")))

# Test original DARs
idx1 <- unique(from(findOverlaps(gr.atac1, dar2)))
idx2 <- unique(from(findOverlaps(gr.atac2, dar2)))

mat <- t(scale(t(adat1[idx1, -c(4,12)])))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
pheatmap(mat, show_rownames = F, angle_col = 45, cellwidth = 20, cluster_cols = F, gaps_col = c(3,4,8,10,14), color=rev(brewer.pal(9, "RdBu")))

mat <- t(scale(t(adat2[idx2, ])))
mat[mat > 2] <- 2
mat[mat < -2] <- -2
pheatmap(mat, show_rownames = F, angle_col = 45, cellwidth = 20, cluster_cols = F, gaps_col = c(5,10,15), color=rev(brewer.pal(9, "RdBu")))

