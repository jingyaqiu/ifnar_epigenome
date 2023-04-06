# Center and resize IDR peaks to fixed-width, non-overlapping peaks
# Count number of insertions in consensus IDR peakset

# Use filtered bam files

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

dir.path <- "/home/jingyaq/Minn/data/R499_CnR/data/run_FINAL/"

metadat <- data.frame(
	condition=rep(c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO"), each=6), 
	rep=rep(c("rep1", "rep2"), 12), 
	histone=rep(rep(c("H3K4me1", "H3K4me3", "H3K27Ac"), each=2), 4))
metadat$sampleID <- paste(metadat$condition, metadat$histone, metadat$rep, sep="_")

peaks.use <- "pooled" # IDR, pooled
bam_type <- "all"
combined_peaks <- "WT" # WT, ALL

##################################
### Generate consensus peakset ###
##################################

### Read in broadPeak file (from macs2), remove peaks that overlap peak with stronger signal ###

marks <- c("H3K27Ac", "H3K4me3", "H3K4me1")

broad.peaks_list <- lapply(1:length(marks), function(x) read.table(paste0(dir.path, "macs2/merged/combined_", combined_peaks, "_", marks[x], "_1e-5_peaks.broadPeak"), sep="\t", header=F, stringsAsFactors=F, fill=T))
sapply(broad.peaks_list, function(x) nrow(x))

gr.peaks_list <- lapply(1:3, function(x) {
	print(marks[x])
	gr.consensus <- makeGRangesFromDataFrame(broad.peaks_list[[x]], seqnames.field="V1", start.field="V2", end.field="V3", ignore.strand=T, keep.extra.columns=T)
	gr.consensus <- keepStandardChromosomes(gr.consensus, species="Mus_musculus", pruning.mode="coarse")
	gr.consensus$idx <- 1:length(gr.consensus)

	# Iteratively remove overlapping peaks (retain only most significant peak in overlapping peakset)
	ol <- findOverlaps(gr.consensus, gr.consensus, minoverlap=250)
	overlapping.peaks <- unique(from(ol)[duplicated(from(ol))])
	nonoverlapping.peaks <- gr.consensus$idx[!gr.consensus$idx %in% overlapping.peaks]

	retain.peaks <- c(nonoverlapping.peaks)
	remove.peaks <- c()
	for(i in 1:length(overlapping.peaks)) {
		if(!overlapping.peaks[i] %in% remove.peaks) {
			ol.idx <- ol[from(ol) == overlapping.peaks[i]]
			gr.ol <- gr.consensus[unique(to(ol.idx))]
			retain.peaks <- c(retain.peaks, gr.ol$idx[gr.ol$V9 == max(gr.ol$V9)])
			remove.peaks <- c(remove.peaks, gr.ol$idx[gr.ol$V9 != max(gr.ol$V9)])
		} else {
			print(paste0(overlapping.peaks[i], " already removed"))
		}
	}
	remove.peaks <- unique(remove.peaks)
	gr.consensus <- gr.consensus[-remove.peaks]
	gr.consensus <- unique(gr.consensus)
	write.table(data.frame(gr.consensus)[,1:3], file=paste0(dir.path, "mat/consensus_peaks_", combined_peaks, "_", marks[x], "_", peaks.use, "1.bed"), sep="\t", quote=F, row.names=F, col.names=F)

	# Write out in fasta format8
	fa.consensus <- getSeq(Mmusculus, gr.consensus)
	names(fa.consensus) <- paste0(as.character(gr.consensus@seqnames), ":", start(gr.consensus), "-", end(gr.consensus))
	writeXStringSet(fa.consensus, file=paste0(dir.path, "mat/consensus_peaks_", combined_peaks, "_", marks[x], "_", peaks.use, "1.fasta"))
	return(gr.consensus)
})
sapply(gr.peaks_list, function(x) length(x))

###########################################
### Count insertions in consensus peaks ###
###########################################

# Using countOverlaps, following Corces et al. 2018 methods

### Convert bam files to GRanges ###

# Load in bam files for each sample
bamReads <- paste0(dir.path, "bowtie2/", metadat$sampleID, ".sorted.filtered.dupMarked.blfilt.bam")
bamFiles_list <- lapply(bamReads, function(x) BamFile(x, asMates=T))

gr.bam_list <- lapply(1:length(bamFiles_list), function(x) {
	print(metadat$sampleID[x])
	if(!file.exists(paste0(dir.path, "mat/gr_", metadat$sampleID[x], "_", bam_type, ".rds"))) {
		bam <- scanBam(bamFiles_list[[x]])
		gr <- GRanges(seqnames=Rle(bam[[1]]$rname), ranges=IRanges(bam[[1]]$pos, end=bam[[1]]$pos+bam[[1]]$qwidth-1), strand=Rle(bam[[1]]$strand))
		gr <- keepStandardChromosomes(gr, pruning.mode="coarse") # Remove non-standard chromosomes
		saveRDS(gr, file=paste0(dir.path, "mat/gr_", metadat$sampleID[x], "_", bam_type, ".rds"))
	} else {
		gr <- readRDS(paste0(dir.path, "mat/gr_", metadat$sampleID[x], "_", bam_type, ".rds"))
	}
	return(gr)
})
names(gr.bam_list) <- metadat$sampleID
saveRDS(gr.bam_list, file=paste0(dir.path, "mat/gr.bam_list_", bam_type, ".rds"))

gr.bam_list <- readRDS(paste0(dir.path, "mat/gr.bam_list_", bam_type, ".rds"))

### Count raw insertion counts at each peak (not single-base resolution) ###

marks <- c("H3K27Ac", "H3K4me3", "H3K4me1")

mat.counts_consensus_list <- lapply(1:3, function(x) {
	print(marks[x])
	file.peaks <- paste0(dir.path, "mat/consensus_peaks_", combined_peaks, "_", marks[x], "_", peaks.use, "1.bed")
	file.mat <- paste0(dir.path, "mat/consensus_mat_insertion_counts_", combined_peaks, "_", marks[x], "_", peaks.use, "1.txt")
	gr <- import.bed(file.peaks)
	bam_list <- gr.bam_list[grep(names(gr.bam_list), pattern=marks[x])]
	insertions_in_peaks_list <- lapply(bam_list, function(x) countOverlaps(gr, x))
	mat.counts_consensus <- do.call(cbind, insertions_in_peaks_list)
	colnames(mat.counts_consensus) <- names(bam_list)
	rownames(mat.counts_consensus) <- paste(as.character(gr@seqnames), start(gr), end(gr), sep="_")
	write.table(mat.counts_consensus, file=file.mat, quote=F, col.names=T, row.names=T, sep="\t")
	return(mat.counts_consensus)
})
sapply(mat.counts_consensus_list, function(x) nrow(x))

##############################
### Normalize count matrix ###
##############################

marks <- c("H3K27Ac", "H3K4me3", "H3K4me1")

lapply(1:3, function(x) {

	file.dds <- paste0(dir.path, "mat/consensus_mat_insertion_counts_", combined_peaks, "_", marks[x], "_", peaks.use, "_dds1.rds")
	file.vst <- paste0(dir.path, "mat/consensus_mat_insertion_counts_", combined_peaks, "_", marks[x], "_", peaks.use, "_vst1.txt")
	file.norm <- paste0(dir.path, "mat/consensus_mat_insertion_counts_", combined_peaks, "_", marks[x], "_", peaks.use, "_norm1.txt")

	mat.counts_consensus <- mat.counts_consensus_list[[x]]

	# Set up data for DESeq2
	coldata <- metadat[match(colnames(mat.counts_consensus), metadat$sampleID), ]
	rownames(coldata) <- colnames(mat.counts_consensus)
	coldata$condition <- factor(coldata$condition, levels=c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO"))

	# DESeq2
	dds <- DESeqDataSetFromMatrix(countData=mat.counts_consensus, colData=coldata, design= ~ condition)
	dds <- DESeq(dds)

	# Count data transformation and normalization
	vst <- vst(dds, blind=FALSE)
	mat.atac_vst <- assay(vst)
	print(nrow(mat.atac_vst))

	mat.norm <- counts(dds, normalized=TRUE) # Normalized to sequencing depth

	# Write out
	saveRDS(dds, file=file.dds)
	write.table(mat.atac_vst, file=file.vst, col.names=T, row.names=T, quote=F, sep="\t")
	write.table(mat.norm, file=, col.names=T, row.names=T, quote=F, sep="\t")

	# rld <- rlog(dds, blind=FALSE)
	# mat.atac_rld <- assay(rld)
	# write.table(mat.atac_rld, file=paste0(dir.path, "mat/consensus_mat_insertion_counts_", peaks.use, "_rlog.txt"), col.names=T, row.names=T, quote=F, sep="\t")
})

####################################
### Differential signal analysis ###
####################################

marks <- c("H3K27Ac", "H3K4me3", "H3K4me1")

remove_samples <- c("R499_SKO_H3K4me1_rep2", "R499_SKO_H3K4me3_rep2", "R499_SKO_H3K27Ac_rep1")

for(i in 1:length(marks)) {

	print(marks[i])

	metadat.sub <- metadat[grep(metadat$sampleID, pattern=marks[i]), ]
	metadat.sub <- metadat.sub[!metadat.sub$sampleID %in% remove_samples, ] # Remove bad samples

	file.peaks <- paste0(dir.path, "mat/consensus_peaks_", combined_peaks, "_", marks[i], "_", peaks.use, "1.bed")
	file.dbDE <- paste0(dir.path, "diffbind/consensus_peaks_", combined_peaks, "_", marks[i], "_", peaks.use, "_Diffbind_DE1.rds")

	### DiffBind ###

	gr <- import.bed(file.peaks)
	names(gr) <- 1:length(gr)

	# Sample sheet for individual replicates
	bamReads <- paste0(dir.path, "bowtie2/", metadat.sub$sampleID, ".sorted.filtered.dupMarked.blfilt.bam")
	Peaks <- paste0(dir.path, "macs2/", metadat.sub$sampleID, "_1e-5_peaks.broadPeak")
	PeakCaller <- rep("bed", nrow(metadat.sub))
	Reps <- as.numeric(sapply(strsplit(metadat.sub$sampleID, split="_"), function(x) substr(x[[4]], 4, 4)))
	ss <- data.frame(Sample=metadat.sub$sampleID, Condition=metadat.sub$condition, Replicate=Reps, bamReads, Peaks, PeakCaller)

	# Load in individual sample sheet
	db <- dba(sampleSheet=ss)

	# Generate normalized counts for each peak in consensus set
	db.counts <- dba.count(db, peaks=gr, score=DBA_SCORE_TMM_READS_FULL, bParallel=TRUE)

	# Differential peak analysis
	contrast.names <- c("R499_WT_v_B16_WT", "B16_WT_v_B16_SKO", "R499_WT_v_R499_SKO", "R499_SKO_v_B16_SKO")

	db.DE_list <- lapply(1:length(contrast.names), function(x) {
		print(contrast.names[x])

		cond1 <- sapply(strsplit(contrast.names[x], split="_"), function(x) paste(x[1:2], collapse="_"))
		cond2 <- sapply(strsplit(contrast.names[x], split="_"), function(x) paste(x[4:5], collapse="_"))

		# Using DiffBind size factors based on full library size
		db.DE <- dba.contrast(db.counts, group1=db.counts$masks[[cond1]], group2=db.counts$masks[[cond2]], name1=cond1, name2=cond2)
		db.DE <- dba.analyze(db.DE, bFullLibrarySize=TRUE, method=c(DBA_EDGER), bSubControl=TRUE, bTagwise=TRUE, bParallel=TRUE)
		db.DE$config$AnalysisMethod <- factor("edgeR")
		return(db.DE)
	})
	saveRDS(db.DE_list, file=file.dbDE)
}
