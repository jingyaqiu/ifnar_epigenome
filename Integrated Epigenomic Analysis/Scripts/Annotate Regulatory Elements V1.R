
### Annotate regulatory elements through integration of multi-omic epigenetic data ###

rm(list=ls())

library(rtracklayer)
library(GenomicRanges)

###############################
### Load in reference files ###
###############################

# Gene annotations
bm <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/mouse_biomart_annotations_2-20-20.txt", sep="\t", header=T, stringsAsFactors=F)
remove_idx <- grep(bm$chromosome_name, pattern="MT|GL|JH|X|Y")
bm <- bm[-remove_idx, ]
bm$strand <- ifelse(bm$strand == 1, "+", "-")
gr.bm <- makeGRangesFromDataFrame(bm, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")

# TSS for protein-coding genes
gr.tss <- makeGRangesFromDataFrame(data.frame(
	chr = as.character(seqnames(gr.bm)), 
	start = gr.bm$tss, 
	end = gr.bm$tss, 
	gene = gr.bm$mgi_symbol, 
	strand = as.character(strand(gr.bm))), keep.extra.columns = T)

##############################
### Load in processed data ###
##############################

marks <- c("H3K4me3", "H3K27Ac", "H3K4me1")
assays <- c("H3K4me3", "H3K27Ac", "H3K4me1", "ATAC")

### CONSENSUS PEAK SETS ###

gr_list <- lapply(1:length(assays), function(x) {
	if(assays[x] == "ATAC") {
		gr.original <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Consensus Peaksets/ATAC_original_consensus_peaks_IDR_fixed750.bed")
		gr.revisions <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Consensus Peaksets/ATAC_revisions_consensus_peaks_IDR_fixed750.bed")
		gr <- GenomicRanges::reduce(c(gr.original, gr.revisions))
	} else {
		gr <- import.bed(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Consensus Peaksets/", assays[x], "_consensus_peaks_WT_pooled.bed"))
	}
	return(gr)
})
names(gr_list) <- assays

### DIFFERENTIAL FEATURES (DIFFBIND) ###

db_list <- lapply(1:length(assays), function(x) {
	if(assays[x] == "ATAC") {
		db.DE <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_revisions_IDR_Diffbind_DE_all.rds")
		db.DE$config$AnalysisMethod <- factor("edgeR")
	} else {
		db.DE <- readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/", assays[x], "_pooled_WT_Diffbind_DE.rds"))
	}
	return(db.DE)
})
names(db_list) <- assays

#####################################
### Annotate histone combinations ###
#####################################

gr.marks_list <- lapply(1:length(marks), function(x) {
	mark <- marks[x]
	gr <- gr_list[[x]]
	idx.other <- setdiff(1:3, x)

	# Annotate if peak overlaps with histone or ATAC signal, or combination
	ol1 <- findOverlaps(gr, gr_list[[idx.other[1]]], minoverlap=100)
	ol2 <- findOverlaps(gr, gr_list[[idx.other[2]]], minoverlap=100)
	idx.ALL <- intersect(unique(from(ol1)), unique(from(ol2)))
	idx1 <- setdiff(from(ol1), from(ol2))
	idx2 <- setdiff(from(ol2), from(ol1))
	idx.only <- c(1:length(gr))[!1:length(gr) %in% unique(c(from(ol1), from(ol2)))]
	gr$anno <- rep(NA, length(gr))
	gr$anno[idx.ALL] <- "ALL"
	gr$anno[idx1] <- paste0(mark, "_", marks[idx.other[1]])
	gr$anno[idx2] <- paste0(mark, "_", marks[idx.other[2]])
	gr$anno[idx.only] <- paste0(mark, "_only")
	ol.atac <- findOverlaps(gr, gr_list[["ATAC"]], minoverlap=100)
	gr$ATAC <- ifelse(1:length(gr) %in% unique(from(ol.atac)), "open", "closed")

	# Does peak have differential histone or ATAC signal? #
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K4me3"]][[1]], name="H3K4me3_R499_WT>B16_WT", contrast=1, fc=0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K4me3"]][[1]], name="H3K4me3_B16_WT>R499_WT", contrast=1, fc=-0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K27Ac"]][[1]], name="H3K27Ac_R499_WT>B16_WT", contrast=1, fc=0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K27Ac"]][[1]], name="H3K27Ac_B16_WT>R499_WT", contrast=1, fc=-0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K4me1"]][[1]], name="H3K4me1_R499_WT>B16_WT", contrast=1, fc=0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["H3K4me1"]][[1]], name="H3K4me1_B16_WT>R499_WT", contrast=1, fc=-0.5, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["ATAC"]], name="ATAC_R499_WT>B16_WT", contrast=2, fc=-0.1, sig=0.05)
	gr <- add_DE_annotation(gr=gr, db.DE=db_list[["ATAC"]], name="ATAC_B16_WT>R499_WT", contrast=2, fc=0.1, sig=0.05)
	saveRDS(gr, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Consensus Peaksets/", marks[x], "_consensus_peaks_WT_pooled_Annotated.rds"))
	return(gr)
})
names(gr.marks_list) <- marks

# gr.marks_list <- lapply(1:length(marks), function(x) readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Consensus Peaksets/", marks[x], "_consensus_peaks_WT_pooled_Annotated.rds")))
# names(gr.marks_list) <- marks

####################################
### Annotate Promoters/Enhancers ###
####################################

# Promoters - must overlap with protein-coding TSS
gr.promoter <- gr.marks_list[["H3K4me3"]]
gr.promoter <- overlap_TSS(gr=gr.promoter, minoverlap=100, gr.tss=gr.tss_window, overlap=TRUE)

# Active enhancer - can't overlap with protein-coding TSS
gr.active_enhancer1 <- gr.marks_list[["H3K27Ac"]]
gr.active_enhancer2 <- gr.marks_list[["H3K4me1"]]
gr.active_enhancer2 <- gr.active_enhancer2[gr.active_enhancer2$anno %in% c("ALL", "H3K4me1_H3K27Ac")]
gr.active_enhancer1 <- overlap_TSS(gr=gr.active_enhancer1, minoverlap=100, gr.tss=gr.tss_window, overlap=FALSE)
gr.active_enhancer2 <- overlap_TSS(gr=gr.active_enhancer2, minoverlap=100, gr.tss=gr.tss_window, overlap=FALSE)

RE.list <- list(gr.promoter, gr.active_enhancer1, gr.active_enhancer2)
names(RE.list) <- c("Promoter", "Active_enhancer_H3K27Ac", "Active_enhancer_H3K4me1")
lapply(1:length(RE.list), function(x) saveRDS(RE.list[[x]], file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/", names(RE.list)[x], "_RE.rds")))

###########################################################
### Identify differential promoter and enhancer regions ###
###########################################################

gr.atac <- gr_list[["ATAC"]]
gr.H3K4me3 <- gr.marks_list[["H3K4me3"]]
gr.H3K27Ac <- gr.marks_list[["H3K27Ac"]]
gr.H3K4me1 <- gr.marks_list[["H3K4me1"]]
gr.H3K4me1$IDX <- 1:length(gr.H3K4me1)

### ACTIVATED PROMOTERS (Protein-coding TSS promoters and other H3K4me3-marked peaks) ###

gained_promoters <- gr.H3K4me3[gr.H3K4me3$"H3K4me3_R499_WT>B16_WT"==1]
gained_TSS_promoters <- overlap_TSS(
	gr = gained_promoters,
	minoverlap = 0,
	gr.tss = gr.tss,
	overlap = TRUE)
gained_TSS_promoters_ATAC <- overlap_ATAC(
	gr = gained_TSS_promoters, 
	gr.atac = gr.atac)

gained_other_promoters <- overlap_TSS(
	gr = gained_promoters,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = FALSE)
gained_other_promoters_ATAC <- overlap_ATAC(
	gr = gained_other_promoters, 
	gr.atac = gr.atac)

### DEACTIVATED PROMOTERS (Protein-coding TSS promoters and other H3K4me3-marked peaks) ###

lost_promoters <- gr.H3K4me3[gr.H3K4me3$"H3K4me3_B16_WT>R499_WT"==1]
lost_TSS_promoters <- overlap_TSS(
	gr = lost_promoters,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = TRUE)
lost_TSS_promoters_ATAC <- overlap_ATAC(
	gr = lost_TSS_promoters, 
	gr.atac = gr.atac)
lost_other_promoters <- overlap_TSS(
	gr = lost_promoters,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = FALSE)
lost_other_promoters_ATAC <- overlap_ATAC(
	gr = lost_other_promoters, 
	gr.atac = gr.atac)

### ACTIVATED ENHANCERS (De novo or poised) ###

# Significant increase in H3K27ac
activated_enhancer_ALL <- gr.H3K4me1[gr.H3K4me1$"H3K27Ac_R499_WT>B16_WT"==1]
activated_enhancer_ALL <- add_DE_annotation(
	gr = activated_enhancer_ALL, 
	db.DE = db_list[["ATAC"]], 
	name = "ATAC_R499_WT>B16_WT_FC1.5", 
	contrast = 2, 
	fc = -1.5, 
	sig = 1e-05)
denovo_enhancer <- activated_enhancer_ALL[activated_enhancer_ALL$B16_ATAC == 0 & activated_enhancer_ALL$R499_ATAC == 1 | activated_enhancer_ALL$"ATAC_R499_WT>B16_WT_FC1.5" == 1]
activated_enhancer_ALL$DENOVO <- ifelse(activated_enhancer_ALL$IDX %in% denovo_enhancer$IDX, 1, 0)
activated_enhancer_ALL <- overlap_TSS(
	gr = activated_enhancer_ALL,
	minoverlap = 0,
	gr.tss = gr.tss,
	overlap = FALSE)
activated_enhancer_ALL_ATAC <- overlap_ATAC(gr=activated_enhancer_ALL, gr.atac=gr.atac)

# Significant increase in H3K27ac OR ATAC
activated_enhancer_ALL_ATAC_inclusive <- gr.H3K4me1[gr.H3K4me1$"H3K27Ac_R499_WT>B16_WT"==1 | gr.H3K4me1$"ATAC_R499_WT>B16_WT"==1]
activated_enhancer_ALL_ATAC_inclusive <- activated_enhancer_ALL_ATAC_inclusive[activated_enhancer_ALL_ATAC_inclusive$anno %in% c("ALL", "H3K4me1_H3K27Ac")]
activated_enhancer_ALL_ATAC_inclusive <- add_DE_annotation(
	gr = activated_enhancer_ALL_ATAC_inclusive, 
	db.DE = db_list[["ATAC"]], 
	name = "ATAC_R499_WT>B16_WT_FC1.5", 
	contrast = 2, 
	fc = -1.5, 
	sig = 1e-05)
denovo_enhancer <- activated_enhancer_ALL_ATAC_inclusive[activated_enhancer_ALL_ATAC_inclusive$B16_ATAC == 0 & activated_enhancer_ALL_ATAC_inclusive$R499_ATAC == 1 | activated_enhancer_ALL_ATAC_inclusive$"ATAC_R499_WT>B16_WT_FC1.5" == 1]
activated_enhancer_ALL_ATAC_inclusive$DENOVO <- ifelse(activated_enhancer_ALL_ATAC_inclusive$IDX %in% denovo_enhancer$IDX, 1, 0)
activated_enhancer_ALL_ATAC_inclusive <- overlap_TSS(
	gr = activated_enhancer_ALL_ATAC_inclusive,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = FALSE)
activated_enhancer_ALL_ATAC_inclusive <- overlap_ATAC(
	gr = activated_enhancer_ALL_ATAC_inclusive, 
	gr.atac = gr.atac)

### DEACTIVATED ENHANCERS ###

deactivated_enhancer_ALL <- gr.H3K4me1[gr.H3K4me1$"H3K27Ac_B16_WT>R499_WT"==1]
deactivated_enhancer_ALL <- overlap_TSS(
	gr = deactivated_enhancer_ALL,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = FALSE)
# deactivated_enhancer_ALL <- overlap_ATAC(gr=deactivated_enhancer_ALL, gr.atac=gr.atac)

# Significant decrease in H3K27ac OR ATAC
deactivated_enhancer_ALL_ATAC_inclusive <- gr.H3K4me1[gr.H3K4me1$"H3K27Ac_B16_WT>R499_WT"==1 | gr.H3K4me1$"ATAC_B16_WT>R499_WT"==1]
deactivated_enhancer_ALL_ATAC_inclusive <- deactivated_enhancer_ALL_ATAC_inclusive[deactivated_enhancer_ALL_ATAC_inclusive$anno %in% c("ALL", "H3K4me1_H3K27Ac")]
deactivated_enhancer_ALL_ATAC_inclusive <- overlap_TSS(
	gr = deactivated_enhancer_ALL_ATAC_inclusive,
	minoverlap = 0,
	gr.tss = gr.tss, 
	overlap = FALSE)
deactivated_enhancer_ALL_ATAC_inclusive <- overlap_ATAC(
	gr = deactivated_enhancer_ALL_ATAC_inclusive, 
	gr.atac = gr.atac)

### UNCHANGED REGULATORY ELEMENTS ###

set.seed(42)
unchanged_REs <- gr.H3K4me1[gr.H3K4me1$"H3K27Ac_R499_WT>B16_WT"==0 & gr.H3K4me1$"H3K27Ac_B16_WT>R499_WT"==0 & gr.H3K4me1$"H3K4me1_R499_WT>B16_WT"==0 & gr.H3K4me1$"H3K4me1_B16_WT>R499_WT"==0 & gr.H3K4me1$"H3K4me3_R499_WT>B16_WT"==0 & gr.H3K4me1$"H3K4me3_B16_WT>R499_WT"==0]
unchanged_REs_ALL <- unchanged_REs
unchanged_REs <- unchanged_REs[sample(1:length(unchanged_REs), 2000)]

# Random
random_peaks <- gr.H3K4me1[sample(1:length(gr.H3K4me1), 2000)]
random_peaks_50k <- gr.H3K4me1[sample(1:length(gr.H3K4me1), 50000)]

#################
### WRITE OUT ###
#################

RE.names <- c(
	"gained_TSS_promoters", "gained_TSS_promoters_ATAC",
	"gained_other_promoters", "gained_other_promoters_ATAC",
	"lost_TSS_promoters", "lost_TSS_promoters_ATAC", 
	"lost_other_promoters", "lost_other_promoters_ATAC", 
	"activated_enhancer", "activated_enhancer_ATAC", "activated_enhancer_ATAC_inclusive", 
	"deactivated_enhancer", "deactivated_enhancer_ATAC_inclusive",
	"unchanged_REs", "random_peaks")
RE_list <- list(gained_TSS_promoters, gained_TSS_promoters_ATAC, 
	gained_other_promoters, gained_other_promoters_ATAC, 
	lost_TSS_promoters, lost_TSS_promoters_ATAC, 
	lost_other_promoters, lost_other_promoters_ATAC, 
	activated_enhancer_ALL, activated_enhancer_ALL_ATAC, activated_enhancer_ALL_ATAC_inclusive,
	deactivated_enhancer_ALL, deactivated_enhancer_ALL_ATAC_inclusive,
	unchanged_REs, random_peaks)
names(RE_list) <- RE.names
sapply(RE_list, function(x) length(x))

lapply(1:length(RE.names), function(x) {
	print(RE.names[x])
	make_HOMER_bed(RE_list[[x]], file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/", RE.names[x], ".bed"))
})

#################
### Functions ###
#################

make_HOMER_bed <- function(gr, file=NULL) {
	bed <- data.frame(gr)[,1:3]
	bed$ID <- paste(bed[,1], bed[,2], bed[,3], sep="_")
	bed$col5 <- NA
	bed$strand <- as.character(strand(gr))
	if(!is.null(file)) write.table(bed, file=file, sep="\t", row.names=F, col.names=F, quote=F)
}

write_fasta <- function(gr, file=NULL) {
    fa <- getSeq(Mmusculus, gr)
    names(fa) <- paste0(as.character(gr@seqnames), ":", start(gr), "-", end(gr))
    writeXStringSet(fa, filepath=file, append=FALSE)
}

# Convert peak IDs (chr_start_end) to GRanges
peakids2GRanges <- function(peakids, delim="_") {
	df <- data.frame(chr=sapply(strsplit(peakids, split=delim), function(x) x[1]), start=as.numeric(sapply(strsplit(peakids, split=delim), function(x) x[2])), end=as.numeric(sapply(strsplit(peakids, split=delim), function(x) x[3])))
	gr <- makeGRangesFromDataFrame(df, starts.in.df.are.0based=FALSE)
	return(gr)
}

# overlap=TRUE: retain only peaks overlapping protein-coding TSS site
# overlap=FALSE: retain only peaks NOT overlapping protein-coding TSS site
overlap_TSS <- function(gr, minoverlap=100, gr.tss, overlap=TRUE) {
	ol <- findOverlaps(gr, gr.tss, minoverlap=minoverlap)
	if(overlap == TRUE) {
		gr <- gr[unique(from(ol))]
	} else {
		gr <- gr[-unique(from(ol))]
	}
	return(gr)
}

# Retain only peaks overlapping accessible peak (gr.atac)
overlap_ATAC <- function(gr, gr.atac) {
	ol <- findOverlaps(gr, gr.atac, minoverlap=100)
	# gr.open <- gr.atac[unique(to(ol))] # Use ATAC peaks (not CUT&RUN peaks)
	gr.open <- gr[unique(from(ol))] # Use peaks from peak set of interest
	print(paste0(length(gr.open), " out of ", length(gr), " peaks selected after ATAC filter"))
	return(gr.open)
}

# Use DiffBind result to annotate peaks with significantly increased or decreased signal
add_DE_annotation <- function(gr, db.DE, name, contrast, fc, sig) {
	gr.DE <- dba.report(db.DE, contrast=contrast, th=sig)
	if(fc > 0) {
		gr.DE <- gr.DE[gr.DE$Fold > fc]
	} else {
		gr.DE <- gr.DE[gr.DE$Fold < fc]
	}
	ol <- findOverlaps(gr, gr.DE, minoverlap=190)
	gr$DE <- ifelse(1:length(gr) %in% unique(from(ol)), 1, 0)
	colnames(gr@elementMetadata)[ncol(gr@elementMetadata)] <- name
	return(gr)
}
