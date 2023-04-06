rm(list=ls())

library(tidyverse)
library(GSVA)
library(rtracklayer)
library(RColorBrewer)
library(ggsci)
library(patchwork)

directory <- "local"

source("~/Dropbox/Minn/resources/useful_protocols/Bulk ATAC processing/ATAC Visualization and Analysis Functions.R")

gr.rmsk <- import("~/Dropbox/Minn/resources/TE/squire/mm10_all.bed")
gr.rmsk$ID <- sapply(strsplit(gr.rmsk$name, split="\\|"), function(x) x[[4]])
gr.rmsk$subF <- sapply(strsplit(gr.rmsk$ID, split="\\:"), function(x) x[[1]])
gr.rmsk$Family <- sapply(strsplit(gr.rmsk$ID, split="\\:"), function(x) x[[2]])
gr.rmsk$Class <- sapply(strsplit(gr.rmsk$ID, split="\\:"), function(x) x[[3]])

identify_TE_overlaps <- function(gr, gr.rmsk) {
	ol <- findOverlaps(gr, gr.rmsk, ignore.strand = T)
	df.te <- lapply(1:length(gr), function(x) {
		idx <- to(ol)[which(from(ol) == x)]
		if(length(idx) > 0) {
			te.id <- paste(unique(gr.rmsk$ID[idx]), collapse=",")
			te.subF <- paste(unique(gr.rmsk$subF[idx]), collapse=",")
			te.class <- paste(unique(gr.rmsk$Family[idx]), collapse=",")
		} else {
			te.id <- NA
			te.subF <- NA
			te.class <- NA
		}
		df <- data.frame(id=te.id, subF=te.subF, class=te.class)
		return(df)
	})
	df.te <- do.call(rbind, df.te)
	gr$ID <- df.te$id
	gr$subF <- df.te$subF
	gr$Class <- df.te$class
	return(gr)
}

#######################################
### StringTie-Assembled Transcripts ###
#######################################

# Res499 clones assembled transcripts
gr.TACO_Transcripts <- import("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/stringtie/TACO_merge/assembly.refcomp.gtf")
gr.TACO_Transcripts$Coordinates <- paste0(as.character(seqnames(gr.TACO_Transcripts)), ":", start(gr.TACO_Transcripts), "-", end(gr.TACO_Transcripts))
print(paste0("Number of transcript ids: ", length(unique(gr.TACO_Transcripts$transcript_id))))

# Res499 clones TE-derived transcripts
exon_overlap_TE_width <- readRDS(paste0("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/Exon_Overlap_TE_Width.rds"))
gr.TE_Transcripts_20 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% exon_overlap_TE_width$transcript_id[exon_overlap_TE_width$te_perc > 20]]
gr.TE_Transcripts_40 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% exon_overlap_TE_width$transcript_id[exon_overlap_TE_width$te_perc > 40]]
gr.TE_Transcripts_50 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% exon_overlap_TE_width$transcript_id[exon_overlap_TE_width$te_perc > 50]]
print(paste0("Number of transcript ids: ", length(unique(gr.TE_Transcripts_20$transcript_id))))
print(paste0("Number of transcript ids: ", length(unique(gr.TE_Transcripts_40$transcript_id))))
print(paste0("Number of transcript ids: ", length(unique(gr.TE_Transcripts_50$transcript_id))))

# Extend TE transcript windows
extend_W <- 2500
gr.TACO_Transcripts_Extended <- extend_regions(gr = gr.TACO_Transcripts, extend_start = extend_W, extend_end = extend_W)
gr.TE_Transcripts_20_Extended <- extend_regions(gr = gr.TE_Transcripts_20, extend_start = extend_W, extend_end = extend_W)
gr.TE_Transcripts_40_Extended <- extend_regions(gr = gr.TE_Transcripts_40, extend_start = extend_W, extend_end = extend_W)
gr.TE_Transcripts_50_Extended <- extend_regions(gr = gr.TE_Transcripts_50, extend_start = extend_W, extend_end = extend_W)

# DE TE transcripts
DEG_list_Res499_PRRKO <- readRDS("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.DEG_list.rds")
DEG_list_B16_R499 <- readRDS("~/Dropbox/Minn/ifnar_epigenome_V2/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.DEG_list.rds")
DEG_list_Res499_Clones <- readRDS("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.DEG_list.rds")

#######################
### Sample metadata ###
#######################

metadat.Res499_Clones <- read.csv("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/metadata.csv", sep=",", header=T)
metadat.Res499_Clones$Class <- factor(metadat.Res499_Clones$Class, levels = c("B16 WT", "R499 WT", "Res499 SCP-Sensitive1", "Res499 SCP-Sensitive2", "Res499 SCP 11", "Res499 SCP-Resistant1", "Res499 SCP-Resistant2"))
metadat.Res499_Clones$Resistant <- factor(metadat.Res499_Clones$Resistant, levels = c("B16 WT", "R499 WT", "Sensitive", "Resistant", "Res499 SCP 11"))

metadat.B16_R499 <- read.csv("~/Dropbox/Minn/ifnar_epigenome_V2/data/RNA/metadata.csv", sep=",", header=T)
metadat.B16_R499$Label <- factor(metadat.B16_R499$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

metadat.Res499_PRRKO <- read.csv("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/metadata.csv", sep=",", header=T)
metadat.Res499_PRRKO$Label <- factor(metadat.Res499_PRRKO$Label, levels = c("B16", "B16y DY", "Res499", "Res499 PRR dKO", "Res499 Cmpd1 RUX 3wk", "Res499 RUX 3wk", "Res499 ITA 3wk"))

######################
### Feature Counts ###
######################

edat.Res499_Clones <- readRDS("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.featureCounts_DESeq2_Normalization.rds")
edat.B16_R499 <- readRDS("~/Dropbox/Minn/ifnar_epigenome_V2/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.featureCounts_DESeq2_Normalization.rds")
edat.Res499_PRRKO <- readRDS("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/FeatureCounts/TACO_Assembly_Res499_Clones.featureCounts_DESeq2_Normalization.rds")

#######################
### Epigenomic data ###
#######################

# STAT/IRF3 KO ATAC-seq normalized counts
mat.ATAC_revisions <- read.table("~/Dropbox/Minn/ifnar_epigenome_V2/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_revisions_consensus_mat_tn5_insertion_counts_IDR_vst.txt", sep="\t", header=T, check.names = F)
colnames(mat.ATAC_revisions) <- gsub("\\.", "_", colnames(mat.ATAC_revisions))
gr.ATAC_revisions <- peakids2GRanges(rownames(mat.ATAC_revisions), delim = "_")
gr.ATAC_revisions_Extended <- extend_regions(gr = gr.ATAC_revisions, extend_start = extend_W, extend_end = extend_W)
metadata.ATAC_revisions <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_ATAC_revisions.csv", header=T)
metadata.ATAC_revisions$Label <- factor(metadata.ATAC_revisions$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO", "Res499 IRF3 KO", "Res499 DKO"))

# H3K4me1 normalized count matrix
mat.H3K4me1 <- read.table("~/Dropbox/Minn/ifnar_epigenome_V2/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/consensus_mat_insertion_counts_WT_H3K4me1_pooled_vst_V2.txt", sep="\t", header=T)
colnames(mat.H3K4me1) <- gsub("_H3K4me1_rep", "_", colnames(mat.H3K4me1))
gr.H3K4me1 <- peakids2GRanges(peakids = rownames(mat.H3K4me1), delim="_")
gr.H3K4me1_Extended <- extend_regions(gr = gr.H3K4me1, extend_start = extend_W, extend_end = extend_W)
metadata.H3K4me1 <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_H3K4me1.csv", header=T)
metadata.H3K4me1$Label <- factor(metadata.H3K4me1$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# H3K4me3 normalized count matrix
mat.H3K4me3 <- read.table("~/Dropbox/Minn/ifnar_epigenome_V2/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/consensus_mat_insertion_counts_WT_H3K4me3_pooled_vst_V2.txt", sep="\t", header=T)
colnames(mat.H3K4me3) <- gsub("_H3K4me3_rep", "_", colnames(mat.H3K4me3))
gr.H3K4me3 <- peakids2GRanges(peakids = rownames(mat.H3K4me3), delim="_")
gr.H3K4me3_Extended <- extend_regions(gr = gr.H3K4me3, extend_start = extend_W, extend_end = extend_W)
metadata.H3K4me3 <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_H3K4me3.csv", header=T)
metadata.H3K4me3$Label <- factor(metadata.H3K4me3$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# H3K27ac normalized count matrix
mat.H3K27ac <- read.table("~/Dropbox/Minn/ifnar_epigenome_V2/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/consensus_mat_insertion_counts_WT_H3K27ac_pooled_vst_V2.txt", sep="\t", header=T)
colnames(mat.H3K27ac) <- gsub("_H3K27Ac_rep", "_", colnames(mat.H3K27ac))
gr.H3K27ac <- peakids2GRanges(peakids = rownames(mat.H3K27ac), delim="_")
gr.H3K27ac_Extended <- extend_regions(gr = gr.H3K27ac, extend_start = extend_W, extend_end = extend_W)
metadata.H3K27ac <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_H3K27ac.csv", header=T)
metadata.H3K27ac$Label <- factor(metadata.H3K27ac$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# PRR KO ATAC-seq normalized counts
metadat.Res499_PRRKO_ATAC <- read.csv("~/Dropbox/Minn/Epi_JAK_ATAC/data/ATAC/processed/metadata.csv", header=T)
metadat.Res499_PRRKO_ATAC <- metadat.Res499_PRRKO_ATAC[!metadat.Res499_PRRKO_ATAC$Label %in% "B16y JB", ]
metadat.Res499_PRRKO_ATAC$Label <- factor(metadat.Res499_PRRKO_ATAC$Label, levels = c("B16", "B16y DY", "Res499", "Res499 takedown 3wk", "Res499 takedown 5wk", "Res499 PRR dKO", "Res499 Cmpd1 RUX 3wk", "Res499 CYT 3wk", "Res499 RUX 3wk", "Res499 RUX 5wk", "Res499 ITA 3wk", "Res499 ITA 5wk", "Res499 AR 3wk", "Res499 AR 5wk"))
adat.Res499_PRRKO <- read.table("~/Dropbox/Minn/Epi_JAK_ATAC/data/ATAC/processed/mat/consensus_mat_tn5_insertion_counts_IDR_vst_all.txt", sep="\t", header=T)
adat.Res499_PRRKO <- adat.Res499_PRRKO[ ,match(metadat.Res499_PRRKO_ATAC$Sample, colnames(adat.Res499_PRRKO))]
gr.Res499_PRRKO <- peakids2GRanges(rownames(adat.Res499_PRRKO), delim = "_")
gr.Res499_PRRKO_Extended <- extend_regions(gr = gr.Res499_PRRKO, extend_start = extend_W, extend_end = extend_W)

###########################
### Regions of interest ###
###########################

extend_W <- 1500

# Clone 11 up regions
gr.Clone11 <- import.bed("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/ATAC/merged_runs/data/processed/diffbind/DAR/DAR_Sensitive3_UP.bed")
gr.Clone11 <- filterAnnotations(
	gr = gr.Clone11,
	gr.annotation = gr.exons_protein_coding,
	minoverlap = 50,
	overlap = FALSE)
gr.Clone11_Extended <- extend_regions(gr = gr.Clone11, extend_start = extend_W, extend_end = extend_W)

# IMDs
gr.IMD <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/IFN Inflammatory Memory Domains ALL.bed")
gr.IMD_Filtered <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/IFN Inflammatory Memory Domains Activated REs.bed")
gr.IMD_Extended <- extend_regions(gr = gr.IMD, extend_start = extend_W, extend_end = extend_W)
gr.IMD_Filtered_Extended <- extend_regions(gr = gr.IMD_Filtered, extend_start = extend_W, extend_end = extend_W)

# R499 activated H3K4me3
gr.R499_Activated_H3K4me3 <- RE_list[["gained_other_promoters"]]
gr.R499_Activated_H3K4me3 <- filterAnnotations(
	gr = gr.R499_Activated_H3K4me3,
	gr.annotation = gr.exons_protein_coding,
	minoverlap = 50,
	overlap = FALSE)
gr.R499_Activated_H3K4me3_Extended <- extend_regions(gr = gr.R499_Activated_H3K4me3, extend_start = extend_W, extend_end = extend_W)

# R499 activated enhancer
gr.R499_Activated_Enhancer <- RE_list[["activated_enhancer_ATAC_inclusive"]]
gr.R499_Activated_Enhancer <- filterAnnotations(
	gr = gr.R499_Activated_Enhancer,
	gr.annotation = gr.exons_protein_coding,
	minoverlap = 50,
	overlap = FALSE)
gr.R499_Activated_Enhancer_Extended <- extend_regions(gr = gr.R499_Activated_Enhancer, extend_start = extend_W, extend_end = extend_W)

# # R499 activated regions
# gr.R499_Activated <- c(gr.R499_Activated_Enhancer, gr.R499_Activated_H3K4me3, gr.IMD)
# gr.R499_Activated <- filterAnnotations(
# 	gr = gr.R499_Activated,
# 	gr.annotation = gr.exons_protein_coding,
# 	minoverlap = 50,
# 	overlap = FALSE)
# gr.R499_Activated_Extended <- extend_regions(gr = gr.R499_Activated, extend_start = extend_W, extend_end = extend_W)

##############################################################
### Overlap TE transcripts with R499-activated Epi Regions ###
##############################################################

# Res 499 SCP 11
ol <- findOverlaps(gr.TE_Transcripts_40_Extended, gr.Clone11_Extended, minoverlap=50, ignore.strand = T)
gr.Clone11_TE_Transcripts <- gr.TE_Transcripts_40_Extended[unique(from(ol))]
print(paste0("# of Res 499 SCP 11 TE Transcripts: ", length(unique(gr.Clone11_TE_Transcripts$transcript_id))))

# B16 vs Res499 (PRR KO dataset)
# res <- DEG_list_Res499_PRRKO[["Res499_vs_B16"]]
res <- DEG_list_Res499_Clones[["Sensitive3_vs_Sensitive1"]]
res <- res[!is.na(res$padj), ]
res.keep1 <- res[which(res$log2FoldChange > 0), ]
res.keep1 <- res.keep1[res.keep1$baseMean > 10, ]
res <- DEG_list_Res499_Clones[["Resistant2_vs_Sensitive1"]]
res <- res[!is.na(res$padj), ]
res.keep2 <- res[which(res$log2FoldChange > 0), ]
res.keep2 <- res.keep2[res.keep2$baseMean > 10, ]
gr.Clone11_TE_Transcripts <- gr.Clone11_TE_Transcripts[gr.Clone11_TE_Transcripts$transcript_id %in% intersect(rownames(res.keep1), rownames(res.keep2))]
gr.Clone11_TE_Transcripts$log2fc_Clone11_vs_Sensitive1 <- res$log2FoldChange[match(gr.Clone11_TE_Transcripts$transcript_id, rownames(res))]
gr.Clone11_TE_Transcripts <- gr.Clone11_TE_Transcripts[sort.int(gr.Clone11_TE_Transcripts$log2fc_Clone11_vs_Sensitive1, decreasing=T, index.return=T)$ix]
gr.Clone11_TE_Transcripts <- identify_TE_overlaps(gr.Clone11_TE_Transcripts, gr.rmsk)
print(paste0("# of Res 499 SCP 11 TE Transcripts: ", length(unique(gr.Clone11_TE_Transcripts$transcript_id))))

# Write out
gr.interest <- gr.Clone11_TE_Transcripts
print(table(gr.interest$type))
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Clone11.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Clone11_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Clone11_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

head(data.frame(gr.Clone11_TE_Transcripts[gr.Clone11_TE_Transcripts$type == "transcript"])[ ,c("Coordinates", "Class", "log2fc_Clone11_vs_Sensitive1", "ref_gene_type")], 25)

##############################################################

# IFN-IMD
ol <- findOverlaps(gr.TE_Transcripts_40_Extended, gr.IMD_Extended, minoverlap=100, ignore.strand = T)
gr.IMD_TE_Transcripts <- gr.TE_Transcripts_40_Extended[unique(from(ol))]
print(paste0("# of IFN-IMD TE Transcripts: ", length(unique(gr.IMD_TE_Transcripts$transcript_id))))

# Write out
gr.interest <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% gr.IMD_TE_Transcripts$transcript_id]
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
print(table(gr.interest$type))
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_IMD.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_IMD_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_IMD_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)
# gr.IMD_TE_Transcripts <- gr.IMD_TE_Transcripts[gr.IMD_TE_Transcripts$transcript_id %in% gr.Res499_vs_B16_50$transcript_id]
# write.table(data.frame(gr.IMD_TE_Transcripts)[,1:3], file = "~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_IMD_Epi.bed", sep="\t", col.names=F, row.names=F, quote=F)

# IFN-IMD peaks
gr.IMD_Extended_Peaks <- gr.IMD_Extended[unique(to(ol))]
print(paste0("# of IFN-IMD Peaks: ", length(gr.IMD_Extended_Peaks)))
write.table(data.frame(gr.IMD_Extended_Peaks)[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_IMD_Peaks.bed", sep="\t", col.names=F, row.names=F, quote=F)

# R499 activated H3K4me3
ol <- findOverlaps(gr.TE_Transcripts_40_Extended, gr.R499_Activated_H3K4me3_Extended, minoverlap=100, ignore.strand = T)
gr.R499_Activated_H3K4me3_TE_Transcripts <- gr.TE_Transcripts_40_Extended[unique(from(ol))]
print(paste0("# of R499_Activated_H3K4me3 TE Transcripts: ", length(unique(gr.R499_Activated_H3K4me3_TE_Transcripts$transcript_id))))

# Write out
gr.interest <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% gr.R499_Activated_H3K4me3_TE_Transcripts$transcript_id]
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
print(table(gr.interest$type))
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_H3K4me3.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_H3K4me3_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_H3K4me3_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

# H3K4me3 peaks
gr.R499_Activated_H3K4me3_H3K4me3 <- gr.R499_Activated_H3K4me3_Extended[unique(to(ol))]
print(paste0("# of R499_Activated_H3K4me3 Peaks: ", length(gr.R499_Activated_H3K4me3_H3K4me3)))
write.table(data.frame(gr.R499_Activated_H3K4me3_H3K4me3)[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_H3K4me3_Peaks.bed", sep="\t", col.names=F, row.names=F, quote=F)

# R499 activated enhancer
ol <- findOverlaps(gr.TE_Transcripts_40_Extended, gr.R499_Activated_Enhancer_Extended, minoverlap=100, ignore.strand = T)
gr.R499_Activated_Enhancer_TE_Transcripts <- gr.TE_Transcripts_40_Extended[unique(from(ol))]
print(paste0("# of R499_Activated_H3K4me3 TE Transcripts: ", length(unique(gr.R499_Activated_Enhancer_TE_Transcripts$transcript_id))))

# Write out
gr.interest <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% gr.R499_Activated_Enhancer_TE_Transcripts$transcript_id]
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
print(table(gr.interest$type))
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_Enhancer.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_Enhancer_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_Enhancer_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

# IFN-IMD peaks
gr.R499_Activated_Enhancer_Peaks <- gr.R499_Activated_Enhancer[unique(to(ol))]
print(paste0("# of R499 Activated Enhancer Peaks: ", length(gr.R499_Activated_Enhancer_Peaks)))
write.table(data.frame(gr.R499_Activated_Enhancer_Peaks)[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_R499_Activated_Enhancer_Peaks.bed", sep="\t", col.names=F, row.names=F, quote=F)

# R499 activated regions
# ol <- findOverlaps(gr.TE_Transcripts_20_Extended, gr.R499_Activated_Extended)
# gr.R499_Activated_TE_Transcripts <- gr.TE_Transcripts_20_Extended[unique(from(ol))]
# gr.R499_Activated_TE_Transcripts <- gr.R499_Activated_TE_Transcripts[gr.R499_Activated_TE_Transcripts$transcript_id %in% gr.Res499_vs_B16_50$transcript_id]
# gr.R499_Activated_TE_Transcripts <- c(gr.IMD_TE_Transcripts, gr.R499_Activated_Enhancer_TE_Transcripts)
# write.table(data.frame(gr.R499_Activated_TE_Transcripts)[,1:3], file = "~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_R499_Activated_Epi.bed", sep="\t", col.names=F, row.names=F, quote=F)

################################
### DE StringTie transcripts ###
################################

# B16 vs Res499 (PRR KO dataset)
res <- DEG_list_Res499_PRRKO[["Res499_vs_B16"]]
res <- res[!is.na(res$padj), ]
res <- res[which(res$pvalue < 0.1), ]
res <- res[res$log2FoldChange > 0, ]
res <- res[sort.int(res$padj, decreasing=F, index.return=T)$ix, ]
res$Coordinates <- gr.TACO_Transcripts$Coordinates[match(rownames(res), gr.TACO_Transcripts$transcript_id)]
gr.Res499_vs_B16_20 <- gr.TE_Transcripts_20[gr.TE_Transcripts_20$transcript_id %in% rownames(res)]
gr.Res499_vs_B16_40 <- gr.TE_Transcripts_40[gr.TE_Transcripts_40$transcript_id %in% rownames(res)]
gr.Res499_vs_B16_50 <- gr.TE_Transcripts_50[gr.TE_Transcripts_50$transcript_id %in% rownames(res)]
print(paste0("# Significant TE Transcripts (20% TE overlap): ", length(unique(gr.Res499_vs_B16_20$transcript_id))))
print(paste0("# Significant TE Transcripts (40% TE overlap): ", length(unique(gr.Res499_vs_B16_40$transcript_id))))
print(paste0("# Significant TE Transcripts (50% TE overlap): ", length(unique(gr.Res499_vs_B16_50$transcript_id))))

gr.Res499_vs_B16_40$log2fc <- res$log2FoldChange[match(gr.Res499_vs_B16_40$transcript_id, rownames(res))]
gr.Res499_vs_B16_40 <- gr.Res499_vs_B16_40[sort.int(gr.Res499_vs_B16_40$log2fc, decreasing=T, index.return=T)$ix]

# Write out
gr.interest <- gr.Res499_vs_B16_40
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

# PRR dKO
res <- DEG_list_Res499_PRRKO[["Res499_vs_Res499_PRR_dKO"]]
res <- res[!is.na(res$padj), ]
res <- res[which(res$pvalue < 0.1), ]
res <- res[res$log2FoldChange > 0, ]
res <- res[sort.int(res$padj, decreasing=F, index.return=T)$ix, ]
res$Coordinates <- gr.TACO_Transcripts$Coordinates[match(rownames(res), gr.TACO_Transcripts$transcript_id)]
gr.Res499_vs_Res499_PRR_dKO_20 <- gr.TE_Transcripts_20[gr.TE_Transcripts_20$transcript_id %in% rownames(res)]
gr.Res499_vs_Res499_PRR_dKO_40 <- gr.TE_Transcripts_40[gr.TE_Transcripts_40$transcript_id %in% rownames(res)]
gr.Res499_vs_Res499_PRR_dKO_50 <- gr.TE_Transcripts_50[gr.TE_Transcripts_50$transcript_id %in% rownames(res)]
print(paste0("# Significant TE Transcripts (20% TE overlap): ", length(unique(gr.Res499_vs_Res499_PRR_dKO_20$transcript_id))))
print(paste0("# Significant TE Transcripts (40% TE overlap): ", length(unique(gr.Res499_vs_Res499_PRR_dKO_40$transcript_id))))
print(paste0("# Significant TE Transcripts (50% TE overlap): ", length(unique(gr.Res499_vs_Res499_PRR_dKO_50$transcript_id))))

# R499 activated, Res499 PRR KO
gr.Res499_PRRdep_40 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% intersect(gr.Res499_vs_B16_40$transcript_id, gr.Res499_vs_Res499_PRR_dKO_40$transcript_id)]
print(paste0("# Res499 PRR-dependent TE Transcripts: ", length(unique(gr.Res499_PRRdep_40$transcript_id))))

gr.Res499_PRRdep_40$log2fc <- res$log2FoldChange[match(gr.Res499_PRRdep_40$transcript_id, rownames(res))]
gr.Res499_PRRdep_40 <- gr.Res499_PRRdep_40[sort.int(gr.Res499_PRRdep_40$log2fc, decreasing=T, index.return=T)$ix]

# Write out
gr.interest <- gr.Res499_PRRdep_40
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

# TBKi
res <- DEG_list_Res499_PRRKO[["Res499_vs_Res499_Cmpd1_RUX_3wk"]]
res <- res[!is.na(res$padj), ]
res <- res[which(res$pvalue < 0.1), ]
res <- res[res$log2FoldChange > 0, ]
res <- res[sort.int(res$padj, decreasing=F, index.return=T)$ix, ]
res$Coordinates <- gr.TACO_Transcripts$Coordinates[match(rownames(res), gr.TACO_Transcripts$transcript_id)]
gr.Res499_vs_Res499_Cmpd1_RUX_3wk_20 <- gr.TE_Transcripts_20[gr.TE_Transcripts_20$transcript_id %in% rownames(res)]
gr.Res499_vs_Res499_Cmpd1_RUX_3wk_40 <- gr.TE_Transcripts_40[gr.TE_Transcripts_40$transcript_id %in% rownames(res)]
gr.Res499_vs_Res499_Cmpd1_RUX_3wk_50 <- gr.TE_Transcripts_50[gr.TE_Transcripts_50$transcript_id %in% rownames(res)]
print(paste0("# Significant TE Transcripts (20% TE overlap): ", length(unique(gr.Res499_vs_Res499_Cmpd1_RUX_3wk_20$transcript_id))))
print(paste0("# Significant TE Transcripts (40% TE overlap): ", length(unique(gr.Res499_vs_Res499_Cmpd1_RUX_3wk_40$transcript_id))))
print(paste0("# Significant TE Transcripts (50% TE overlap): ", length(unique(gr.Res499_vs_Res499_Cmpd1_RUX_3wk_50$transcript_id))))

# R499 activated, Res499 PRR KO
gr.Res499_Cmpd1_RUX_40 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% intersect(gr.Res499_vs_B16_40$transcript_id, gr.Res499_vs_Res499_Cmpd1_RUX_3wk_40$transcript_id)]
print(paste0("# Res499 Cmpd1 RUX TE Transcripts: ", length(unique(gr.Res499_Cmpd1_RUX_40$transcript_id))))

gr.Res499_Cmpd1_RUX_40$log2fc <- res$log2FoldChange[match(gr.Res499_Cmpd1_RUX_40$transcript_id, rownames(res))]
gr.Res499_Cmpd1_RUX_40 <- gr.Res499_Cmpd1_RUX_40[sort.int(gr.Res499_Cmpd1_RUX_40$log2fc, decreasing=T, index.return=T)$ix]

# Write out
gr.interest <- gr.Res499_Cmpd1_RUX_40
gr.interest <- identify_TE_overlaps(gr.interest, gr.rmsk)
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX_Transcripts.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX_Exons.bed", sep="\t", col.names=F, row.names=F, quote=F)

# See differential expression results
res <- DEG_list_Res499_PRRKO[["Res499_vs_B16"]]
res <- res[!is.na(res$padj), ]
res <- res[which(res$pvalue < 0.1), ]
res <- res[res$log2FoldChange > 0, ]
res <- res[sort.int(res$padj, decreasing=F, index.return=T)$ix, ]

gr.Res499_PRRdep_40$log2fc <- res$log2FoldChange[match(gr.Res499_PRRdep_40$transcript_id, rownames(res))]
gr.IMD_TE_Transcripts$log2fc <- res$log2FoldChange[match(gr.IMD_TE_Transcripts$transcript_id, rownames(res))]

test <- gr.Res499_PRRdep_40
test <- test[!is.na(test$log2fc)]
test <- test[sort.int(test$log2fc, decreasing=T, index.return=T)$ix]
head(data.frame(test), 25)

##################################
### Alternative B16/Res499 TEs ###
##################################

# B16 vs Res499 (ifnar_epigenome dataset)
res <- DEG_list_B16_R499[["R499_WT_vs_B16_WT"]]
res <- res[!is.na(res$padj), ]
res <- res[which(res$pvalue < 0.1), ]
res <- res[res$log2FoldChange > 0, ]
res <- res[sort.int(res$padj, decreasing=F, index.return=T)$ix, ]
res$Coordinates <- gr.TACO_Transcripts$Coordinates[match(rownames(res), gr.TACO_Transcripts$transcript_id)]
gr.Res499_vs_B16_40 <- gr.TE_Transcripts_40[gr.TE_Transcripts_40$transcript_id %in% rownames(res)]
print(paste0("# Significant TE Transcripts (40% TE overlap): ", length(unique(gr.Res499_vs_B16_40$transcript_id))))

# Write out
gr.interest <- gr.Res499_vs_B16_40
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16_V1.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16_Transcripts_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_vs_B16_Exons_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)

# R499 activated, Res499 PRR KO
gr.Res499_PRRdep_40 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% intersect(gr.Res499_vs_B16_40$transcript_id, gr.Res499_vs_Res499_PRR_dKO_40$transcript_id)]
print(paste0("# Res499 PRR-dependent TE Transcripts: ", length(unique(gr.Res499_PRRdep_40$transcript_id))))

# Write out
gr.interest <- gr.Res499_PRRdep_40
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep_V1.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep_Transcripts_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Res499_PRRdep_Exons_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)

# R499 activated, Res499 PRR KO
gr.Res499_Cmpd1_RUX_40 <- gr.TACO_Transcripts[gr.TACO_Transcripts$transcript_id %in% intersect(gr.Res499_vs_B16_40$transcript_id, gr.Res499_vs_Res499_Cmpd1_RUX_3wk_40$transcript_id)]
print(paste0("# Res499 Cmpd1 RUX TE Transcripts: ", length(unique(gr.Res499_Cmpd1_RUX_40$transcript_id))))

# Write out
gr.interest <- gr.Res499_Cmpd1_RUX_40
export(gr.interest, "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX_V1.gtf")
write.table(data.frame(gr.interest[gr.interest$type == "transcript"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX_Transcripts_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)
write.table(data.frame(gr.interest[gr.interest$type == "exon"])[,1:3], file = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Interest/TE_Transcripts_Cmpd1_RUX_Exons_V1.bed", sep="\t", col.names=F, row.names=F, quote=F)

#################################
### DE TE-derived transcripts ###
#################################

# gr.Res499_vs_B16 <- import("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_20_Res499_vs_B16.gtf")
# gr.Res499_vs_Res499_PRR_dKO <- import("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_Res499_vs_Res499_PRR_dKO.gtf")
# gr.Res499_vs_Res499_Cmpd1_RUX_3wk <- import("~/Dropbox/Minn/Epi_JAK_ATAC/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_Res499_vs_Res499_Cmpd1_RUX_3wk.gtf")

#################
### Visualize ###
#################

gr.interest <- gr.Res499_vs_B16_noProteinCoding_Clone11
gr <- gr.TACO_Transcripts_Extended[gr.TACO_Transcripts_Extended$transcript_id %in% gr.interest$transcript_id]
gr <- gr[gr$type == "transcript"]

mat.plot <- edat.Res499_Clones
use.metadat <- metadat.Res499_Clones
use.pal <- c("navyblue", "firebrick", brewer.pal(9, "Blues")[c(3,4,6)], pal_jama()(7)[5], brewer.pal(9, "Reds")[c(3,5)])
names(use.pal) <- c("B16 WT", "R499 WT", "B16 SCP", "Res499 SCP-Sensitive1", "Res499 SCP-Sensitive2", "Res499 SCP 11", "Res499 SCP-Resistant1", "Res499 SCP-Resistant2")

# Res 499 clones
p1 <- plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(unique(gr$transcript_id)),
    geneset_names = paste0("TE-Derived Transcripts\n(n=", length(unique(gr$transcript_id)), ")"),
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Class",
    colorBy = "Class",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    point_size = 2.5,
    # stroke = 1.5,
    plot.legend = F)

# PRR KO
mat.plot <- edat.Res499_PRRKO
use.metadat <- metadat.Res499_PRRKO
use.pal <- pal_jama()(7)
names(use.pal) <- levels(metadat.Res499_PRRKO)
p2 <- plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(unique(gr$transcript_id)),
    geneset_names = paste0("TE-Derived Transcripts\n(n=", length(unique(gr$transcript_id)), ")"),
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = FALSE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    point_size = 3.5,
    # stroke = 1.5,
    plot.legend = F)

# B16 Res499
mat.plot <- edat.B16_R499
use.metadat <- metadat.B16_R499
use.pal <- pal_nejm()(4)
names(use.pal) <- levels(metadat.B16_R499$Label)
p3 <- plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(unique(gr$transcript_id)),
    geneset_names = paste0("TE-Derived Transcripts\n(n=", length(unique(gr$transcript_id)), ")"),
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    point_size = 3.5,
    # stroke = 1.5,
    plot.legend = F)

p1 + p2 + p3

# Random control
set.seed(3)
plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(rownames(mat.plot)[sample(1:nrow(mat.plot), 156)]),
    geneset_names = "Random Transcripts",
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Class",
    colorBy = "Class",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    # point_size = 3.5,
    # stroke = 1.5,
    plot.legend = F)

######################################################
### Epigenetic profiling at TE-derived transcripts ###
######################################################

# H3K4me3
ol <- findOverlaps(gr.H3K4me3, gr, ignore.strand = T)
print(all(colnames(mat.H3K4me3) == metadata.H3K4me3$Sample))
p.H3K4me3 <- plot_summary(
	mat.plot = mat.H3K4me3,
	geneset_list = list(gr.H3K4me3$ID[unique(from(ol))]),
	geneset_names = paste0("H3K4me3 with Transcripts\n(n = ", length(gr.H3K4me3$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
	metadat = metadata.H3K4me3,
	groupBy = "Label",
	colorBy = "Label",
	use.pal = pal_nejm()(6),
	show.xaxis.labels = FALSE,
	plot.legend = FALSE)

# H3K27ac
ol <- findOverlaps(gr.H3K27ac, gr, ignore.strand = T)
print(all(colnames(mat.H3K27ac) == metadata.H3K27ac$Sample))
p.H3K27ac <- plot_summary(
	mat.plot = mat.H3K27ac,
	geneset_list = list(gr.H3K27ac$ID[unique(from(ol))]),
	geneset_names = paste0("H3K27ac with Transcripts\n(n = ", length(gr.H3K27ac$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
	metadat = metadata.H3K27ac,
	groupBy = "Label",
	colorBy = "Label",
	use.pal = pal_nejm()(6),
	show.xaxis.labels = FALSE,
	plot.legend = FALSE)

# H3K4me1
ol <- findOverlaps(gr.H3K4me1, gr, ignore.strand = T)
print(all(colnames(mat.H3K4me1) == metadata.H3K4me1$Sample))
print(length(gr.H3K4me1$ID[unique(from(ol))]))
p.H3K4me1 <- plot_summary(
    mat.plot = mat.H3K4me1,
    geneset_list = list(gr.H3K4me1$ID[unique(from(ol))]),
    geneset_names = paste0("H3K4me1 with Transcripts\n(n = ", length(gr.H3K4me1$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
    metadat = metadata.H3K4me1,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(6),
	show.xaxis.labels = FALSE,
	plot.legend = FALSE)

# ATAC
ol <- findOverlaps(gr.ATAC_revisions, gr, ignore.strand = T)
print(all(colnames(mat.ATAC_revisions) == metadata.ATAC_revisions$Sample))
print(length(gr.ATAC_revisions$ID[unique(from(ol))]))
p.ATAC_revision <- plot_summary(
	mat.plot = mat.ATAC_revisions,
	geneset_list = list(gr.ATAC_revisions$ID[unique(from(ol))]),
	geneset_names = paste0("ATAC with Transcripts\n(n = ", length(gr.ATAC_revisions$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
	metadat = metadata.ATAC_revisions,
	groupBy = "Label",
	colorBy = "Label",
	use.pal = pal_nejm()(6),
    method = "gsva",
	show.xaxis.labels = FALSE,
	plot.legend = FALSE)

# PRR KO ATAC
ol <- findOverlaps(gr.Res499_PRRKO, gr, ignore.strand = T)
print(all(colnames(adat.Res499_PRRKO) == metadat.Res499_PRRKO_ATAC$Sample))
p.ATAC_PRRKO <- plot_summary(
	mat.plot = adat.Res499_PRRKO,
	geneset_list = list(gr.Res499_PRRKO$ID[unique(from(ol))]),
	geneset_names = paste0("ATAC with Transcripts\n(n = ", length(gr.Res499_PRRKO$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
	metadat = metadat.Res499_PRRKO_ATAC,
	groupBy = "Label",
	colorBy = "Label",
	use.pal = colorRampPalette(pal_jama()(7))(length(levels(metadat.Res499_PRRKO_ATAC$Label))),
    method = "gsva",
    aspect.ratio = 0.8,
	show.xaxis.labels = FALSE,
	plot.legend = FALSE)

p.H3K4me3 + p.H3K27ac + p.H3K4me1 + p.ATAC_revision + p.ATAC_PRRKO + plot_layout(nrow = 1)

#####################################
### Subfamily enrichment (GIGGLE) ###
#####################################

gr.rmsk_TE <- import.bed("~/Dropbox/Minn/resources/TE/squire/mm10_all.bed")
gr.rmsk_TE$ID <- sapply(strsplit(gr.rmsk_TE$name, split="\\|"), function(x) x[[4]])
gr.rmsk_TE$subF <- sapply(strsplit(gr.rmsk_TE$ID, split="\\:"), function(x) x[[1]])
gr.rmsk_TE$Family <- sapply(strsplit(gr.rmsk_TE$ID, split="\\:"), function(x) x[[2]])

TE.families <- names(table(gr.rmsk_TE$Family)[table(gr.rmsk_TE$Family) > 8])
TE.families <- TE.families[!grepl(TE.families, pattern="\\?|Other")]

load_GIGGLE_Results <- function(
	filename,
	pattern,
	sig = 0.5,
	combo_score = 0) {
	results <- read.table(filename, sep="\t", header=F)
	colnames(results) <- c("file_name", "file_size", "overlaps", "odds_ratio", "fishers_two_tail", "fishers_left_tail", "fishers_right_tail", "combo_score")
	results <- results[sort.int(results$combo_score, decreasing=T, index.return=T)$ix, ]
	results$Name <- gsub("giggle_index/regions_bed/sorted/|.bed.gz", "", results$file_name)
	results$subF <- results$Name
	results$ID <- gr.rmsk_TE$ID[match(results$subF, gr.rmsk_TE$subF)]
	results$ID[is.na(results$ID)] <- results$Name[is.na(results$ID)]
	results$Family <- gr.rmsk_TE$Family[match(results$subF, gr.rmsk_TE$subF)]
	results$Family[is.na(results$Family)] <- results$Name[is.na(results$Family)]
	print(head(results[grep(results$file_name, pattern=pattern), ], 25))
	results.sig <- results[which(results$fishers_two_tail < sig & results$combo_score > combo_score), ]
	return(results.sig)
}

GIGGLE.IFN_IMD <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_IMD_Exons.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.Res499_vs_B16 <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_Res499_vs_B16_Exons.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

# GIGGLE.B16y_3wk_vs_B16 <- load_GIGGLE_Results(
# 	filename = "~/Dropbox/Minn/B16_B16y_R499/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_20_B16y_3wk_vs_B16.txt",
# 	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
# 	sig = 0.1,
# 	combo_score = 5)

GIGGLE.Res499_Activated_PRR_Dependent <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_Res499_PRRdep_Exons.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.Res499_Activated_TBKi_Dependent <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_Cmpd1_RUX_Exons.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.KPy_vs_KP <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_KPy_Transcripts.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.ResKP2_vs_KP <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_ResKP2_Transcripts.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.ResResKP_vs_KP <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_ResResKP_Transcripts.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.LLC1y_vs_LLC1 <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/LLC1_LLC1y_ResLLC1/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_LLC1y_Transcripts.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.ResResLLC1_vs_LLC1 <- load_GIGGLE_Results(
	filename = "~/Dropbox/Minn/LLC1_LLC1y_ResLLC1/data/RNA/GRCm38_Build/GIGGLE/Results/TE_Transcripts_ResResLLC1_Transcripts.txt",
	pattern = "MMERVK|RLTR4|B2_Mm2|IAPEY3|IAPEY4|IAP|B2|B4",
	sig = 0.1,
	combo_score = 5)

GIGGLE.list <- list(
	IFN_IMD = GIGGLE.IFN_IMD,
	Res499_vs_B16 = GIGGLE.Res499_vs_B16,
	# B16y_3wk_vs_B16 = GIGGLE.B16y_3wk_vs_B16,
	Res499_Activated_PRR_Dependent = GIGGLE.Res499_Activated_PRR_Dependent,
	Res499_Activated_TBKi_Dependent = GIGGLE.Res499_Activated_TBKi_Dependent,
	KPy_vs_KP = GIGGLE.KPy_vs_KP,
	ResKP2_vs_KP = GIGGLE.ResKP2_vs_KP,
	ResResKP_vs_KP = GIGGLE.ResResKP_vs_KP,
	# LLC1y_vs_LLC1 = GIGGLE.LLC1y_vs_LLC1,
	# ResLLC1A_vs_LLC1 = GIGGLE.ResLLC1A_vs_LLC1,
	# ResLLC1B_vs_LLC1 = GIGGLE.ResLLC1B_vs_LLC1,
	ResResLLC1_vs_LLC1 = GIGGLE.ResResLLC1_vs_LLC1)
sapply(GIGGLE.list, function(x) nrow(x))

# GIGGLE.list_Filtered <- lapply(1:length(GIGGLE.list), function(x) {
# 	giggle <- GIGGLE.list[[x]]
# 	if(nrow(giggle) >= 25) {
# 		giggle <- giggle[1:25, ]
# 	}
# 	return(giggle)
# })
# names(GIGGLE.list_Filtered) <- names(GIGGLE.list)
# GIGGLE.list <- GIGGLE.list_Filtered

te.plot <- names(sort(table(unlist(sapply(GIGGLE.list, function(x) x$Name)))))[sort(table(unlist(sapply(GIGGLE.list, function(x) x$Name)))) > 2]

mat.plot <- lapply(1:length(GIGGLE.list), function(x) {
	print(x)
	giggle <- GIGGLE.list[[x]]
	giggle <- giggle[match(te.plot, giggle$Name), c("ID", "subF", "Family", "fishers_two_tail", "combo_score", "file_size", "overlaps", "odds_ratio")]
	giggle$Name <- te.plot

	giggle$subF <- giggle$Name
	giggle$ID <- gr.rmsk_TE$ID[match(giggle$subF, gr.rmsk_TE$subF)]
	giggle$ID[is.na(giggle$ID)] <- giggle$Name[is.na(giggle$ID)]
	giggle$Family <- gr.rmsk_TE$Family[match(giggle$subF, gr.rmsk_TE$subF)]
	giggle$Family[is.na(giggle$Family)] <- giggle$Name[is.na(giggle$Family)]
	giggle$Regions <- names(GIGGLE.list)[x]
	return(giggle)
})
mat.plot <- do.call(rbind, mat.plot)
mat.plot$Name <- factor(mat.plot$Name, levels = te.plot)
mat.plot$Regions <- factor(mat.plot$Regions, levels = names(GIGGLE.list))

mat.plot1 <- mat.plot
# mat.plot1 <- mat.plot[mat.plot$Family == "ERVK", ]
# mat.plot1 <- mat.plot[mat.plot$Name %in% TE.families, ]

ggplot(mat.plot1, aes(x = Regions, y = ID, color = combo_score, size = combo_score)) +
    # geom_point(alpha=0.8, size = 5) +
    geom_point(alpha=0.8) +
    scale_color_gradientn(
        colors=brewer.pal(9, "Reds")[3:9], 
        # limits = c(0,3), 
        oob = scales::squish) +
    labs(color="log2FoldChange", y = "", x = "", title = "Resistance") +
    guides(size = "none") +
    scale_size(range = c(2,8)) +
    theme_bw(base_size=12) +
    theme(
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size=12, color = "black", hjust = 1, angle = 45),
        axis.ticks = element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

####################
### Repeat pairs ###
####################

# Repeat pairs
repeat_pairs_list <- readRDS("~/Dropbox/Minn/scATAC_B16_R499_Stat1KO/resistant_clones/data/RNA/GRCm38_Build/TE_Derived_Transcripts/StringTie_TACO_Assembly_NoProteinCoding_Repeat_Pairs.rds")

scores_list <- sapply(1:length(repeat_pairs_list), function(x) {
    names <- names(repeat_pairs_list[[x]])
    scores <- as.numeric(sapply(strsplit(names, split="_"), function(x) x[[length(x)]]))
    return(scores)
})
names(scores_list) <- names(repeat_pairs_list)

lapply(scores_list, function(x) quantile(x, probs=seq(0,1,0.1)))

alignment_threshold <- 50
repeat_pairs_list_Filtered <- lapply(1:length(repeat_pairs_list), function(x) {
    rp <- repeat_pairs_list[[x]]
    rp <- rp[scores_list[[x]] > alignment_threshold]
    return(rp)
})
names(repeat_pairs_list_Filtered) <- names(repeat_pairs_list)
sort(sapply(repeat_pairs_list_Filtered, function(x) length(x)))
repeat_pairs_list_Filtered <- repeat_pairs_list_Filtered[as.numeric(sapply(repeat_pairs_list_Filtered, function(x) length(x))) != 0]
repeat_pairs_list_Filtered <- lapply(repeat_pairs_list_Filtered, function(x) GRangesList(x))

# Check assembled transcripts for repeat pairs
gr.interest <- gr.Res499_vs_B16_20
gr <- gr.TE_Transcripts_20_Extended[gr.TE_Transcripts_20_Extended$transcript_id %in% gr.interest$transcript_id]
gr <- gr[gr$type == "transcript"]

gr.rp <- lapply(1:length(gr), function(x) {
	gr.t <- gr[x]
	rp.overlap <- lapply(1:length(repeat_pairs_list_Filtered), function(y) {
		rp <- repeat_pairs_list_Filtered[[y]]
		check_rp_overlaps <- overlapsAny(rp, gr.t)
		return(rp[check_rp_overlaps])
	})
	rp.overlap <- rp.overlap[sapply(rp.overlap, function(x) length(x)) > 0]
	return(rp.overlap)
})
names(gr.rp) <- gr$transcript_id

rp.ids_Family <- sapply(1:length(gr.rp), function(x) {
	gr.list <- gr.rp[[x]]
	if(length(gr.list) > 0) {
		gr.list <- do.call(c, lapply(1:length(gr.list), function(y) {
			tmp <- gr.list[[y]]
			names(tmp) <- NULL
			return(tmp)
		}))
		rp.ids <- paste(sapply(gr.list, function(x) unique(x$Family)), collapse="|")
	} else {
		rp.ids <- NA
	}
	return(rp.ids)
})

rp.ids_Subfamily <- sapply(1:length(gr.rp), function(x) {
	gr.list <- gr.rp[[x]]
	if(length(gr.list) > 0) {
		gr.list <- do.call(c, lapply(1:length(gr.list), function(y) {
			tmp <- gr.list[[y]]
			names(tmp) <- NULL
			return(tmp)
		}))
		rp.ids <- paste(sapply(gr.list, function(x) paste0(unique(x$subF), collapse=",")), collapse="|")
	} else {
		rp.ids <- NA
	}
	return(rp.ids)
})

gr$repeat_pairs_Family <- rp.ids_Family
gr$repeat_pairs_Subfamily <- rp.ids_Subfamily

data.frame(gr[grepl(gr$repeat_pairs_Family, pattern="B2")])

#################
### Visualize ###
#################

gr <- gr.ResResLLC1_vs_LLC1_40

# RNA
mat.plot <- edat
use.metadat <- metadat
use.pal <- colorRampPalette(brewer.pal(9, "Greens")[2:8])(5)
names(use.pal) <- levels(metadat$Label)
p.rna <- plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(unique(gr$transcript_id)),
    geneset_names = paste0("TE-Derived Transcripts\n(n=", length(unique(gr$transcript_id)), ")"),
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    point_size = 2.5,
    # stroke = 1.5,
    plot.legend = T)

# ATAC
gr.Extended <- extend_regions(gr, extend_start = 2500, extend_end = 2500)

mat.plot <- adat
use.metadat <- metadat_ATAC
use.pal <- colorRampPalette(brewer.pal(9, "Purples")[2:8])(5)
names(use.pal) <- levels(metadat_ATAC$Label)

ol <- findOverlaps(gr.atac, gr.Extended, ignore.strand = T)
print(all(colnames(adat) == use.metadat$Sample))

p.atac <- plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(unique(gr.atac$ID[unique(from(ol))])),
    geneset_names = paste0("ATAC with Transcripts\n(n = ", length(gr.atac$ID[unique(from(ol))]), "/", length(unique(gr$transcript_id)), ")"),
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    point_size = 2.5,
    # stroke = 1.5,
    plot.legend = T)

p.rna + p.atac

# Random
set.seed(2)
plot_summary(
    mat.plot = mat.plot,
    geneset_list = list(rownames(mat.plot)[sample(1:nrow(mat.plot), 156)]),
    geneset_names = "Random Transcripts",
    metadat = use.metadat[match(colnames(mat.plot), use.metadat$Sample), ],
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal,
    method = "gsva",
    plot.boxplot = TRUE,
    show.xaxis.labels = FALSE,
    font_size = 14,
    aspect.ratio = 1,
    # point_size = 3.5,
    # stroke = 1.5,
    plot.legend = T)

##############################
### TE-derived transcripts ###
##############################

# gr <- import("~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_ResResKP_vs_KP.gtf")
# gr <- import("~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_ResKP2_vs_KP.gtf")
# gr <- import("~/Dropbox/Minn/KP_KPy_ResKP/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_KPy_vs_KP.gtf")
# gr <- import("~/Dropbox/Minn/LLC1_LLC1y_ResLLC1/data/RNA/GRCm38_Build/TE_Derived_Transcripts/TE_Transcripts_50_ResResLLC1_vs_LLC1.gtf")
