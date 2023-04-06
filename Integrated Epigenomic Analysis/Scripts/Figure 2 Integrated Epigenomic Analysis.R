rm(list=ls())

library(GenomicRanges)
library(rtracklayer)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(matrixStats)
library(DiffBind)
library(ggrastr)

directory <- "local"

source("~/Dropbox/Minn/resources/useful_protocols/Bulk ATAC processing/ATAC Visualization and Analysis Functions.R")
source("~/Dropbox/Minn/ifnar_epigenome/scripts/Integrated Epigenomic Analysis/Integrated Epigenomic Analysis Functions.R")

# ###############################
# ### Load in reference files ###
# ###############################

# # Gene annotations
# bm <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/mouse_biomart_annotations_2-20-20.txt", sep="\t", header=T, stringsAsFactors=F)
# remove_idx <- grep(bm$chromosome_name, pattern="MT|GL|JH|X|Y")
# bm <- bm[-remove_idx, ]
# bm$strand <- ifelse(bm$strand == 1, "+", "-")
# gr.bm <- makeGRangesFromDataFrame(bm, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")

# # TSS for protein-coding genes
# gr.tss <- makeGRangesFromDataFrame(data.frame(
# 	chr=as.character(seqnames(gr.bm)), 
# 	start=gr.bm$transcription_start_site, 
# 	end=gr.bm$transcription_start_site, 
# 	gene=gr.bm$mgi_symbol, 
# 	strand = as.character(strand(gr.bm))), keep.extra.columns = T)
# gr.tss_window <- gr.tss
# start(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) <- start(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) - 2500
# end(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) <- end(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) + 1000
# start(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) <- start(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) - 1000
# end(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) <- end(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) + 2500

# # Gene sets
# gs.names <- c("IFN.I", "ISG.RS")
# gs_list <- lapply(1:length(gs.names), function(x) {
#   gs <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/", gs.names[x], "_mouse.txt"), sep="\t", header=F)
#   gs <- unique(as.character(gs$V1))
#   gs <- as.character(bm$ensembl_gene_id[na.omit(match(gs, bm$mgi_symbol))])
#   return(gs)
# })
# names(gs_list) <- gs.names

##################################
### Load processed data files ####
##################################

marks <- c("H3K4me3", "H3K27Ac", "H3K4me1")
assays <- c(marks, "ATAC")
RE.names <- c(
	"gained_TSS_promoters", "gained_other_promoters", "gained_other_promoters_ATAC", 
	"lost_TSS_promoters", "lost_other_promoters", "lost_other_promoters_ATAC", 
	"activated_enhancer", "activated_enhancer_ATAC", "activated_enhancer_ATAC_inclusive", 
	"deactivated_enhancer", "deactivated_enhancer_ATAC_inclusive",
	"unchanged_REs", "random_peaks")

### ANNOTATED REGULATORY ELEMENT ANNOTATIONS ###

print("Load annotated RE sets")

RE_list <- lapply(1:length(RE.names), function(x) import.bed(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/", RE.names[x], ".bed")))
names(RE_list) <- RE.names

RE.repertoire.names <- c("Promoter", "Active_enhancer_H3K27Ac", "Active_enhancer_H3K4me1")
RE_repertoire_list <- lapply(1:length(RE.repertoire.names), function(x) readRDS(file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/", RE.repertoire.names[x], "_RE.rds")))

### PAIRED RNA/ATAC DATA ###

print("Load paired data")

# RNA data
rna.dat <- read.table(file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/RNA_counts_rlog.txt"), header=T, sep="\t", stringsAsFactors=F)
rna.dat <- rna.dat[!is.na(match(rownames(rna.dat), bm$ensembl_gene_id)), ]

# ATAC data
atac.dat_idr <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_original_consensus_mat_tn5_insertion_counts_IDR_rlog.txt"), sep="\t", stringsAsFactors=F, header=T)

# Match ATAC and RNA samples
edat <- rna.dat
adat <- atac.dat_idr[ ,intersect(colnames(rna.dat), colnames(atac.dat_idr))]

# Remove B16_SKO_3 poor quality sample (low FRiP, not consistent with other replicates)
remove_samples <- c("B16_SKO_3")
edat <- edat[ ,!colnames(edat) %in% remove_samples]
adat <- adat[ ,!colnames(adat) %in% remove_samples]
colnames(edat) <- gsub("cas", "WT", colnames(edat))
colnames(adat) <- gsub("cas", "WT", colnames(adat))

# Original paired RNA/ATAC metadata
metadat <- data.frame(
	Sample = colnames(edat),
	Condition = factor(sapply(strsplit(colnames(edat), split="_"), function(x) paste(x[1:2], collapse="_")), levels=c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO")))
levels(metadat$Condition) <- c("B16 WT", "B16 SKO", "Res499 WT", "Res499 SKO")

### ATAC/H3K27ac/H3K4me1 COUNT MATRICES ###

print("Load normalized count matrices")

mat_list <- lapply(1:length(assays), function(x) {
	if(assays[x] == "ATAC") {
		mat <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_original_consensus_mat_tn5_insertion_counts_IDR_rlog.txt", sep="\t", header=T, stringsAsFactors=F)
		# mat <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/consensus_mat_tn5_insertion_counts_IDR_vst.txt", sep="\t", header=T, stringsAsFactors=F)
		colnames(mat) <- gsub("cas", "WT", colnames(mat))
	} else {
		mat <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/", assays[x], "_consensus_mat_insertion_counts_WT_pooled_vst.txt"), sep="\t", header=T, stringsAsFactors=F)
	}
	return(mat)
})
names(mat_list) <- assays
mat_list[["RNA"]] <- edat

# Revisions ATAC data
mat_list1 <- mat_list
mat_list1[["ATAC"]] <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_revisions_consensus_mat_tn5_insertion_counts_IDR_vst.txt"), sep="\t", stringsAsFactors=F, header=T)
colnames(mat_list1[["ATAC"]]) <- gsub("\\.", "_", colnames(mat_list1[["ATAC"]]))

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

### Differential peaks ###

db.DE <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/ATAC_revisions_IDR_Diffbind_DE_all.rds")
db.DE$config$AnalysisMethod <- factor("edgeR")

#################################
### Fig 2B: Summary PCA plots ###
#################################

p.rna <- plot_PCA(
	dat = edat, 
	title = "RNA", 
	n_var_features = nrow(edat)/3, 
	annotation = metadat$Condition, 
	color.pal = pal_nejm()(4))

p.atac <- plot_PCA(
	dat = adat, 
	title = "ATAC", 
	n_var_features = nrow(adat)/3, 
	annotation = metadat$Condition, 
	color.pal = pal_nejm()(4))

p.combined <- ggarrange(p.rna, p.atac, common.legend = TRUE, legend="right", ncol=2)
ggsave(p.combined, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2B PCA plots.pdf", device="pdf", height=8, width=12)

# # Write out plot data
# dat <- rbind(p.rna$data[,c("PC1", "PC2", "sample", "anno")], p.atac$data[,c("PC1", "PC2", "sample", "anno")])
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2B.csv"), sep=",", quote=F, col.names=T, row.names=F)

###############################################
### Fig 2C: Differential ATAC peaks MA plot ###
###############################################

sig <- 0.05

# Use all peaks differentially accessible between B16 and R499
db.DE_B16_v_R499 <- data.frame(dba.report(db.DE, contrast=2, bUsePval=TRUE, th=1))

f <- ggplot(db.DE_B16_v_R499, aes(Conc, Fold)) +
    geom_bin2d(bins = 300) +
    scale_fill_gradientn(colours = "grey75") +
    geom_hline(yintercept=0, linetype="dashed", color="black", size=0.6) +
    geom_hline(yintercept=1, linetype="dashed", color="black", size=0.4) +
    geom_hline(yintercept=-1, linetype="dashed", color="black", size=0.4) +
    geom_point(data=db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR < sig), ], col="orangered", size=0.5) + 
    ylim(c(-5,5)) +
    labs(x = "Log2 accessibility", y = "Log2 fold change",
    	title = paste0("Differentially accessible regions \n between ", strsplit("B16_v_R499", split="_")[[1]][1], " and ", strsplit("B16_v_R499", split="_")[[1]][3], " (", nrow(db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR < sig), ]), " FDR < ", sig, ")")) +
    theme_bw() +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=18),
          legend.position="none",
          text=element_text(size=20),
          legend.text=element_text(size=12), legend.title=element_text(size=12),
          axis.text=element_text(size=18, color="black"),
          plot.margin=unit(c(0.5,0.5,0.75,0.5), "cm"),
          aspect.ratio=1,
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
ggsave(f, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2C MA plot B16 v R499 ATAC.pdf", width=6, height=6)

# # Write out plot data
# dat <- f$data[ ,c("Conc", "Fold", "FDR")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2C.csv"), sep=",", quote=F, col.names=T, row.names=F)

################################################
### Fig 2D: Rank genes by variance explained ###
################################################

mRFARobj_ALL <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/mRF/mRFARobj_ALL_genes.rds")

Rsq <- sapply(mRFARobj_ALL, function(x) x$Rsq_rf)
df.Rsq <- data.frame(gene=sapply(mRFARobj_ALL, function(x) x$geneSymb), Rsq=Rsq)
df.Rsq <- na.omit(df.Rsq)
df.Rsq <- df.Rsq[sort.int(df.Rsq$Rsq, decreasing=T, index.return=T)$ix, ]
df.Rsq$rank <- 1:nrow(df.Rsq)
df.Rsq$label <- ifelse(df.Rsq$gene %in% c("Oas1a", "Oas1g"), "firebrick", "black")
df.Rsq$size <- ifelse(df.Rsq$gene %in% c("Oas1a", "Oas1g"), 6, 1)
fig <- ggplot(df.Rsq, aes(x=rank, Rsq)) +
    geom_point_rast(fill=df.Rsq$label, col=df.Rsq$label, size=df.Rsq$size) +
    geom_label_repel(data=df.Rsq[grep(df.Rsq$gene, pattern="Oas1"), ], aes(label=gene), point.padding=0.2, segment.size=1, segment.color="firebrick", col="firebrick", fontface="bold", size=5) +
    ylim(c(0,1)) +
    xlab("Rank") +
    ylab("Variance explained") +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 12),
          aspect.ratio = 1,
          text=element_text(size=14))
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2D All genes variance explained ranked.pdf", width=4, height=3)

# # Write out plot data
# dat <- fig$data[ ,c("rank", "gene", "Rsq")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2D.csv"), sep=",", quote=F, col.names=T, row.names=F)

##################################################
### Fig 2F: Summary profiles for annotated REs ###
##################################################

# Run computeMatrix (DeepTools) to calculate ATAC/H3K27ac signals at annotated RE regions
# ~/Dropbox/Minn/ifnar_epigenome/scripts/make_summary_profile_plots.sh

library(data.table)

cond1 <- "B16_WT"
cond2 <- "R499_WT"

pal.colors <- pal_nejm()(4)
names(pal.colors) <- c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO")

p_theme <- theme_minimal(base_size=11) +
	theme(aspect.ratio = 1,
		axis.text = element_text(size=11, color="black"),
		axis.title = element_text(size=11, color="black"),
		panel.border = element_rect(color = "black", fill=NA, size=1),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

RE.name <- "activated_enhancer"

mat.ATAC <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/summary_profiles/", RE.name, ".mat_ATAC_", cond1, "v", cond2, ".gz"), skip=1)
mat.ATAC.signal <- data.frame(t(mat.ATAC[ ,7:ncol(mat.ATAC)]))
colnames(mat.ATAC.signal) <- mat.ATAC$V4

mat.H3K27Ac <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/summary_profiles/", RE.name, ".mat_H3K27Ac_", cond1, "v", cond2, ".gz"), skip=1)
mat.H3K27Ac.signal <- data.frame(t(mat.H3K27Ac[ ,7:ncol(mat.H3K27Ac)]))
colnames(mat.H3K27Ac.signal) <- mat.H3K27Ac$V4

# mat.H3K4me1 <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, ".mat_H3K4me1_", cond1, "v", cond2, ".gz"), skip=1)
# mat.H3K4me1.signal <- data.frame(t(mat.H3K4me1[ ,7:ncol(mat.H3K4me1)]))
# colnames(mat.H3K4me1.signal) <- mat.H3K4me1$V4

# mat.H3K4me3 <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, ".mat_H3K4me3_", cond1, "v", cond2, ".gz"), skip=1)
# mat.H3K4me3.signal <- data.frame(t(mat.H3K4me3[ ,7:ncol(mat.H3K4me3)]))
# colnames(mat.H3K4me3.signal) <- mat.H3K4me3$V4

df.Peaks <- data.frame(
	bp = c(seq(-5000,(5000-1),100), seq(-5000,(5000-1),100)),
	Condition = c(rep(cond1, 5000*2 / 100), rep(cond2, 5000*2 / 100)),
	Signal_H3K27Ac = rowMeans(mat.H3K27Ac.signal))
RE1 <- "H3K27Ac"

df.Peaks.ATAC <- data.frame(
	bp = c(seq(-2500,(2500-1),10), seq(-2500,(2500-1),10)),
	Condition = c(rep(cond1, 2500*2 / 10), rep(cond2, 2500*2 / 10)),
	Signal_ATAC = rowMeans(mat.ATAC.signal))
df.Peaks.ATAC <- df.Peaks.ATAC[df.Peaks.ATAC$bp %in% -1000:1000, ]

p.H3K27ac <- ggplot(df.Peaks, aes(x = bp, y = Signal_H3K27Ac, color = Condition)) +
	geom_line(size = 0.5) +
	scale_color_manual(values = pal.colors[match(c(cond1, cond2), names(pal.colors))]) +
	facet_wrap(~Condition) +
	# ylim(c(0.2,2)) +
	labs(title = paste0(gsub("_", " ", RE.name), "\nH3K27Ac"), 
		x = "Distance from TSS (bp)", 
		y = "H3K27Ac Signal") +
	p_theme
p.ATAC <- ggplot(df.Peaks.ATAC, aes(x = bp, y = Signal_ATAC, color = Condition)) +
	geom_line(size = 0.5) +
	scale_color_manual(values = pal.colors[match(c(cond1, cond2), names(pal.colors))]) +
	facet_wrap(~Condition) +
	# ylim(c(0.2,2)) +
	labs(title = "ATAC", 
		x = "Distance from TSS (bp)", 
		y = "ATAC Signal") +
	p_theme

ggsave(p.H3K27ac, file=paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, "_summary_H3K27Ac.pdf"), width=5, height=3)
ggsave(p.ATAC, file=paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, "_summary_ATAC.pdf"), width=5, height=3)

# # Write out plot data
# dat1 <- p.H3K27ac$data[ ,c("bp", "Condition", "Signal_H3K27Ac")]
# dat2 <- p.ATAC$data[ ,c("bp", "Condition", "Signal_ATAC")]
# colnames(dat1) <- c("bp", "Condition", "Signal")
# dat1$Assay <- "H3K27ac"
# colnames(dat2) <- c("bp", "Condition", "Signal")
# dat2$Assay <- "ATAC"
# dat <- rbind(dat1, dat2)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2F.Summary.csv"), sep=",", quote=F, col.names=T, row.names=F)

### HEATMAPS ###

mat.plot <- t(mat.H3K27Ac.signal)
mat.plot <- mat.plot[sort.int(rowMeans(mat.plot), decreasing=T, index.return=T)$ix, ]
colnames(mat.plot) <- df.Peaks$bp
col_max <- as.numeric(quantile(mat.plot, probs=seq(0,1,0.01))[100])
hm.H3K27Ac <- Heatmap(mat.plot,
	name = "H3K27Ac",
	use_raster = TRUE,
	raster_quality = 1,
	cluster_columns = F,
	cluster_rows = F,
	show_column_names = F,
	show_row_names = F,
	col = colorRamp2(seq(0,col_max,length=12), colorRampPalette(brewer.pal(9, "Reds"))(12)),
	column_split = factor(c(rep(gsub("_", " ", cond1), ncol(mat.plot)/2), rep(gsub("_", " ", cond2), ncol(mat.plot)/2)), levels=c(gsub("_", " ", cond1), gsub("_", " ", cond2))),
	border = TRUE)
# dat <- mat.plot
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2F.H3K27ac.csv"), sep=",", quote=F, col.names=T, row.names=F)

mat.plot <- t(scale(mat.ATAC.signal))
mat.plot <- mat.plot[sort.int(rowMeans(mat.plot), decreasing=T, index.return=T)$ix, ]
colnames(mat.plot) <- df.Peaks.ATAC$bp
col_max <- as.numeric(quantile(mat.plot, probs=seq(0,1,0.01))[100])
hm.ATAC <- Heatmap(mat.plot,
	name = "ATAC",
	use_raster = TRUE,
	raster_quality = 1,
	cluster_columns = F,
	cluster_rows = F,
	show_column_names = F,
	show_row_names = F,
	col = colorRamp2(seq(-5,5,length=12), colorRampPalette(rev(brewer.pal(9, "RdBu")))(12)),
	column_split = factor(c(rep(gsub("_", " ", cond1), ncol(mat.plot)/2), rep(gsub("_", " ", cond2), ncol(mat.plot)/2)), levels=c(gsub("_", " ", cond1), gsub("_", " ", cond2))),
	border = TRUE)
# dat <- mat.plot
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2F.ATAC.csv"), sep=",", quote=F, col.names=T, row.names=F)

pdf(paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, " H3K27Ac Heatmap.pdf"), width=3, height=5)
draw(hm.H3K27Ac)
dev.off()
pdf(paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, " ATAC Heatmap.pdf"), width=3, height=5)
draw(hm.ATAC)
dev.off()

###############################################################
### Fig 2H: Gene expression in cis-regulatory window of REs ###
###############################################################

# Genes within 50-500 Kb window of activated enhancers have increased gene expression compared to expected
# GO on genes within 50-500 Kb window of activated enhancers -> IFN pathway genes enriched

### IDENTIFY CIS-GENES AND PLOT EXPRESSION ###

names <- c("gained_TSS_promoters", "activated_enhancer_ATAC_inclusive", "deactivated_enhancer_ATAC_inclusive", "unchanged_REs")
plot.conditions <- c("B16_WT", "R499_WT")

mat_avg_list <- lapply(1:length(names), function(x) {
	RE.idx <- names[x]
	if(grepl(RE.idx, pattern="promoter")) {
		W <- 5000
	} else {
		W <- 50000
	}
	gr.RE <- RE_list[[RE.idx]]

	# Get genes near annotated REs
	gr.RE_window <- gr.RE
	start(gr.RE_window) <- start(gr.RE_window) - W
	end(gr.RE_window) <- end(gr.RE_window) + W
	ol <- findOverlaps(gr.RE_window, gr.bm)
	genes.ol <- unique(gr.bm[unique(to(ol))]$ensembl_gene_id)

	# Expressed genes only
	genes.ol <- genes.ol[genes.ol %in% rownames(edat)]

	# GSVA
	gsva <- gsva(
		expr=as.matrix(edat[,grep(colnames(edat), pattern=paste(plot.conditions, collapse="|"))]), 
		gset.idx.list=list(genes.ol), 
		kcdf="Gaussian", 
		min.sz=10, 
		max.sz=20000)
	mat <- data.frame(Score=as.numeric(gsva))
	mat$Sample <- factor(colnames(gsva), levels=colnames(gsva))
	mat$Condition <- sapply(strsplit(as.character(mat$Sample), split="_"), function(x) paste(x[1:2], collapse="_"))
	mat$Condition <- factor(gsub("_", " ", mat$Condition), levels=gsub("_", " ", plot.conditions))
	mat$RE <- RE.idx

	# Write out cis-genes for GO analysis
	genes.ol_mgi <- gr.bm$mgi_symbol[match(genes.ol, gr.bm$ensembl_gene_id)]
	# write.table(genes.ol_mgi, file=paste0("~/Dropbox/Minn/R499_CnR/data/run_FINAL/cis_genes/", RE.idx, "_cis_genes2.txt"), quote=F, row.names=F, col.names=F)
	return(mat)
})
mat <- do.call(rbind, mat_avg_list)
mat$RE <- gsub("_ATAC_inclusive", "", mat$RE)
mat$RE <- gsub("_", " ", mat$RE)

p <- ggplot(mat, aes(x=Condition, y=Score, color=Condition)) +
	geom_point(size=4, shape=1, stroke=1.5) +
	stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
	facet_wrap(~RE, ncol=2) +
	scale_color_manual(values=pal_nejm()(4)[c(1,3)]) +
	labs(y="Expression score", x="") + 
	ylim(c(min(mat$Score)-0.02, max(mat$Score)+0.02)) +
	theme_bw() +
	theme(aspect.ratio=0.8,
		axis.text.y = element_text(size=12, color="black"),
		axis.text.x = element_text(size=14, color="black", hjust=1, angle=45),
		text=element_text(size=14),
		strip.text.x = element_text(margin = margin(6,0,6,0, "pt"), size=14),
		strip.background=element_rect(color="black", fill="white", size=0.8),
		legend.position="bottom")
ggsave(p, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2H Annotated REs Cis Genes Expression.pdf"), width=5.5, height=5)

# dat <- mat
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2H.csv"), sep=",", quote=F, col.names=T, row.names=F)

### Calculate significance values ###

t.test(mat_avg_list[[1]]$Score[1:3], mat_avg_list[[1]]$Score[4:6])
t.test(mat_avg_list[[2]]$Score[1:3], mat_avg_list[[2]]$Score[4:6])
t.test(mat_avg_list[[3]]$Score[1:3], mat_avg_list[[3]]$Score[4:6])
t.test(mat_avg_list[[4]]$Score[1:3], mat_avg_list[[4]]$Score[4:6])

###############################################################
### Extended Data Figure 3A: RNA/ATAC Pairwise Correlations ###
###############################################################

plot_pairwise_correlations <- function(dat, anno, n_features=2500, title) {

	# Get variable features
	rv <- rowVars(as.matrix(dat))
	select <- order(rv, decreasing=TRUE)[seq_len(min(n_features, length(rv)))]

	# Calculate spearman correlation
	cor.mat <- cor(dat[select, ], method="spearman")
	df.anno <- data.frame(anno=anno)
	rownames(df.anno) <- colnames(dat)
	anno_colors <- list(Anno=pal_npg()(length(levels(anno))))
	names(anno_colors$Anno) <- levels(anno)
	hm <- pheatmap(cor.mat, color=brewer.pal(9, "YlOrRd"), annotation_row=df.anno, annotation_colors=anno_colors, main=title, cellwidth=15, cellheight=15, angle="45", legend=T)

	# # Write out plot data
	# write.table(hm@matrix, file = paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 2A ", title, ".csv"), sep=",", quote=F, col.names=T, row.names=T)
	return(as.ggplot(hm))
}

hm1 <- plot_pairwise_correlations(dat=edat, anno=metadat$Condition, title="RNA", n_features=2500)
hm2 <- plot_pairwise_correlations(dat=adat, anno=metadat$Condition, title="ATAC", n_features=2500)
fig <- hm1 + hm2
ggsave(fig, file=paste0(dir.path, "results/summary/pairwise_correlations.pdf"), width=15, height=7)

######################################################################################
### Extended Data Figure 3H: GO enrichment on Res 499 activated enhancer cis-genes ###
######################################################################################

# Activated enhancers
res.go <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/activated_enhancer_ATAC_inclusive_GO_PANTHER.txt"), sep="\t", header=T, stringsAsFactors=F, skip=11)
colnames(res.go) <- gsub("upload_1..", "", colnames(res.go))
res.go <- res.go[-grep(res.go$PANTHER.Pathways, pattern="Unclassified"), ]
res.go$PANTHER.Pathways <- sapply(strsplit(res.go$PANTHER.Pathways, split="\\("), function(x) substr(x[[1]], 1, nchar(x[[1]])-1))
res.go$PANTHER.Pathways <- factor(res.go$PANTHER.Pathways, levels=rev(res.go$PANTHER.Pathways))
res.go <- res.go[res.go$FDR. < 0.01, ]

p <- ggplot(res.go, aes(x=PANTHER.Pathways, y=fold.Enrichment.)) +
	geom_bar(stat="identity", width=0.75, size=1, fill=pal_jama()(1), color=pal_jama()(1)) +
	coord_flip() +
	labs(x="", y="Fold enrichment") +
	theme_bw() +
	theme(aspect.ratio=0.8,
		legend.position="none",
		axis.text=element_text(size=14, color="black"),
		axis.title.x=element_text(size=14, color="black"),
		panel.grid.major=element_blank(),
		panel.grid.minor=element_blank(),
		plot.margin = margin(1,1,1,1, "cm"))
ggsave(p, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2H Activated Enhancer GO Enrichment PANTHER.pdf"), width=7.5, height=7.5)

# # Write out plot data
# dat <- res.go[,c("PANTHER.Pathways", "fold.Enrichment.", "FDR.")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 2H.csv"), sep=",", quote=F, col.names=T, row.names=F)

#####################################################################
### Fig. 2I: Acute vs. Chronic IFNG ATAC at Res499 activated loci ###
#####################################################################

# Reviewer 4 Chronic B16y Analysis.R

dir.path_ATAC <- "~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/"

metadat_ATAC <- read.csv("~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/metadata.csv", sep=",", header=T, stringsAsFactors=F)
metadat_ATAC$Label <- factor(metadat_ATAC$Condition, levels=c("B16", "B16y_6hr", "B16y_3wk", "Res499A"))
levels(metadat_ATAC$Label) <- c("B16", "B16y 6hr Acute IFN", "B16y 3.5wk Chronic IFN", "Res499")
use.pal_ATAC <- colorRampPalette(pal_nejm()(4)[c(1,3)])(4)
names(use.pal_ATAC) <- levels(metadat_ATAC$Label)

# Normalized counts
adat <- read.table("~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/mat/consensus_mat_tn5_insertion_counts_IDR_vst_all.txt", sep="\t", header=T, stringsAsFactors=F)
DARs_list <- readRDS("~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/diffbind/DARs_list_DESeq2_Treatment.rds")

adat <- adat[ ,match(metadat_ATAC$Sample, colnames(adat))]
gr.atac <- peakids2GRanges(rownames(adat), delim="_")
start(gr.atac) <- start(gr.atac) - 1
gr.atac$ID <- paste(as.character(seqnames(gr.atac)), start(gr.atac), end(gr.atac), sep="_")
rownames(adat) <- gr.atac$ID
print(all(rownames(adat) == gr.atac$ID))

# Motif deviations
mat.motifs <- readRDS("~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/chromvar/mat.motifs_Vierstra_annotations.rds")
mat.motifs_collapsed <- readRDS("~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/chromvar/mat.motifs_Vierstra_annotations_collapsed_average.rds")
mat.motifs <- mat.motifs[ ,match(metadat_ATAC$Sample, colnames(mat.motifs))]
mat.motifs_collapsed <- mat.motifs_collapsed[ ,match(metadat_ATAC$Sample, colnames(mat.motifs_collapsed))]

# Chronic IFNG Res499 enhancers
gr <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/Chronic IFNG Res499 Enhancers.rds")

# Summary
features_list <- list(gr$ID)
names(features_list) <- "DARs"
fig <- plot_summary(
    mat.plot = adat,
    geneset_list = features_list,
    geneset_names = "Res499 Chronic IFNG\nEnhancers (n = 900)",
    metadat = metadat_ATAC,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal_ATAC,
    method = "gsva",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 15,
    font_size_strip = 14,
    point_size = 5,
    stroke = 1.5,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")
ggsave(fig, file = "~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2I Chronic IFNG Res499 Enhancers Summary.pdf", width = 4.5, height = 5)

# Signficance test
df <- plot_summary(
    mat.plot = adat,
    geneset_list = features_list,
    geneset_names = "Res499 Chronic IFNG\nEnhancers (n = 900)",
    metadat = metadat_ATAC,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal_ATAC,
    method = "gsva",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 15,
    font_size_strip = 14,
    point_size = 5,
    stroke = 1.5,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "B16y 6hr Acute IFN"])
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "B16y 3.5wk Chronic IFN"])
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "Res499"])

# # Write out plot data
# dat <- df[ ,c("Sample", "Group", "Signal")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Figure 2I.csv"), sep=",", quote=F, col.names=T, row.names=F)

#################################################################
### Extended Data Fig. 3D: Summary profiles for annotated REs ###
#################################################################

# Run computeMatrix (DeepTools) to calculate ATAC/H3K27ac signals at annotated RE regions
# ~/Dropbox/Minn/ifnar_epigenome/scripts/make_summary_profile_plots.sh

plot.REs <- c("activated_enhancer", "deactivated_enhancer", "gained_TSS_promoters", "lost_TSS_promoters", "unchanged_REs")

cond1 <- "B16_WT"
cond2 <- "R499_WT"

pal.colors <- pal_nejm()(4)
names(pal.colors) <- c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO")

p_theme <- theme_minimal(base_size=11) +
	theme(aspect.ratio = 1,
		axis.text = element_text(size=11, color="black"),
		axis.title = element_text(size=11, color="black"),
		panel.border = element_rect(color = "black", fill=NA, size=1),
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

for(i in 1:length(plot.REs)) {
	RE.name <- plot.REs[i]

	mat.ATAC <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/summary_profiles/", RE.name, ".mat_ATAC_", cond1, "v", cond2, ".gz"), skip=1)
	mat.ATAC.signal <- data.frame(t(mat.ATAC[ ,7:ncol(mat.ATAC)]))
	colnames(mat.ATAC.signal) <- mat.ATAC$V4

	mat.H3K27Ac <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/summary_profiles/", RE.name, ".mat_H3K27Ac_", cond1, "v", cond2, ".gz"), skip=1)
	mat.H3K27Ac.signal <- data.frame(t(mat.H3K27Ac[ ,7:ncol(mat.H3K27Ac)]))
	colnames(mat.H3K27Ac.signal) <- mat.H3K27Ac$V4

	mat.H3K4me3 <- fread(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/summary_profiles/", RE.name, ".mat_H3K4me3_", cond1, "v", cond2, ".gz"), skip=1)
	mat.H3K4me3.signal <- data.frame(t(mat.H3K4me3[ ,7:ncol(mat.H3K4me3)]))
	colnames(mat.H3K4me3.signal) <- mat.H3K4me3$V4

	df.Peaks <- data.frame(
		bp = c(seq(-5000,(5000-1),100), seq(-5000,(5000-1),100)),
		Condition = c(rep(cond1, 5000*2 / 100), rep(cond2, 5000*2 / 100)))
	if(grepl(RE.name, pattern="promoters")) {
		df.Peaks$Signal <- rowMeans(mat.H3K4me3.signal)
		} else {
			df.Peaks$Signal <- rowMeans(mat.H3K27Ac.signal)
		}
	RE1 <- ifelse(grepl(RE.name, pattern="promoters"), "H3K4me3", "H3K27Ac")

	df.Peaks.ATAC <- data.frame(
		bp = c(seq(-2500,(2500-1),10), seq(-2500,(2500-1),10)),
		Condition = c(rep(cond1, 2500*2 / 10), rep(cond2, 2500*2 / 10)),
		Signal_ATAC = rowMeans(mat.ATAC.signal))
	df.Peaks.ATAC <- df.Peaks.ATAC[df.Peaks.ATAC$bp %in% -1000:1000, ]

	### SUMMARY PLOTS ###

	p1 <- ggplot(df.Peaks, aes_string(x = "bp", y = paste0("Signal_", RE1), color = "Condition")) +
		geom_line(size = 0.5) +
		scale_color_manual(values = pal.colors[match(c(cond1, cond2), names(pal.colors))]) +
		# ylim(c(0.2,2)) +
		labs(title = paste0(gsub("_", " ", RE.name), "\n", RE1), 
			x = "Distance from TSS (bp)", 
			y = "H3K27Ac Signal") +
		p_theme
	p2 <- ggplot(df.Peaks.ATAC, aes(x = bp, y = Signal_ATAC, color = Condition)) +
		geom_line(size = 0.5) +
		scale_color_manual(values = pal.colors[match(c(cond1, cond2), names(pal.colors))]) +
		# ylim(c(0.2,2)) +
		labs(title = "ATAC", 
			x = "Distance from TSS (bp)", 
			y = "ATAC Signal") +
		p_theme
	fig <- p1 + p2 + plot_layout(nrow=1, guides="collect")
	# ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome/results/summary_profiles/", RE.name, "_summary.pdf"), width=10, height=5)

	dat1 <- p1$data
	colnames(dat1) <- c("bp", "Condition", "Signal")
	dat2 <- p2$data
	colnames(dat2) <- c("bp", "Condition", "Signal")
	dat <- rbind(dat1, dat2)
	dat$Assay <- c(rep(RE1, nrow(p1$data)), rep("ATAC", nrow(p2$data)))
	write.csv(dat, file = paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 3D ", RE.name, ".csv"), quote=F, col.names=T, row.names=F)
}

#####################################################################################
### Extended Data Figure 3I: Acute vs. Chronic IFNG ATAC at Res499 activated loci ###
#####################################################################################

mat.plot <- adat[match(gr$ID, rownames(adat)), ]
mat.plot <- t(scale(t(mat.plot)))
hm <- Heatmap(
    mat.plot,
    name = "Res499\nActivated Enhancers",
    cluster_columns = F,
    show_column_names=F, 
    show_row_names=F,
    column_names_gp = gpar(fontsize = 10),
    col=colorRamp2(seq(-2.5,2.5, length=12), colorRampPalette(rev(brewer.pal(9, "RdBu")))(12)),
    column_split = factor(metadat_ATAC$Label, levels=names(use.pal_ATAC)),
    bottom_annotation = ha_col_ATAC,
    # row_split = 2,
    border = TRUE,
    use_raster=TRUE, 
    raster_quality=1,
    column_title = NULL,
    row_title = NULL)
ht <- draw(hm)
pdf("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/Res499 Activated Enhancers Chronic IFNG Filtered Heatmap.pdf", width = 6, height = 4)
ht
dev.off()

# # Write out plot data
# write.table(hm@matrix, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 2I.csv", sep=",", quote=F, col.names=T, row.names=T)

##############################################################
### Extended Data Figure 3J: Chronic IFNG Motif Enrichment ###
##############################################################

library(universalmotif)

dir.path <- "~/Dropbox/Minn/B16_B16y_R499/data/ATAC/processed/Motif/"
names <- c("Res499_Activated_Enhancers_Filtered_fixed1500_V2_noMask", "Res499_Activated_Enhancers_Other_fixed1500_V2_noMask")
labels <- c("Res499\nChronic IFNG Enhancers", "Res499\nOther Enhancers")

Vierstra.Motifs <- read.csv("~/Dropbox/Minn/resources/Motifs/Vierstra/V2/Vierstra_Archetype_Motifs_metadata.tsv", header=T, sep="\t")

# De novo
motif_list <- lapply(1:length(names), function(x) {
    dir.motif <- paste0(dir.path, names[x])
    denovo <- load_denovo_motifs(
        dir.motif = dir.motif, 
        name = names[x],
        label = labels[x])
    denovo <- denovo[ ,c("Motif_Name", "Similar_motif", "Archetype_Motif", "Family", "pval", "OE", "Perc_target", "Perc_background", "Label")]
    return(denovo)
})
mat.combined <- do.call(rbind, motif_list)
mat.combined <- mat.combined[sort.int(mat.combined$pval, decreasing=F, index.return=T)$ix, ]
print(mat.combined[grep(tolower(mat.combined$Motif_Name), pattern="irf|stat|nfkb"), ])

# Known
motif_list <- lapply(1:length(names), function(x) {
    dir.motif <- paste0(dir.path, names[x])
    print(dir.motif)
    known <- load_known_motifs(
        dir.motif = dir.motif, 
        name = names[x],
        label = labels[x])
    known <- known[ ,c("Motif_Name", "Similar_motif", "Archetype_Motif", "pval", "OE", "Perc_target", "Perc_background", "Label")]
    return(known)
})
mat.combined <- do.call(rbind, motif_list)
mat.combined <- mat.combined[sort.int(mat.combined$pval, decreasing=F, index.return=T)$ix, ]
print(mat.combined[grep(tolower(mat.combined$Motif_Name), pattern="irf|stat|nfkb"), ])

pval_threshold <- 1e-03
OE_threshold <- 1

mat.plot <- mat.combined[which(mat.combined$pval < pval_threshold & mat.combined$OE > OE_threshold), ]

pattern.interest <- "irf|stat|ets|creb|atf|fos|jun|ap1|maf"
mat.plot <- mat.plot %>%
    group_by(Label, Motif_Name) %>%
    slice_min(pval, n=1)
mat.plot <- mat.plot[sort.int(mat.plot$pval, decreasing=F, index.return=T)$ix, ]
mat.plot$Label <- str_wrap(gsub("_", " ", mat.plot$Label), width = 8)
mat.plot$Label <- factor(mat.plot$Label, levels=unique(mat.plot$Label))
mat.plot$Highlight <- ifelse(grepl(tolower(mat.plot$Motif_Name), pattern=pattern.interest), T, F)
mat.plot$Motif_ID <- factor(mat.plot$Motif_Name, levels=rev(unique(mat.plot$Motif_Name)))
# mat.plot$Archetype_Motif <- factor(mat.plot$Archetype_Motif, levels=rev(unique(mat.plot$Archetype_Motif)))
# head(data.frame(mat.plot[!mat.plot$Motif_ID %in% Vierstra.Motifs$motif_cluster, ]), 50)

fig <- ggplot(mat.plot, aes(x = Label, y = Motif_ID, color = -log10(pval), size = -log10(pval))) +
    # geom_point(alpha=0.8, size = 5) +
    geom_point(alpha=0.8) +
    scale_color_gradientn(
        colors=brewer.pal(9, "Greens")[3:9], 
        limits = c(ifelse(min(-log10(mat.plot$pval)) < 5, 5, min(-log10(mat.plot$pval))), ifelse(max(-log10(mat.plot$pval)) > 25, 25, max(-log10(mat.plot$pval)))), 
        oob = scales::squish) +
    labs(color="Fold enrichment\nover background\n(-log10 pval)", y = "", x = "") +
    guides(size = "none") +
    scale_size(range = c(2,8)) +
    theme_bw(base_size=12) +
    theme(
        axis.text.y = element_text(
            size = 11, 
            color = ifelse(grepl(tolower(levels(mat.plot$Motif_ID)), pattern=pattern.interest), "firebrick", "black"),
            face = ifelse(grepl(tolower(levels(mat.plot$Motif_ID)), pattern=pattern.interest), "bold", "plain")),
        axis.text.x = element_text(size=12, color = "black"),
        axis.ticks = element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))


# Write out plot data
dat <- data.frame(fig$data)
dat$Label <- gsub("\n", " ", dat$Label)
dat <- dat[ ,c("Archetype_Motif", "pval", "OE", "Label")]
write.table(dat, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 2J.csv", sep=",", quote=F, col.names=T, row.names=F)


