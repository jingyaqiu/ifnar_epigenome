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

source("~/Dropbox/Minn/ifnar_epigenome/scripts/Integrated Epigenomic Analysis/Integrated Epigenomic Analysis Functions.R")

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
	chr=as.character(seqnames(gr.bm)), 
	start=gr.bm$transcription_start_site, 
	end=gr.bm$transcription_start_site, 
	gene=gr.bm$mgi_symbol, 
	strand = as.character(strand(gr.bm))), keep.extra.columns = T)
gr.tss_window <- gr.tss
start(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) <- start(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) - 2500
end(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) <- end(gr.tss_window[as.character(strand(gr.tss_window)) == "+"]) + 1000
start(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) <- start(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) - 1000
end(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) <- end(gr.tss_window[as.character(strand(gr.tss_window)) == "-"]) + 2500

# Gene sets
gs.names <- c("IFN.I", "ISG.RS")
gs_list <- lapply(1:length(gs.names), function(x) {
  gs <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/", gs.names[x], "_mouse.txt"), sep="\t", header=F)
  gs <- unique(as.character(gs$V1))
  gs <- as.character(bm$ensembl_gene_id[na.omit(match(gs, bm$mgi_symbol))])
  return(gs)
})
names(gs_list) <- gs.names

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

### GO ENRICHMENT ON CIS-GENES ###

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
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2H.GO.csv"), sep=",", quote=F, col.names=T, row.names=F)

###########################################################
### Fig 3A: Plot ATAC, histone signals at annotated REs ###
###########################################################

### Activated enhancer signal density plot ###

p1 <- plot_signals_density(
	gr.peaks = RE_list[["activated_enhancer_ATAC"]], 
	assay.plot = c("H3K27Ac", "H3K4me1", "ATAC"), 
	plot.conditions = c("B16_WT", "B16_SKO"), 
	mat_list = mat_list, 
	colors = pal_nejm()(4)[1:2])
p2 <- plot_signals_density(
	gr.peaks = RE_list[["activated_enhancer_ATAC"]], 
	assay.plot = c("H3K27Ac", "H3K4me1", "ATAC"), 
	plot.conditions = c("R499_WT", "R499_SKO"), 
	mat_list = mat_list, 
	colors = pal_nejm()(4)[3:4])
fig <- ggarrange(p1, p2, nrow=2, common.legend = F, legend = "right")
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3A Activated Enhancer STAT1 dependency density plot.pdf"), width=7.5, height=6.5)

# # Write out plot data
# dat <- rbind(p1$data, p2$data)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3A.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
assay.plot <- c("H3K27Ac", "H3K4me1", "ATAC")
for(i in 1:length(assay.plot)) {
	assay <- assay.plot[i]
	mat.plot <- extract_signal_matrix(
		gr.peaks = RE_list[["activated_enhancer_ATAC"]],
		mat.assay = mat_list[[assay]],
		gr.assay = peakids2GRanges(rownames(mat_list[[assay]])))
	df <- data.frame(
		Sample = colnames(mat.plot),
		Condition = sapply(strsplit(colnames(mat.plot), split="_"), function(x) paste(x[1:2], collapse="_")),
		Signal = as.numeric(colMeans(mat.plot)))
	res.anova <- aov(Signal ~ Condition, data = df)
	print(assay)
	print(summary(res.anova))
	print(TukeyHSD(res.anova)) # Pairwise comparisons
}

# Activated Enhancers (Fig 3A)

# H3K27ac (pval = 1.8e-05):
# B16 WT vs R499 WT padj = 2e-05
# B16 WT vs B16 SKO padj = 7.6e-03
# R499 WT vs R499 SKO padj = 1.5e-03

# H3K4me1 (pval = 1.3e-03):
# B16 WT vs R499 WT padj = 2.6e-03
# B16 WT vs B16 SKO padj = 6.3e-01
# R499 WT vs R499 SKO padj = 9.8e-01

# ATAC (pval = 1.7e-11):
# B16 WT vs R499 WT padj = 0 (***)
# B16 WT vs B16 SKO padj = 2.8e-02
# R499 WT vs R499 SKO padj = 7.2e-01

### ~85% activated enhancers are pre-existing/primed in B16
### ~15% activated enhancers are de novo acquired?

###############################
### Fig 3E: ISG epigenetics ###
###############################

# Combine annotated promoter and enhancer REs
gr.promoter <- RE_repertoire_list[[1]]
gr.promoter$RE <- "promoter"
gr.enhancer <- RE_repertoire_list[[2]]
gr.enhancer$RE <- "enhancer"
ol <- findOverlaps(RE_list[["activated_enhancer_ATAC_inclusive"]], gr.enhancer)
gr.enhancer <- gr.enhancer[unique(to(ol))]
RE.combined <- c(gr.promoter, gr.enhancer)
RE.combined$idx <- 1:length(RE.combined)

gs.name <- "ISG.RS" # ISG.RS, IFN.I
gs <- gs_list[[gs.name]]
gs <- gs[gs %in% rownames(edat)]
gs <- gr.bm$mgi_symbol[match(gs, gr.bm$ensembl_gene_id)]

##############################################################

W <- 92000

# Link REs to target genes
# Promoter - overlap with TSS
# Enhancer - ATAC correlation with RNA with r > 0.1
gs_linked_REs <- identify_cis_REs_ATAC(
	gr.peaks = RE.combined,
	gs = gs, 
	W = W, 
	gr.bm = gr.bm, 
	gr.tss = gr.tss_window,
	adat = mat_list[["ATAC"]], 
	edat = edat,
	cor_cutoff = 0.2)
gs_linked_REs_filt <- sapply(gs_linked_REs, function(x) na.omit(unique(c(x[[1]], x[[2]]))))

# Remove genes with no associated REs
remove_idx <- unlist(sapply(1:length(gs_linked_REs_filt), function(x) {
	if(all(is.na(gs_linked_REs_filt[[x]])) | length(gs_linked_REs_filt[[x]])==0) return(x)
}))
if(length(remove_idx) > 0) gs_linked_REs_filt <- gs_linked_REs_filt[-remove_idx]

### Summary plots ###

assay.plot <- c("RNA", "H3K27Ac", "H3K4me1")

p_list <- calculate_cis_genes_average_signal(
	gs_linked_REs = gs_linked_REs_filt, 
	gr.peaks = RE.combined,
	mat_list = mat_list,
	assay.plot = assay.plot,
	plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO"),
	use.pal = pal_nejm()(4),
	method = "both",
	summarize_genes = TRUE,
	return_values = FALSE,
	plot.boxplot = FALSE,
	point_size = 4,
	stroke = 1,
	title = "",
	ylim = c(-2,2))
p.atac <- calculate_cis_genes_average_signal(
	gs_linked_REs = gs_linked_REs_filt, 
	gr.peaks = RE.combined,
	mat_list = mat_list,
	assay.plot = "ATAC",
	plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO"),
	use.pal = pal_nejm()(4),
	method = "gsva",
	summarize_genes = TRUE,
	return_values = FALSE,
	plot.boxplot = TRUE,
	point_size = 2,
	stroke = 1,
	title = "",
	ylim = c(-2,2))
p_list <- c(p_list, p.atac)
p <- p_list[[1]] + p_list[[2]] + p_list[[3]] + p_list[[4]] + plot_layout(nrow=1)
ggsave(p, file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", gs.name, " signals summary.pdf"), width=10, height=6)

# Write out plot data
dat <- rbind(p_list[[1]]$data, p_list[[2]]$data, p_list[[3]]$data, p_list[[4]]$data)
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Figure 3E.csv"), sep=",", quote=F, col.names=T, row.names=F)


# ISG.RS DKO
p <- calculate_cis_genes_average_signal(
	gs_linked_REs = gs_linked_REs_filt, 
	gr.peaks = RE.combined,
	mat_list = mat_list1,
	assay.plot = "ATAC",
	plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO"),
	use.pal = pal_nejm()(6),
	method = "gsva",
	summarize_genes = TRUE,
	return_values = FALSE,
	plot.boxplot = FALSE,
	point_size = 5,
	stroke = 1.2,
	title = "ISG.RS cis-REs",
	ylim = c(-2,2))
ggsave(p[[1]], file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", gs.name, " signals summary DKO.pdf"), width=5, height=4)

df <- calculate_cis_genes_average_signal(
	gs_linked_REs = gs_linked_REs_filt, 
	gr.peaks = RE.combined,
	mat_list = mat_list1,
	assay.plot = "ATAC",
	plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO"),
	use.pal = pal_nejm()(6),
	method = "gsva",
	summarize_genes = TRUE,
	return_values = TRUE,
	plot.boxplot = FALSE,
	point_size = 5,
	stroke = 1.2,
	title = "ISG.RS cis-REs",
	ylim = c(-2,2))
df <- df[[1]]
t.test(df$Signal[df$Condition == "B16_WT"], df$Signal[df$Condition == "R499_WT"])
t.test(df$Signal[df$Condition == "B16_WT"], df$Signal[df$Condition == "B16_SKO"])
t.test(df$Signal[df$Condition == "B16_SKO"], df$Signal[df$Condition == "R499_SKO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_SKO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_IRF3KO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_DKO"])

# # Write out plot data
# dat <- df[ ,c("Sample", "Condition", "Signal")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3F.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Heatmap
metadat_revisions <- data.frame(
    Sample = colnames(mat_list1[["ATAC"]]),
    Condition = sapply(strsplit(colnames(mat_list1[["ATAC"]]), split="_|-"), function(x) paste(x[1:2], collapse=" ")))
metadat_revisions$Label <- factor(metadat_revisions$Condition, levels = c("B16 WT", "B16 SKO", "R499 WT", "R499 SKO", "R499 IRF3KO", "R499 DKO"))
use.pal_revisions <- pal_nejm()(6)
names(use.pal_revisions) <- levels(metadat_revisions$Label)
ha_col_revisions <- HeatmapAnnotation(
    Condition=factor(metadat_revisions$Label, levels=names(use.pal_revisions)), 
    col=list(Condition=use.pal_revisions), 
    simple_anno_size=unit(0.7, "cm"),
    border = TRUE)

set.seed(12)
gr.atac_revisions <- peakids2GRanges(rownames(mat_list1[["ATAC"]]), delim="_")
gr.atac_revisions$ID <- rownames(mat_list1[["ATAC"]])

ol <- findOverlaps(gr.atac_revisions, RE.combined[as.numeric(unlist(gs_linked_REs_filt))])
mat.plot <- mat_list1[["ATAC"]][unique(from(ol)), ]
mat.plot <- t(scale(t(mat.plot)))
hm <- Heatmap(
    mat.plot,
    name = "ATAC",
    use_raster = TRUE, 
    raster_quality = 1, 
    cluster_columns = F, 
    cluster_rows = T, 
    show_column_names = F, 
    show_row_names = F, 
    row_names_gp = gpar(fontsize = 7),
    col = colorRamp2(seq(-2,2,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)), 
    column_split = metadat_revisions$Label[match(colnames(mat.plot), metadat_revisions$Sample)],
    bottom_annotation = ha_col_revisions,
    column_title = NULL,
    row_title = NULL,
    border = TRUE)
ht <- draw(hm)

pdf(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/ISG.RS Memory Domains DKO Dependence Heatmap.pdf"), width = 4.5, height = 3.25)
ht
dev.off()

# ggsave(p, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3E ISG signals summary.pdf", width=6, height=10)

# # Write out plot data
# dat <- rbind(p_list[[1]]$data, p_list[[2]]$data, p_list[[3]]$data, p_list[[4]]$data)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3E.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
for(i in 1:length(p_list)) {
	df <- p_list[[i]]$data
	res.anova <- aov(Signal ~ Condition, data = df)
	print(assay.plot[i])
	print(summary(res.anova))
	print(TukeyHSD(res.anova)) # Pairwise comparisons
}

# ISG features (Fig 3E)

# RNA (pval = 1.28-05):
# B16 WT vs R499 WT padj = 1.1e-02
# B16 WT vs B16 SKO padj = 1.5e-03
# R499 WT vs R499 SKO padj = 2.1e-05
# B16 SKO vs R499 SKO padj = 4.9e-01

# H3K27ac (pval = 1.5e-02):
# B16 WT vs R499 WT padj = 7.6e-02
# B16 WT vs B16 SKO padj = 1.9e-01
# R499 WT vs R499 SKO padj = 3.5e-02
# B16 SKO vs R499 SKO padj = 4.9e-01

# H3K4me1 (pval = 1.4e-02)
# B16 WT vs R499 WT padj = 8.4e-02
# B16 WT vs B16 SKO padj = 1.3e-01
# R499 WT vs R499 SKO padj = 2.7e-01
# B16 SKO vs R499 SKO padj = 4.6e-02

# ATAC (pval = 1.2e-06)
# B16 WT vs R499 WT padj = 5.3e-02
# B16 WT vs B16 SKO padj = 7.0e-05
# R499 WT vs R499 SKO padj = 5.7e-03
# B16 SKO vs R499 SKO padj =5.9e-04

###########################################################################
### Fig 3I: Activated enhancers coordinately mediated by STAT1 and IRF3 ###
###########################################################################

gr.atac_revisions <- peakids2GRanges(rownames(mat_list1[["ATAC"]]), delim="_")
gr.atac_revisions$name <- rownames(mat_list1[["ATAC"]])

gr.RE <- RE_list[["activated_enhancer"]]
ol <- findOverlaps(gr.atac_revisions, gr.RE)
gr.RE_ATAC <- gr.atac_revisions[unique(from(ol))]

### Density plot of ATAC signals (Fig. 3I) ###

p1 <- plot_signals_density(
	gr.peaks = gr.RE_ATAC, 
	assay.plot = "ATAC", 
	plot.conditions = c("B16_WT", "R499_WT"), 
	mat_list = mat_list1,
	colors = pal_nejm()(6)[c(1,3)])
p2 <- plot_signals_density(
	gr.peaks = gr.RE_ATAC, 
	assay.plot = "ATAC", 
	mat_list = mat_list1, 
	plot.conditions = c("R499_WT", "R499_SKO"), 
	colors = pal_nejm()(6)[c(3,4)])
p3 <- plot_signals_density(
	gr.peaks = gr.RE_ATAC, 
	assay.plot = "ATAC", 
	mat_list = mat_list1, 
	plot.conditions = c("R499_WT", "R499_IRF3KO"), 
	colors = pal_nejm()(6)[c(3,5)])
p4 <- plot_signals_density(
	gr.peaks = gr.RE_ATAC, 
	assay.plot = "ATAC", 
	mat_list = mat_list1,
	plot.conditions = c("R499_WT", "R499_DKO"), 
	colors = pal_nejm()(6)[c(1,6)])
fig <- p1 + p2 + p3 + p4 + plot_layout(nrow=1)
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3I Activated Enhancers DKO Effect density plot.pdf", width=10, height=6)

# # Write out plot data
# dat <- rbind(p1$data, p2$data, p3$data, p4$data)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3I.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
mat.plot <- extract_signal_matrix(
	gr.peaks = gr.RE_ATAC,
	mat.assay = mat_list1[["ATAC"]],
	gr.assay = peakids2GRanges(rownames(mat_list1[["ATAC"]])))
df <- data.frame(
	Sample = colnames(mat.plot),
	Condition = sapply(strsplit(colnames(mat.plot), split="_"), function(x) paste(x[1:2], collapse="_")),
	Signal = as.numeric(colMeans(mat.plot)))
res.anova <- aov(Signal ~ Condition, data = df)
print(assay)
print(summary(res.anova))
print(TukeyHSD(res.anova)) # Pairwise comparisons

# Activated enhancers (Fig. 3I)

# ATAC (pval = 4.25e-11)
# B16 WT vs R499 WT padj = 0 (***)
# B16 WT vs B16 SKO padj = 4.4e-01
# R499 WT vs R499 SKO padj = 4.0e-02
# R499 WT vs R499 IRF3 KO padj = 4.8e-03
# R499 WT vs R499 DKO padj = 0 (***)

##################################################
### Fig 3J: Activated Enhancers Memory Regions ###
##################################################

# library(DiffBind)

# # Differential features
# gr.B16SKOvsR499SKO_ATAC <- dba.report(db_list[["ATAC"]], contrast=7, th=1)
# gr.R499SKOvsB16SKO_H3K4me1 <- dba.report(db_list[["H3K4me1"]][[4]], contrast=1, th=1)

# # Persistent chromatin R499 SKO > B16 SKO
# ol <- findOverlaps(RE_list[["activated_enhancer"]], gr.B16SKOvsR499SKO_ATAC[gr.B16SKOvsR499SKO_ATAC$Fold < 0])
# gr.memory <- RE_list[["activated_enhancer"]][unique(from(ol))]
# ol <- findOverlaps(gr.memory, gr.R499SKOvsB16SKO_H3K4me1[gr.R499SKOvsB16SKO_H3K4me1$Fold > 0])
# gr.memory <- gr.memory[unique(from(ol))]

# # Must overlap with H3K4me1 called peak
# ol <- findOverlaps(gr.memory, gr_list[["H3K4me1"]])
# gr.memory <- gr_list[["H3K4me1"]][unique(to(ol))]
# make_HOMER_bed(gr.memory, file="~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/activated_enhancers_memory_domains_H3K4me1Centered.bed")

gr.memory <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/activated_enhancers_memory_domains_H3K4me1Centered.bed")

##################################################
##################################################

### Memory Regions ATAC Signal ###

ol <- findOverlaps(gr.memory, peakids2GRanges(rownames(mat_list1[["ATAC"]])))
mat.plot <- mat_list1[["ATAC"]][unique(to(ol)), ]
mat.plot <- t(scale(t(mat.plot)))

# Summary
df <- data.frame(
    Sample = colnames(mat.plot),
    Signal = colMeans(mat.plot))
df$Condition <- sapply(strsplit(as.character(df$Sample), split="_"), function(x) paste(x[1:2], collapse="_"))
df$Label <- factor(gsub("_", " ", df$Condition),
    levels = c("B16 WT", "B16 SKO", "R499 WT", "R499 SKO", "R499 IRF3KO", "R499 DKO"))
p1 <- ggplot(df, aes(x = Label, y = Signal, color = Label)) +
    geom_point(shape = 1, size = 3, stroke = 1.5) +
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
    scale_color_nejm() +
    labs(x = "", y = "ATAC Signal", color = "Condition") +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))
# # Write out plot data
# dat <- p1$data[ ,c("Sample", "Label", "Signal")]
# write.table(dat, file="~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3J.ATAC.summary.csv", sep=",", quote=F, col.names=T, row.names=F)

# Individual loci
hm <- Heatmap(
    mat.plot, 
    name="Persistent Enhancers", 
    use_raster=TRUE, 
    raster_quality=1, 
    cluster_columns=F, 
    cluster_rows=T, 
    show_column_names=F, 
    show_row_names=F, 
    row_names_gp = gpar(fontsize = 6),
    col=colorRamp2(seq(-2.5,2.5,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)), 
    column_split=droplevels(factor(sapply(strsplit(colnames(mat.plot), split="_"), function(x) paste(x[1:2], collapse="_")), levels=c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO"))),
    column_title=NULL, 
    width=ncol(mat.plot)*unit(3, "mm"), 
    column_names_rot=45, 
    column_names_gp = gpar(fontsize = 10),
    border = TRUE)
pdf("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3J Activated Enhancers Memory Regions Heatmap ATAC.pdf", width=6.5, height=2.75)
draw(hm)
dev.off()
# # Write out plot data
# dat <- mat.plot
# write.table(dat, file="~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3J.ATAC.heatmap.csv", sep=",", quote=F, col.names=T, row.names=F)

### Memory Regions H3K4me1 Signal ###

ol <- findOverlaps(gr.memory, peakids2GRanges(rownames(mat_list1[["H3K4me1"]])))
mat.plot <- mat_list1[["H3K4me1"]][unique(to(ol)), ]
mat.plot <- t(scale(t(mat.plot)))

# Summary
df <- data.frame(
    Sample = colnames(mat.plot),
    Signal = colMeans(mat.plot))
df$Condition <- sapply(strsplit(as.character(df$Sample), split="_"), function(x) paste(x[1:2], collapse="_"))
df$Label <- factor(gsub("_", " ", df$Condition),
    levels = c("B16 WT", "B16 SKO", "R499 WT", "R499 SKO"))
p2 <- ggplot(df, aes(x = Label, y = Signal, color = Label)) +
    geom_point(shape = 1, size = 3, stroke = 1.5) +
    stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
    scale_color_nejm() +
    labs(x = "", y = "H3K4me1 Signal", color = "Condition") +
    theme_classic(base_size = 14) +
    theme(aspect.ratio = 1,
        axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black", size = 12),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "grey90"), 
        panel.grid.minor = element_line(size = 0.25, linetype = 'solid', colour = "grey90"))
# # Write out plot data
# dat <- p2$data[ ,c("Sample", "Label", "Signal")]
# write.table(dat, file="~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3J.H3K4me1.summary.csv", sep=",", quote=F, col.names=T, row.names=F)

# Individual loci
hm <- Heatmap(
    mat.plot, 
    name="Persistent Enhancers", 
    use_raster=TRUE, 
    raster_quality=1, 
    cluster_columns=F, 
    cluster_rows=T, 
    show_column_names=F, 
    show_row_names=F, 
    row_names_gp = gpar(fontsize = 6),
    col=colorRamp2(seq(-2.5,2.5,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)), 
    column_split=droplevels(factor(sapply(strsplit(colnames(mat.plot), split="_"), function(x) paste(x[1:2], collapse="_")), levels=c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO"))),
    column_title=NULL, 
    width=ncol(mat.plot)*unit(5, "mm"), 
    column_names_rot=45, 
    column_names_gp = gpar(fontsize = 10),
    border = TRUE)
pdf("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3J Activated Enhancers Memory Regions Heatmap H3K4me1.pdf", width=4.5, height=2.75)
draw(hm)
dev.off()
# # Write out plot data
# dat <- mat.plot
# write.table(dat, file="~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3J.H3K4me1.heatmap.csv", sep=",", quote=F, col.names=T, row.names=F)

fig <- p2 + p1 + plot_layout(guides = "collect")
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3J Activated Enhancers Memory Regions Summary.pdf", width=7, height=4.5)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
df <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3J.H3K4me1.summary.csv", sep=",", header=T)
res.anova <- aov(Signal ~ Label, data = df)
print(assay)
print(summary(res.anova))
print(TukeyHSD(res.anova)) # Pairwise comparisons

# Persistent enhancers (Fig. 3J)

# ATAC (pval = 7.25e-12)
# B16 WT vs R499 WT padj = 0 (***)
# B16 WT vs B16 SKO padj = 4.3e-02
# R499 WT vs R499 SKO padj = 9.0e-01
# R499 WT vs R499 IRF3 KO padj = 1.0e-04
# R499 WT vs R499 DKO padj = 0 (***)

# H3K4me1 (pval = 1.0e-03)
# B16 WT vs R499 WT padj = 3.9e-03
# B16 WT vs B16 SKO padj = 9.9e-01
# R499 WT vs R499 SKO padj = 6.6e-01
# B16 SKO vs R499 SKO padj = 2.1e-03

