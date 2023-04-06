
### Inflammatory memory signature ###

rm(list=ls())

library(tidyverse)
library(GSVA)
library(rtracklayer)
library(Biostrings)
library(BSgenome.Mmusculus.UCSC.mm10)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggsci)
library(ggridges)
library(patchwork)
library(ggpubr)
library(ggrepel)

directory <- "local"

source("~/Dropbox/Minn/resources/useful_protocols/Bulk ATAC processing/ATAC Visualization and Analysis Functions.R")
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

############################################
### Identify inflammatory memory domains ###
############################################

# Criteria:
# Condition-specific STAT1 KO effect, STAT1 KO is stronger in B16 compared to Res499; in other words, persistent chromatin accessibility with STAT1 loss in Res499 compared to B16
# Res499 activated enhancer?
# Center around H3K4me1 peak call?

Interaction Testing ATAC Combined.R

###################################
### Inflammatory memory domains ###
###################################

# Preferential STAT1-dependency (interaction term)

# Memory domains
gr.memory <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/IFN Inflammatory Memory Domains ALL.bed")
start(gr.memory) <- start(gr.memory) - 1

# Memory domains activated REs
gr.memory_filtered <- import.bed(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/IFN Inflammatory Memory Domains Activated REs.bed"))
start(gr.memory_filtered) <- start(gr.memory_filtered) - 1

# Memory genes
cis_genes <- read.delim("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/Inflammatory Memory Genes ALL V1.txt", header=F)$V1

#################################################
### H3K27ac/H3K4me1 normalized count matrices ###
#################################################

# H3K27ac normalized count matrix
mat.H3K27ac <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/H3K27Ac_consensus_mat_insertion_counts_WT_pooled_vst.txt", sep="\t", header=T)
colnames(mat.H3K27ac) <- gsub("_H3K27Ac_rep", "_", colnames(mat.H3K27ac))
gr.H3K27ac <- peakids2GRanges(peakids = rownames(mat.H3K27ac), delim="_")
metadat.H3K27ac <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_H3K27ac.csv", header=T)
metadat.H3K27ac$Label <- factor(metadat.H3K27ac$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# H3K4me1 normalized count matrix
mat.H3K4me1 <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/H3K4me1_consensus_mat_insertion_counts_WT_pooled_vst.txt", sep="\t", header=T)
colnames(mat.H3K4me1) <- gsub("_H3K4me1_rep", "_", colnames(mat.H3K4me1))
gr.H3K4me1 <- peakids2GRanges(peakids = rownames(mat.H3K4me1), delim="_")
metadat.H3K4me1 <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_H3K4me1.csv", header=T)
metadat.H3K4me1$Label <- factor(metadat.H3K4me1$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# ATAC-seq normalized counts
mat.ATAC_original <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_original_consensus_mat_tn5_insertion_counts_IDR_rlog.txt", sep="\t", header=T, check.names = F)
gr.ATAC_original <- peakids2GRanges(rownames(mat.ATAC_original), delim = "_")
metadat.ATAC_original <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_ATAC_original.csv", header=T)
metadat.ATAC_original$Label <- factor(metadat.ATAC_original$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))

# RNA-seq normalized counts
mat.RNA_original <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/RNA_counts_rlog.txt", sep="\t", header=T, check.names = F)
mat.RNA_original <- mat.RNA_original[!is.na(match(rownames(mat.RNA_original), bm$ensembl_gene_id)), ]
metadat.RNA_original <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_RNA_original.csv", header=T)
metadat.RNA_original$Label <- factor(metadat.RNA_original$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO"))
DEGs_list <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Diffbind/DEG_list.rds")

# Remove B16_SKO_3 poor quality sample (low FRiP, not consistent with other replicates)
remove_samples <- c("B16_SKO_3")
mat.RNA_original <- mat.RNA_original[ ,!colnames(mat.RNA_original) %in% remove_samples]
mat.ATAC_original <- mat.ATAC_original[ ,!colnames(mat.ATAC_original) %in% remove_samples]
metadat.RNA_original <- metadat.RNA_original[match(colnames(mat.RNA_original), metadat.RNA_original$Sample), ]
metadat.ATAC_original <- metadat.ATAC_original[match(colnames(mat.ATAC_original), metadat.ATAC_original$Sample), ]
colnames(mat.RNA_original) <- gsub("cas", "WT", colnames(mat.RNA_original))
colnames(mat.ATAC_original) <- gsub("cas", "WT", colnames(mat.ATAC_original))
metadat.ATAC_original$Sample <- metadat.ATAC_original$ID
metadat.RNA_original$Sample <- metadat.RNA_original$ID

use.pal_original <- pal_nejm()(4)
names(use.pal_original) <- levels(metadat.RNA_original$Label)
ha_col_original <- HeatmapAnnotation(
    Condition=factor(metadat.RNA_original$Label, levels=names(use.pal_original)), 
    col=list(Condition=use.pal_original), 
    simple_anno_size=unit(0.7, "cm"),
    border = TRUE)

##################################################################################
### Fig. 3B: Plot ATAC, H3K27ac, H3K4me1 signal at inflammatory memory regions ###
##################################################################################

# H3K27ac
ol <- findOverlaps(gr.memory, gr.H3K27ac)
gs <- gr.H3K27ac$ID[unique(to(ol))]
p.H3K27ac <- plot_summary(
    mat.plot = mat.H3K27ac,
    geneset_list = list(gs),
    geneset_names = paste0("H3K27ac"),
    metadat = metadat.H3K27ac,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

# H3K4me1
ol <- findOverlaps(gr.memory, gr.H3K4me1)
gs <- gr.H3K4me1$ID[unique(to(ol))]
p.H3K4me1 <- plot_summary(
    mat.plot = mat.H3K4me1,
    geneset_list = list(gs),
    geneset_names = paste0("H3K4me1"),
    metadat = metadat.H3K4me1,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

# ATAC
ol <- findOverlaps(gr.memory, gr.ATAC_original)
gs <- gr.ATAC_original$ID[unique(to(ol))]
p.ATAC <- plot_summary(
    mat.plot = mat.ATAC_original,
    geneset_list = list(gs),
    geneset_names = paste0("ATAC"),
    metadat = metadat.ATAC_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

fig <- p.H3K27ac + p.H3K4me1 + p.ATAC
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", name, " Domains Epigenetic Profile Summaries.pdf"), width = 8, height = 4)

### Significance tests ###

# H3K27ac
ol <- findOverlaps(gr.memory, gr.H3K27ac)
gs <- gr.H3K27ac$ID[unique(to(ol))]
df.H3K27ac <- plot_summary(
    mat.plot = mat.H3K27ac,
    geneset_list = list(gs),
    geneset_names = paste0("H3K27ac"),
    metadat = metadat.H3K27ac,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 WT"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "B16 STAT1 KO"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 STAT1 KO"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 STAT1 KO"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 STAT1 KO"])

# H3K4me1
ol <- findOverlaps(gr.memory, gr.H3K4me1)
gs <- gr.H3K4me1$ID[unique(to(ol))]
df.H3K4me1 <- plot_summary(
    mat.plot = mat.H3K4me1,
    geneset_list = list(gs),
    # geneset_names = paste0("H3K4me1\n(n = ", length(gs), ")"),
    geneset_names = paste0("H3K4me1"),
    metadat = metadat.H3K4me1,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 WT"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "B16 STAT1 KO"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 STAT1 KO"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 STAT1 KO"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 STAT1 KO"])

# ATAC
ol <- findOverlaps(gr.memory, gr.ATAC_original)
gs <- gr.ATAC_original$ID[unique(to(ol))]
df.ATAC <- plot_summary(
    mat.plot = mat.ATAC_original,
    geneset_list = list(gs),
    # geneset_names = paste0("ATAC\n(n = ", length(gr.interest), ")"),
    geneset_names = paste0("ATAC"),
    metadat = metadat.ATAC_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 WT"], df.ATAC$Signal[df.ATAC$Group == "Res499 WT"])
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 WT"], df.ATAC$Signal[df.ATAC$Group == "B16 STAT1 KO"])
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 STAT1 KO"], df.ATAC$Signal[df.ATAC$Group == "Res499 STAT1 KO"])
t.test(df.ATAC$Signal[df.ATAC$Group == "Res499 WT"], df.ATAC$Signal[df.ATAC$Group == "Res499 STAT1 KO"])

# * = 0.05
# ** = 0.01
# *** = 0.001

# # Write out plot data
# dat <- rbind(df.H3K27ac, df.H3K4me1, df.ATAC)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3B.csv"), sep=",", quote=F, col.names=T, row.names=F)

############################################
### Fig. 3D: Genomic regions annotations ###
############################################

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

peakAnno <- annotatePeak(resize(gr.memory, fix = "center", width = 750),
	tssRegion = c(-500,500),
	TxDb = txdb,
	annoDb = "org.Mm.eg.db")
# p.anno <- upsetplot(peakAnno)

# Visualize

df <- peakAnno@anno
df$Annotation1 <- sapply(strsplit(df$annotation, split="\\("), function(x) x[[1]])
df$Annotation1 <- gsub(" ", "", df$Annotation1)
df <- data.frame(
    Annotation = factor(names(table(df$Annotation1))),
    Count = as.numeric(table(df$Annotation1)))
df$Percent <- df$Count / sum(df$Count)
df <- df[sort.int(df$Count, decreasing=T, index.return=T)$ix, ]
df$Annotation <- gsub("DistalIntergenic", "Distal Intergenic", df$Annotation)
# df$Annotation <- factor(df$Annotation, levels = c("Promoter", "Intron", "Exon", "Downstream", "DistalIntergenic", "5'UTR", "3'UTR"))
df$Annotation <- factor(df$Annotation, levels = df$Annotation)
df$ymax <- cumsum(df$Percent)
df$ymin <- c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$Label <- paste0(round(df$Percent*100, 1), "%")

p.anno <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Annotation)) +
    geom_rect(alpha=0.8, color="black") +
    geom_text_repel(data=df, aes(x=4, label = Label, y = labelPosition), color="black", size=4, nudge_x = .3, segment.size = .6) +
    labs(fill = "") +
    scale_fill_manual(values=colorRampPalette(brewer.pal(9, "Set1"))(length(levels(df$Annotation)))) +
    coord_polar(theta="y") +
    xlim(c(1,4)) +
    theme_void(base_size = 14) +
    theme(
        plot.margin=margin(1,1,1,1, "cm"))

### Percentage overlapping Res499 activated regulatory elements ###

df <- data.frame(
    Annotation = factor(c("Res499 Activated REs", "Other"), levels = c("Res499 Activated REs", "Other")),
    Percentage = c(length(gr.memory_filtered)/length(gr.memory),
        1 - (length(gr.memory_filtered)/length(gr.memory))))
df$ymax <- cumsum(df$Percentage)
df$ymin <- c(0, head(df$ymax, n=-1))
df$labelPosition <- (df$ymax + df$ymin) / 2
df$Label <- paste0(round(df$Percent*100, 1), "%")

p.REs <- ggplot(df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Annotation)) +
    geom_rect(alpha=0.8, color="black") +
    geom_text_repel(data=df, aes(x=4, label = Label, y = labelPosition), color="black", size=4, nudge_x = .3, segment.size = .6) +
    labs(fill = "") +
    scale_fill_manual(values=pal_nejm()(2)) +
    coord_polar(theta="y") +
    xlim(c(1,4)) +
    theme_void(base_size = 14) +
    theme(
        plot.margin=margin(1,1,1,1, "cm"))

fig <- p.anno + p.REs
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", name, " Domains Annotations.pdf"), width=10, height=5)

# # Write out plot data
# dat <- df[ ,c("Annotation", "Percent")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3D_left.csv"), sep=",", quote=F, col.names=T, row.names=F)

# # Write out plot data
# dat <- df[ ,c("Annotation", "Percentage")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3D_right.csv"), sep=",", quote=F, col.names=T, row.names=F)

#############################################
### Fig. 3E: Motif Enrichment in IFN-IMDs ###
#############################################



############################################################################
### Fig. 3F: Inflammatory memory gene expression and associated pathways ###
############################################################################

### Link memory domains to gene pathways ###

# Summary of inflammatory memory genes
gs <- unique(gr.bm$ensembl_gene_id[match(cis_genes, gr.bm$mgi_symbol)])
p <- plot_summary(
    mat.plot = mat.RNA_original,
    geneset_list = list(gs),
    geneset_names = paste0("Inflammatory Memory\nGenes (n = ", length(gs), ")"),
    metadat = metadat.RNA_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")
ggsave(p, file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/Inflammatory Memory Genes Summary ", memory_subset, ".pdf"), width = 4.5, height = 4.5)

# Significance test
gs <- unique(gr.bm$ensembl_gene_id[match(cis_genes, gr.bm$mgi_symbol)])
df <- plot_summary(
    mat.plot = mat.RNA_original,
    geneset_list = list(gs),
    geneset_names = paste0("Inflammatory Memory\nGenes (n = ", length(gs), ")"),
    metadat = metadat.RNA_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "Res499 WT"])
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "B16 STAT1 KO"])
t.test(df$Signal[df$Group == "B16 STAT1 KO"], df$Signal[df$Group == "Res499 STAT1 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 STAT1 KO"])

# Write out plot data
dat <- df[ ,c("Sample", "Group", "Signal")]
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3F_left.csv"), sep=",", quote=F, col.names=T, row.names=F)

### Pathway enrichment in inflammatory memory genes ###

library(clusterProfiler)
library(org.Mm.eg.db)
library(GO.db)

# Run clusterProfiler GO enrichment
go <- enrichGO(
    gene = cis_genes,
    OrgDb = org.Mm.eg.db,
    ont = "BP", 
    keyType = "SYMBOL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.1)

go.terms <- c(
    "response to interferon-beta", "response to interferon-alpha",
    "response to virus", "negative regulation of viral process", 
    "peptide antigen assembly with MHC class II protein complex", 
    "regulation of innate immune response",
    "response to interferon-gamma")
redundant.terms <- c(
    "cellular response to interferon-beta", "cellular response to interferon-alpha",
    "viral process", "defense response to virus", "defense response to symbiont", "viral life cycle", "regulation of viral life cycle", "negative regulation of viral genome replication", "viral entry into host cell", "regulation of viral process", "regulation of response to biotic stimulus", "biological process involved in symbiotic interaction", "biological process involved in interaction with host", "entry into host", "movement in host environment",
    "antigen processing and presentation", "antigen processing and presentation of exogenous antigen", "antigen processing and presentation of peptide antigen", "MHC class II protein complex assembly", "peptide antigen assembly with MHC protein complex")

# Plot GO pathway enrichment plot
go.plot <- data.frame(go)
go.plot <- go.plot[1:ifelse(nrow(go.plot) > 25, 25, nrow(go.plot)), ]
go.plot <- go.plot[!go.plot$Description %in% redundant.terms, ]
go.plot$Description <- factor(go.plot$Description, levels = rev(go.plot$Description))
p.go <- ggplot(go.plot, aes(x=Description, y=-log10(p.adjust))) +
    geom_bar(stat="identity", width=0.75, size=1, fill=brewer.pal(9, "Reds")[4], color=brewer.pal(9, "Reds")[4]) +
    coord_flip() +
    labs(x="", y="-log10(padj)", title = "GO BP Enrichment") +
    theme_bw() +
    theme(
        legend.position="none",
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))
ggsave(p.go, file = paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", name, " Genes GO BP Enrichment ", memory_subset, ".pdf"), width = 8, height = 3.5)

# # Write out plot data
# dat <- go.plot[ ,c("Description", "p.adjust")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3F_right.csv"), sep=",", quote=F, col.names=T, row.names=F)

#############################################
### Fig. 3G: Many ISG.RS are memory genes ###
#############################################

go_name <- "ISG.RS"
fontsize <- 11
gs <- gs_list[["ISG.RS"]]

res <- DEGs_list[["B16_v_R499"]]
res.gs <- res[na.omit(match(gs, rownames(res))), ]
res.gs <- res.gs[sort.int(res.gs$log2FoldChange, decreasing=T, index.return=T)$ix, ]
res.gs$MGI <- gr.bm$mgi_symbol[match(rownames(res.gs), gr.bm$ensembl_gene_id)]
res.gs$MGI <- factor(res.gs$MGI, levels = rev(res.gs$MGI))
res.gs$Memory <- factor(ifelse(res.gs$MGI %in% cis_genes, "Memory Gene", "Non-memory Gene"), levels = c("Memory Gene", "Non-memory Gene"))
# res.gs$log2FoldChange[res.gs$log2FoldChange > 3] <- 3

p <- ggplot(res.gs, aes(x = log2FoldChange, y = MGI, fill = Memory)) +
    geom_bar(stat="identity", width=0.75, size=1) +
    scale_fill_manual(values = c(brewer.pal(9, "Reds")[4], pal_jama()(1))) +
    labs(x="log2FoldChange\n(Res499 vs B16)", y="", title = paste0(str_wrap(go_name, width = 50)), fill = "") +
    xlim(c(-2.5,3)) +
    theme_bw(base_size = 11) +
    theme(
        axis.text.y=element_text(size=10, color="black"),
        axis.text.x = element_text(size = fontsize, color = "black", hjust = 1, angle = 45),
        axis.title.x = element_text(size=12, color="black"),
        axis.title.y = element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))
ggsave(p, file = paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", go_name, ".pdf"), width = 5, height = 7.5)

# # Write out plot data
# dat <- res.gs[ ,c("MGI", "log2FoldChange")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3G.csv"), sep=",", quote=F, col.names=T, row.names=F)

#######################################################
### Fig 3H: ISG expression and epigenetic profiling ###
#######################################################

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
mat_list[["RNA"]] <- mat.RNA_original

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
gs <- gs[gs %in% rownames(mat.RNA_original)]
gs <- gr.bm$mgi_symbol[match(gs, gr.bm$ensembl_gene_id)]

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
    adat = mat.ATAC_original, 
    edat = mat.RNA_original,
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

# # Write out plot data
# dat <- rbind(p_list[[1]]$data, p_list[[2]]$data, p_list[[3]]$data, p_list[[4]]$data)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3E.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
for(i in 1:length(p_list)) {
    # df <- p_list[[i]]$data
    df <- dat[dat$Assay == assay.plot[i], ]
    res.anova <- aov(Signal ~ Condition, data = df)
    print(assay.plot[i])
    print(summary(res.anova))
    print(TukeyHSD(res.anova)) # Pairwise comparisons
}

############################################################################
### Extended Data Figure 4A: Plot ATAC, histone signals at annotated REs ###
############################################################################

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

###############################################################################
### Extended Data Figure 4B: GO pathway enrichment on Resolved Domain Genes ###
###############################################################################

cis_genes <- read.table("~/Dropbox/Xu Qiu IFN paper/Revsion 2 manuscript and figures/Revision/Revision drafts with figures v17/Supplementary Files/IFN Resolved Domains Genes.txt", header=F)[,1]

# Run GO
go <- enrichGO(
  gene = cis_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP", 
  keyType = "SYMBOL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.1)

# Plot GO pathway enrichment plot
go.plot <- data.frame(go)
go.plot <- go.plot[1:ifelse(nrow(go.plot) > 25, 25, nrow(go.plot)), ]
go.plot <- go.plot[!go.plot$Description %in% redundant.terms, ]
go.plot$Description <- factor(go.plot$Description, levels = rev(go.plot$Description))
p.go <- ggplot(go.plot, aes(x=Description, y=-log10(p.adjust))) +
    geom_bar(stat="identity", width=0.75, size=1, fill=brewer.pal(9, "Reds")[4], color=brewer.pal(9, "Reds")[4]) +
    coord_flip() +
    labs(x="", y="-log10(padj)", title = "GO BP Enrichment") +
    theme_bw() +
    theme(
        legend.position="none",
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.title.x=element_text(size=12, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))
p.go

# # Write out plot data
# dat <- data.frame(p.go$data)
# dat <- dat[ ,c("ID", "Description", "p.adjust")]
# write.table(dat, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 4B.csv", sep=",", quote=F, col.names=T, row.names=F)

##############################################################################
### Extended Data Figure 4C: Summary enrichment scores at Resolved Domains ###
##############################################################################

# Resolved
gr.resolved <- import.bed("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/IFN Resolved Domains ALL.bed")
start(gr.resolved) <- start(gr.resolved) - 1

# H3K27ac
ol <- findOverlaps(gr.resolved, gr.H3K27ac)
gs <- gr.H3K27ac$ID[unique(to(ol))]
p.H3K27ac <- plot_summary(
    mat.plot = mat.H3K27ac,
    geneset_list = list(gs),
    geneset_names = paste0("H3K27ac"),
    metadat = metadat.H3K27ac,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

# H3K4me1
ol <- findOverlaps(gr.resolved, gr.H3K4me1)
gs <- gr.H3K4me1$ID[unique(to(ol))]
p.H3K4me1 <- plot_summary(
    mat.plot = mat.H3K4me1,
    geneset_list = list(gs),
    geneset_names = paste0("H3K4me1"),
    metadat = metadat.H3K4me1,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

# ATAC
ol <- findOverlaps(gr.resolved, gr.ATAC_original)
gs <- gr.ATAC_original$ID[unique(to(ol))]
p.ATAC <- plot_summary(
    mat.plot = mat.ATAC_original,
    geneset_list = list(gs),
    geneset_names = paste0("ATAC"),
    metadat = metadat.ATAC_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

fig <- p.H3K27ac + p.H3K4me1 + p.ATAC
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/Resolved Domains Epigenetic Profile Summaries.pdf"), width = 8, height = 4)

### Significance tests ###

# H3K27ac
ol <- findOverlaps(gr.resolved, gr.H3K27ac)
gs <- gr.H3K27ac$ID[unique(to(ol))]
df.H3K27ac <- plot_summary(
    mat.plot = mat.H3K27ac,
    geneset_list = list(gs),
    geneset_names = paste0("H3K27ac"),
    metadat = metadat.H3K27ac,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 WT"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "B16 STAT1 KO"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "B16 STAT1 KO"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 STAT1 KO"])
t.test(df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 WT"], df.H3K27ac$Signal[df.H3K27ac$Group == "Res499 STAT1 KO"])

# H3K4me1
ol <- findOverlaps(gr.resolved, gr.H3K4me1)
gs <- gr.H3K4me1$ID[unique(to(ol))]
df.H3K4me1 <- plot_summary(
    mat.plot = mat.H3K4me1,
    geneset_list = list(gs),
    # geneset_names = paste0("H3K4me1\n(n = ", length(gs), ")"),
    geneset_names = paste0("H3K4me1"),
    metadat = metadat.H3K4me1,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 WT"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "B16 STAT1 KO"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "B16 STAT1 KO"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 STAT1 KO"])
t.test(df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 WT"], df.H3K4me1$Signal[df.H3K4me1$Group == "Res499 STAT1 KO"])

# ATAC
ol <- findOverlaps(gr.resolved, gr.ATAC_original)
gs <- gr.ATAC_original$ID[unique(to(ol))]
df.ATAC <- plot_summary(
    mat.plot = mat.ATAC_original,
    geneset_list = list(gs),
    # geneset_names = paste0("ATAC\n(n = ", length(gr.interest), ")"),
    geneset_names = paste0("ATAC"),
    metadat = metadat.ATAC_original,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(4),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 WT"], df.ATAC$Signal[df.ATAC$Group == "Res499 WT"])
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 WT"], df.ATAC$Signal[df.ATAC$Group == "B16 STAT1 KO"])
t.test(df.ATAC$Signal[df.ATAC$Group == "B16 STAT1 KO"], df.ATAC$Signal[df.ATAC$Group == "Res499 STAT1 KO"])
t.test(df.ATAC$Signal[df.ATAC$Group == "Res499 WT"], df.ATAC$Signal[df.ATAC$Group == "Res499 STAT1 KO"])

# # Write out plot data
# dat <- rbind(df.H3K27ac, df.H3K4me1, df.ATAC)
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 4C.csv"), sep=",", quote=F, col.names=T, row.names=F)

##########################################################################
### Extended Data Figure 4E: Chronic IFNG on inflammatory memory genes ###
##########################################################################

# RNA-seq normalized counts
edat_chronicIFNG <- read.table("~/Dropbox/Minn/B16_B16y_R499/data/RNA/GRCm38_Build/mat/RNA_counts_vst.txt", sep="\t", header=T, check.names = F)

# Metadata
metadat_chronicIFNG <- read.csv("~/Dropbox/Minn/B16_B16y_R499/data/RNA/metadata.csv", sep=",", header=T, stringsAsFactors=F)
metadat_chronicIFNG$Condition <- gsub("A", "", metadat_chronicIFNG$Condition)
metadat_chronicIFNG$Label <- factor(metadat_chronicIFNG$Condition, levels=c("B16", "B16y_6hr", "B16y_3wk", "Res499"))
levels(metadat_chronicIFNG$Label) <- c("B16", "B16y 6hr Acute IFNG", "B16y 3.5wk Chronic IFNG", "Res499")

# Memory genes
cis_genes <- read.delim("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Annotated REs/Inflammatory Memory Genes ALL V1.txt", header=F)$V1

use.pal_chronicIFNG <- colorRampPalette(pal_nejm()(4)[c(1,3)])(4)
names(use.pal_chronicIFNG) <- levels(metadat_chronicIFNG$Label)

p <- plot_summary(
    mat.plot = edat_chronicIFNG,
    geneset_list = list(na.omit(gr.bm$ensembl_gene_id[match(cis_genes, gr.bm$mgi_symbol)])),
    geneset_names = paste0("Inflammatory Memory\nGenes (n = ", length(cis_genes), ")"),
    metadat = metadat_chronicIFNG,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal_chronicIFNG,
    method = "gsva",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")
ggsave(p, file = "~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/Inflammatory Memory Genes Chronic IFNG Summary.pdf", width = 5, height = 5)

df <- plot_summary(
    mat.plot = edat_chronicIFNG,
    geneset_list = list(na.omit(gr.bm$ensembl_gene_id[match(cis_genes, gr.bm$mgi_symbol)])),
    geneset_names = paste0("Inflammatory Memory\nGenes (n = ", length(cis_genes), ")"),
    metadat = metadat_chronicIFNG,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = use.pal_chronicIFNG,
    method = "gsva",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "B16y 6hr Acute IFNG"])
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "B16y 3.5wk Chronic IFNG"])
t.test(df$Signal[df$Group == "B16"], df$Signal[df$Group == "Res499"])

# Write out plot data
dat <- df[ ,c("Sample", "Group", "Signal")]
dat$Sample <- gsub("A", "", dat$Sample)
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 4E.csv"), sep=",", quote=F, col.names=T, row.names=F)

#########################################################################################
### Extended Data Figure 4F: Res499 vs. B16 log2 fold change for pathways of interest ###
#########################################################################################

# response to interferon-beta (GO:0035456)
go_term <- "GO:0035456"
go_name <- "response to interferon-beta"
plot.genes <- "all"
fontsize <- 9

results <- AnnotationDbi::select(
    org.Mm.eg.db, 
    keys=c(go_term), 
    columns = c('SYMBOL'), 
    keytype = "GOALL")
gs <- unique(results$SYMBOL)
gs <- na.omit(gr.bm$ensembl_gene_id[match(gs, gr.bm$mgi_symbol)])

res <- DEGs_list[["B16_v_R499"]]
res.gs <- res[na.omit(match(gs, rownames(res))), ]
res.gs <- res.gs[sort.int(res.gs$log2FoldChange, decreasing=T, index.return=T)$ix, ]
res.gs$MGI <- gr.bm$mgi_symbol[match(rownames(res.gs), gr.bm$ensembl_gene_id)]
res.gs$MGI <- factor(res.gs$MGI, levels = res.gs$MGI)
res.gs$Memory <- factor(ifelse(res.gs$MGI %in% cis_genes, "Memory Gene", "Non-memory Gene"), levels = c("Memory Gene", "Non-memory Gene"))

if(plot.genes == "memory") {
    res.gs <- res.gs[res.gs$Memory == "Memory Gene", ]
}
p <- ggplot(res.gs, aes(x = MGI, y = log2FoldChange, fill = Memory)) +
    # geom_bar(stat="identity", width=0.75, size=1, fill=pal_jama()(1), color=pal_jama()(1)) +
    geom_bar(stat="identity", width=0.75, size=1) +
    scale_fill_manual(values = c(brewer.pal(9, "Reds")[4], pal_jama()(1))) +
    labs(x="", y="log2FoldChange\n(Res499 vs B16)", title = paste0(str_wrap(go_name, width = 50), " (", go_term, ")"), fill = "") +
    theme_bw(base_size = 11) +
    theme(
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x = element_text(size = fontsize, color = "black", hjust = 1, angle = 45),
        axis.title.x = element_text(size=10, color="black"),
        axis.title.y = element_text(size=10, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))
ggsave(p, file = paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/response to interferon-beta (GO:0035456).pdf"), width = 12, height = 4)

# Write out plot data
dat <- p$data
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 4F GO:0035456.csv"), sep=",", quote=F, col.names=T, row.names=F)

# response to virus (GO:0009615)
go_term <- "GO:0009615"
go_name <- "response to virus"
plot.genes <- "memory"
fontsize <- 11

results <- AnnotationDbi::select(
    org.Mm.eg.db, 
    keys=c(go_term), 
    columns = c('SYMBOL'), 
    keytype = "GOALL")
gs <- unique(results$SYMBOL)
gs <- na.omit(gr.bm$ensembl_gene_id[match(gs, gr.bm$mgi_symbol)])

res <- DEGs_list[["B16_v_R499"]]
res.gs <- res[na.omit(match(gs, rownames(res))), ]
res.gs <- res.gs[sort.int(res.gs$log2FoldChange, decreasing=T, index.return=T)$ix, ]
res.gs$MGI <- gr.bm$mgi_symbol[match(rownames(res.gs), gr.bm$ensembl_gene_id)]
res.gs$MGI <- factor(res.gs$MGI, levels = res.gs$MGI)
res.gs$Memory <- factor(ifelse(res.gs$MGI %in% cis_genes, "Memory Gene", "Non-memory Gene"), levels = c("Memory Gene", "Non-memory Gene"))

if(plot.genes == "memory") {
    res.gs <- res.gs[res.gs$Memory == "Memory Gene", ]
}
p <- ggplot(res.gs, aes(x = MGI, y = log2FoldChange, fill = Memory)) +
    # geom_bar(stat="identity", width=0.75, size=1, fill=pal_jama()(1), color=pal_jama()(1)) +
    geom_bar(stat="identity", width=0.75, size=1) +
    scale_fill_manual(values = c(brewer.pal(9, "Reds")[4], pal_jama()(1))) +
    labs(x="", y="log2FoldChange\n(Res499 vs B16)", title = paste0(str_wrap(go_name, width = 50), " (", go_term, ")"), fill = "") +
    theme_bw(base_size = 11) +
    theme(
        axis.text.y=element_text(size=12, color="black"),
        axis.text.x = element_text(size = fontsize, color = "black", hjust = 1, angle = 45),
        axis.title.x = element_text(size=10, color="black"),
        axis.title.y = element_text(size=10, color="black"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.margin = margin(1,1,1,1, "cm"))
ggsave(p, file = paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/response to virus (GO:0009615).pdf"), width = 12, height = 4)

# # Write out plot data
# dat <- p$data
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 4F GO:0009615.csv"), sep=",", quote=F, col.names=T, row.names=F)
