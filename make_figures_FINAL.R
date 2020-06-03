rm(list=ls())

library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggfortify)
library(ggpubr)
library(DiffBind)
library(ggplotify)
library(reshape)
library(rtracklayer)
library(grid)
library(ggsci)
library(ggpubr)
library(dplyr)

dir.path <- "~/Dropbox/Minn/"

source(paste(dir.path, "ATAC_RNA_integration/scripts/mRFAR_v7/v1/ATAC_reanalysis_v1_functions.R", sep=""))
source(paste(dir.path, "ATAC_RNA_integration/scripts/mRFAR_v7/mRFAR_v7_functions.R", sep=""))

##########################
### Plotting functions ###
##########################

draw_colnames_45 <- function (coln, gaps, ...) {
    coord = pheatmap:::find_coordinates(length(coln), gaps)
    x = coord$coord - 0.5 * coord$size
    res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
    return(res)}

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
ns=asNamespace("pheatmap"))

save_pheatmap_png <- function(x, filename, width=1200, height=1200, res = 200) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

p_theme <- theme_bw(base_size = 16) +
    theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 16),
        axis.title.x=element_blank(),
        legend.position="none")

plot <- FALSE

################################
### Load in ATAC-seq samples ###
################################

# Load in peak coordinates df
atac_IDR_counts <- read.table(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/consensus_peaksets/v1/counts/all/consensus_peakset_IDR_tn5_insertion_count.txt", sep=""), header=F, sep="\t")
colnames(atac_IDR_counts) <- c("seqnames", "start", "end", paste("V", 4:23, sep=""))

# Load in DESeq2 normalized counts
atac.dat_idr <- read.table(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/consensus_peaksets/v1/DESeq2/all/consensus_peakset_IDR_tn5_insertion_counts_DESeq2Norm.txt", sep=""), header=T, sep="\t")

# ATAC peak coordinates
gr.atac_IDR <- makeGRangesFromDataFrame(atac_IDR_counts[ ,1:3], starts.in.df.are.0based=FALSE)

###############################
### Load in RNA-seq samples ###
###############################

samples_rna <- c("B16_cas_1", "B16_cas_4", "B16_cas_5", paste("B16_SKO_", 1:5, sep=""), "R499_cas_1", "R499_cas_2", "R499_cas_4", paste("R499_SKO_", 1:5, sep=""))
cond_rna <- data.frame(condition = factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))

dds <- readRDS(file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/dds_mRFAR_v5.rds", sep=""))
rna.dat <- read.table(file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/rna.dat_mRFAR_v5.txt", sep=""), header=T, sep="\t")
rna.anno <- read.table(file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/rna.anno_mRFAR_v5.txt", sep=""), header=T, sep="\t")
chr_names <- paste("chr", rna.anno$chromosome_name, sep="")
rna.anno$chromosome_name <- chr_names
gr.rna <- GRanges(seqnames=chr_names, IRanges(start=rna.anno$start_position, end=rna.anno$end_position))

gs.names <- c("IFN.I", "ISG", "EMT", "TGFB")

gs_list <- lapply(1:length(gs.names), function(x) {
  gs <- read.table(paste(dir.path, "resources/", gs.names[x], "_mouse.txt", sep=""), sep="\t", header=F)
  gs <- unique(as.character(gs$V1))
  gs <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs, rna.anno$mgi_symbol))])
  return(gs)
})

##################################
### Match ATAC and RNA samples ###
##################################

edat <- rna.dat
adat <- atac.dat_idr[ ,intersect(colnames(rna.dat), colnames(atac.dat_idr))]

#########################
### Summary PCA plots ###
#########################

color_palette <- brewer.pal(9, "YlOrRd")
anno_color <- list(condition=brewer.pal(11, "Spectral")[c(10,9,2,4)])
names(anno_color$condition) <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
anno_df <- data.frame(condition=factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))
rownames(anno_df) <- colnames(adat)

# Set up ATAC pca
adat_t <- t(adat)
pca_atac <- prcomp(adat_t, scale=T, center=T)
pca_atac_df <- data.frame(pca_atac$x, id=colnames(adat), condition=factor(c(rep("B16 cas", 3), rep("B16 SKO", 5), rep("R499 cas", 3), rep("R499 SKO", 5))))
levels(pca_atac_df$condition) <- c("B16 cas", "B16 SKO", "R499 cas", "R499 SKO")

# Extract variance explained
eigs <- (pca_atac$sdev)^2
var_exp_pc1 <- round(eigs[1]*100/sum(eigs), 2)
var_exp_pc2 <- round(eigs[2]*100/sum(eigs), 2)

fig.pca_atac <- ggplot(pca_atac_df, aes(x=PC1, y=PC2, col=condition)) +
    geom_point(size=10, alpha=1, shape=18) +
    scale_color_manual(values=c(as.character(anno_color[[1]][1]), as.character(anno_color[[1]][2]), as.character(anno_color[[1]][3]), as.character(anno_color[[1]][4]))) +
    theme_classic() + 
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    xlab(paste("PC1 [", var_exp_pc1, "%]", sep="")) +
    ylab(paste("PC2 [", var_exp_pc2, "%]", sep="")) +
    ggtitle("ATAC") +
    theme(axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA, size=0.8),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          text = element_text(size=22),
          aspect.ratio=1,
          plot.title = element_text("Helvetica"))

# Set up RNA pca
edat_t <- t(edat)
pca_rna <- prcomp(edat_t, scale=T, center=T)
pca_rna_df <- data.frame(pca_rna$x, id=colnames(edat), condition=factor(c(rep("B16 cas", 3), rep("B16 SKO", 5), rep("R499 cas", 3), rep("R499 SKO", 5))))
levels(pca_rna_df$condition) <- c("B16 cas", "B16 SKO", "R499 cas", "R499 SKO")

# Extract variance explained
eigs <- (pca_rna$sdev)^2
var_exp_pc1 <- round(eigs[1]*100/sum(eigs), 2)
var_exp_pc2 <- round(eigs[2]*100/sum(eigs), 2)

fig.pca_rna <- ggplot(pca_rna_df, aes(x=PC1, y=PC2, col=condition)) +
    geom_point(size=10, alpha=1, shape=18) +
    scale_color_manual(values=c(as.character(anno_color[[1]][1]), as.character(anno_color[[1]][2]), as.character(anno_color[[1]][3]), as.character(anno_color[[1]][4]))) +
    theme_classic() + 
    geom_vline(xintercept=0, linetype="dashed") +
    geom_hline(yintercept=0, linetype="dashed") +
    xlab(paste("PC1 [", var_exp_pc1, "%]", sep="")) +
    ylab(paste("PC2 [", var_exp_pc2, "%]", sep="")) +
    ggtitle("RNA") +
    theme(axis.line=element_line(colour="black"),
          panel.border=element_rect(colour="black", fill=NA, size=0.8),
          legend.title=element_text(size=15),
          legend.text=element_text(size=15),
          text = element_text(size=22),
          aspect.ratio=1,
          plot.title = element_text("Helvetica"))

fig.pca_combined <- ggarrange(fig.pca_rna, fig.pca_atac, common.legend = TRUE, legend="right")

if(plot == TRUE) {
  ggsave(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/summary/PCA_plots_combined.pdf", sep=""), device="pdf", plot=fig.pca_combined, height=6, width=13)
}

################
### MA plots ###
################

# MA plot to show fold changes in accessibility between B16 and R499
db.DE <- readRDS(paste(dir.path, "ATAC_RNA_integration/diffbind/v7_tn5/DE/consensus_peakset_DE_IDR.rds", sep=""))

# Use all peaks differentially accessible between B16 and R499
db.DE_B16_v_R499 <- data.frame(dba.report(db.DE, contrast=2, bUsePval=TRUE, th=1))

f1 <- ggplot(db.DE_B16_v_R499, aes(Conc, Fold)) +
    geom_bin2d(bins = 175) +
    scale_fill_gradientn(colours = brewer.pal(9, "YlOrRd")[3:9]) +
    theme_classic() + 
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    xlab("Log2 accessibility") +
    ylab("log2 fold change") +
    ggtitle(paste("Differential accessibility: ", strsplit("B16_v_R499", split="_")[[1]][1], " vs ", strsplit("B16_v_R499", split="_")[[1]][3], "\n (", nrow(db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR <= 0.05), ]), " FDR < 0.05)", sep="")) +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=18),
          text=element_text(size=22),
          legend.text=element_text(size=12), legend.title=element_text(size=12),
          axis.text=element_text(size=18),
          plot.margin=unit(c(0.5,0.5,0.75,0.5), "cm"))

f2 <- ggplot(db.DE_B16_v_R499, aes(Conc, Fold)) +
    geom_bin2d(bins = 300) +
    scale_fill_gradientn(colours = "grey75") +
    theme_classic() + 
    geom_hline(yintercept=0, linetype="dashed", color="black") +
    geom_point(data=db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR <= 0.05), ], col="orangered", size=0.5) + 
    xlab("Log2 accessibility") +
    ylab("log2 fold change") +
    ggtitle(paste("Differentially accessible regions \n between ", strsplit("B16_v_R499", split="_")[[1]][1], " and ", strsplit("B16_v_R499", split="_")[[1]][3], " (", nrow(db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR <= 0.05), ]), " FDR < 0.05)", sep="")) +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=18),
          legend.position="none",
          text=element_text(size=20),
          legend.text=element_text(size=12), legend.title=element_text(size=12),
          axis.text=element_text(size=18),
          plot.margin=unit(c(0.5,0.5,0.75,0.5), "cm"))

if(plot == TRUE) {
    ggsave(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/summary/B16_v_R499_differential_peaks_MAplot1.pdf", sep=""), device="pdf", width=8, height=8, plot=f1)
    ggsave(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/summary/B16_v_R499_differential_peaks_MAplot2.pdf", sep=""), device="pdf", width=8, height=8, plot=f2)
}

#################################################################
### Make bar plot of which condition has higher accessibility ###
#################################################################

db.DE_B16_v_R499_sig <- db.DE_B16_v_R499[which(db.DE_B16_v_R499$FDR <= 0.05), ]

colors <- brewer.pal(11, "Spectral")[c(10,9,2,4)]

df_bar <- data.frame(condition=c(strsplit("B16_v_R499", split="_")[[1]][1], strsplit("B16_v_R499", split="_")[[1]][3]), num_up=c(length(which(db.DE_B16_v_R499_sig$Fold > 0)), length(which(db.DE_B16_v_R499_sig$Fold < 0))))

f3 <- ggplot(df_bar, aes(x=condition, y=num_up, fill=condition)) +
    geom_bar(stat="identity", width=0.65) +
    scale_fill_manual(values=colors[c(1,3)]) +
    theme_classic() +
    geom_text(aes(label=num_up), vjust=-0.3, size=6, fontface="bold") +
    xlab("") +
    ylab("# upregulated accessible regions") +
    theme(plot.title=element_text(hjust=0.5, face="bold", size=18),
          text=element_text(size=18),
          legend.text=element_text(size=12), legend.title=element_text(size=12),
          axis.text.x=element_text(size=18, face="bold"),
          legend.position="none")

if(plot == TRUE) {
    f3_gt <- ggplot_gtable(ggplot_build(f3))
    f3_gt$layout$clip[f3_gt$layout$name == "panel"] <- "off"
    ggsave(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/summary/db.DE_B16_v_R499_differential_peaks_barPlot.pdf", sep=""), device="pdf", width=5, height=5, plot=grid.draw(f3_gt))
}

#########################
### Plot RNA heatmaps ###
#########################

# Plot RNA expression for gene set of interest
cellheights <- c(1.65, 4.9, 0.8, 2.8)
for(i in 1:length(gs_list)) {
  gs <- gs_list[[i]]

  edat.gs <- edat[unique(match(gs, rownames(edat))), ]

  # Plot heatmap
  anno_color <- list(condition=brewer.pal(11, "Spectral")[c(10,9,2,4)])
  names(anno_color$condition) <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
  anno_df <- data.frame(condition=factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))
  rownames(anno_df) <- colnames(adat)

  color_palette <- rev(brewer.pal(11, "RdBu"))
  hm.rna <- pheatmap(t(scale(t(edat.gs))), cluster_cols=FALSE, show_colnames = F, show_rownames=F, annotation_col=anno_df, annotation_colors=anno_color, color=color_palette, border_color=NA, legend=T, cellheight=cellheights[i], cellwidth=8, fontsize=6, gaps_col=c(3,8,11), treeheight_row=0, annotation_legend=F, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.heatmaps/RNA_", gs.names[i], "_heatmap.pdf"))
}

### Boxplot ###

for(i in 1:length(gs_list)) {
  gs <- gs_list[[i]]

  edat.gs <- edat[unique(match(gs, rownames(edat))), ]

  df.rna_meta <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(edat.gs), rna.meta=colMeans(t(scale(t(edat.gs)))))
  df.rna_meta_summary <- df.rna_meta %>% 
    group_by(condition) %>% 
    summarize(rna.mean=mean(rna.meta))

  fig.rna <- ggplot(df.rna_meta, aes(x=condition, y=rna.meta, color=factor(condition))) + 
    # geom_boxplot() +
    geom_jitter(position=position_jitter(0.1), size=4) +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6) +
    scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
    ylab("Average IFN.I expression") +
    theme_bw() +
    theme(text=element_text(size=18),
      axis.text.x=element_text(size=18),
      legend.position="none")
  ggsave(fig.rna, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/RNA_", gs.names[i], "_boxplot.pdf"), width=3, height=2.5)
}

### Get p-values ###
# * 0.1
# ** 0.01
# *** 0.001
wilcox.test(df.rna_meta[which(df.rna_meta$condition == "B16"), 3], df.rna_meta[which(df.rna_meta$condition == "B16_SKO"), 3])
wilcox.test(df.rna_meta[which(df.rna_meta$condition == "B16"), 3], df.rna_meta[which(df.rna_meta$condition == "R499"), 3], alternative="less")
wilcox.test(df.rna_meta[which(df.rna_meta$condition == "R499"), 3], df.rna_meta[which(df.rna_meta$condition == "R499_SKO"), 3])
wilcox.test(df.rna_meta[which(df.rna_meta$condition == "B16_SKO"), 3], df.rna_meta[which(df.rna_meta$condition == "R499_SKO"), 3])

### Oas1 plot ###

gs.Oas1 <- as.character(rna.anno$ensembl_gene_id)[grep(rna.anno$mgi_symbol, pattern="Oas1")]

edat.gs <- edat[unique(match(gs.Oas1, rownames(edat))), ]
rownames(edat.gs) <- as.character(rna.anno$mgi_symbol)[match(gs.Oas1, as.character(rna.anno$ensembl_gene_id))]

# Plot heatmap
anno_color <- list(condition=brewer.pal(11, "Spectral")[c(10,9,2,4)])
names(anno_color$condition) <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
anno_df <- data.frame(condition=factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))
rownames(anno_df) <- colnames(adat)

color_palette <- rev(brewer.pal(11, "RdBu"))
hm.rna <- pheatmap(t(scale(t(edat.gs))), cluster_cols=FALSE, show_colnames = F, show_rownames=T, annotation_col=anno_df, annotation_colors=anno_color, color=color_palette, border_color=NA, legend=T, cellheight=8, cellwidth=8, fontsize=8, gaps_col=c(3,8,11), treeheight_row=0, annotation_legend=F, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.heatmaps/RNA_Oas1_heatmap.pdf"))

### Random genes ###

gs.rand_names <- c("rand1", "rand2", "rand3")
gs.rand <- lapply(1:3, function(x) {
  gs <- read.table(file=paste0(dir.path, "resources/random_gene_set_", x, "_mouse.txt"), sep="\t", header=F, stringsAsFactors=F)
  gs <- as.character(gs$V1)
  return(gs)
})

anno_color <- list(condition=brewer.pal(11, "Spectral")[c(10,9,2,4)])
names(anno_color$condition) <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
anno_df <- data.frame(condition=factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))
rownames(anno_df) <- colnames(adat)

color_palette <- rev(brewer.pal(11, "RdBu"))

lapply(1:3, function(x) {
  edat.gs <- edat[unique(match(gs.rand[[x]], rownames(edat))), ]

  # Plot heatmap
  hm.rna <- pheatmap(t(scale(t(edat.gs))), cluster_cols=FALSE, show_colnames = F, show_rownames=F, annotation_col=anno_df, annotation_colors=anno_color, color=color_palette, border_color=NA, legend=T, cellheight=0.8, cellwidth=8, fontsize=6, gaps_col=c(3,8,11), treeheight_row=0, annotation_legend=F, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.heatmaps/RNA_", gs.rand_names[x], "_heatmap.pdf"))

  # Boxplot
  df.rna_meta <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(edat.gs), rna.meta=colMeans(t(scale(t(edat.gs)))))
  df.rna_meta_summary <- df.rna_meta %>% 
    group_by(condition) %>% 
    summarize(rna.mean=mean(rna.meta))

  # Get pvals
  t.b16_v_b16sko <- wilcox.test(df.rna_meta$rna.meta[which(df.rna_meta$condition == "B16")], df.rna_meta$rna.meta[which(df.rna_meta$condition == "B16_SKO")])
  p.b16_v_b16sko <- t.b16_v_b16sko$p.value

  t.b16_v_r499 <- wilcox.test(df.rna_meta$rna.meta[which(df.rna_meta$condition == "B16")], df.rna_meta$rna.meta[which(df.rna_meta$condition == "R499")])
  p.b16_v_r499 <- t.b16_v_r499$p.value

  t.r499_v_r499sko <- wilcox.test(df.rna_meta$rna.meta[which(df.rna_meta$condition == "R499")], df.rna_meta$rna.meta[which(df.rna_meta$condition == "R499_SKO")])
  p.r499_v_r499sko <- t.r499_v_r499sko$p.value

  t.b16sko_v_r499sko <- wilcox.test(df.rna_meta$rna.meta[which(df.rna_meta$condition == "B16_SKO")], df.rna_meta$rna.meta[which(df.rna_meta$condition == "R499_SKO")])
  p.b16sko_v_r499sko <- t.b16sko_v_r499sko$p.value

  p.annotate <- paste0("B16 v B16 SKO: ", round(p.b16_v_b16sko, 2), "\n", "B16 v R499: ", round(p.b16_v_r499, 2), "\n", "R499 v R499 SKO: ", round(p.r499_v_r499sko, 2), "\n", "B16 SKO v R499 SKO: ", round(p.b16sko_v_r499sko, 2))

  fig.rna <- ggplot(df.rna_meta, aes(x=condition, y=rna.meta, color=factor(condition))) + 
    geom_jitter(position=position_jitter(0.1), size=4) +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6) +
    scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
    ylab("Average random genes expression") +
    theme_bw() +
    theme(text=element_text(size=18),
      axis.text.x=element_text(size=18),
      legend.position="none")
  ggsave(fig.rna, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/RNA_", gs.rand_names[x], "_boxplot.pdf"), width=3, height=2.5)
})

############################################
### Annotate peaks as promoter or distal ###
############################################

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(org.Mm.eg.db)

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

# Annotate peaks
# Downstream is defined as the downstream of gene end
peakAnno <- annotatePeak(gr.atac_IDR, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Mm.eg.db")

# # Visualize annotations
# plotAnnoBar(peakAnno)
# plotDistToTSS(peakAnno, title="Distance to TSS")

# gs.IFN.I <- read.table(paste(dir.path, "resources/IFN.I_mouse.txt", sep=""), sep="\t", header=F)
# gs.IFN.I <- unique(as.character(gs.IFN.I$V1))

# Get peaks that are annotated as promoters
peakAnno.promoter <- peakAnno@anno[grep(peakAnno@anno$annotation, pattern="Promoter")]
peakAnno.promoter$peakid <- grep(peakAnno@anno$annotation, pattern="Promoter")

#####################################

conditions <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")

name <- peaks <- "rand2"

gs.df_peak <- read.table(file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/peaks_of_interest/v1/", name, "/peaks/peaks_", name, "_ALL.txt"), sep="\t", header=T, stringsAsFactors=F)

# gs <- as.character(rna.anno$mgi_symbol)[match(gs_list[[which(gs.names == name)]], as.character(rna.anno$ensembl_gene_id))]
gs <- as.character(rna.anno$mgi_symbol)[match(gs.rand[[which(gs.rand_names == name)]], as.character(rna.anno$ensembl_gene_id))]

# Annotate promoter peaks
promoter_ids <- lapply(1:length(gs), function(x) {
  gene <- gs[x]
  idx <- match(gene, peakAnno.promoter$SYMBOL)
  if(!is.na(idx)) {
    gr <- peakAnno.promoter[idx]
    return(gr$peakid)
  } else {
    return(NA)
  }
})
names(promoter_ids) <- gs
promoter_ids <- unlist(promoter_ids)

# Distinguish top distal elements from other distal elements
gs.df_peak_filt <- gs.df_peak
distal1_ids <- lapply(1:length(gs), function(x) {
  if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
    df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
    df <- df[!df$peak_ID %in% promoter_ids, ]
    return(as.numeric(df$peak_ID[which(df$VIMP_rf == max(df$VIMP_rf))]))
  } else {
    return(NA)
  }  
})
names(distal1_ids) <- gs
distal1_ids <- unlist(distal1_ids)

df.anno <- data.frame(promoter=promoter_ids, distal1=distal1_ids)

distal2_ids <- lapply(1:length(gs), function(x) {
  ids <- as.numeric(df.anno[which(rownames(df.anno) == gs[x]), ])
  
  if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
    df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
    df <- df[!df$peak_ID  %in% na.omit(ids), ]
    distal2 <- df$peak_ID[which(df$VIMP_rf == max(df$VIMP_rf))]
    return(distal2)
  } else {
    return(NA)
  } 
})
names(distal2_ids) <- gs
distal2_ids <- unlist(distal2_ids)

df.anno$distal2 <- distal2_ids

vimp_cutoff <- 2.5
distal_ids <- lapply(1:length(gs), function(x) {
  
  if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
    df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
    df <- df[!df$peak_ID %in% promoter_ids, ]
    df <- df[which(df$VIMP_rf > vimp_cutoff), ]
    distal <- df$peak_ID
    return(distal)
  } else {
    return(NA)
  } 
})
names(distal_ids) <- gs

gs.df_peak_filt$name <- paste0(gs.df_peak_filt$gene, "_", gs.df_peak_filt$peak_ID)

# Subset by regulatory element, order by VIMP
promoter.names <- paste0(rownames(df.anno), "_", df.anno[,1])
df.promoter <- gs.df_peak_filt[gs.df_peak_filt$name %in% promoter.names, ]
df.promoter <- df.promoter[sort.int(df.promoter$VIMP_rf, decreasing=T, index.return=T)$ix, ]

distal1.names <- paste0(rownames(df.anno), "_", df.anno[,2])
df.distal1 <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal1.names, ]
df.distal1 <- df.distal1[sort.int(df.distal1$VIMP_rf, decreasing=T, index.return=T)$ix, ]

distal2.names <- paste0(rownames(df.anno), "_", df.anno[,3])
df.distal2 <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal2.names, ]
df.distal2 <- df.distal2[sort.int(df.distal2$VIMP_rf, decreasing=T, index.return=T)$ix, ]

df.distal <- lapply(1:length(gs), function(x) {
  if(length(distal_ids[[x]]) != 0 && !is.na(distal_ids[[x]])) {
    distal.names <- paste0(gs[x], "_", distal_ids[[x]])
    df <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal.names, ]
    return(df)
  }
})
df.distal <- do.call(rbind, df.distal)
df.distal <- df.distal[sort.int(df.distal$VIMP_rf, decreasing=T, index.return=T)$ix, ]

##########################
### Plot ATAC heatmaps ###
##########################

### Plot promoter and distal peaks separately ###

cellheight <- 0.5

elements <- c("PROMOTER", "DISTAL")

for(i in 1:length(elements)) {
  element <- elements[i]
  if(element == "PROMOTER") {
    adat.gs <- adat[df.promoter$peak_ID, ]
  #     cellheight <- 5
  } else if(element == "DISTAL1") {
      adat.gs <- adat[df.distal1$peak_ID, ]
  #     cellheight <- 2.25
  } else if(element == "DISTAL") {
      adat.gs <- adat[unique(df.distal$peak_ID), ]
  #     cellheight <- 2
  } else if(element == "combined") {
    adat.gs <- adat[unique(c(df.promoter$peak_ID, df.distal$peak_ID)), ]
  }
                                        
  # Plot ATAC accessibility with Stat1 and Irf motif annotation
  anno_color <- list(condition=brewer.pal(11, "Spectral")[c(10,9,2,4)], Stat1.peak=c("black", "white"), Irf3.peak=c("black", "white"))
  names(anno_color$condition) <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")
  # names(anno_color$Stat1.peak) <- c("present", "no")
  # names(anno_color$Irf3.peak) <- c("present", "no")
  anno_df <- data.frame(condition=factor(c(rep("B16_cas", 3), rep("B16_SKO", 5), rep("R499_cas", 3), rep("R499_SKO", 5))))
  # anno_df_row <- data.frame(Stat1.peak=Stat1.peak, Irf3.peak=Irf3.peak)
  rownames(anno_df) <- colnames(adat)
  # rownames(anno_df_row) <- rownames(adat.gs)

  # Plot
  # color_palette <- rev(brewer.pal(11, "RdBu"))
  # color_palette <- rev(brewer.pal(11, "RdGy"))
  color_palette <- rev(brewer.pal(11, "Spectral"))

  # hm <- pheatmap(t(scale(t(adat.gs))), cluster_cols=FALSE, cluster_rows=T, show_colnames = F, show_rownames=F, annotation_col=anno_df, annotation_colors=anno_color, color=color_palette, border_color=NA, legend=F, cellheight=cellheight, cellwidth=8, fontsize=6, gaps_col=c(3,8,11), treeheight_row=0)

  # Promoter or distal separately
  hm <- pheatmap(t(scale(t(adat.gs))), cluster_cols=FALSE, cluster_rows=T, show_colnames = F, show_rownames=F, annotation_col=anno_df, annotation_colors=anno_color, color=color_palette, border_color=NA, legend=T, cellheight=cellheight, cellwidth=8, fontsize=6, gaps_col=c(3,8,11), treeheight_row=0, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.heatmaps/ATAC_impt_RF_", element, "_", name, "_heatmap.pdf"))
}

### Boxplot ###

# Promoter peaks
adat.gs_promoter <- adat[df.promoter$peak_ID, ]
df.atac_meta_promoter <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(adat.gs_promoter), atac.meta=colMeans(t(scale(t(adat.gs_promoter)))))

fig.atac_promoter <- ggplot(df.atac_meta_promoter, aes(x=condition, y=atac.meta, color=factor(condition))) + 
  # geom_boxplot() +
  geom_jitter(position=position_jitter(0.1), size=4) +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6, size=1) +
  scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
  ylab("Average IFN.I expression") +
  theme_bw() +
  theme(text=element_text(size=12),
    axis.text.x=element_text(size=12),
    legend.position="none")

# Distal peaks
adat.gs_distal <- adat[unique(df.distal$peak_ID), ]
df.atac_meta_distal <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(adat.gs_distal), atac.meta=colMeans(t(scale(t(adat.gs_distal)))))

fig.atac_distal <- ggplot(df.atac_meta_distal, aes(x=condition, y=atac.meta, color=factor(condition))) + 
  # geom_boxplot() +
  geom_jitter(position=position_jitter(0.1), size=4) +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6, size=1) +
    scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
    ylab("Average IFN.I expression") +
    theme_bw() +
    theme(text=element_text(size=12),
      axis.text.x=element_text(size=12),
      legend.position="none")

fig.atac <- ggarrange(fig.atac_promoter, fig.atac_distal, ncol=1)

ggsave(fig.atac, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/ATAC_impt_RF_", name, "_boxplot.pdf"), width=4, height=3)

#############################

adat.gs_combined <- rbind(adat.gs_promoter, adat.gs_distal)

df.atac_meta_combined <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(adat.gs_combined), atac.meta=colMeans(t(scale(t(adat.gs_combined)))))

fig.atac_combined <- ggplot(df.atac_meta_combined, aes(x=condition, y=atac.meta, color=factor(condition))) + 
  # geom_boxplot() +
  geom_jitter(position=position_jitter(0.1), size=4) +
    stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6, size=1) +
    scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
    ylab("Average IFN.I expression") +
    theme_bw() +
    theme(text=element_text(size=12),
      axis.text.x=element_text(size=12),
      legend.position="none")
ggsave(fig.atac_combined, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/ATAC_impt_RF_combined_", name, "_boxplot.pdf"), width=3, height=2.5)

# Get p-values
wilcox.test(df.atac_meta_combined[which(df.atac_meta_combined$condition == "B16"), 3], df.atac_meta_combined[which(df.atac_meta_combined$condition == "B16_SKO"), 3])
wilcox.test(df.atac_meta_combined[which(df.atac_meta_combined$condition == "B16"), 3], df.atac_meta_combined[which(df.atac_meta_combined$condition == "R499"), 3])
wilcox.test(df.atac_meta_combined[which(df.atac_meta_combined$condition == "R499"), 3], df.atac_meta_combined[which(df.atac_meta_combined$condition == "R499_SKO"), 3])
wilcox.test(df.atac_meta_combined[which(df.atac_meta_combined$condition == "B16_SKO"), 3], df.atac_meta_combined[which(df.atac_meta_combined$condition == "R499_SKO"), 3])

###############################
### Random gene sets (n=20) ###
###############################

### Take average of random gene sets ###

### RNA ###

n <- 20
gs.random_RNA <- lapply(1:n, function(x) {
  gs <- read.table(file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/random_gene_sets/random_gene_set_", x, "_mouse.txt"), sep="\t", header=F, stringsAsFactors=F)
  gs <- as.character(gs$V1)
  edat.gs <- edat[unique(match(gs, rownames(edat))), ]
  gs.meta <- as.numeric(colMeans(t(scale(t(edat.gs)))))
  df.meta <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(edat.gs), rna.meta=gs.meta)
  return(df.meta)
})
gs.random_RNA <- do.call(rbind, gs.random_RNA)

# Boxplot
gs.random_RNA_meta <- gs.random_RNA %>%
  group_by(sample) %>%
  summarize(rna.meta=mean(rna.meta))
gs.random_RNA_meta$condition <- factor(sapply(strsplit(as.character(gs.random_RNA_meta$sample), split="_"), function(x) paste(x[1:2], collapse="_")))

fig.rna <- ggplot(gs.random_RNA_meta, aes(x=condition, y=rna.meta, color=factor(condition))) +
  geom_jitter(position=position_jitter(0.1), size=4) +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6) +
  scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
  ylab("RNA expression (n=20 random gene sets)") +
  theme_bw() +
  theme(text=element_text(size=12),
    axis.text.x=element_text(size=12),
    legend.position="none")
ggsave(fig.rna, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/random/RNA_20_random_gene_sets_boxplot.pdf"), width=3, height=2.5)

### ATAC ###

conditions <- c("B16_cas", "B16_SKO", "R499_cas", "R499_SKO")

gs.random_ATAC <- lapply(1:n, function(x) {
  name <- paste0("rand", x)
  print(name)

  gs.df_peak <- read.table(file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/peaks_of_interest/v1/random/", name, "/peaks/peaks_", name, "_ALL.txt"), sep="\t", header=T, stringsAsFactors=F)
  gs <- read.table(file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/random_gene_sets/random_gene_set_", x, "_mouse.txt"), sep="\t", header=F, stringsAsFactors=F)
  gs <- as.character(gs$V1)
  gs <- as.character(rna.anno$mgi_symbol)[match(gs, as.character(rna.anno$ensembl_gene_id))]

  # Annotate promoter peaks
  promoter_ids <- lapply(1:length(gs), function(x) {
    gene <- gs[x]
    idx <- match(gene, peakAnno.promoter$SYMBOL)
    if(!is.na(idx)) {
      gr <- peakAnno.promoter[idx]
      return(gr$peakid)
    } else {
      return(NA)
    }
  })
  names(promoter_ids) <- gs
  promoter_ids <- unlist(promoter_ids)

  # Distinguish top distal elements from other distal elements
  gs.df_peak_filt <- gs.df_peak
  distal1_ids <- lapply(1:length(gs), function(x) {
    if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
      df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
      df <- df[!df$peak_ID %in% promoter_ids, ]
      return(as.numeric(df$peak_ID[which(df$VIMP_rf == max(df$VIMP_rf))]))
    } else {
      return(NA)
    }  
  })
  names(distal1_ids) <- gs
  distal1_ids <- unlist(distal1_ids)

  df.anno <- data.frame(promoter=promoter_ids, distal1=distal1_ids)

  distal2_ids <- lapply(1:length(gs), function(x) {
    ids <- as.numeric(df.anno[which(rownames(df.anno) == gs[x]), ])
    
    if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
      df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
      df <- df[!df$peak_ID  %in% na.omit(ids), ]
      if(nrow(df) == 0) {
        distal2 <- NA
      } else {
        distal2 <- df$peak_ID[which(df$VIMP_rf == max(df$VIMP_rf))]
      }
      return(distal2)
    } else {
      return(NA)
    } 
  })
  names(distal2_ids) <- gs
  distal2_ids <- unlist(distal2_ids)

  df.anno$distal2 <- distal2_ids

  vimp_cutoff <- 2.5
  distal_ids <- lapply(1:length(gs), function(x) {
    
    if(length(which(gs.df_peak_filt$gene == gs[x])) > 0) {
      df <- gs.df_peak_filt[which(gs.df_peak_filt$gene == gs[x]), ]
      df <- df[!df$peak_ID %in% promoter_ids, ]
      df <- df[which(df$VIMP_rf > vimp_cutoff), ]
      distal <- df$peak_ID
      return(distal)
    } else {
      return(NA)
    } 
  })
  names(distal_ids) <- gs

  gs.df_peak_filt$name <- paste0(gs.df_peak_filt$gene, "_", gs.df_peak_filt$peak_ID)

  # Subset by regulatory element, order by VIMP
  promoter.names <- paste0(rownames(df.anno), "_", df.anno[,1])
  df.promoter <- gs.df_peak_filt[gs.df_peak_filt$name %in% promoter.names, ]
  df.promoter <- df.promoter[sort.int(df.promoter$VIMP_rf, decreasing=T, index.return=T)$ix, ]

  distal1.names <- paste0(rownames(df.anno), "_", df.anno[,2])
  df.distal1 <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal1.names, ]
  df.distal1 <- df.distal1[sort.int(df.distal1$VIMP_rf, decreasing=T, index.return=T)$ix, ]

  distal2.names <- paste0(rownames(df.anno), "_", df.anno[,3])
  df.distal2 <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal2.names, ]
  df.distal2 <- df.distal2[sort.int(df.distal2$VIMP_rf, decreasing=T, index.return=T)$ix, ]

  df.distal <- lapply(1:length(gs), function(x) {
    if(length(distal_ids[[x]]) != 0 && !is.na(distal_ids[[x]])) {
      distal.names <- paste0(gs[x], "_", distal_ids[[x]])
      df <- gs.df_peak_filt[gs.df_peak_filt$name %in% distal.names, ]
      return(df)
    }
  })
  df.distal <- do.call(rbind, df.distal)
  df.distal <- df.distal[sort.int(df.distal$VIMP_rf, decreasing=T, index.return=T)$ix, ]

  ### Plot combined peaks ###
  adat.gs_promoter <- adat[df.promoter$peak_ID, ]
  adat.gs_distal <- adat[unique(df.distal$peak_ID), ]
  adat.gs_combined <- rbind(adat.gs_promoter, adat.gs_distal)

  df.atac_meta_combined <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(adat.gs_combined), atac.meta=colMeans(t(scale(t(adat.gs_combined)))))
  return(df.atac_meta_combined)
})
gs.random_ATAC <- do.call(rbind, gs.random_ATAC)

# Boxplot
gs.random_ATAC_meta <- gs.random_ATAC %>%
  group_by(sample) %>%
  summarize(atac.meta=mean(atac.meta))
gs.random_ATAC_meta$condition <- factor(sapply(strsplit(as.character(gs.random_ATAC_meta$sample), split="_"), function(x) paste(x[1:2], collapse="_")))

fig.atac_combined <- ggplot(gs.random_ATAC_meta, aes(x=condition, y=atac.meta, color=factor(condition))) + 
  geom_jitter(position=position_jitter(0.1), size=4) +
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6, size=1) +
  scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
  ylab("ATAC accessibility (n=200 random gene sets)") +
  theme_bw() +
  theme(text=element_text(size=12),
    axis.text.x=element_text(size=12),
    legend.position="none")
ggsave(fig.atac_combined, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.boxplots/random/ATAC_20_random_gene_sets_boxplot.pdf"), width=3, height=2.5)

######################
### Motif analysis ###
######################

bm <- read.table(paste0(dir.path, "resources/mouse_biomart_annotations_2-20-20.txt"), sep="\t", header=T, stringsAsFactors=F)

color_palette <- brewer.pal(n=9, name="Blues")[c(1,9)]

name <- "EMT"

fimo <- load_fimo(fimo_dir.path=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/fimo/v1/peaks_of_interest/"), name=paste0(name, "_ALL"), mouse2human=mouse2human, sig=0.05, score_cutoff=5)
gr.motifs <- get_motif_occurrences(fimo, motif_min=2)
gr.motifs <- remove_unexpressedTFs(gr.motifs, rna.anno, species="mouse")

consensus_TFs <- names(gr.motifs)
consensus_TFs_ensembl <- bm$ensembl_gene_id[match(consensus_TFs, bm$mgi_symbol)]

### Filter by correlation ###

TF.cor <- sapply(1:length(consensus_TFs), function(x) {
  gr <- gr.motifs[[x]]
  atac <- as.numeric(colMeans(t(scale(t(adat[unique(gr$peak_id), ])))))
  rna <- as.numeric(t(scale(t(edat[match(consensus_TFs_ensembl[x], rownames(edat)), ]))))
  return(cor(rna, atac))
})
consensus_TFs <- consensus_TFs[which(TF.cor > 0.1)]
# consensus_TFs[which(abs(TF.cor) > 0.1)]

### Plot TF motifs in important ATAC peaks ###

### Promoter peaks ###

promoter_ids <- df.promoter$peak_ID

mat.promoter_TFs <- lapply(1:length(promoter_ids), function(x) {
  TF.present <- ifelse(sapply(1:length(consensus_TFs), function(y) promoter_ids[x] %in% gr.motifs[[which(names(gr.motifs) == consensus_TFs[y])]]$peak_id), 1, 0)
  return(TF.present)
})
mat.promoter_TFs <- do.call(rbind, mat.promoter_TFs)
colnames(mat.promoter_TFs) <- consensus_TFs

mat <- mat.promoter_TFs[,sort.int(colSums(mat.promoter_TFs), decreasing=T, index.return=T)$ix]
remove_idx <- which(colSums(mat) < 3)
mat <- mat[,-remove_idx]
# pheatmap(t(mat), cellwidth=7, cellheight=7, cluster_cols = T, cluster_rows=F, color=color_palette, fontsize_col=6, fontsize_row=7, show_colnames = T, show_rownames=T, legend = F, clustering_method = "ward.D2")

### Distal peaks ###

distal_ids <- df.distal$peak_ID

mat.distal_TFs <- lapply(1:length(distal_ids), function(x) {
  TF.present <- ifelse(sapply(1:length(consensus_TFs), function(y) distal_ids[x] %in% gr.motifs[[which(names(gr.motifs) == consensus_TFs[y])]]$peak_id), 1, 0)
  return(TF.present)
})
mat.distal_TFs <- do.call(rbind, mat.distal_TFs)
colnames(mat.distal_TFs) <- consensus_TFs

mat <- mat.distal_TFs[,sort.int(colSums(mat.distal_TFs), decreasing=T, index.return=T)$ix]
remove_idx <- which(colSums(mat) < 3)
mat <- mat[,-remove_idx]
# pheatmap(t(mat), cellwidth=1.5, cellheight=7, cluster_cols = T, cluster_rows=F, color=color_palette, fontsize_col=6, fontsize_row=7.5, show_colnames = T, show_rownames=T, legend = F, clustering_method = "ward.D2")

### Combined peaks ###

combined_ids <- unique(c(promoter_ids, distal_ids))

mat.combined_TFs <- lapply(1:length(combined_ids), function(x) {
  TF.present <- ifelse(sapply(1:length(consensus_TFs), function(y) combined_ids[x] %in% gr.motifs[[which(names(gr.motifs) == consensus_TFs[y])]]$peak_id), 1, 0)
  return(TF.present)
})
mat.combined_TFs <- do.call(rbind, mat.combined_TFs)
colnames(mat.combined_TFs) <- consensus_TFs

mat <- mat.combined_TFs[,sort.int(colSums(mat.combined_TFs), decreasing=T, index.return=T)$ix]
remove_idx <- which(colSums(mat) < 3)
mat <- mat[,-remove_idx]
# pheatmap(t(mat), cellwidth=1, cellheight=9, cluster_cols = T, cluster_rows=F, color=color_palette, fontsize_col=6, fontsize_row=7.5, show_colnames = T, show_rownames=T, legend = F, clustering_method = "ward.D2")

remove_idx <- which(rowSums(mat) < 1)
mat <- mat[-remove_idx, ]

# Stat1_idx <- which(colnames(mat) == "Stat1")
# Stat1 <- ifelse(mat[,which(colnames(mat) == "Stat1")] == 1, 2, 0)

# mat <- mat[ ,-Stat1_idx]
# mat <- cbind(Stat1, mat)

color_palette <- c(brewer.pal(n=9, name="Blues")[1], brewer.pal(n=9, name="Blues")[9], "firebrick")
# color_palette <- c("white", "black", "firebrick")

if(name == "IFN.I") {
  cellwidth <- 1
} else if(name == "ISG") {
  cellwidth <- 2
}
pheatmap(t(mat), cellwidth=cellwidth, cellheight=8.5, cluster_cols = T, cluster_rows=F, color=color_palette, fontsize_col=6, fontsize_row=8, show_colnames = T, show_rownames=T, legend = F, treeheight_row = 15, treeheight_col = 15, filename=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/gs.motifs/motifs in impt ", name, " ATAC peaks.pdf"))

####################################################################
### Predict gene expression with putative cis-regulatory regions ###
####################################################################

fit.rsq <- lapply(1:length(gs), function(x) {
  gs.ensembl <- as.character(rna.anno$ensembl_gene_id)[which(as.character(rna.anno$mgi_symbol) == gs[x])]

  adat.p <- adat[df.promoter$peak_ID[which(df.promoter$gene == gs[x])], ]
  adat.d <- adat[df.distal$peak_ID[which(df.distal$gene == gs[x])], ]

  edat.x <- edat[which(as.character(rna.anno$ensembl_gene_id) == gs.ensembl), ]
  rna.x <- as.numeric(scale(t(edat.x)))

  if(nrow(adat.p) == 0) {
    atac.p <- NA
  } else if(nrow(adat.p) == 1) {
    atac.p <- as.numeric(scale(t(adat.p)))
  } else if(nrow(adat.p) > 1) {
    atac.p <- as.numeric(colMeans(t(scale(t(adat.p)))))
  }

  if(nrow(adat.d) == 0) {
    atac.d <- NA
  } else if(nrow(adat.d) == 1) {
    atac.d <- as.numeric(scale(t(adat.d)))
  } else if(nrow(adat.d) > 1) {
    atac.d <- as.numeric(colMeans(t(scale(t(adat.d)))))
  }

  df <- data.frame(rna=rna.x, atac.p=atac.p, atac.d=atac.d)
  # df$atac.p.d <- rowMeans(df[,2:3])

  if(!is.na(atac.p[1]) && !is.na(atac.d[1])) {
    fit1 <- lm(formula = rna ~ atac.p, data=df)
    rsq1 <- summary(fit1)$adj.r.squared

    fit2 <- lm(formula = rna ~ atac.d, data=df)
    rsq2 <- summary(fit2)$adj.r.squared

    fit3 <- lm(formula = rna ~ atac.p + atac.d, data=df)
    rsq3 <- summary(fit3)$adj.r.squared
  } else {
    rsq1 <- NA
    rsq2 <- NA
    rsq3 <- NA
  }

  return(c(rsq1, rsq2, rsq3))
})

# p = 2.3e-4 (IFN.I promoter vs promoter+distal)
df.fit.rsq <- do.call(rbind, fit.rsq)
colnames(df.fit.rsq) <- c("promoter", "distal", "promoter + distal")
rownames(df.fit.rsq) <- gs
df.fit.rsq <- na.omit(df.fit.rsq)

df.fit.rsq.m <- melt(df.fit.rsq)
df.fit.rsq.m$X2 <- factor(df.fit.rsq.m$X2, levels=c("promoter", "distal", "promoter + distal"))

fig <- ggplot(df.fit.rsq.m[df.fit.rsq.m$X2 %in% c("promoter", "promoter + distal"), ], aes(x=X2, y=value, color=X2)) +
  geom_boxplot(size=1) +
  geom_jitter(shape=1, size=2.5, , stroke=1.25, position=position_jitter(0.15)) +
  scale_color_jco() +
  ylab("Adjusted R-squared") +
  theme_bw() +
  p_theme
ggsave(fig, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/promoter distal peak models Rsq.pdf"), width=4, height=4)



# ###############################################################
# ### Plot accessibility at regions with TF motif of interest ###
# ###############################################################

# TF <- "Sp2"

# consensus_peakids_TF <- gr.motifs[which(names(gr.motifs) == TF)][[1]]$peak_id

# adat.gs_TF <- adat[unique(consensus_peakids_TF), ]
# df.atac_meta_TF <- data.frame(condition=c(rep("B16", 3), rep("B16_SKO", 5), rep("R499", 3), rep("R499_SKO", 5)), sample=colnames(adat.gs_TF), atac.meta=colMeans(t(scale(t(adat.gs_TF)))))

# fig.atac_TF <- ggplot(df.atac_meta_TF, aes(x=condition, y=atac.meta, color=factor(condition))) + 
#   # geom_boxplot() +
#   geom_jitter(position=position_jitter(0.1), size=4) +
#     stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.6, size=1) +
#     scale_color_manual(values=brewer.pal(11, "Spectral")[c(10,9,2,4)]) +
#     ylab("Average IFN.I expression") +
#     theme_bw() +
#     theme(text=element_text(size=12),
#       axis.text.x=element_text(size=12),
#       legend.position="none")
# fig.atac_TF

# # # Get motif matches in gs peaks of interest
# # if(peak_type == "full_peak") {
# #     fimo <- load_fimo(fimo_dir.path=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/fimo/v1/peaks_of_interest/"), name=paste0(name, "_ALL"), mouse2human=mouse2human, sig=0.05, score_cutoff=5)
# #     gr.motifs <- get_motif_occurrences(fimo, motif_min=5)
# #     gr.motifs <- remove_unexpressedTFs(gr.motifs, rna.anno, species="mouse")
# # } else {
# #     gr.motifs_list <- lapply(1:length(conditions), function(x) {
# # #         fimo <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/fimo/v1/peaks_of_interest_fp/ALL_genes_ALL/fimo/", conditions[x], "_ALL_genes_ALL.rds"))
# #         fimo <- load_fimo(fimo_dir.path=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/fimo/v1/peaks_of_interest_fp/"), name=paste0(conditions[x], "_IFN.I_ALL"), mouse2human=mouse2human, sig=0.05, score_cutoff=5)
# #         gr.motifs <- get_motif_occurrences(fimo, motif_min=8)
# #         gr.motifs <- remove_unexpressedTFs(gr.motifs, rna.anno, species="mouse")
# #     })
# # }

# # # Get peaks with TF motif (Stat1, Irf)
# # if(peak_type == "full_peak") {    
# #     consensus_peakids_Stat1 <- gr.motifs[which(names(gr.motifs) == "Stat1")][[1]]$peak_id
# #     consensus_peakids_Irf3 <- gr.motifs[which(names(gr.motifs) == "Irf3")][[1]]$peak_id
# # } else {
# #     consensus_peakids_Stat1 <- unique(unlist(sapply(gr.motifs_list, function(x) x[[which(names(x) == "Stat1")]]$peak_id)))
# #     consensus_peakids_Irf3 <- unique(unlist(sapply(gr.motifs_list, function(x) x[[which(names(x) == "Irf3")]]$peak_id)))
# # }

# # Stat1.peak <- ifelse(as.numeric(sapply(rownames(adat.gs), function(x) as.numeric(substr(x, start=5, stop=nchar(x))))) %in% consensus_peakids_Stat1, "present", "no")
# # paste0("Percentage of peaks with Stat1 motif: ", round(length(which(Stat1.peak == "present"))/length(Stat1.peak), 2))

# # Irf3.peak <- ifelse(as.numeric(sapply(rownames(adat.gs), function(x) as.numeric(substr(x, start=5, stop=nchar(x))))) %in% consensus_peakids_Irf3, "present", "no")
# # paste0("Percentage of peaks with Irf3 motif: ", round(length(which(Irf3.peak == "present"))/length(Irf3.peak), 2))

########################################
### Rank genes by variance explained ###
########################################

library(ggrepel)
library(ggrastr)

mRFARobj_ALL <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/v1/mRFARobj_ALL_genes_unscaledVIMPs.rds"))

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
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          plot.title = element_text(size = 10))
ggsave(fig, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/final/v3/All_genes_variance_explained_ranked.pdf"), width=4, height=3)
