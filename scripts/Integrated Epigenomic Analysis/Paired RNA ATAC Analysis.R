### Paired RNA/ATAC analysis (original) ###

rm(list=ls())

library(GenomicRanges)
library(tidyverse)
library(matrixStats)
library(DiffBind)
library(ggsci)
library(ggpubr)
library(ggrastr)
library(ggrepel)

pal <- pal_nejm()(4)

###############################
### Load in reference files ###
###############################

# Gene annotations
bm <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/mouse_biomart_annotations_2-20-20.txt", sep="\t", header=T, stringsAsFactors=F)
remove_idx <- grep(bm$chromosome_name, pattern="MT|GL|JH|X|Y")
bm <- bm[-remove_idx, ]
bm$strand <- ifelse(bm$strand == 1, "+", "-")
gr.bm <- makeGRangesFromDataFrame(bm, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")

##################################
### Load processed data files ####
##################################

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

### DIFFERENTIAL FEATURES (DIFFBIND) ###

print("Load differential peaks")

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
	color.pal = pal)

p.atac <- plot_PCA(
	dat = adat, 
	title = "ATAC", 
	n_var_features = nrow(adat)/3, 
	annotation = metadat$Condition, 
	color.pal = pal)

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


