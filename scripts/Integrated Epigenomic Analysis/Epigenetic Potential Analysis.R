### Epigenetic potential analysis ###

rm(list=ls())

library(rtracklayer)
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(limma)
library(data.table)
library(clusterProfiler)
library(ggplot2)
library(ggsci)
library(viridis)
library(patchwork)

source("~/Dropbox/Minn/ifnar_epigenome/scripts/Final Scripts/Integrated Epigenomic Analysis/Integrated Epigenomic Analysis Functions.R")

txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

###############################
### Load in reference files ###
###############################

# Gene annotations
bm <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Resource Files/mouse_biomart_annotations_2-20-20.txt", sep="\t", header=T, stringsAsFactors=F)
remove_idx <- grep(bm$chromosome_name, pattern="MT|GL|JH|X|Y")
bm <- bm[-remove_idx, ]
bm$strand <- ifelse(bm$strand == 1, "+", "-")
gr.bm <- makeGRangesFromDataFrame(bm, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")

##############################
### Load in processed data ###
##############################

conditions <- c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO")

# RNA data
rna.dat_raw <- read.table(file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/RNA_counts_raw.txt"), header=T, sep="\t", stringsAsFactors=F)
colnames(rna.dat_raw) <- paste(colnames(rna.dat_raw), "RNA", sep="_")
colnames(rna.dat_raw) <- gsub("cas", "WT", colnames(rna.dat_raw))

# ATAC data
atac.dat_raw <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_original_consensus_mat_tn5_insertion_counts_IDR_raw.txt", sep="\t", header=T, stringsAsFactors=F)
colnames(atac.dat_raw) <- paste0(colnames(atac.dat_raw), "_ATAC")
gr.atac_IDR <- peakids2GRanges(rownames(atac.dat_raw), delim="_")

# H3K4me1 data 
H3K4me1.vst <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/H3K4me1_consensus_mat_insertion_counts_WT_pooled_vst.txt"), sep="\t", header=T, stringsAsFactors=F)
colnames(H3K4me1.vst) <- paste(rep(conditions, each=2), 1:2, "H3K4me1", sep="_")
gr.H3K4me1 <- peakids2GRanges(rownames(H3K4me1.vst), delim="_")

####################################
### Annotate ATAC promoter peaks ###
####################################

gr.atac_IDR$peak_ID <- 1:length(gr.atac_IDR)
peakAnno <- annotatePeak(gr.atac_IDR, tssRegion=c(-3000,3000), TxDb=txdb, annoDb="org.Mm.eg.db")
peakAnno.promoter <- peakAnno@anno[grep(peakAnno@anno$annotation, pattern="Promoter")]

# Genes with promoter peaks
gs <- intersect(rownames(rna.dat_raw), peakAnno.promoter$ENSEMBL)

# Select ATAC peaks
peakids <- peakAnno.promoter$peak_ID[na.omit(match(gs, peakAnno.promoter$ENSEMBL))]
names(peakids) <- gs
atac <- atac.dat_raw[as.numeric(peakids), ]
rownames(atac) <- names(peakids)

###################################
### Batch-correct RNA/ATAC data ###
###################################

# Combine RNA and ATAC count matrices
all(rownames(atac) %in% rownames(rna.dat_raw)) # Check
mat.combined <- cbind(rna.dat_raw[match(rownames(atac), rownames(rna.dat_raw)), ], atac[match(rownames(atac), rownames(atac)), ])

# Sample metadata
sample <- sapply(strsplit(colnames(mat.combined), split="_"), function(x) paste(x[1:3], collapse="_"))
cond <- sapply(strsplit(colnames(mat.combined), split="_"), function(x) paste(x[1:2], collapse="_"))
assay <- sapply(strsplit(colnames(mat.combined), split="_"), function(x) x[4])

ann <- data.frame(
	ID=colnames(mat.combined), 
	sample=sample, 
	condition=cond, 
	assay=factor(assay, levels=c("RNA", "ATAC")), row.names=colnames(mat.combined))
all(rownames(ann) == colnames(mat.combined))

# log2 transform and normalize count data with voom
design <- model.matrix(data = ann, ~ 0 + condition / assay)
v <- voom(counts=mat.combined, design=design, normalize.method="quantile")

# https://stat.ethz.ch/pipermail/bioconductor/2014-August/060895.html
# When calling removeBatchEffect, you should use the same design that you used for limma, but with with batch effect term removed from the design. Then you would pass the batch effect factor as the batch argument instead. So, if the design matrix that you used for limma was constructed as: model.matrix(~Condition + Batch), then for removeBatchEffect, you would use design=model.matrix(~Condition), and batch=Batch. In other words, you take the batch effect out of your model design and pass it as the batch argument instead.
v$E <- removeBatchEffect(
	v$E, 
	covariates=model.matrix(data=ann, ~ assay)[,2, drop=F], 
	design=model.matrix(data=ann, ~0 + condition))
# v$E <- removeBatchEffect(v$E, batch=ann$assay, design=model.matrix(data=ann, ~0 + condition))
fit <- lmFit(v, design)
fit <- eBayes(fit)

mat.batchCorrected <- v$E

###########################################################
### Identify genes with epigenetic potential (ATAC>RNA) ###
###########################################################

# Get coefficients
coefss <- grep(":assayATAC", colnames(coef(fit)), value=T)
res <- data.table()
for(coefx in coefss) {
	message(coefx)
	res <- rbind(res, data.table(topTable(fit, coef=coefx, number=nrow(mat.combined)), keep.rownames=T, condition=coefx))
}
res$mgi_symbol <- gr.bm$mgi_symbol[match(res$rn, gr.bm$ensembl_gene_id)]
res.sig <- res[AveExpr > 0 & logFC > 0.5 & adj.P.Val < 0.1]
table(res.sig$condition)

#################
### Visualize ###
#################

### Number of epigenetic potential genes ###

res.sig$condition <- factor(res.sig$condition)
levels(res.sig$condition) <- c("B16_SKO", "B16_WT", "R499_SKO", "R499_WT")

mat.count <- data.frame(
	condition = names(table(res.sig$condition)), 
	count = as.numeric(table(res.sig$condition)))
mat.count$condition <- factor(mat.count$condition, levels=conditions)
p1 <- ggplot(data=mat.count, aes(x=condition, y=count, fill=condition)) +
    geom_bar(stat="identity", width=0.7) +
    scale_fill_nejm() +
    labs(title="Number of genes with epigenetic potential", x="", y="Count", fill="Condition") +
    theme_bw(base_size=14) +
    theme(axis.text.x = element_text(size=14))
# # Write out plot data
# dat <- p1$data
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3C.left.csv"), sep=",", quote=F, col.names=T, row.names=F)

### GO analysis on genes with epigenetic potential ###

go_list <- lapply(1:length(conditions), function(x) {
	print(conditions[x])
	potential.genes <- res.sig[grepl(res.sig$condition, pattern=conditions[x]), ]$rn
	potential.genes_mgi <- res.sig[grepl(res.sig$condition, pattern=conditions[x]), ]$mgi_symbol
	go <- enrichGO(
		gene = unique(potential.genes), 
		universe = unique(gr.bm$ensembl_gene_id), 
		OrgDb = org.Mm.eg.db, 
		ont = "BP", 
		keyType = 'ENSEMBL', 
		pAdjustMethod = "BH", 
		pvalueCutoff  = 0.1, 
		qvalueCutoff  = 0.1, 
		readable = TRUE)
	return(go)
})

# Plot GO enrichment
sig <- 0.1
res.GO <- lapply(1:4, function(x) {
	print(x)
	res <- go_list[[x]]@result
	mat <- data.frame(Description=res$Description, p.adjust=res$p.adjust, qvalue=res$qvalue, count=res$Count, condition=conditions[x])
	return(mat)
})
res.GO <- do.call(rbind, res.GO)
res.GO_sig <- res.GO[res.GO$qvalue <= sig, ]
print(nrow(res.GO_sig))

res.GO_sorted <- res.GO[sort.int(res.GO$qvalue, decreasing=F, index.return=T)$ix, ]
prune <- c("regulation of interferon-alpha production", "interferon-alpha production", "regulation of neurotransmitter levels", "facial nerve morphogenesis", "type I interferon secretion", "interferon-alpha secretion", "positive regulation of interferon-alpha secretion", "regulation of type I interferon production", "type I interferon production")
res.GO_sorted <- res.GO_sorted[!res.GO_sorted$Description %in% prune, ]
plot.terms <- as.character(res.GO_sorted[1:10, ]$Description)

mat.go <- res.GO_sorted[res.GO_sorted$Description %in% plot.terms, ]
mat.go$Description <- factor(mat.go$Description, levels=rev(plot.terms))
mat.go$condition <- factor(mat.go$condition, levels=conditions)

p2 <- ggplot(mat.go, aes(x=condition, y=Description)) + 
	geom_point(mapping=aes(color=qvalue, size=-log10(qvalue))) +
	scale_color_gradientn(colours=viridis(16), trans = 'reverse') +
	labs(title="Epi potential genes", x ="", y = "", color="qvalue") +
	guides(size=FALSE) +
	# coord_fixed() +
	scale_size(range = c(2,8)) +
	theme_bw(base_size=14) +
	theme(axis.text.x = element_text(size=14, color="black", angle=30, hjust=1),
		axis.text.y=element_text(size=14, color="black"))
# # Write out plot data
# dat <- p2$data[ ,c("condition", "Description", "qvalue")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3C.right.csv"), sep=",", quote=F, col.names=T, row.names=F)

fig <- p1 + p2 +
	plot_layout(widths=c(0.8, 1))
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Fig 3C Epigenetic Potential Genes.pdf"), width=11, height=4.5)

##################################################
### H3K4me1 signal at genes with epi potential ###
##################################################

gr.atac_epiPotential <- peakAnno.promoter[match(unique(res.sig$rn[grepl(res.sig$condition, pattern="R499_SKO")]), peakAnno.promoter$ENSEMBL)]

ol <- findOverlaps(resize(gr.H3K4me1, fix = "center", width = 5000), gr.atac_epiPotential)
mat.plot <- H3K4me1.vst[unique(from(ol)), ]
mat.plot <- t(scale(t(mat.plot)))
df <- data.frame(
	Sample = colnames(mat.plot),
	Signal = colMeans(mat.plot), stringsAsFactors=F)
df$Condition <- sapply(strsplit(df$Sample, split="_"), function(x) paste(x[1:2], collapse="_"))
df$Label <- factor(gsub("_", " ", df$Condition),
	levels = c("B16 WT", "B16 SKO", "R499 WT", "R499 SKO"))

p <- ggplot(df, aes(x = Label, y = Signal, color = Label)) +
	geom_point(shape = 1, size = 3, stroke = 1.5) +
	stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
	scale_color_nejm() +
	labs(x = "", y = "H3K4me1 Signal", color = "Condition") +
	theme_bw(base_size = 14) +
	theme(aspect.ratio = 1,
		panel.border = element_rect(color = "black", fill=NA, size=1),
		axis.text.x = element_text(color = "black", angle = 45, hjust = 1, size = 12),
		axis.text.y = element_text(color = "black", size = 12))
ggsave(p, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 3/Figure 3D H3K4me1 Epigenetic Potential Genes H3K4me1.pdf"), width=5, height=4.5)

# Write out plot data
dat <- p$data[ ,c("Signal", "Label")]
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3D.left.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Significance test (one-way ANOVA, Tukey HSD pairwise comparisons)
df <- read.csv(paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 3D.left.csv"), sep=",", header=T)
res.anova <- aov(Signal ~ Label, data = df)
print(assay)
print(summary(res.anova))
print(TukeyHSD(res.anova)) # Pairwise comparisons

# H3K4me1 at epi potential genes (Fig 3D)

# H3K4me1 (pval = 3.1e-03)
# B16 WT vs R499 WT padj = 8.1e-03
# B16 WT vs B16 SKO padj = 8.0e-01
# R499 WT vs R499 SKO padj = 6.9e-01
# B16 SKO vs R499 SKO padj = 9.1e-03
