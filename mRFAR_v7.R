# mRFAR with Tn5 insertion counts, normalized by DESeq2
# Run with mRFAR_v7_functions.R
# Both RF and glmnet models
# window size = 150000

# Run on peaks of interest (Oas1, ISG, IFNG, ICB.RS, IFN.I)
# Run on all genes
# Run on DE genes (B16cas_v_R499cas, B16cas_v_R499cas_UP, B16cas_v_B16sko, R499cas_v_R499sko)

# Run with R version 3.5.1
# BSgenome.Mmusculus.UCSC.mm10_1.4.0
# DESeq2_1.22.1
# rtracklayer_1.42.1
# caret_6.0-81
# GenomicRanges_1.34.0
# ggplot2_3.1.0
# randomForestSRC_2.9.3

rm(list=ls())

library(GenomicRanges)
library(DESeq2)
library(caret)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(randomForestSRC)

dir.path <- "~/Dropbox/Minn/"

source(paste(dir.path, "ATAC_RNA_integration/scripts/mRFAR_v7/v1/final_scripts/mRFAR_v7_functions.R", sep=""))

### Mouse RNA annotation data ###

# # Extract TSS
# # For genes with multiple TSS annotations, take TSS that overlaps with start position
# get_tss <- function(geneID, bm) {
# 	df <- bm[which(bm$ensembl_gene_id == geneID), ]
# 	if(nrow(df) == 1) {
# 		tss <- df$transcription_start_site
# 	} else {
# 		if(df$strand[1] == 1) {
# 			tss <- df$transcription_start_site[which.min(abs(df[,7] - df[,4]))]
# 		} else {
# 			tss <- df$transcription_start_site[which.min(abs(df[,7] - df[,5]))]
# 		}
# 	}
# 	return(tss)
# }

# # Get biomart annotations for mouse genes
# ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
# bm <- getBM(attributes=c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position", "strand", "transcription_start_site"), mart=ensembl)
# remove_idx <- which(sapply(strsplit(bm$chromosome_name, split="_"), function(x) length(x)) > 1)
# bm <- bm[-remove_idx, ]
# bm$chromosome_name <- paste0("chr", bm$chromosome_name)
# bm$mgi_symbol[which(bm$mgi_symbol == "")] <- NA
# bm$tss <- sapply(bm$ensembl_gene_id, function(x) get_tss(geneID=x, bm=bm))
# write.table(bm, paste0(dir.path, "resources/mouse_biomart_annotations_2-20-20.txt"), quote=F, col.names=T, row.names=F, sep="\t")

# Get gene annotation info
bm <- read.table(paste0(dir.path, "resources/mouse_biomart_annotations_2-20-20.txt"), sep="\t", header=T, stringsAsFactors=F)
gr.bm <- makeGRangesFromDataFrame(bm, seqnames.field="chromosome_name", start.field="start_position", end.field="end_position", ignore.strand=T)
gr.bm$gene <- bm$mgi_symbol

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

##################################
### Match ATAC and RNA samples ###
##################################

edat <- rna.dat
adat <- atac.dat_idr[ ,intersect(colnames(rna.dat), colnames(atac.dat_idr))]

# ###########################################
# ### Select genes of interest (DE genes) ###
# ###########################################

# res_B16cas_v_R499cas <- read.table(paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/DESeq2/res_B16cas_v_R499cas.txt", sep=""), sep="\t", header=T)
# res_B16cas_v_B16sko <- read.table(paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/DESeq2/res_B16cas_v_B16sko.txt", sep=""), sep="\t", header=T)
# res_R499cas_v_R499sko <- read.table(paste(dir.path, "ATAC_RNA_integration/mRFAR_v5/DESeq2/res_R499cas_v_R499sko.txt", sep=""), sep="\t", header=T)

# res_B16cas_v_R499cas_UP <- res_B16cas_v_R499cas[which(res_B16cas_v_R499cas$log2FoldChange > 0), ]
# res_B16cas_v_B16sko_DOWN <- res_B16cas_v_B16sko[which(res_B16cas_v_B16sko$log2FoldChange < 0), ]
# res_R499cas_v_R499sko_DOWN <- res_R499cas_v_R499sko[which(res_R499cas_v_R499sko$log2FoldChange < 0), ]

# # Lists of DE genes (Ensembl IDs)
# gs.res_0.05 <- rownames(res_B16cas_v_R499cas)[which(res_B16cas_v_R499cas$padj <= 0.05)]
# gs.res_0.1 <- rownames(res_B16cas_v_R499cas)[which(res_B16cas_v_R499cas$padj <= 0.1)]
# gs.res_UP_0.05 <- rownames(res_B16cas_v_R499cas_UP)[which(res_B16cas_v_R499cas_UP$padj <= 0.05)]
# gs.res_UP_0.1 <- rownames(res_B16cas_v_R499cas_UP)[which(res_B16cas_v_R499cas_UP$padj <= 0.1)]

# gs.b16cas_v_b16sko <- rownames(res_B16cas_v_B16sko)[which(res_B16cas_v_B16sko$padj <= 0.1)]
# gs.b16cas_v_b16sko_DOWN <- rownames(res_B16cas_v_B16sko_DOWN)[which(res_B16cas_v_B16sko_DOWN$padj <= 0.1)]

# gs.r499cas_v_r499sko <- rownames(res_R499cas_v_R499sko)[which(res_R499cas_v_R499sko$padj <= 0.1)]
# gs.r499cas_v_r499sko_DOWN <- rownames(res_R499cas_v_R499sko_DOWN)[which(res_R499cas_v_R499sko_DOWN$padj < 0.1)]

# gs.DE_list <- list(gs.res_0.05, gs.res_UP_0.05, gs.res_0.1, gs.res_UP_0.1, gs.b16cas_v_b16sko, gs.b16cas_v_b16sko_DOWN, gs.r499cas_v_r499sko, gs.r499cas_v_r499sko_DOWN)

#####################################
### Load in gene sets of interest ###
#####################################

# All genes
gs.all <- rownames(edat)

# Make random gene sets #

# # v1 (made representative heatmaps from this)
# set.seed(1)
# gs.rand1 <- rownames(edat)[sample(1:nrow(edat), 200)]
# gs.rand2 <- rownames(edat)[sample(1:nrow(edat), 200)]
# gs.rand3 <- rownames(edat)[sample(1:nrow(edat), 200)]
# gs.rand <- list(gs.rand1, gs.rand2, gs.rand3)
# lapply(1:3, function(x) write.table(gs.rand[[x]], file=paste0(dir.path, "resources/random_gene_set_", x, "_mouse.txt"), sep="\t", row.names=F, col.names=F, quote=F))

# v2 - iterate 100x
n <- 100
set.seed(1)
lapply(1:n, function(x) {
	print(x)
	gs.rand <- rownames(edat)[sample(1:nrow(edat), 200)]
	write.table(gs.rand, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/random_gene_sets/random_gene_set_", x, "_mouse.txt"), row.names=F, col.names=F, quote=F)
})

n <- 100
gs.rand <- lapply(1:n, function(x) {
	# gs <- read.table(file=paste0(dir.path, "resources/random_gene_set_", x, "_mouse.txt"), sep="\t", header=F, stringsAsFactors=F)
	gs <- read.table(file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/random_gene_sets/random_gene_set_", x, "_mouse.txt"), sep="\t", header=F, stringsAsFactors=F)
	gs <- as.character(gs$V1)
	return(gs)
})

# Oas1 family genes
gs.Oas1 <- as.character(rna.anno$ensembl_gene_id[grep(rna.anno$mgi_symbol, pattern="Oas1")])

# ISGs
gs.ISG <- read.table(paste(dir.path, "resources/ISG_mouse.txt", sep=""), sep="\t", header=F)
gs.ISG <- as.character(gs.ISG$V1)
gs.ISG <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.ISG, rna.anno$mgi_symbol))])

gs.IFNG <- read.table(paste(dir.path, "resources/IFNG_mouse.txt", sep=""), sep="\t", header=F)
gs.IFNG <- as.character(gs.IFNG$V1)
gs.IFNG <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.IFNG, rna.anno$mgi_symbol))])

gs.IFN.I <- read.table(paste(dir.path, "resources/IFN.I_mouse.txt", sep=""), sep="\t", header=F)
gs.IFN.I <- unique(as.character(gs.IFN.I$V1))
gs.IFN.I <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.IFN.I, rna.anno$mgi_symbol))])

gs.Twyman <- read.table(paste(dir.path, "resources/Twyman_up_mouse.txt", sep=""), sep="\t", header=F)
gs.Twyman <- unique(as.character(gs.Twyman$V1))
gs.Twyman <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.Twyman, rna.anno$mgi_symbol))])

# Negative control gene sets
gs.EMT <- read.table(paste(dir.path, "resources/EMT_mouse.txt", sep=""), sep="\t", header=F)
gs.EMT <- unique(as.character(gs.EMT$V1))
gs.EMT <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.EMT, rna.anno$mgi_symbol))])

gs.TGFB <- read.table(paste(dir.path, "resources/TGFB_mouse.txt", sep=""), sep="\t", header=F)
gs.TGFB <- unique(as.character(gs.TGFB$V1))
gs.TGFB <- as.character(rna.anno$ensembl_gene_id[na.omit(match(gs.TGFB, rna.anno$mgi_symbol))])

# TCIR ligands
gs.tcir <- as.character(rna.anno$ensembl_gene_id[match(c("Cd274", "Tnfrsf14", "Lgals9"), as.character(rna.anno$mgi_symbol))])

#################
### Run mRFAR ###
#################

# v1 - univariate RF implemented by caret
# v2 - multivariate RF implemented by randomForestSRC

version <- "v1"

# Window size from Rao et al.
W <- 92500

# names <- c("Oas1", "ISG", "IFNG", "IFN.I", "Twyman", "EMT", "TGFB")
# gs_list <- list(gs.Oas1, gs.ISG, gs.IFNG, gs.IFN.I, gs.Twyman, gs.EMT, gs.TGFB)

names <- paste0("rand", 1:n)
gs_list <- gs.rand

# x <- 1
# gs <- gs_list[[x]]
# gr_rna <- gr.rna
# rna_anno <- rna.anno
# rna_dat <- edat
# gr_atac <- gr.atac_IDR
# atac_dat <- adat
# scaled <- F
# seed <- T
# version <- "v2"

mRFARobj_list <- lapply(1:length(names), function(x) {
	mRFARobj <- mRFAR(gs=gs_list[[x]], W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_IDR, atac_dat=adat, scaled=F, seed=T, version=version)
	saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/mRFARobj_", names[x], "_unscaledVIMPs.rds", sep=""))
	# saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/random/mRFARobj_", names[x], "_unscaledVIMPs.rds", sep=""))
})

# Scale importance scores
mRFARobj_list <- lapply(1:length(names), function(x) {
	mRFARobj <- mRFAR(gs=gs_list[[x]], W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_IDR, atac_dat=adat, scaled=T, seed=T, version=version)
	saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/mRFARobj_", names[x], "_scaledVIMPs.rds", sep=""))
})

# Run mRFAR on all genes
gs <- gs.all

mRFARobj_ALL <- mRFAR(gs=gs, W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_IDR, atac_dat=adat, scaled=F, seed=T, version=version)
saveRDS(mRFARobj_ALL, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/mRFARobj_ALL_genes_unscaledVIMPs.rds", sep=""))

mRFARobj_ALL <- mRFAR(gs=gs, W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_IDR, atac_dat=adat, scaled=T, seed=T, version=version)
saveRDS(mRFARobj_ALL, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/mRFARobj_ALL_genes_scaledVIMPs.rds", sep=""))

##############################
### Extract mRFARobj peaks ###
##############################

names <- c("Oas1", "ISG", "IFNG", "IFN.I", "Twyman", "EMT", "TGFB")
names <- c("ALL_genes")

# Load in mRFAR objects for gene sets of interest
mRFARobj_list <- lapply(1:length(names), function(x) readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/mRFARobj_", names[x], "_unscaledVIMPs.rds", sep="")))

# Extract peaks for genes sets of interest
sapply(1:length(names), function(x) extract_mRFARobj_peaks(mRFARobj=mRFARobj_list[[x]], gr.atac=gr.atac_IDR, dir_path=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/peaks_of_interest/", version, "/", names[x], sep=""), peak_type=names[x], trimN=500, rna.anno=rna.anno, version=version))


### Random ###

names <- paste0("rand", 1:n)

mRFARobj_list <- lapply(1:length(names), function(x) readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/random/mRFARobj_", names[x], "_unscaledVIMPs.rds", sep="")))

# Extract peaks for genes sets of interest
sapply(1:length(names), function(x) extract_mRFARobj_peaks(mRFARobj=mRFARobj_list[[x]], gr.atac=gr.atac_IDR, dir_path=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/peaks_of_interest/", version, "/random/", names[x], sep=""), peak_type=names[x], trimN=500, rna.anno=rna.anno, version=version))

########################################
### Variability of importance scores ###
########################################

W <- 92500

mRFARobj_repeats <- lapply(1:10, function(x) mRFAR(gs=gs.Oas1, W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_IDR, atac_dat=adat, scaled=F, seed=F))

Rsq_rf <- sapply(mRFARobj_repeats, function(x) x[[1]]$Rsq_rf)
Rsq_glmnet <- sapply(mRFARobj_repeats, function(x) x[[1]]$Rsq_glmnet)

error_rf <- sapply(mRFARobj_repeats, function(x) x[[1]]$error_rf)
error_glmnet <- sapply(mRFARobj_repeats, function(x) x[[1]]$error_glmnet)

vimp_scores <- lapply(mRFARobj_repeats, function(x) x[[1]]$vimps_rf$Overall)
vimps <- do.call(cbind, vimp_scores)
rownames(vimps) <- paste("p", 1:length(mRFARobj_repeats[[1]][[1]]$atacIndices), sep="")
# rownames(vimps) <- 1:length(mRFARobj_repeats[[1]][[1]]$atacIndices)
colnames(vimps) <- paste("rep_", 1:10, sep="")

vimps_m <- melt(vimps)

fig <- ggplot(vimps_m, aes(x=Var1, y=value)) +
	geom_boxplot(outlier.shape=NA) +
	geom_point() +
	xlab("Peak ID") +
	ylab("VIMP score") +
	ggtitle("VIMP score variance (Oas1 cis-regulatory peaks)") +
	theme(axis.text.x=element_text(size=12))
ggsave(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/mRFAR/VIMP_variance_boxplot_10repeats_Oas1_peaks.pdf", sep=""), device="pdf", width=8, height=6, plot=fig)

# Plot average VIMP score vs variance
pdf(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/reanalysis/mRFAR/VIMP_vs_variance_scatter_10repeats_Oas1_peaks.pdf", sep=""), width=6, height=6)
par(mar=c(5,5,2.5,2.5))
df <- data.frame(mean=as.numeric(apply(vimps, 1, mean)), var=as.numeric(apply(vimps, 1, var)))
fit <- lm(var~mean, data=df)
plot(df$mean, df$var, pch=18, cex=1, xlab="VIMP (mean)", ylab="Variance", cex.lab=1.5)
abline(fit, col="red", lty="dotted", lwd=2.5)
text(0.4, 0.25, labels="F=9.3, p=8.6e-3", cex=1.1)
dev.off()

#################################
#################################

mRFARobj_ALL <- readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_ALL_genes_unscaledVIMPs.rds", sep=""))
geneids <- sapply(mRFARobj_ALL, function(x) x$geneid)
genenames <- sapply(mRFARobj_ALL, function(x) x$geneSymb)

# Get mRFAR objects for DE gene sets
names_DE <- c("B16cas_v_R499cas_0.05", "B16cas_v_R499cas_UP_0.05", "B16cas_v_R499cas_0.1", "B16cas_v_R499cas_UP_0.1", "B16cas_v_B16sko_0.1", "B16cas_v_B16sko_DOWN_0.1", "R499cas_v_R499sko_0.1", "R499cas_v_R499sko_DOWN_0.1")
mRFARobj_DE_list <- lapply(gs.DE_list, function(x) mRFARobj_ALL[match(x, geneids)])
names(mRFARobj_DE_list) <- names_DE

# Get mRFAR objects for peaks of interest
mRFARobj_list <- lapply(1:length(names), function(x) readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_", names[x], "_unscaledVIMPs.rds", sep="")))
mRFARobj_list <- lapply(1:length(names), function(x) readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_", names[x], ".rds", sep="")))

########################

mRFARobj_ALL <- readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_ALL_genes_unscaledVIMPs.rds", sep=""))
geneids <- sapply(mRFARobj_ALL, function(x) x$geneid)
genenames <- sapply(mRFARobj_ALL, function(x) x$geneSymb)
mRFARobj <- mRFARobj_ALL[match(gs.Twyman, geneids)]
saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_Twyman_unscaledVIMPs.rds", sep=""))
process_toFasta(mRFARobj, gr.atac_IDR, dir_name="Twyman", file_name="Twyman", trimN=NA)

## Explore ##

# Variance explained of all genes
Rsq_rf <- sapply(mRFARobj_ALL, function(x) x$Rsq_rf)
Rsq_glmnet <- sapply(mRFARobj_ALL, function(x) x$Rsq_glmnet)
quantile(Rsq_rf, na.rm=T)
quantile(Rsq_glmnet, na.rm=T)
hist(Rsq_rf, breaks=50, col="yellow")
hist(Rsq_glmnet, breaks=50, col="yellow")

# DE genes between B16 and R499 explained significantly more by chromatin accessibility than background
# Glmnet seems to have better variance explained models
hist(sapply(mRFARobj_DE_list[[1]], function(x) x$Rsq_rf), breaks=50, col="yellow")
hist(sapply(mRFARobj_DE_list[[1]], function(x) x$Rsq_glmnet), breaks=50, col="yellow")
quantile(sapply(mRFARobj_DE_list[[1]], function(x) x$Rsq_rf), na.rm=T)
quantile(sapply(mRFARobj_DE_list[[1]], function(x) x$Rsq_glmnet), na.rm=T)

# Oas1, ICB.RS well explained by accessibility. ISGs not as much and IFNG even less.
quantile(sapply(mRFARobj_list[[1]], function(x) x$Rsq_rf), na.rm=T)
quantile(sapply(mRFARobj_list[[1]], function(x) x$Rsq_glmnet), na.rm=T)
quantile(sapply(mRFARobj_list[[2]], function(x) x$Rsq_rf), na.rm=T)
quantile(sapply(mRFARobj_list[[2]], function(x) x$Rsq_glmnet), na.rm=T)
quantile(sapply(mRFARobj_list[[3]], function(x) x$Rsq_rf), na.rm=T)
quantile(sapply(mRFARobj_list[[3]], function(x) x$Rsq_glmnet), na.rm=T)
quantile(sapply(mRFARobj_list[[4]], function(x) x$Rsq_rf), na.rm=T)
quantile(sapply(mRFARobj_list[[4]], function(x) x$Rsq_glmnet), na.rm=T)

# Errors
Error_rf <- sapply(mRFARobj_ALL, function(x) x$error_rf)
Error_glmnet <- sapply(mRFARobj_ALL, function(x) x$error_glmnet)
hist(Error_rf, breaks=50, col="yellow")
hist(Error_glmnet, breaks=50, col="yellow")
quantile(Error_rf, na.rm=T)
quantile(Error_glmnet, na.rm=T)

# Genes that are well explained - not great overlap between two methods! Complementary? Or something else..
genes_explained_rf <- genenames[which(Rsq_rf > 0.75)]
genes_explained_glmnet <- genenames[which(Rsq_glmnet > 0.75)]
c(genes_explained_rf, genes_explained_glmnet)
intersect(genes_explained_rf, genes_explained_glmnet)

genes_explained_rf <- genenames[which(Rsq_rf > 0.6)]
genes_explained_glmnet <- genenames[which(Rsq_glmnet > 0.6)]
# c(genes_explained_rf, genes_explained_glmnet)
intersect(genes_explained_rf, genes_explained_glmnet)

genes_explained_rf <- genenames[which(Rsq_rf > 0.5)]
genes_explained_glmnet <- genenames[which(Rsq_glmnet > 0.5)]
# c(genes_explained_rf, genes_explained_glmnet)
intersect(genes_explained_rf, genes_explained_glmnet)

# Sort genes that are well explained by RF and glmnet by Rsq
genes_explained_int_Rsq_glmnet <- Rsq_glmnet[match(intersect(genes_explained_rf, genes_explained_glmnet), genenames)]
genes_explained_int_Rsq_rf <- Rsq_rf[match(intersect(genes_explained_rf, genes_explained_glmnet), genenames)]
genes_explained_int_Rsq_glmnet_sorted <- sort.int(genes_explained_int_Rsq_glmnet, decreasing=T, index.return=T)$ix
genes_explained_int_Rsq_rf_sorted <- sort.int(genes_explained_int_Rsq_rf, decreasing=T, index.return=T)$ix
intersect(genes_explained_rf, genes_explained_glmnet)[genes_explained_int_Rsq_glmnet_sorted]
intersect(genes_explained_rf, genes_explained_glmnet)[genes_explained_int_Rsq_rf_sorted]

#############
mRFARobj_ALL <- readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_ALL_genes_unscaledVIMPs.rds", sep=""))
process_toFasta(mRFARobj_ALL, gr.atac_IDR, dir_name="ALL", file_name="ALL", trimN=NA)

############

mRFARobj <- readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_TCIR_unscaledVIMPs.rds", sep=""))
process_toFasta(mRFARobj, gr.atac_IDR, dir_name="TCIR", file_name="TCIR", trimN=NA)

#########################################
### Make Oas1 mRFARobj with all peaks ###
#########################################

gr.atac_Oas1 <- gr.atac_IDR[49231:49238]
adat_Oas1 <- adat[49231:49238,]

W <- 150000
mRFARobj <- mRFAR(gs=gs.Oas1, W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_Oas1, atac_dat=adat_Oas1)

saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_Oas1_ALL_peaks.rds", sep=""))
process_toFasta(mRFARobj, gr.atac_Oas1, dir_name="Oas1_ALL_peaks", file_name="Oas1_ALL_peaks", trimN=NA)

########################################
### Make Oas mRFARobj with all peaks ###
########################################

gs.Oas <- grep(rna.anno$mgi_symbol, pattern="Oas1|Oas2|Oas3", value=T)
gs.Oas <- as.character(rna.anno$ensembl_gene_id[match(gs.Oas, rna.anno$mgi_symbol)])

gr.atac_Oas <- gr.atac_IDR[49227:49238]
adat_Oas <- adat[49227:49238,]

W <- 500000
mRFARobj <- mRFAR(gs=gs.Oas, W=W, gr_rna=gr.rna, rna_anno=rna.anno, rna_dat=edat, gr_atac=gr.atac_Oas, atac_dat=adat_Oas)

saveRDS(mRFARobj, file=paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_Oas_genes_ALL_peaks.rds", sep=""))
process_toFasta(mRFARobj, gr.atac_Oas, dir_name="Oas_genes_ALL_peaks", file_name="Oas1_genes_ALL_peaks", trimN=NA)

########################################
### Rank genes by variance explained ###
########################################

mRFARobj_ALL <- readRDS(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/mRFARobj_ALL_genes_unscaledVIMPs.rds", sep=""))

genes <- sapply(mRFARobj_ALL, function(x) x$geneSymb)
Rsq_rf <- sapply(mRFARobj_ALL, function(x) x$Rsq_rf)

genes_filt <- genes[!is.na(Rsq_rf)]
Rsq_rf_filt <- Rsq_rf[!is.na(Rsq_rf)]

Rsq_rf.sorted <- Rsq_rf_filt[sort.int(Rsq_rf_filt, decreasing=TRUE, index.return=T)$ix]
genes.sorted <- genes_filt[sort.int(Rsq_rf_filt, decreasing=TRUE, index.return=T)$ix]

pdf(paste(dir.path, "ATAC_RNA_integration/mRFAR_v7/figures/v4/Rsq_all_genes_ranked.pdf", sep=""), width=6, height=5)
plot(1:length(Rsq_rf.sorted), Rsq_rf.sorted, pch=18, cex=0.4, xlab="Rank", ylab="Variance explained (R squared)")
points(grep(genes.sorted, pattern="Oas1a|Oas1g"), Rsq_rf.sorted[grep(genes.sorted, pattern="Oas1a|Oas1g")], pch=16, cex=2.3, col="red")
dev.off()



