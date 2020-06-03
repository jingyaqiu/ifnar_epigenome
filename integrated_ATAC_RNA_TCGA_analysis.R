rm(list=ls())

library(TCGAutils)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(umap)
library(BuenColors)
library(betareg)
library(ggpubr)

dir.path <- "~/Dropbox/Minn/"

cancer_types <- c("SKCM", "ACCx", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBMx", "HNSC", "KIRC", "KIRP", "LGGx", "LIHC", "LUAD", "LUSC", "MESO", "PCPG", "PRAD", "STAD", "TGCT", "THCA", "UCEC")

isg.types <- c("ISG.RS", "IFNG.GS")

# Get gene annotation info
bm <- read.table(paste0(dir.path, "resources/human_biomart_annotations_2-20-20.txt"), header=T, stringsAsFactors=F)

# Load in ISG.RS and IFNG.GS genes
gs.ISG.RS_hgnc <- read.delim(paste(dir.path, "resources/Andy_Dropbox/ISG.RS genes.txt", sep=""), header=F, stringsAsFactors=F)
gs.ISG.RS_hgnc <- gs.ISG.RS_hgnc$V1
gs.ISG.RS <- na.omit(bm$ensembl_gene_id[match(c(gs.ISG.RS_hgnc), bm$hgnc_symbol)])

gs.IFNG.GS_hgnc  <- read.delim(paste(dir.path, "resources/Andy_Dropbox/IFNG.GS genes.txt", sep=""), header=F, stringsAsFactors=F)
gs.IFNG.GS_hgnc  <- gs.IFNG.GS_hgnc $V1
gs.IFNG.GS <- na.omit(bm$ensembl_gene_id[match(c(gs.IFNG.GS_hgnc ), bm$hgnc_symbol)])

plot <- FALSE

#################
### Functions ###
#################

### Merge technical replicates ###

merge_technical_replicates <- function(dat) {
	unique_samples <- unique(colnames(dat))
	dat_merged <- lapply(1:length(unique_samples), function(x) {
		dat_sub <- as.matrix(dat[ ,which(colnames(dat) == unique_samples[x])])
		if(ncol(dat_sub) > 1) {
			return(as.numeric(rowMeans(dat_sub)))
		} else {
			return(dat_sub[,1])
		}
	})
	dat_merged <- do.call(cbind, dat_merged)
	colnames(dat_merged) <- unique_samples
	rownames(dat_merged) <- rownames(dat)
	return(dat_merged)
}

### Get important ISG ATAC peaks ###

get_impt_ISG_peaks <- function(var_exp_cutoff, vimp_cutoff, cancer_type, isg.type, version, method, scaled=F) {

	scale <- ifelse(scaled == TRUE, "scaled", "unscaled")
	
	### Load in mRFAR objects ###
	# mRFARobj <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/mRFARobj/mRFARobj_", isg.type, "_unscaledVIMPs_", cancer_type, "_RF.rds"))
	mRFARobj <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/mRFARobj/", version, "/Corces/ISG/mRFARobj_", isg.type, "_", scale, "VIMPs_", cancer_type, "_", method, ".rds"))

	### Get important peaks for each gene ###

	# Remove genes not explained well by ATAC
	remove_idx <- which(sapply(mRFARobj, function(x) x$Rsq_rf) < var_exp_cutoff)
	if(length(remove_idx) > 0) {
		mRFARobj_filt <- mRFARobj[-remove_idx]
	} else {
		mRFARobj_filt <- mRFARobj
	}

	if(method == "rfSRC") {
		impt_peaks <- lapply(mRFARobj_filt, function(x) names(x$vimps_rf)[which(as.numeric(x$vimps_rf) > vimp_cutoff)])
	} else {
		# impt_peaks <- lapply(mRFARobj_filt, function(x) rownames(x$vimps_glmnet)[which(x$vimps_glmnet[,1] > vimp_cutoff)])
		impt_peaks <- lapply(mRFARobj_filt, function(x) rownames(x$vimps_rf)[which(as.numeric(x$vimps_rf[,1]) > vimp_cutoff)])
	}
	impt_peaks <- unique(do.call(c, impt_peaks))
	return(impt_peaks)
}

### Extract estimate, 95% confidence intervals for beta regression model ###
get_betareg_stats <- function(fit, name) {
	summary <- confint(fit)
	coefs <- summary(fit)$coefficients$mean

	est <- as.numeric(coef(fit)[which(names(coef(fit)) == name)])
	se <- coefs[rownames(coefs) == name, 2]
	lower <- summary[which(rownames(summary) == name), 1]
	upper <- summary[which(rownames(summary) == name), 2]
	return(c(est, se, lower, upper))
}

### Make RNA metagenes ###

make_ISG_metagenes <- function(cancer_type, gs1, gs2, rna.dat, cibersort.dat, rna.meta) {
	caseUUIDs <- as.character(rna.meta$case_UUID)[which(as.character(rna.meta$cancer_type) == cancer_type)]

    edat <- rna.dat[ ,match(caseUUIDs, colnames(rna.dat))]
    ciber.dat <- cibersort.dat[ ,match(caseUUIDs, colnames(cibersort.dat))]

  	# RNA data
  	edat_gs1 <- edat[na.omit(match(gs1, rownames(edat))), ]
  	edat_gs2 <- edat[na.omit(match(gs2, rownames(edat))), ]

  	### Make RNA metagenes ###

  	gs1_meta_RNA <- colMeans(t(scale(t(edat_gs1))))
  	gs2_meta_RNA  <- colMeans(t(scale(t(edat_gs2))))

  	### Merge ISG.RS/IFNG.GS accessibility and immune cell infiltration ###

  	df <- data.frame(gs1_RNA=gs1_meta_RNA, gs2_RNA=gs2_meta_RNA, T.cells.CD8=as.numeric(ciber.dat[which(rownames(ciber.dat) == "T.cells.CD8"), ]), Macrophages.M0=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M0"), ]), Macrophages.M1=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M1"), ]), Macrophages.M2=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M2"), ]))
  	return(df)
}

### Dimensionality reduction on ISG features (RNA or ATAC) ###

plot_ISG_umap <- function(dat, labels) {
	umap <- umap(dat, n_components=10, random_state=42)
	df.umap <- data.frame(umap$layout)
	colnames(df.umap) <- paste0("PC", 1:10)
	df.umap$label <- labels

	if(length(unique(labels)) > 2) {
		pal <- colorRampPalette(brewer.pal(11, "Spectral"))
		cols <- pal(length(unique(labels)))
	} else {
		cols <- c(brewer.pal(11, "Spectral")[10], brewer.pal(11, "Spectral")[2])
	}
	fig <- ggplot(shuf(df.umap), aes(x=PC1, y=PC2, col=label)) +
		geom_point(size=3.5, alpha=1, shape=1, stroke=2) +
		scale_color_manual(values=cols) +
		geom_vline(xintercept=0, linetype="dashed") +
		geom_hline(yintercept=0, linetype="dashed") +
		xlab("Component 1") +
		ylab("Component 2") +
		theme_classic() +
		theme(axis.line=element_line(colour="black"), 
			panel.border=element_rect(colour="black", fill=NA, size=0.8), 
			legend.title=element_text(size=15),
			legend.text=element_text(size=15),
			text = element_text(size=22),
			aspect.ratio=1,
			plot.title = element_text("Helvetica"),
			plot.margin=unit(c(0.25,1,0.5,0.25), "cm"))
	return(fig)
}

####################################
### Make processed data matrices ###
####################################

##########################################
### Make consensus RNA-seq data matrix ###
##########################################

### 8738 samples total, 23162 consensus genes ###

# Load in manifests downloaded from GDC Data Portal 
# https://portal.gdc.cancer.gov/repository
manifest_files <- list.files(paste(dir.path, "ATAC_RNA_integration/resources/TCGA/", sep=""), full.names=T, pattern="gdc_manifest.")

# Load in RNA-seq data for each cancer type
rna.dat_list <- lapply(1:length(cancer_types), function(x) {
	cancer_type <- cancer_types[x]
	print(cancer_type)

	# Load in normalized RNA-seq data
	rna.dat <- read.table(file=paste(dir.path, "ATAC_RNA_integration/resources/TCGA/TCGA_", cancer_type, "/TCGA_", cancer_type, "_RNAseq_normalized_counts_DESeq2_vst.txt", sep=""), sep="\t", header=T, check.names=F)

	# Load in manifest
    manifest <- read.table(manifest_files[grep(manifest_files, pattern=cancer_type)], sep="\t", header=T, stringsAsFactors=F)

    # Convert manifest file_UUID to case_UUID
    df_fileUUIDtoCaseUUID <- UUIDtoUUID(manifest$id)

    # Convert file UUIDs (colnames) to case UUIDs with manifest
    rna.dat_fileUUID <- manifest$id[match(colnames(rna.dat), sapply(strsplit(manifest$filename, split="\\."), function(x) x[[1]]))]
    rna.dat_caseUUID <- df_fileUUIDtoCaseUUID$cases.case_id[match(rna.dat_fileUUID, df_fileUUIDtoCaseUUID$file_id)]
    
    # RNA metadata
    rna.meta <- data.frame(cancer_type=cancer_type, rep=1:ncol(rna.dat), RNA_ID=colnames(rna.dat), file_UUID=rna.dat_fileUUID, case_UUID=rna.dat_caseUUID)

    colnames(rna.dat) <- rna.dat_caseUUID
    rna.dat <- merge_technical_replicates(rna.dat)

	return(list(rna.dat, rna.meta))
})

# Consensus genes
consensus_genes <- Reduce(intersect, lapply(rna.dat_list, function(x) rownames(x[[1]])))

# Subset only consensus genes
rna.dat_consensus <- lapply(1:length(cancer_types), function(x) {
	rna.dat <- rna.dat_list[[x]][[1]][match(consensus_genes, rownames(rna.dat_list[[x]][[1]])), ]
	return(rna.dat)
})
rna.dat_consensus <- do.call(cbind, rna.dat_consensus)

# Get sample information
rna.dat_meta <- lapply(1:length(cancer_types), function(x) {
	cancer_type <- cancer_types[x]
	meta <- rna.dat_list[[x]][[2]]
	return(meta)
})
rna.dat_meta <- do.call(rbind, rna.dat_meta)

saveRDS(rna.dat_consensus, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts.rds"))
saveRDS(rna.dat_meta, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_metadata.rds"))

###########################################
### Make consensus ATAC-seq data matrix ###
###########################################

pan_cancer_atac_counts_log2norm <- readRDS(paste(dir.path, "ATAC_RNA_integration/resources/Corces_2018/TCGA-ATAC_PanCan_Log2Norm_Counts.rds", sep=""))

# Convert ATAC names to case UUIDs
TCGA_identifier <- read.table(paste(dir.path, "ATAC_RNA_integration/resources/Corces_2018/TCGA_identifier_mapping.txt", sep=""), header=T, sep="\t", stringsAsFactors=F)

# Convert stanford UUIDs (colnames) to case UUIDs with TCGA_identifier
ATAC_ids <- colnames(pan_cancer_atac_counts_log2norm)[8:ncol(pan_cancer_atac_counts_log2norm)]
atac.dat_stanfordUUID <- sapply(strsplit(ATAC_ids, split="_"), function(x) paste(x[[2]], "-", x[[3]], "-", x[[4]], "-", x[[5]], "-", x[[6]], sep=""))
atac.dat_caseUUID <- as.character(TCGA_identifier$Case_UUID[match(atac.dat_stanfordUUID, TCGA_identifier$stanfordUUID)])

# Final ATAC dataset
atac.dat <- as.matrix(pan_cancer_atac_counts_log2norm[ ,8:ncol(pan_cancer_atac_counts_log2norm)])
colnames(atac.dat) <- atac.dat_caseUUID

atac.dat_merged <- merge_technical_replicates(atac.dat)

saveRDS(atac.dat_merged, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-ATAC_PanCan_Log2Norm_Counts_merged.rds"))

#######################################
### Make consensus CIBERSORT matrix ###
#######################################

cibersort_list <- lapply(1:length(cancer_types), function(x) {
	cancer_type <- cancer_types[x]
	print(cancer_type)

	immune.cibersort <- read.table(paste(dir.path, "ATAC_RNA_integration/resources/TCGA/CIBERSORT_RNA_matrices/RESULTS/CIBERSORT.Output_", cancer_type, ".txt", sep=""), sep="\t", header=T, stringsAsFactors=F)
	rownames(immune.cibersort) <- immune.cibersort$Input.Sample
	immune.cibersort$Input.Sample <- NULL
	print(nrow(immune.cibersort))
	return(t(immune.cibersort))
})
cibersort_combined <- do.call(cbind, cibersort_list)
colnames(cibersort_combined) <- as.character(rna.dat_meta$case_UUID[match(colnames(cibersort_combined), as.character(rna.dat_meta$RNA_ID))])

# cibersort_combined <- merge_technical_replicates(cibersort_combined)

saveRDS(cibersort_combined, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan.rds"))

# # Take only unique case UUIDs for RNA metadata (don't need RNA IDs to match cibersort data anymore)
# unique_ids <- colnames(rna.dat_consensus)
# rna.dat_meta <- rna.dat_meta[match(unique_ids, rna.dat_meta$case_UUID), ]
# saveRDS(rna.dat_meta, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_metadata.rds"))

####################
### Pair samples ###
####################

rna.dat_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts.rds"))
rna.dat_meta <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_metadata.rds"))

atac.dat_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-ATAC_PanCan_Log2Norm_Counts_merged.rds"))

cibersort_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan.rds"))

paired.samples <- Reduce(intersect, list(colnames(rna.dat_consensus), colnames(atac.dat_consensus), colnames(cibersort_consensus)))

df.meta_paired <- rna.dat_meta[match(paired.samples, rna.dat_meta$case_UUID), ]

### Remove cancer types with <10 paired samples ###
cancer_types <- names(table(df.meta_paired$cancer_type))[which(as.numeric(table(df.meta_paired$cancer_type)) >= 10)]

# Remove LGG and GBM #
cancer_types <- cancer_types[!cancer_types %in% c("GBMx", "LGGx")]

df.meta_paired <- df.meta_paired[df.meta_paired$cancer_type %in% cancer_types, ]

rna.dat_consensus_paired <- rna.dat_consensus[ ,match(as.character(df.meta_paired$case_UUID), colnames(rna.dat_consensus))]
atac.dat_consensus_paired <- atac.dat_consensus[ ,match(as.character(df.meta_paired$case_UUID), colnames(atac.dat_consensus))]
cibersort_consensus_paired <- cibersort_consensus[ ,match(as.character(df.meta_paired$case_UUID), colnames(cibersort_consensus))]

saveRDS(rna.dat_consensus_paired, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts_PAIRED.rds"))
saveRDS(atac.dat_consensus_paired, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-ATAC_PanCan_Log2Norm_Counts_merged_PAIRED.rds"))
saveRDS(cibersort_consensus_paired, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan_PAIRED.rds"))
saveRDS(df.meta_paired, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-PAIRED_metadata.rds"))

##################
### START HERE ###
##################

##################################################
### Load in TCGA RNA, ATAC, CIBERSORT datasets ###
##################################################

rna.dat_consensus_paired <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts_PAIRED.rds"))
atac.dat_consensus_paired <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-ATAC_PanCan_Log2Norm_Counts_merged_PAIRED.rds"))
cibersort_consensus_paired <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan_PAIRED.rds"))
df.meta_paired <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-PAIRED_metadata.rds"))

cancer_types <- as.character(unique(df.meta_paired$cancer_type))

#######################################################                
### Use immune peak-to-gene links for IFNG.GS peaks ###
#######################################################

immune.links <- read.table(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/resources/ATAC_immune_infiltration_links.txt"), sep="\t", header=T)
                           
# # Take peaks in top quadrant
# plot(immune.links$peak_CorrelationToCytolyticActivity, immune.links$maximum_immuneL2FC, cex=0.6)
# abline(v=0.38, col="red")
# abline(h=0, col="red")

immune.links_interest <- subset(immune.links, peak_CorrelationToCytolyticActivity > 0 & maximum_immuneL2FC > 0.38)

# Remove progenitor peaks
immune.links_interest <- immune.links_interest[immune.links_interest$enriched_cell_type %in% c("mDC", "pDC", "Mono", "CD8", "CD4", "NK", "Bcell"), ]

# # Retain only peaks linked to IFNG.GS genes
immune.links_IFNG.GS <- immune.links_interest[as.character(immune.links_interest$gene) %in% gs.IFNG.GS_hgnc, ]

#####################################################
### Identify important peaks for each cancer type ###
#####################################################

version <- "v1"
method <- "rf2"
scaled <- FALSE

# var_exp_cutoff <- c(0.3, 0.2, 0.4)[3]
# vimp_cutoff <- c(2.5, 1, 5)[3]

# v1
var_exp_cutoff <- 0.5
vimp_cutoff <- 3.5

impt_peaks.list <- lapply(1:length(cancer_types), function(i) {
	impt_peaks_isg <- lapply(1:length(isg.types), function(j) {
		cancer_type <- cancer_types[i]
		isg.type <- isg.types[j]
		impt_peaks <- get_impt_ISG_peaks(var_exp_cutoff=var_exp_cutoff, vimp_cutoff=vimp_cutoff, cancer_type=cancer_type, isg.type=isg.type, version=version, method=method, scaled=scaled)
		return(impt_peaks)
	})
    return(impt_peaks_isg)
})
length(unique(unlist(lapply(impt_peaks.list, function(x) x[[1]]))))

# sort(table(unlist(lapply(1:length(cancer_types), function(x) {
#     genes <- impt_peaks.list[[x]][[2]]
#     return(genes)
# }))))

# ############################
# ### Correlation heatmaps ###
# ############################

# color_palette <- brewer.pal(9, "YlOrRd")
# pal <- colorRampPalette(brewer.pal(11, "Spectral"))
# anno_color <- list(cancer_type=pal(length(cancer_types)))
# names(anno_color$cancer_type) <- cancer_types
# anno_df <- data.frame(cancer_type=df.meta_paired$cancer_type)
# rownames(anno_df) <- colnames(rna.dat_consensus_paired)

# rna.dat_ISG <- rna.dat_consensus_paired[match(c(gs.ISG.RS, gs.IFNG.GS), rownames(rna.dat_consensus_paired)), ]
# rna.dat_ISG <- t(scale(t(rna.dat_ISG)))
# anno_df_row <- data.frame(ISG_type=c(rep("ISG.RS", length(gs.ISG.RS)), rep("IFNG.GS", length(gs.IFNG.GS))))
# rownames(anno_df_row) <- rownames(rna.dat_ISG)

# mat.cor <- cor(t(rna.dat_ISG))
# hm1 <- pheatmap(as.matrix(mat.cor), cluster_cols=TRUE, cluster_rows=TRUE, show_rownames=F, show_colnames=F, annotation_col=anno_df, annotation_colors=anno_color,color=color_palette)

# mat.cor <- cor(rna.dat_ISG)
# hm2 <- pheatmap(as.matrix(mat.cor), cluster_cols=TRUE, cluster_rows=TRUE, show_rownames=F, show_colnames=F, annotation_col=anno_df_row, annotation_colors=anno_color,color=color_palette)

############################
### UMAP plot (RNA ISGs) ###
############################

# umap (version 0.2.3.1)

rna.dat_ISG <- rna.dat_consensus_paired[match(c(gs.ISG.RS, gs.IFNG.GS), rownames(rna.dat_consensus_paired)), ]
rna.dat_ISG <- t(scale(t(rna.dat_ISG)))

p1 <- plot_ISG_umap(dat=rna.dat_ISG, labels=c(rep("ISG.RS", length(gs.ISG.RS)), rep("IFNG.GS", length(gs.IFNG.GS))))
if(plot) ggsave(p1, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_umap/umap_ISGs_RNA.pdf"), width=7.5, height=6)

p2 <- plot_ISG_umap(dat=t(rna.dat_ISG), labels=df.meta_paired$cancer_type)
if(plot) ggsave(p2, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_umap/umap_cancer_types_RNA.pdf"), width=7.5, height=6)

############################
### Make ISG ATAC matrix ###
############################

# # Retain only peaks linked to IFNG.GS genes
# immune.links_IFNG.GS1 <- immune.links_IFNG.GS[as.character(immune.links_IFNG.GS$gene) %in% gs.IFNG.GS_hgnc, ]

consensus_peaks_ISG.RS <- unique(unlist(lapply(impt_peaks.list, function(x) x[[1]])))
consensus_peaks_IFNG.GS <- unique(as.character(immune.links_IFNG.GS$Peak_Name))	# Use Corces et al. immune peaks
consensus_peaks_immune_interest <- unique(as.character(immune.links_interest$Peak_Name))

# Remove peaks from ISG.RS that are in Corces immune link set
remove_idx <- na.omit(match(consensus_peaks_immune_interest, consensus_peaks_ISG.RS))
if(length(remove_idx) > 0) {consensus_peaks_ISG.RS <- consensus_peaks_ISG.RS[-remove_idx]}

# ISG ATAC matrix
adat_ISG_peaks <- atac.dat_consensus_paired[match(c(consensus_peaks_ISG.RS, consensus_peaks_IFNG.GS), rownames(atac.dat_consensus_paired)), ]

p1 <- plot_ISG_umap(dat=adat_ISG_peaks, labels=c(rep("ISG.RS", length(consensus_peaks_ISG.RS)), rep("IFNG.GS", length(consensus_peaks_IFNG.GS))))
if(plot) ggsave(p1, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_umap/umap_ISGs_ATAC.pdf"), width=7.5, height=6)

p2 <- plot_ISG_umap(dat=t(adat_ISG_peaks), labels=df.meta_paired$cancer_type)
if(plot) ggsave(p2, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_umap/umap_cancer_types_ATAC.pdf"), width=7.5, height=6)

#############################################################
### Make ATAC and RNA metafeatures + CIBERSORT data frame ###
#############################################################

version <- "v1"
method <- "rf2"
scaled <- F

### Identify important ATAC peaks ###
# var_exp_cutoff <- c(0.3, 0.2, 0.4)[2]
# vimp_cutoff <- c(2.5, 1, 5)[2]

# v1 - unscaled
var_exp_cutoff <- 0.2
vimp_cutoff <- 1

impt_peaks.list <- lapply(1:length(cancer_types), function(i) {
	impt_peaks_isg <- lapply(1:length(isg.types), function(j) {
		cancer_type <- cancer_types[i]
		isg.type <- isg.types[j]
		impt_peaks <- get_impt_ISG_peaks(var_exp_cutoff=var_exp_cutoff, vimp_cutoff=vimp_cutoff, cancer_type=cancer_type, isg.type=isg.type, version=version, method=method, scaled=scaled)
		return(impt_peaks)
	})
    return(impt_peaks_isg)
})
length(unique(unlist(impt_peaks.list)))

### Make ATAC metafeatures ###

df.meta_list <- lapply(1:length(cancer_types), function(i) {
	cancer_type <- cancer_types[i]
	print(cancer_type)

	caseUUIDs <- as.character(df.meta_paired$case_UUID)[which(as.character(df.meta_paired$cancer_type) == cancer_type)]

    adat <- atac.dat_consensus_paired[ ,match(caseUUIDs, colnames(atac.dat_consensus_paired))]
    edat <- rna.dat_consensus_paired[ ,match(caseUUIDs, colnames(rna.dat_consensus_paired))]
    ciber.dat <- cibersort_consensus_paired[ ,match(caseUUIDs, colnames(cibersort_consensus_paired))]

    # ATAC data with important peaks from mRFAR # 
    adat_ISG.RS <- adat[match(impt_peaks.list[[i]][[1]], rownames(adat)), ]

    # Remove peaks from ISG.RS that are in Corces immune link set
    # remove_idx <- na.omit(match(as.character(unique(immune.links_IFNG.GS$Peak_Name)), rownames(adat_ISG.RS)))
    remove_idx <- na.omit(match(unique(as.character(immune.links_interest$Peak_Name)), rownames(adat_ISG.RS)))
    if(length(remove_idx) > 0) {adat_ISG.RS <- adat_ISG.RS[-remove_idx, ]}
  
  	# ATAC data with immune peaks
  	adat_IFNG.GS_Corces <- adat[match(unique(as.character(immune.links_IFNG.GS$Peak_Name)), rownames(adat)), ]

  	# RNA data
  	edat_ISG.RS <- edat[na.omit(match(gs.ISG.RS, rownames(edat))), ]
  	edat_IFNG.GS <- edat[na.omit(match(gs.IFNG.GS, rownames(edat))), ]

  	### Make ATAC and RNA metagenes ###

  	# ATAC #
  	ISG.RS_meta_ATAC <- colMeans(t(scale(t(adat_ISG.RS))))
  	IFNG.GS_meta_ATAC  <- colMeans(t(scale(t(adat_IFNG.GS_Corces))))
  	# ISG.RS_meta_ATAC <- colMeans(adat_ISG.RS)
  	# IFNG.GS_meta_ATAC  <- colMeans(adat_IFNG.GS_Corces)

  	ISG.ratio_ATAC <- IFNG.GS_meta_ATAC / ISG.RS_meta_ATAC 
  	ISG.binary_ATAC <- ifelse(IFNG.GS_meta_ATAC > ISG.RS_meta_ATAC, 1, 0)

  	# RNA #
  	ISG.RS_meta_RNA <- colMeans(t(scale(t(edat_ISG.RS))))
  	IFNG.GS_meta_RNA  <- colMeans(t(scale(t(edat_IFNG.GS))))

  	ISG.ratio_RNA <- IFNG.GS_meta_RNA / ISG.RS_meta_RNA
  	ISG.binary_RNA <- ifelse(IFNG.GS_meta_RNA > ISG.RS_meta_RNA, 1, 0)

  	### Merge ISG.RS/IFNG.GS accessibility and immune cell infiltration ###

  	df <- data.frame(ISG.RS_ATAC=ISG.RS_meta_ATAC, IFNG.GS_ATAC=IFNG.GS_meta_ATAC, ISG.binary_ATAC=ISG.binary_ATAC, ISG.RS_RNA=ISG.RS_meta_RNA, IFNG.GS_RNA=IFNG.GS_meta_RNA, ISG.binary_RNA=ISG.binary_RNA, T.cells.CD8=as.numeric(ciber.dat[which(rownames(ciber.dat) == "T.cells.CD8"), ]), Macrophages.M0=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M0"), ]), Macrophages.M1=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M1"), ]), Macrophages.M2=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M2"), ]))
  	return(df)
})

###################################################
### Correlation of ISG.RS and IFNG.GS gene sets ###
###################################################

df.meta_pancancer <- do.call(rbind, df.meta_list)

# Decide on CD8 class cutoffs
quartiles <- quantile(df.meta_pancancer$T.cells.CD8, probs=seq(0,1,0.25))
quantiles <- quantile(df.meta_pancancer$T.cells.CD8, probs=seq(0,1,0.1))

# cd8_cutoff <- 0.2
cd8_upper_cutoff <- quantiles[10]
cd8_lower_cutoff <- quantiles[2]

plot(1:length(df.meta_pancancer$T.cells.CD8), sort(df.meta_pancancer$T.cells.CD8))
# abline(h=cd8_cutoff, col="red")
# length(which(df.meta_pancancer$T.cells.CD8 > cd8_cutoff))
abline(h=cd8_upper_cutoff, col="red")
abline(h=cd8_lower_cutoff, col="red")

# df.meta_pancancer$cd8.class <- ifelse(df.meta_pancancer$T.cells.CD8 > cd8_cutoff, "high", "low")
df.meta_pancancer$cd8.class <- ifelse(df.meta_pancancer$T.cells.CD8 > cd8_upper_cutoff, "high", ifelse(df.meta_pancancer$T.cells.CD8 < cd8_lower_cutoff, "low", "int"))

# ggplot(df.meta_pancancer, aes(ISG.RS_RNA, IFNG.GS_RNA)) +
ggplot(df.meta_pancancer[df.meta_pancancer$cd8.class %in% c("high", "low"), ], aes(ISG.RS_RNA, IFNG.GS_RNA)) +
	# geom_point(aes(color=cd8.class, size=T.cells.CD8), shape=1, stroke=0.7) +
	geom_point(aes(color=cd8.class), size=4, shape=1, stroke=1) +
	scale_color_manual(values=c("red", "blue")) +
	# scale_color_manual(values=c("red", "yellow", "blue")) +
	geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
	# xlim(c(-2,2)) +
	# ylim(c(-2,2)) +
	theme_bw()

# ggplot(df.meta_pancancer, aes(ISG.RS_ATAC, IFNG.GS_ATAC)) +
ggplot(df.meta_pancancer[df.meta_pancancer$cd8.class %in% c("high", "low"), ], aes(ISG.RS_ATAC, IFNG.GS_ATAC)) +
	# geom_point(aes(color=cd8.class, size=T.cells.CD8), shape=1, stroke=0.7) +
	geom_point(aes(color=cd8.class), size=4, shape=1, stroke=1) +
	scale_color_manual(values=c("red", "blue")) +
	# scale_color_manual(values=c("red", "yellow", "blue")) +
	geom_abline(intercept=0, slope=1, col="red", linetype="dashed") +
	xlim(c(-1.2,1.2)) +
	ylim(c(-1.2,1.2)) +
	theme_bw()

plot(df.meta_pancancer$ISG.RS_RNA, df.meta_pancancer$IFNG.GS_RNA, xlab="ISG.RS", ylab="IFNG.GS", xlim=c(-2,2), ylim=c(-2,2), main=paste0("RNA (cor=", round(cor(df.meta_pancancer$ISG.RS_RNA, df.meta_pancancer$IFNG.GS_RNA), 2), ")"))
abline(a=0, b=1, col="red")

plot(df.meta_pancancer$ISG.RS_ATAC, df.meta_pancancer$IFNG.GS_ATAC, xlab="ISG.RS", ylab="IFNG.GS", xlim=c(-1.5,1.5), ylim=c(-1.5,1.5), main=paste0("ATAC (cor=", round(cor(df.meta_pancancer$ISG.RS_ATAC, df.meta_pancancer$IFNG.GS_ATAC), 2), ")"))
abline(a=0, b=1, col="red")

#####################################################

ggplot(df.meta_pancancer, aes(ISG.RS_RNA, T.cells.CD8)) +
# ggplot(df.meta_pancancer[df.meta_pancancer$cd8.class %in% c("high", "low"), ], aes(ISG.RS_RNA, IFNG.GS_RNA)) +
	# geom_point(aes(color=cd8.class, size=T.cells.CD8), shape=1, stroke=0.7) +
	geom_point(aes(color=cd8.class), size=4, shape=1, stroke=1) +
	# scale_color_manual(values=c("red", "blue")) +
	scale_color_manual(values=c("red", "yellow", "blue")) +
	# xlim(c(-2,2)) +
	# ylim(c(-2,2)) +
	theme_bw()

par(mfrow=c(1,2))
fit <- lm(df.meta_pancancer$T.cells.CD8~df.meta_pancancer$IFNG.GS_RNA)
plot(df.meta_pancancer$IFNG.GS_RNA, df.meta_pancancer$T.cells.CD8, main=paste0("IFNG.GS RNA (cor=", round(cor(df.meta_pancancer$IFNG.GS_RNA, df.meta_pancancer$T.cells.CD8), 2), ")"))
abline(a=fit$coefficients[1], b=fit$coefficients[2], col="red", lty="dotted")

fit <- lm(df.meta_pancancer$T.cells.CD8~df.meta_pancancer$ISG.RS_RNA)
plot(df.meta_pancancer$ISG.RS_RNA, df.meta_pancancer$T.cells.CD8, main=paste0("ISG.RS RNA (cor=", round(cor(df.meta_pancancer$ISG.RS_RNA, df.meta_pancancer$T.cells.CD8), 2), ")"))
abline(a=fit$coefficients[1], b=fit$coefficients[2], col="red", lty="dotted")

par(mfrow=c(1,2))
fit <- lm(df.meta_pancancer$T.cells.CD8~df.meta_pancancer$IFNG.GS_ATAC)
plot(df.meta_pancancer$IFNG.GS_ATAC, df.meta_pancancer$T.cells.CD8, main=paste0("IFNG.GS ATAC (cor=", round(cor(df.meta_pancancer$IFNG.GS_ATAC, df.meta_pancancer$T.cells.CD8), 2), ")"))
abline(a=fit$coefficients[1], b=fit$coefficients[2], col="red", lty="dotted")

fit <- lm(df.meta_pancancer$T.cells.CD8~df.meta_pancancer$ISG.RS_ATAC)
plot(df.meta_pancancer$ISG.RS_ATAC, df.meta_pancancer$T.cells.CD8, main=paste0("ISG.RS ATAC (cor=", round(cor(df.meta_pancancer$ISG.RS_ATAC, df.meta_pancancer$T.cells.CD8), 2), ")"))
abline(a=fit$coefficients[1], b=fit$coefficients[2], col="red", lty="dotted")

# ## Boxplots ##

# ggplot(df.meta_pancancer, aes(x=factor(ISG.binary_RNA), y=T.cells.CD8, color=factor(ISG.binary_RNA))) +
# # ggplot(df.meta_pancancer[df.meta_pancancer$cd8.class %in% c("high", "low"), ], aes(ISG.RS_RNA, IFNG.GS_RNA)) +
# 	# geom_point(aes(color=cd8.class, size=T.cells.CD8), shape=1, stroke=0.7) +
# 	# geom_point(aes(color=cd8.class), size=4, shape=1, stroke=1) +
# 	geom_boxplot() +
# 	geom_jitter(position=position_jitter(0.1), size=0.2) +
# 	# scale_color_manual(values=c("red", "blue")) +
# 	scale_color_manual(values=c("blue", "red")) +
# 	# xlim(c(-2,2)) +
# 	# ylim(c(-2,2)) +
# 	theme_bw()

# ggplot(df.meta_pancancer, aes(x=factor(ISG.binary_ATAC), y=T.cells.CD8, color=factor(ISG.binary_ATAC))) +
# # ggplot(df.meta_pancancer[df.meta_pancancer$cd8.class %in% c("high", "low"), ], aes(ISG.RS_RNA, IFNG.GS_RNA)) +
# 	# geom_point(aes(color=cd8.class, size=T.cells.CD8), shape=1, stroke=0.7) +
# 	# geom_point(aes(color=cd8.class), size=4, shape=1, stroke=1) +
# 	geom_boxplot() +
# 	geom_jitter(position=position_jitter(0.1), size=0.2) +
# 	# scale_color_manual(values=c("red", "blue")) +
# 	scale_color_manual(values=c("blue", "red")) +
# 	# xlim(c(-2,2)) +
# 	# ylim(c(-2,2)) +
# 	theme_bw()

################################
### Scatter plots ISG vs CD8 ###
################################

assay <- c("RNA", "ATAC")
isgs <- c("IFNG.GS", "ISG.RS")
pseudocount <- 0.0001

cancer <- "LIHC"

df <- df.meta_list[[which(cancer_types == cancer)]]
df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

par(mfrow=c(2,2))
for(i in 1:2) {
	for(j in 1:2) {
		name <- paste0(isgs[i], "_", assay[j])

		df$isg <- df[,which(colnames(df) == name)]
		fit <- betareg(T.cells.CD8 ~ isg, data = df)

		plot(df$isg, df$T.cells.CD8, main=paste0(name, " (cor=", round(cor(df$isg, df$T.cells.CD8), 2), ")"), xlab=isgs[i], ylab="CD8+ infiltration")
		lines(seq(min(df$isg), max(df$isg), 0.05), predict(fit, newdata=data.frame(isg=seq(min(df$isg), max(df$isg), 0.05))), col="red", lwd=2, lty="dotted")
	}	
}

########################################

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/RNA_IFNG.GS_v_CD8_correlation PAIRED ", method, ".pdf"), width=12, height=4)
par(mfrow=c(2,7))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount
	
	fit <- betareg(T.cells.CD8 ~ IFNG.GS_RNA, data = df)

	plot(df$IFNG.GS_RNA, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$IFNG.GS_RNA, df$T.cells.CD8), 2), ")"))
	lines(seq(min(df$IFNG.GS_RNA), max(df$IFNG.GS_RNA), 0.05), predict(fit, newdata=data.frame(IFNG.GS_RNA=seq(min(df$IFNG.GS_RNA), max(df$IFNG.GS_RNA), 0.05))), col="red", lwd=2, lty="dotted")
}
dev.off()

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/RNA_ISG.RS_v_CD8_correlation PAIRED ", method, ".pdf"), width=12, height=4)
par(mfrow=c(2,7))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	fit <- betareg(T.cells.CD8 ~ ISG.RS_RNA, data = df)

	plot(df$ISG.RS_RNA, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$ISG.RS_RNA, df$T.cells.CD8), 2), ")"))
	lines(seq(min(df$ISG.RS_RNA), max(df$ISG.RS_RNA), 0.05), predict(fit, newdata=data.frame(ISG.RS_RNA=seq(min(df$ISG.RS_RNA), max(df$ISG.RS_RNA), 0.05))), col="red", lwd=2, lty="dotted")
}
dev.off()

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/ATAC_IFNG.GS_v_CD8_correlation PAIRED ", method, ".pdf"), width=12, height=4)
par(mfrow=c(2,7))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	fit <- betareg(T.cells.CD8 ~ IFNG.GS_ATAC, data = df)

	plot(df$IFNG.GS_ATAC, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$IFNG.GS_ATAC, df$T.cells.CD8), 2), ")"))
	lines(seq(min(df$IFNG.GS_ATAC), max(df$IFNG.GS_ATAC), 0.01), predict(fit, newdata=data.frame(IFNG.GS_ATAC=seq(min(df$IFNG.GS_ATAC), max(df$IFNG.GS_ATAC), 0.01))), col="red", lwd=2, lty="dotted")
}
dev.off()

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/ATAC_ISG.RS_v_CD8_correlation PAIRED ", method, ".pdf"), width=12, height=4)
par(mfrow=c(2,7))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	if(length(na.omit(df$ISG.RS_ATAC)) > 0) {
		fit <- betareg(T.cells.CD8 ~ ISG.RS_ATAC, data = df)

		plot(df$ISG.RS_ATAC, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$ISG.RS_ATAC, df$T.cells.CD8), 2), ")"))
		lines(seq(min(df$ISG.RS_ATAC), max(df$ISG.RS_ATAC), 0.01), predict(fit, newdata=data.frame(ISG.RS_ATAC=seq(min(df$ISG.RS_ATAC), max(df$ISG.RS_ATAC), 0.01))), col="red", lwd=2, lty="dotted")
	}
}
dev.off()

###################################################################
### Predict immune infiltrate with ISG accessibility/expression ###
###################################################################

### Fit model for all cancers (PAIRED data) ###
assay <- "ATAC"
immune.cells <- c("CD8", "M1", "M2")[1]
pseudocount <- 1e-8

fits_list <- lapply(1:length(cancer_types), function(x) {
	print(cancer_types[x])

	df <- df.meta_list[[x]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	# ATAC #
	if(length(na.omit(df$ISG.RS_ATAC)) == 0) {
		f1.atac <- NA
	} else {
		f1.atac <- betareg(T.cells.CD8 ~ ISG.RS_ATAC + IFNG.GS_ATAC, data = df)
	}

	# RNA #
	f1.rna <- betareg(T.cells.CD8 ~ ISG.RS_RNA + IFNG.GS_RNA, data = df)
	return(list(f1.atac, f1.rna))
})

### Get model statistics ###

fit.stats <- lapply(1:length(cancer_types), function(x) {
	cancer_type <- cancer_types[x]
    print(cancer_type)

    if(assay == "ATAC") {
		fit <- fits_list[[x]][[1]]
	} else {
		fit <- fits_list[[x]][[2]]
	}

	if(!is.na(fit)) {
		ISG.RS_summary <- get_betareg_stats(fit, name=paste0("ISG.RS_", assay))
		IFNG.GS_summary <- get_betareg_stats(fit, name=paste0("IFNG.GS_", assay))
	} else {
		ISG.RS_summary <- rep(NA, 4)
		IFNG.GS_summary <- rep(NA, 4)
	}

	ISG.summary <- list(ISG.RS_summary, IFNG.GS_summary)
	df.stats <- data.frame(cancer_type=cancer_type, immune.cell="CD8", gs=c("ISG.RS", "IFNG.GS"), est=sapply(ISG.summary, function(x) x[1]), se=sapply(ISG.summary, function(x) x[2]), lower=sapply(ISG.summary, function(x) x[3]), upper=sapply(ISG.summary, function(x) x[4]))
	return(df.stats)
})
fit.stats <- do.call(rbind, fit.stats)

### Forest plot ###

sample_sizes <- sapply(1:length(cancer_types), function(x) nrow(df.meta_list[[x]]))

# ISG.RS #
df.ISG.RS <- fit.stats[which(as.character(fit.stats$gs) == "ISG.RS"), ]
df.ISG.RS$name <- paste0(df.ISG.RS$cancer_type, " (n=", sample_sizes, ")")
df.ISG.RS.sorted <- df.ISG.RS[sort.int(df.ISG.RS$est, decreasing = TRUE, index.return = TRUE)$ix, ]
df.ISG.RS.sorted$name <- factor(df.ISG.RS.sorted$name, levels=rev(as.character(df.ISG.RS.sorted$name)))

# IFNG.GS #
df.IFNG.GS <- fit.stats[which(as.character(fit.stats$gs) == "IFNG.GS"), ]
df.IFNG.GS$name <- paste0(df.IFNG.GS$cancer_type, " (n=", sample_sizes, ")")
df.IFNG.GS.sorted <- df.IFNG.GS[sort.int(df.IFNG.GS$est, decreasing = TRUE, index.return = TRUE)$ix, ]
df.IFNG.GS.sorted$name <- factor(df.IFNG.GS.sorted$name, levels=rev(as.character(df.IFNG.GS.sorted$name)))

p_coeff.isg.rs <- ggplot(data=df.ISG.RS.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
    geom_pointrange(color = "firebrick", size = 0.7) +
    geom_hline(yintercept=0, lty=2) +
    coord_flip() +
    xlab("Cancer type") +
    ylab("Effect size (95% CI)") +
    ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.ISG.RS.sorted$upper) + 0.01) +
    # ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
    ggtitle("ISG.RS") +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 9.8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))

p_coeff.ifng.gs <- ggplot(data=df.IFNG.GS.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
    geom_pointrange(color = "royalblue", size = 0.7) +
    geom_hline(yintercept=0, lty=2) +
    coord_flip() +
    xlab("Cancer type") +
    ylab("Effect size (95% CI)") +
    ylim(min(df.IFNG.GS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
    # ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
    ggtitle("IFNG.GS") +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 9.8),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))

# ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, ncol=1)

##########################################
### Fit linear model for combined data ###
##########################################

# ALL PAIRED SAMPLES #
df <- do.call(rbind, df.meta_list)
df$T.cells.CD8 <- df$T.cells.CD8 + pseudocount

# ATAC #
f1.atac <- betareg(T.cells.CD8 ~ ISG.RS_ATAC + IFNG.GS_ATAC, data = df)

# RNA #
f1.rna <- betareg(T.cells.CD8 ~ ISG.RS_RNA + IFNG.GS_RNA, data = df)

assay <- "ATAC"
if(assay == "ATAC") {
	fit <- f1.atac
} else {
	fit <- f1.rna
}
ISG.RS_summary <- get_betareg_stats(fit, name=paste0("ISG.RS_", assay))
IFNG.GS_summary <- get_betareg_stats(fit, name=paste0("IFNG.GS_", assay))
ISG.summary <- list(ISG.RS_summary, IFNG.GS_summary)

df.stats <- data.frame(cancer_type="ALL", immune.cell="CD8", gs=c("ISG.RS", "IFNG.GS"), est=sapply(ISG.summary, function(x) x[1]), se=sapply(ISG.summary, function(x) x[2]), lower=sapply(ISG.summary, function(x) x[3]), upper=sapply(ISG.summary, function(x) x[4]))
df.stats$gs <- factor(df.stats$gs, levels=c("ISG.RS", "IFNG.GS"))

p_coeff.comb <- ggplot(data=df.stats, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
	geom_pointrange(aes(color = gs), size = 0.75) +
	geom_hline(yintercept=0, lty = 2) +
	scale_color_manual(values = c("firebrick", "royalblue")) +
	coord_flip() +
	xlab("Cancer type") +
	ylab("Effect size (95% CI)") +
	ylim(min(df.stats$lower) - 0.01, max(df.stats$upper) + 0.01) +
	# ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
	ggtitle(paste0("ALL Cancers ", assay, " (n=", nrow(df), ")")) +
	theme_bw() +
	theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 9.8),
        axis.title.x=element_text(size=10),
        axis.title.y=element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# egg::ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, p_coeff.comb, ncol = 1, heights = c(1, 1, 0.2))

if (plot) pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression ATAC PAIRED rf2.pdf"), height = 6.8, width = 3.5)
egg::ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, p_coeff.comb, ncol = 1, heights = c(1, 1, 0.2))
if (plot) dev.off()

# if (plot) pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression ATAC PAIRED rf2 centered.pdf"), height = 5, width = 3)
# egg::ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, p_coeff.comb, ncol = 1, heights = c(1, 1, 0.2))
# if (plot) dev.off()

###########################################
### Forest plot for all RNA-seq samples ###
###########################################

rna.dat_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts.rds"))
rna.dat_meta <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_metadata.rds"))

cibersort_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan.rds"))

cancer_types <- as.character(unique(rna.dat_meta$cancer_type))

### Make RNA metagenes ###

# df.meta_list <- lapply(1:length(cancer_types), function(i) {
# 	cancer_type <- cancer_types[i]
# 	print(cancer_type)
# 	df <- make_ISG_metagenes(cancer_type=cancer_type, gs1=gs.IFNG.GS, gs2=gs.IFNG.GS, rna.dat=rna.dat_consensus, cibersort.dat=cibersort_consensus, rna.meta=rna.dat_meta)
#   	return(df)
# })

df.meta_list <- lapply(1:length(cancer_types), function(i) {
	cancer_type <- cancer_types[i]
	print(cancer_type)

	caseUUIDs <- unique(as.character(rna.dat_meta$case_UUID)[which(as.character(rna.dat_meta$cancer_type) == cancer_type)])

    edat <- rna.dat_consensus[ ,match(caseUUIDs, colnames(rna.dat_consensus))]
    ciber.dat <- cibersort_consensus[ ,match(caseUUIDs, colnames(cibersort_consensus))]

  	# RNA data
  	edat_ISG.RS <- edat[na.omit(match(gs.ISG.RS, rownames(edat))), ]
  	edat_IFNG.GS <- edat[na.omit(match(gs.IFNG.GS, rownames(edat))), ]

  	### Make RNA metagenes ###

  	ISG.RS_meta_RNA <- colMeans(t(scale(t(edat_ISG.RS))))
  	IFNG.GS_meta_RNA  <- colMeans(t(scale(t(edat_IFNG.GS))))

  	ISG.ratio_RNA <- IFNG.GS_meta_RNA / ISG.RS_meta_RNA
  	ISG.binary_RNA <- ifelse(IFNG.GS_meta_RNA > ISG.RS_meta_RNA, 1, 0)

  	### Merge ISG.RS/IFNG.GS accessibility and immune cell infiltration ###

  	df <- data.frame(ISG.RS_RNA=ISG.RS_meta_RNA, IFNG.GS_RNA=IFNG.GS_meta_RNA, ISG.binary_RNA=ISG.binary_RNA, T.cells.CD8=as.numeric(ciber.dat[which(rownames(ciber.dat) == "T.cells.CD8"), ]), Macrophages.M0=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M0"), ]), Macrophages.M1=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M1"), ]), Macrophages.M2=as.numeric(ciber.dat[which(rownames(ciber.dat) == "Macrophages.M2"), ]))
  	return(df)
})

### ISG vs CD8 scatterplot ###

pseudocount <- 0.01

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/RNA_IFNG.GS_v_CD8_correlation_ALL.pdf"), width=12, height=10.5)
par(mfrow=c(4,6))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount
	
	fit <- betareg(T.cells.CD8 ~ IFNG.GS_RNA, data = df)

	plot(df$IFNG.GS_RNA, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$IFNG.GS_RNA, df$T.cells.CD8), 2), ")"), cex=0.4)
	lines(seq(min(df$IFNG.GS_RNA), max(df$IFNG.GS_RNA), 0.05), predict(fit, newdata=data.frame(IFNG.GS_RNA=seq(min(df$IFNG.GS_RNA), max(df$IFNG.GS_RNA), 0.05))), col="red", lwd=2, lty="dotted")
}
dev.off()

pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/ISG_v_CD8_correlation/RNA_ISG.RS_v_CD8_correlation_ALL.pdf"), width=12, height=10.5)
par(mfrow=c(4,6))
for(i in 1:length(cancer_types)) {
	df <- df.meta_list[[i]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	fit <- betareg(T.cells.CD8 ~ ISG.RS_RNA, data = df)

	plot(df$ISG.RS_RNA, df$T.cells.CD8, main=paste0(cancer_types[i], " (cor=", round(cor(df$ISG.RS_RNA, df$T.cells.CD8), 2), ")"), cex=0.4)
	lines(seq(min(df$ISG.RS_RNA), max(df$ISG.RS_RNA), 0.05), predict(fit, newdata=data.frame(ISG.RS_RNA=seq(min(df$ISG.RS_RNA), max(df$ISG.RS_RNA), 0.05))), col="red", lwd=2, lty="dotted")
}
dev.off()

### Fit beta regression model for all cancers ###

immune.cells <- c("CD8", "M1", "M2")[1]
pseudocount <- 0.01

fits_list <- lapply(1:length(cancer_types), function(x) {
	print(cancer_types[x])

	df <- df.meta_list[[x]]
	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

	# RNA #
	f1.rna <- betareg(T.cells.CD8 ~ ISG.RS_RNA + IFNG.GS_RNA, data = df)
	return(f1.rna)
})

### Get model statistics ###

fit.stats <- lapply(1:length(cancer_types), function(x) {
	cancer_type <- cancer_types[x]
    print(cancer_type)

	fit <- fits_list[[x]]

	ISG.RS_summary <- get_betareg_stats(fit, name=paste0("ISG.RS_RNA"))
	IFNG.GS_summary <- get_betareg_stats(fit, name=paste0("IFNG.GS_RNA"))

	ISG.summary <- list(ISG.RS_summary, IFNG.GS_summary)
	df.stats <- data.frame(cancer_type=cancer_type, immune.cell="CD8", gs=c("ISG.RS", "IFNG.GS"), est=sapply(ISG.summary, function(x) x[1]), se=sapply(ISG.summary, function(x) x[2]), lower=sapply(ISG.summary, function(x) x[3]), upper=sapply(ISG.summary, function(x) x[4]))
	return(df.stats)
})
fit.stats <- do.call(rbind, fit.stats)

### Forest plot ###

sample_sizes <- sapply(1:length(cancer_types), function(x) nrow(df.meta_list[[x]]))

# ISG.RS #
df.ISG.RS <- fit.stats[which(as.character(fit.stats$gs) == "ISG.RS"), ]
df.ISG.RS$name <- paste0(df.ISG.RS$cancer_type, " (n=", sample_sizes, ")")
df.ISG.RS.sorted <- df.ISG.RS[sort.int(df.ISG.RS$est, decreasing = TRUE, index.return = TRUE)$ix, ]
df.ISG.RS.sorted$name <- factor(df.ISG.RS.sorted$name, levels=rev(as.character(df.ISG.RS.sorted$name)))

# IFNG.GS #
df.IFNG.GS <- fit.stats[which(as.character(fit.stats$gs) == "IFNG.GS"), ]
df.IFNG.GS$name <- paste0(df.IFNG.GS$cancer_type, " (n=", sample_sizes, ")")
df.IFNG.GS.sorted <- df.IFNG.GS[sort.int(df.IFNG.GS$est, decreasing = TRUE, index.return = TRUE)$ix, ]
df.IFNG.GS.sorted$name <- factor(df.IFNG.GS.sorted$name, levels=rev(as.character(df.IFNG.GS.sorted$name)))

p_coeff.isg.rs <- ggplot(data=df.ISG.RS.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
    geom_pointrange(color = "firebrick", size = 0.4) +
    geom_hline(yintercept=0, lty=2) +
    coord_flip() +
    xlab("Cancer type") +
    ylab("Effect size (95% CI)") +
    ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
    ggtitle("ISG.RS") +
    theme_classic() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 7.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))

p_coeff.ifng.gs <- ggplot(data=df.IFNG.GS.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
    geom_pointrange(color = "royalblue", size = 0.4) +
    geom_hline(yintercept=0, lty=2) +
    coord_flip() +
    xlab("Cancer type") +
    ylab("Effect size (95% CI)") +
    ylim(min(df.ISG.RS.sorted$lower) - 0.01, max(df.IFNG.GS.sorted$upper) + 0.01) +
    ggtitle("IFNG.GS") +
    theme_classic() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 7.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10))

### Combine samples ###

df <- do.call(rbind, df.meta_list)
df$T.cells.CD8 <- df$T.cells.CD8 + pseudocount

# RNA #
fit <- betareg(T.cells.CD8 ~ ISG.RS_RNA + IFNG.GS_RNA, data = df)

ISG.RS_summary <- get_betareg_stats(fit, name=paste0("ISG.RS_RNA"))
IFNG.GS_summary <- get_betareg_stats(fit, name=paste0("IFNG.GS_RNA"))
ISG.summary <- list(ISG.RS_summary, IFNG.GS_summary)

df.stats <- data.frame(cancer_type="ALL", immune.cell="CD8", gs=c("ISG.RS", "IFNG.GS"), est=sapply(ISG.summary, function(x) x[1]), se=sapply(ISG.summary, function(x) x[2]), lower=sapply(ISG.summary, function(x) x[3]), upper=sapply(ISG.summary, function(x) x[4]))
df.stats$gs <- factor(df.stats$gs, levels=c("ISG.RS", "IFNG.GS"))

p_coeff.comb <- ggplot(data=df.stats, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
	geom_pointrange(aes(color = gs), size = 0.4) +
	geom_hline(yintercept=0, lty = 2) +
	scale_color_manual(values = c("firebrick", "royalblue")) +
	coord_flip() +
	xlab("Cancer type") +
	ylab("Effect size (95% CI)") +
	ylim(min(df.stats$lower) - 0.01, max(df.stats$upper) + 0.01) +
	ggtitle(paste0("ALL Cancers RNA (n=", nrow(df), ")")) +
	theme_classic() +
	theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 7.5),
        axis.title.x=element_text(size=10),
        axis.title.y=element_blank(),
        plot.title = element_text(size = 10),
        legend.position = "none")
# egg::ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, p_coeff.comb, ncol = 1, heights = c(1, 1, 0.2))

if (plot) pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression RNA ALL.pdf"), height = 7.5, width = 3.5)
# p_coeff.ifng.gs / p_coeff.isg.rs / p_coeff.comb
egg::ggarrange(p_coeff.ifng.gs, p_coeff.isg.rs, p_coeff.comb, ncol = 1, heights = c(1, 1, 0.2))
if (plot) dev.off()

################################################################
### Test random gene sets association with CD8+ infiltration ###
################################################################

rna.dat_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_DESeq2_vst_Counts.rds"))
rna.dat_meta <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-RNA_PanCan_metadata.rds"))

cibersort_consensus <- readRDS(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/data/final/TCGA-CIBERSORT_PanCan.rds"))

cancer_types <- as.character(unique(rna.dat_meta$cancer_type))

####################################################################

n <- 100
pseudocount <- 0.01

num_reps <- 1

set.seed(1)

random_beta_est <- lapply(1:num_reps, function(x) {
	df.stats_random <- lapply(1:n, function(i) {
		print(i)

		gs.rand1 <- rownames(rna.dat_consensus)[sample(1:nrow(rna.dat_consensus), length(gs.IFNG.GS))]
		# gs.rand1 <- gs.IFNG.GS
		gs.rand2 <- rownames(rna.dat_consensus)[sample(1:nrow(rna.dat_consensus), length(gs.ISG.RS))]
		# gs.rand2 <- gs.ISG.RS

		### Make RNA metagenes ###

		df.meta_list <- lapply(1:length(cancer_types), function(x) {
			cancer_type <- cancer_types[x]
			df <- make_ISG_metagenes(cancer_type=cancer_type, gs1=gs.rand1, gs2=gs.rand2, rna.dat=rna.dat_consensus, cibersort.dat=cibersort_consensus, rna.meta=rna.dat_meta)
		  	return(df)
		})

		### Combine samples ###

		df <- do.call(rbind, df.meta_list)
		df$T.cells.CD8 <- df$T.cells.CD8 + pseudocount

		# RNA #
		fit <- betareg(T.cells.CD8 ~ gs1_RNA + gs2_RNA, data = df)

		gs1_summary <- get_betareg_stats(fit, name=paste0("gs1_RNA"))
		gs2_summary <- get_betareg_stats(fit, name=paste0("gs2_RNA"))
		gs.summary <- list(gs1_summary, gs2_summary)

		df.stats <- data.frame(rep=paste0("rep", i), immune.cell="CD8", gs=c("gs.random1", "gs.random2"), est=sapply(gs.summary, function(x) x[1]), se=sapply(gs.summary, function(x) x[2]), lower=sapply(gs.summary, function(x) x[3]), upper=sapply(gs.summary, function(x) x[4]))
		df.stats$gs <- factor(df.stats$gs, levels=c("gs.random1", "gs.random2"))
		return(df.stats)
	})
	df.stats_random <- do.call(rbind, df.stats_random)

	# Compute mean, CI for coefficient estimate
	df.gs1 <- df.stats_random[which(as.character(df.stats_random$gs) == "gs.random1"), ]
	df.gs2 <- df.stats_random[which(as.character(df.stats_random$gs) == "gs.random2"), ]

	t1 <- t.test(df.gs1$est)
	t2 <- t.test(df.gs2$est)

	df.stats <- data.frame(gs=c("gs.random1", "gs.random2"), est=c(t1$estimate, t2$estimate), lower=c(t1$conf.int[1], t2$conf.int[1]), upper=c(t1$conf.int[2], t2$conf.int[2]), rep=x)
	print(df.stats)
	return(df.stats)
})
random_beta_est <- do.call(rbind, random_beta_est)
random_beta_est$gs <- factor(random_beta_est$gs, levels=c("gs.random2", "gs.random1"))

fig.rand <- ggplot(data=random_beta_est, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
    geom_pointrange(aes(color=gs), size = 1.2) +
    geom_hline(yintercept=0, lty=2) +
    scale_color_manual(values = c("firebrick", "royalblue")) +
    coord_flip() +
    xlab("") +
    ylab("Effect size (95% CI)") +
    ylim(-0.25, 0.25) +
    ggtitle(paste0("Random gene sets (n=", n, ")")) +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 16),
          legend.position = "none")

if(plot) ggsave(fig.rand, file=paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression RNA random gene sets.pdf"), height=2.5, width=4)

# random_beta_est1 <- random_beta_est[which(random_beta_est$gs == "gs.random1"), ]
# random_beta_est2 <- random_beta_est[which(random_beta_est$gs == "gs.random2"), ]

# p_rand1 <- ggplot(data=random_beta_est1, aes(x=rep, y=est, ymin=lower, ymax=upper)) +
#     geom_pointrange(color = "royalblue", size = 0.4) +
#     geom_hline(yintercept=0, lty=2) +
#     coord_flip() +
#     xlab("Rep") +
#     ylab("Effect size (95% CI)") +
#     ylim(-0.5, 0.5) +
#     ggtitle("gs.random1 (n=50)") +
#     theme_bw() +
#     theme(axis.ticks = element_line(color = "black"),
#           axis.text = element_text(color = "black", size = 8),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           plot.title = element_text(size = 10))

# p_rand2 <- ggplot(data=random_beta_est2, aes(x=rep, y=est, ymin=lower, ymax=upper)) +
#     geom_pointrange(color = "firebrick", size = 0.4) +
#     geom_hline(yintercept=0, lty=2) +
#     coord_flip() +
#     xlab("Rep") +
#     ylab("Effect size (95% CI)") +
#     ylim(-0.5, 0.5) +
#     ggtitle("gs.random2 (n=50)") +
#     theme_bw() +
#     theme(axis.ticks = element_line(color = "black"),
#           axis.text = element_text(color = "black", size = 8),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           plot.title = element_text(size = 10))

# if (plot) pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression RNA random gene sets.pdf"), height = 6.8, width = 3.5)
# ggarrange(p_rand1, p_rand2, ncol=1)
# dev.off()

# ######################################################
# ######################################################

# set.seed(1)
# gs.rand1 <- rownames(rna.dat_consensus)[sample(1:nrow(rna.dat_consensus), length(gs.IFNG.GS))]
# gs.rand2 <- rownames(rna.dat_consensus)[sample(1:nrow(rna.dat_consensus), length(gs.ISG.RS))]

# ### Make RNA metagenes ###

# df.meta_list <- lapply(1:length(cancer_types), function(x) {
# 	cancer_type <- cancer_types[x]
# 	df <- make_ISG_metagenes(cancer_type=cancer_type, gs1=gs.rand1, gs2=gs.rand2, rna.dat=rna.dat_consensus, cibersort.dat=cibersort_consensus, rna.meta=rna.dat_meta)
# 	return(df)
# })

# ### Fit beta regression model for all cancers ###

# fits_list <- lapply(1:length(cancer_types), function(x) {
# 	print(cancer_types[x])

# 	df <- df.meta_list[[x]]
# 	df$T.cells.CD8 <- df$T.cells.CD8+pseudocount

# 	# RNA #
# 	f1.rna <- betareg(T.cells.CD8 ~ gs1_RNA + gs2_RNA, data = df)
# 	return(f1.rna)
# })

# ### Get model statistics ###

# fit.stats <- lapply(1:length(cancer_types), function(x) {
# 	cancer_type <- cancer_types[x]
# 	   print(cancer_type)

# 	fit <- fits_list[[x]]

# 	gs1_summary <- get_betareg_stats(fit, name=paste0("gs1_RNA"))
# 	gs2_summary <- get_betareg_stats(fit, name=paste0("gs2_RNA"))

# 	gs.summary <- list(gs1_summary, gs2_summary)
# 	df.stats <- data.frame(cancer_type=cancer_type, immune.cell="CD8", gs=c("gs1", "gs2"), est=sapply(gs.summary, function(x) x[1]), se=sapply(gs.summary, function(x) x[2]), lower=sapply(gs.summary, function(x) x[3]), upper=sapply(gs.summary, function(x) x[4]))
# 		return(df.stats)
# })
# fit.stats <- do.call(rbind, fit.stats)

# ### Forest plot ###

# sample_sizes <- sapply(1:length(cancer_types), function(x) nrow(df.meta_list[[x]]))

# # gs1 #
# df.gs1 <- fit.stats[which(as.character(fit.stats$gs) == "gs1"), ]
# df.gs1$name <- paste0(df.gs1$cancer_type, " (n=", sample_sizes, ")")
# df.gs1.sorted <- df.gs1[sort.int(df.gs1$est, decreasing = TRUE, index.return = TRUE)$ix, ]
# df.gs1.sorted$name <- factor(df.gs1.sorted$name, levels=rev(as.character(df.gs1.sorted$name)))

# # gs2 #
# df.gs2 <- fit.stats[which(as.character(fit.stats$gs) == "gs2"), ]
# df.gs2$name <- paste0(df.gs2$cancer_type, " (n=", sample_sizes, ")")
# df.gs2.sorted <- df.gs2[sort.int(df.gs2$est, decreasing = TRUE, index.return = TRUE)$ix, ]
# df.gs2.sorted$name <- factor(df.gs2.sorted$name, levels=rev(as.character(df.gs2.sorted$name)))

# p_coeff.gs1 <- ggplot(data=df.gs1.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
#     geom_pointrange(color = "black", size = 0.4) +
#     geom_hline(yintercept=0, lty=2) +
#     coord_flip() +
#     xlab("Cancer type") +
#     ylab("Effect size (95% CI)") +
#     ylim(min(c(df.gs1.sorted$lower, df.gs2.sorted$lower)) - 0.01, max(c(df.gs1.sorted$upper, df.gs2.sorted$upper)) + 0.01) +
#     ggtitle("gs1") +
#     theme_bw() +
#     theme(axis.ticks = element_line(color = "black"),
#           axis.text = element_text(color = "black", size = 8),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           plot.title = element_text(size = 10))

# p_coeff.gs2 <- ggplot(data=df.gs2.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
#     geom_pointrange(color = "black", size = 0.4) +
#     geom_hline(yintercept=0, lty=2) +
#     coord_flip() +
#     xlab("Cancer type") +
#     ylab("Effect size (95% CI)") +
#     ylim(min(c(df.gs1.sorted$lower, df.gs2.sorted$lower)) - 0.01, max(c(df.gs1.sorted$upper, df.gs2.sorted$upper)) + 0.01) +
#     ggtitle("gs2") +
#     theme_bw() +
#     theme(axis.ticks = element_line(color = "black"),
#           axis.text = element_text(color = "black", size = 8),
#           axis.title.x = element_blank(),
#           axis.title.y = element_blank(),
#           plot.title = element_text(size = 10))

# ### Combine samples ###

# df <- do.call(rbind, df.meta_list)
# df$T.cells.CD8 <- df$T.cells.CD8 + pseudocount

# # RNA #
# fit <- betareg(T.cells.CD8 ~ gs1_RNA + gs2_RNA, data = df)

# gs1_summary <- get_betareg_stats(fit, name=paste0("gs1_RNA"))
# gs2_summary <- get_betareg_stats(fit, name=paste0("gs2_RNA"))
# gs.summary <- list(gs1_summary, gs2_summary)

# df.stats <- data.frame(rep=paste0("rep", i), immune.cell="CD8", gs=c("gs1", "gs2"), est=sapply(gs.summary, function(x) x[1]), lower=sapply(gs.summary, function(x) x[2]), upper=sapply(gs.summary, function(x) x[3]))
# df.stats$gs <- factor(df.stats$gs, levels=c("gs2", "gs1"))

# p_coeff.rand.comb <- ggplot(data=df.stats, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
# 		geom_pointrange(aes(color = gs), size = 0.8) +
# 		geom_hline(yintercept=0, lty = 2) +
# 		scale_color_manual(values = c("black", "black")) +
# 		coord_flip() +
# 		xlab("Cancer type") +
# 		ylab("Effect size (95% CI)") +
# 		ylim(min(c(df.gs1.sorted$lower, df.gs2.sorted$lower)) - 0.01, max(c(df.gs1.sorted$upper, df.gs2.sorted$upper)) + 0.01) +
# 		ggtitle(paste0("ALL Cancers RNA (n=", nrow(df), ")")) +
# 		theme_bw() +
# 		theme(axis.ticks = element_line(color = "black"),
# 	        axis.text = element_text(color = "black", size = 8),
# 	        axis.title.x=element_text(size=10),
# 	        axis.title.y=element_blank(),
# 	        plot.title = element_text(size = 10),
# 	        legend.position = "none")

# dev.off()
# egg::ggarrange(p_coeff.gs1, p_coeff.gs2, p_coeff.rand.comb, ncol = 1, heights = c(1, 1, 0.2))

# # if (plot) pdf(paste0(dir.path, "ATAC_RNA_integration/mRFAR_v7/TCGA_RNA_ATAC/figures/forest/BetaRegression RNA random gene sets.pdf"), height = 6.8, width = 3.5)
# # # p_coeff.ifng.gs / p_coeff.isg.rs / p_coeff.comb
# # egg::ggarrange(p_coeff.gs1, p_coeff.gs2, p_coeff.rand.comb, ncol = 1, heights = c(1, 1, 0.2))
# # if (plot) dev.off()






