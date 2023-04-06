### FIGURE 1B AND 1C ###

rm(list=ls())

library(data.table)
library(TCGAutils)
library(tidyverse)
library(GenomicRanges)
library(GSVA)
library(broom)
library(RColorBrewer)
library(patchwork)
library(ggrepel)

setwd("ifnar_epigenome/TCGA")
source("ifnar_epigenome/TCGA/Scripts/TCGA Analysis Functions.R")

#############################################################################
### Prep TCGA Pan-Cancer Atlas RNA-seq data (batch-corrected, normalized) ###
#############################################################################

# Download Batch-Corrected RNA-seq Count Matrix from Xena
# https://dev.xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443
# Sample metadata: https://dev.xenabrowser.net/datapages/?dataset=TCGA_phenotype_denseDataOnlyDownload.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# Load in batch corrected pancancer data (20,531 genes x 11,069 samples)
rna.batchCorrected <- fread("Resource Files/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz", sep="\t", header=T, stringsAsFactors=F)
colnames(rna.batchCorrected)[1] <- "Gene"
rna.batchCorrected <- rna.batchCorrected[!duplicated(rna.batchCorrected$Gene), ] # Remove duplicate genes
rna.batchCorrected <- as.matrix(rna.batchCorrected, rownames=1)

# Get caseUUIDs associated with TCGA barcodes
df.barcodeToUUID <- barcodeToUUID(barcodes = TCGAbarcode(colnames(rna.batchCorrected)))
colnames(df.barcodeToUUID) <- c("TCGA_barcode", "caseUUID")
table(TCGAbarcode(colnames(rna.batchCorrected)) %in% df.barcodeToUUID$TCGA_barcode)

# Load in sample metadata
metadat.rna <- read.table("Resource Files/TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", header=T, stringsAsFactors=F)
colnames(metadat.rna) <- c("ID", "sample_type_id", "sample_type", "Cancer")
metadat.rna <- metadat.rna[match(colnames(rna.batchCorrected), metadat.rna$ID), ]
metadat.rna$TCGA_barcode_short <- TCGAbarcode(colnames(rna.batchCorrected))
metadat.rna$caseUUID <- df.barcodeToUUID$caseUUID[match(metadat.rna$TCGA_barcode_short, df.barcodeToUUID$TCGA_barcode)]
all(metadat.rna$ID == colnames(rna.batchCorrected))

# Recode cancer types
metadat.rna$cancer_type <- recode(metadat.rna$Cancer,
	'breast invasive carcinoma' = "BRCA",
	'kidney clear cell carcinoma' = "KIRC",
	'lung adenocarcinoma' = "LUAD",
	'thyroid carcinoma' = "THCA",
	'uterine corpus endometrioid carcinoma' = "UCEC",
	'head & neck squamous cell carcinoma' = "HNSC",
	'lung squamous cell carcinoma' = "LUSC",
	'prostate adenocarcinoma' = "PRAD",
	'brain lower grade glioma' = "LGGx",
	'colon adenocarcinoma' = "COAD",
	'skin cutaneous melanoma' = "SKCM",
	'stomach adenocarcinoma' = "STAD",
	'bladder urothelial carcinoma' = "BLCA",
	'liver hepatocellular carcinoma' = "LIHC",
	'kidney papillary cell carcinoma' = "KIRP",
	'cervical & endocervical cancer' = "CESC",
	'ovarian serous cystadenocarcinoma' = "OV",
	'sarcoma' = "SARC",
	'esophageal carcinoma' = "ESCA",
	'pheochromocytoma & paraganglioma' = "PCPG",
	'pancreatic adenocarcinoma' = "PAAD",
	'glioblastoma multiforme' = "GBMx",
	'acute myeloid leukemia' = "LAML",
	'rectum adenocarcinoma' = "READ",
	'testicular germ cell tumor' = "TGCT",
	'thymoma' = "THYM",
	'kidney chromophobe' = "KICH",
	'mesothelioma' = "MESO",
	'uveal melanoma' = "UVM",
	'adrenocortical cancer' = "ACCx",
	'uterine carcinosarcoma' = "UCS",
	'diffuse large B-cell lymphoma' = "DLBC",
	'cholangiocarcinoma' = "CHOL")

# Cytolytic score (sqrt(PRF1 x GZMA))
all(c("PRF1", "GZMA") %in% rownames(rna.batchCorrected))
metadat.rna$cytolytic_score <- as.numeric(sqrt(rna.batchCorrected["PRF1", ] * rna.batchCorrected["GZMA", ])) # cytolytic score

# # Write out data to use for downstream analysis
# saveRDS(rna.batchCorrected, file = "Processed Data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena_Reformat.rds")
# saveRDS(metadat.rna, file = "Processed Data/TCGA_RNA_Pancancer_Metadata.rds")

######################
### Load TCGA data ###
######################

# # Pancancer RNA
# rna.batchCorrected <- readRDS("Processed Data/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena_Reformat.rds")
# metadat.rna <- readRDS("Processed Data/TCGA_RNA_Pancancer_Metadata.rds")
# print(all(metadat.rna$ID == colnames(rna.batchCorrected))) # Check

############################
### Load resources files ###
############################

### CANCER TYPES ###

cancer_types <- c("SKCM", "ACCx", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBMx", "HNSC", "KIRC", "KIRP", "LGGx", "LIHC", "LUAD", "LUSC", "MESO", "PCPG", "PRAD", "STAD", "TGCT", "THCA", "UCEC")

# Remove cancer types with fewer than <8 samples, remove brain cancers with low immune infiltration
cancer_types.retain <- cancer_types[!cancer_types %in% c("CESC", "CHOL", "MESO", "GBMx", "LGGx")] 

# Cancer types with >15 paired samples
cancer_types.paired <- c("BRCA", "COAD", "KIRP", "PRAD", "LUAD",
						 "LIHC", "STAD", "LUSC", "KIRC", "ESCA")

### ISG GENE SETS ###

gs.IFNG <- read.delim("Resource Files/IFNG.txt", stringsAsFactors=F, header=F)[,1]
gs.ISG.RS <- read.delim("Resource Files/ISG.RS.txt", stringsAsFactors=F, header=F)[,1]
gs.IFN.I <- read.delim("Resource Files/IFN.I.txt", stringsAsFactors=F, header=F)[,1]
gs.IFNG <- setdiff(gs.IFNG, gs.ISG.RS)

############################################
### Fig 1B: IFNG.GS vs ISG.RS expression ###
############################################

# ISG.RS AND IFNG.GS EXPRESSION ARE CORRELATED, BUT IFNG.GS > ISG.RS RATIO TRENDS TOWARDS HIGHER CD8 ACTIVITY (SIMILAR TO HUGO ANALYSIS IN BENCI 2019 PAPER)

# ISG summary scores
gsva <- gsva(
	expr = as.matrix(rna.batchCorrected),
	gset.idx.list = list(ISG.RS=gs.ISG.RS, IFNG.GS=gs.IFNG),
	method = "zscore",
	kcdf = "Gaussian")
gsva <- t(scale(t(gsva)))
metadat.rna$RNA_ISG.RS <- gsva[1, match(metadat.rna$ID, colnames(gsva))]
metadat.rna$RNA_IFNG.GS <- gsva[2, match(metadat.rna$ID, colnames(gsva))]

### With matched normals ###

# r = 0.825, p < 2.2e-16
cor.test(metadat.rna$RNA_IFNG.GS, metadat.rna$RNA_ISG.RS)

r <- round(cor(metadat.rna$RNA_IFNG.GS, metadat.rna$RNA_ISG.RS), 2)
p1 <- ggplot(metadat.rna, aes(x = RNA_IFNG.GS, y = RNA_ISG.RS, color = cytolytic_score)) +
	geom_point(shape=1, size=0.8, stroke=0.2) +
	geom_abline(intercept = 0, slope = 1, size = 0.5, color = "firebrick") +
	scale_color_gradientn(name = "Cytolytic \n Score", colors = rev(brewer.pal(9, "RdYlBu"))) +
	labs(x = "IFNG.GS expression", y = "ISG.RS expression", title = "All samples") +
	theme_bw() +
	theme(axis.text.y = element_text(size=12, color="black"),
		axis.text.x = element_text(size=12, color="black"),
		text=element_text(size=14),
		aspect.ratio=1)

### Without matched normals ###

metadat.normal <- metadat.rna[metadat.rna$sample_type == "Solid Tissue Normal", ]
rna.batchCorrected_noNormal <- rna.batchCorrected[ ,!colnames(rna.batchCorrected) %in% metadat.normal$ID]
metadat.rna_noNormal <- metadat.rna[match(colnames(rna.batchCorrected_noNormal), metadat.rna$ID), ]

# r = 0.82, p < 2.2e-16
cor.test(metadat.rna_noNormal$RNA_IFNG.GS, metadat.rna_noNormal$RNA_ISG.RS)

r <- round(cor(metadat.rna_noNormal$RNA_IFNG.GS, metadat.rna_noNormal$RNA_ISG.RS), 2)
p2 <- ggplot(metadat.rna_noNormal, aes(x = RNA_IFNG.GS, y = RNA_ISG.RS, color = cytolytic_score)) +
	geom_point(shape=1, size=0.8, stroke=0.2) +
	geom_abline(intercept = 0, slope = 1, size = 0.5, color = "firebrick") +
	scale_color_gradientn(name = "Cytolytic \n Score", colors = rev(brewer.pal(9, "RdYlBu"))) +
	labs(x = "IFNG.GS expression", y = "ISG.RS expression", title = "Matched normals excluded") +
	theme_bw() +
	theme(axis.text.y = element_text(size=12, color="black"),
		axis.text.x = element_text(size=12, color="black"),
		text=element_text(size=14),
		aspect.ratio=1)

fig <- p1 + p2

########################################################
### Fig 1C: Cancer/Immune ISG effect on CD8 activity ###
########################################################

# CANCER/IMMUNE ISGs HAVE OPPOSING EFFECTS ON CD8 CYTOLYTIC ACTIVITY (ALL TCGA TUMORS)

### No matched normals ###

assay <- "RNA"
use.metric <- "cytolytic_score"

# Regression model
fits_list <- lapply(1:length(cancer_types), function(x) {
	print(cancer_types[x])
	df <- metadat.rna_noNormal[metadat.rna_noNormal$cancer_type == cancer_types[x], ]
	df.stats <- fit_lm_regression(
		df = df, 
		use.metric = use.metric, 
		predictors = c(paste0(assay, "_ISG.RS"), paste0(assay, "_IFNG.GS")),
		cancer_type = cancer_types[x])
	return(df.stats)
})
fit.stats <- do.call(rbind, fits_list)

# Forest plot
p <- visualize_forest_plot(
	fit.stats = fit.stats, 
	metadat = metadat.rna_noNormal,
	predictors = c("RNA_IFNG.GS", "RNA_ISG.RS"),
	colors = c("navyblue", "firebrick"),
	labels = c("Immune ISGs", "Cancer ISGs"),
	plot.cancers = cancer_types.retain)
