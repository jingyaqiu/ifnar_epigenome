### PREP TCGA RNA AND ATAC DATA ###

rm(list=ls())

library(data.table)
library(TCGAutils)
library(tidyverse)
library(GenomicRanges)

########################################################################
### TCGA Pan-Cancer Atlas RNA-seq data (batch-corrected, normalized) ###
########################################################################

# Download from Xena: https://dev.xenabrowser.net/datapages/?dataset=EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

# Load in batch corrected pancancer data (20,531 genes x 11,069 samples)
rna.batchCorrected <- fread("~/Dropbox/Minn/resources/TCGA/RNA/pancancer_merged/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz", sep="\t", header=T, stringsAsFactors=F)
colnames(rna.batchCorrected)[1] <- "Gene"
rna.batchCorrected <- rna.batchCorrected[!duplicated(rna.batchCorrected$Gene), ] # Remove duplicate genes
rna.batchCorrected <- as.matrix(rna.batchCorrected, rownames=1)

# Get caseUUIDs associated with TCGA barcodes
df.barcodeToUUID <- barcodeToUUID(barcodes = TCGAbarcode(colnames(rna.batchCorrected)))
colnames(df.barcodeToUUID) <- c("TCGA_barcode", "caseUUID")
table(TCGAbarcode(colnames(rna.batchCorrected)) %in% df.barcodeToUUID$TCGA_barcode)

# Load in phenotype metadata
metadat.rna <- read.table("~/Dropbox/Minn/resources/TCGA/RNA/pancancer_merged/TCGA_phenotype_denseDataOnlyDownload.tsv.gz", sep="\t", header=T, stringsAsFactors=F)
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

# Write out data to use for downstream analysis
saveRDS(rna.batchCorrected, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena_Reformat.rds")
saveRDS(metadat.rna, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_RNA_Pancancer_Metadata.rds")

################################################################
### TCGA Pan-Cancer Atlas ATAC-seq data (Corces et al. 2018) ###
################################################################

# Download from NCI GDC: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG

# Load in normalized ATAC-seq insertion counts (562,709 pan-cancer peaks x 796 samples)
pan_cancer_atac_counts_log2norm <- readRDS("~/Dropbox/Minn/resources/TCGA/ATAC/TCGA-ATAC_PanCan_Log2Norm_Counts.rds")
gr.atac <- makeGRangesFromDataFrame(pan_cancer_atac_counts_log2norm[,1:7], keep.extra.columns=T)
atac.dat <- pan_cancer_atac_counts_log2norm[ ,8:ncol(pan_cancer_atac_counts_log2norm)]
colnames(atac.dat) <- gsub("_", "-", colnames(atac.dat))

# Load in TCGA ATAC-seq metadata, re-format
metadat.atac <- read.table("~/Dropbox/Minn/resources/TCGA/ATAC/TCGA_identifier_mapping.txt", header=T, sep="\t", stringsAsFactors=F)
colnames(metadat.atac) <- c("ID", "stanfordUUID", "aliquot_id", "caseUUID", "TCGA_barcode")
metadat.atac <- metadat.atac[match(colnames(atac.dat), metadat.atac$ID), ]
metadat.atac$TCGA_barcode_short <- TCGAbarcode(metadat.atac$TCGA_barcode)
metadat.atac$cancer_type <- sapply(strsplit(metadat.atac$ID, split="-"), function(x) x[[1]])
print(all(colnames(atac.dat) == metadat.atac$ID))
table(metadat.atac$TCGA_barcode_short %in% metadat.rna$TCGA_barcode_short)
table(metadat.atac$caseUUID %in% metadat.rna$caseUUID)

# Write out data to use for downstream analysis
saveRDS(atac.dat, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Pancancer_Log2Norm_Counts.rds")
saveRDS(metadat.atac, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Pancancer_Metadata.rds")

############################
### Paired RNA/ATAC data ###
############################

# Paired samples
int.samples <- intersect(metadat.atac$TCGA_barcode_short, metadat.rna$TCGA_barcode_short)

# Paired samples metadata
metadat.int <- data.frame(
	TCGA_barcode_short = int.samples,
	caseUUID = metadat.rna$caseUUID[match(int.samples, metadat.rna$TCGA_barcode_short)],
	cancer_type = metadat.rna$cancer_type[match(int.samples, metadat.rna$TCGA_barcode_short)],
	RNA_ID = sapply(int.samples, function(x) paste(metadat.rna$ID[metadat.rna$TCGA_barcode_short %in% x], collapse = ", ")),
	ATAC_ID = sapply(int.samples, function(x) paste(metadat.atac$ID[metadat.atac$TCGA_barcode_short %in% x], collapse = ", ")), stringsAsFactors=F)

### Merge replicates ###

int.sample_RNA <- strsplit(metadat.int$RNA_ID, split=", ")
names(int.sample_RNA) <- metadat.int$TCGA_barcode_short
rna.dat.int <- merge_replicates(
	dat = rna.batchCorrected,
	ids = int.sample_RNA)

int.sample_ATAC <- strsplit(metadat.int$ATAC_ID, split=", ")
names(int.sample_ATAC) <- metadat.int$TCGA_barcode_short
atac.dat.int <- merge_replicates(
	dat = atac.dat,
	ids = int.sample_ATAC)

# Check
print(all(colnames(rna.dat.int) == metadat.int$TCGA_barcode_short))
print(all(colnames(atac.dat.int) == metadat.int$TCGA_barcode_short))

# Cytolytic score
all(c("PRF1", "GZMA") %in% rownames(rna.dat.int))
metadat.int$cytolytic_score <- as.numeric(sqrt(rna.dat.int["PRF1", ] * rna.dat.int["GZMA", ])) # cytolytic score

# Write out data to use for downstream analysis
saveRDS(rna.dat.int, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_RNA_Paired_Samples_Merged_Data.rds")
saveRDS(atac.dat.int, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Paired_Samples_Merged_Data.rds")
saveRDS(metadat.int, file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_Paired_Samples_Metadata.rds")

###############################################
### WRITE OUT CD8 IMMUNE PEAK-TO-GENE LINKS ###
###############################################

# Supplementary Data File S9
# Download from NCI GDC: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG

# Infiltrating immune cells can contribute to ATAC-seq data through actions on tumor cells and through increased chromatin accessbility at known immune-specific regulatory elements
# Compare accessibility of each linked peak in immune cell types vs. bulk cancer samples - peaks more accessible in immune cells might be generated from immune cells associated with tumor tissue
immune.links <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Corces2018_peak-to-gene_links_immune_response.txt", sep="\t", stringsAsFactors=F, header=T)
gr.immune.links <- makeGRangesFromDataFrame(immune.links, keep.extra.columns=T)
gr.immune.links_CD8 <- gr.immune.links[which(gr.immune.links$enriched_cell_type == "CD8")]
# gr.immune.links_CD8 <- gr.immune.links[gr.immune.links$L2FC_CD8 > 0]

# Write out immune links
for(i in 1:length(cancer_types)) {
	cancer_type <- cancer_types[i]
	print(cancer_type)
	if(cancer_type %in% c("CESC", "CHOL", "GBMx")) next
	metadat <- metadat.int[metadat.int$cancer_type == cancer_type, ]

	atac <- atac.dat.int[match(gr.immune.links_CD8$Peak_Name, rownames(atac.dat.int)), match(metadat$TCGA_barcode_short, colnames(atac.dat.int))]
	cor.cytolytic_score <- apply(atac, 1, cor, metadat$cytolytic_score)
	gr.immune <- gr.immune.links_CD8[gr.immune.links_CD8$Peak_Name %in% names(cor.cytolytic_score)[cor.cytolytic_score > 0.2]]
	gr.immune <- unique(gr.immune)
	print(paste0("Number of CD8 peaks: ", length(gr.immune)))

	# Write out immune peaks
	make_HOMER_bed(gr.immune,
		file = paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Immune_CD8_Cancer_Type-specific/Immune_CD8_", cancer_type, ".bed"))
	saveRDS(gr.immune, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Immune_CD8_Cancer_Type-specific/Immune_CD8_", cancer_type, ".rds"))
}

# All CD8 peaks
cd8_peaks <- unique(unlist(sapply(Immune_CD8_peaks_list, function(x) x$Peak_Name)))
gr.cd8 <- gr.immune.links_CD8[match(cd8_peaks, gr.immune.links_CD8$Peak_Name)]
gr.cd8 <- gr.cd8[sort.int(gr.cd8$L2FC_CD8, decreasing=T, index.return=T)$ix]
make_HOMER_bed(gr.cd8, 
	file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Motif Enrichment/Immune_CD8_ALL_peaks.bed")

#################
### Functions ###
#################

merge_replicates <- function(dat, ids) {
	dat_merged <- lapply(1:length(ids), function(x) {
		dat_sub <- as.matrix(dat[ ,which(colnames(dat) %in% ids[[x]])])
		if(ncol(dat_sub) > 1) {
			return(as.numeric(rowMeans(dat_sub)))
		} else {
			return(dat_sub[,1])
		}
	})
	dat_merged <- do.call(cbind, dat_merged)
	colnames(dat_merged) <- names(ids)
	rownames(dat_merged) <- rownames(dat)
	return(dat_merged)
}
