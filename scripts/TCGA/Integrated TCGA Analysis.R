### Integrated TCGA analysis, plot figures (Figure 1) ###

library(tidyverse)
library(broom)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(ggrepel)

######################
### Load TCGA data ###
######################

# Pancancer RNA
rna.batchCorrected <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena_Reformat.rds")
metadat.rna <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_RNA_Pancancer_Metadata.rds")

# Pancancer ATAC
atac.dat <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Pancancer_Log2Norm_Counts.rds")
metadat.atac <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Pancancer_Metadata.rds")

# Paired RNA/ATAC
rna.dat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_RNA_Paired_Samples_Merged_Data.rds")
atac.dat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Paired_Samples_Merged_Data.rds")
metadat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_Paired_Samples_Metadata.rds")

# Checks
print(all(metadat.rna$ID == colnames(rna.batchCorrected)))
print(all(metadat.atac$ID == colnames(atac.dat)))

############################
### Load resources files ###
############################

### CANCER TYPES ###

cancer_types <- c("SKCM", "ACCx", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBMx", "HNSC", "KIRC", "KIRP", "LGGx", "LIHC", "LUAD", "LUSC", "MESO", "PCPG", "PRAD", "STAD", "TGCT", "THCA", "UCEC")

# Remove cancer types with fewer than <8 samples, remove brain cancers with low immune infiltration
cancer_types.retain <- cancer_types[!cancer_types %in% c("CESC", "CHOL", "MESO", "GBMx", "LGGx")] 

# Keep cancer types with >15 paired samples
cancer_types.paired <- c("BRCA", "COAD", "KIRP", "PRAD", "LUAD",
						 "LIHC", "STAD", "LUSC", "KIRC", "ESCA")

### ISG GENE SETS ###

gs.IFNG <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Resource Files/IFNG.txt", stringsAsFactors=F)
gs.ISG.RS <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Resource Files/ISG.RS.txt", stringsAsFactors=F)
gs.IFN.I <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Resource Files/IFN.I.txt", stringsAsFactors=F)
gs.IFNG <- gs.IFNG$V1
gs.ISG.RS <- gs.ISG.RS$V1
gs.IFN.I <- gs.IFN.I$V1
gs.IFNG <- setdiff(gs.IFNG, gs.ISG.RS)

## ISG.RS mRF peaks ###

name <- "ISG.RS"
ISG.RS_peaks_list <- lapply(1:length(cancer_types.paired), function(x) {
	cancer_type <- cancer_types.paired[x]
	gr <- readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRF_peaks_Significant_", cancer_type, ".rds"))
	return(gr$name)
})
names(ISG.RS_peaks_list) <- cancer_types.paired
sapply(ISG.RS_peaks_list, function(x) length(x))

### Immune CD8 peak-to-gene links ###

Immune_CD8_peaks_list <- lapply(1:length(cancer_types.retain), function(x) readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Immune_CD8_Cancer_Type-specific/Immune_CD8_", cancer_types.retain[x], ".rds")))
names(Immune_CD8_peaks_list) <- cancer_types.retain
sapply(Immune_CD8_peaks_list, function(x) length(x))

######################################################
######################################################

### MAIN FIGURES ###

############################################
### Fig 1B: IFNG.GS vs ISG.RS expression ###
############################################

# ISG.RS AND IFNG.GS EXPRESSION ARE CORRELATED, BUT IFNG.GS > ISG.RS RATIO TRENDS TOWARDS HIGHER CD8 ACTIVITY (SIMILAR TO HUGO ANALYSIS IN BENCI 2019 PAPER)

metadat.rna$RNA_ISG.RS <- scale(colMeans(na.omit(rna.batchCorrected[match(gs.ISG.RS, rownames(rna.batchCorrected)), ])))
metadat.rna$RNA_IFNG.GS <- scale(colMeans(na.omit(rna.batchCorrected[na.omit(match(gs.IFNG, rownames(rna.batchCorrected))), ])))

# r = 0.82, p < 2.2e-16
cor.test(metadat.rna$RNA_IFNG.GS, metadat.rna$RNA_ISG.RS)

r <- round(cor(metadat.rna$RNA_IFNG.GS, metadat.rna$RNA_ISG.RS), 2)
p <- ggplot(metadat.rna, aes(x = RNA_IFNG.GS, y = RNA_ISG.RS, color = cytolytic_score)) +
	geom_point(shape=1, size=0.8, stroke=0.2) +
	geom_abline(intercept = 0, slope = 1, size = 0.5, color = "firebrick") +
	scale_color_gradientn(name = "Cytolytic \n Score", colors = rev(brewer.pal(9, "RdYlBu"))) +
	annotate("text", x=-1.8,  y=3, label=paste0("r = ", r), size=5) +
	labs(x = "IFNG.GS expression", y = "ISG.RS expression") +
	theme_bw() +
	theme(axis.text.y = element_text(size=12, color="black"),
		axis.text.x = element_text(size=12, color="black"),
		text=element_text(size=14),
		aspect.ratio=1)
ggsave(p, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1B IFNG.GS vs ISG.RS.pdf", width=5, height=5)

# Write out plot data
# dat <- p$data[ ,c("IFNG.GS", "ISG.RS", "cytolytic_score")]
# write.table(dat, file=paste0(dir.path, "ifnar_epigenome/Figures Data/Source/Figure 1B.csv"), sep=",", quote=F, col.names=T, row.names=T)

########################################################
### Fig 1C: Cancer/Immune ISG effect on CD8 activity ###
########################################################

# CANCER/IMMUNE ISGs HAVE OPPOSING EFFECTS ON CD8 CYTOLYTIC ACTIVITY (ALL TCGA TUMORS)

assay <- "RNA"
use.metric <- "cytolytic_score"

# Regression model
fits_list <- lapply(1:length(cancer_types), function(x) {
	print(cancer_types[x])
	mat.plot <- metadat.rna[metadat.rna$cancer_type == cancer_types[x], ]
	df.stats <- fit_lm_regression(
		mat.plot = mat.plot, 
		use.metric = use.metric, 
		predictors = c(paste0(assay, "_ISG.RS"), paste0(assay, "_IFNG.GS")),
		name = cancer_types[x])
	return(df.stats)
})
fit.stats <- do.call(rbind, fits_list)

# Forest plot
p.isg <- visualize_forest_plot(
	fit.stats = fit.stats, 
	metadat = metadat.rna,
	predictors = c("RNA_IFNG.GS", "RNA_ISG.RS"),
	colors = c("navyblue", "firebrick"),
	names = c("Immune ISGs", "Cancer ISGs"),
	plot.cancers = cancer_types.retain)
ggsave(p.isg, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1C RNA pancancer ISGs forest plot.pdf", width=4.5, height=7)

# Write out plot data
# dat <- fit.stats[ ,c("est", "upper", "lower", "cancer_type", "gs")]
# write.table(dat, file="~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 1C.csv"), sep=",", quote=F, col.names=T, row.names=F)

########################################################
### Fig. 1D: TSS distance vs. Regression coefficient ###
########################################################

mRF_peaks <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/ISG.RS_mRF_peaks_ALL_KIRP.txt", sep="\t", header=T, stringsAsFactors=F)
mRF_peaks <- mRF_peaks[mRF_peaks$Peak_Name != "(Intercept)", ]
gr.mRF_peaks <- gr.atac[match(mRF_peaks$Peak_Name, gr.atac$name)]
mRF_peaks$DistanceFromTSS <- start(resize(gr.mRF_peaks, fix = "center", width=1)) - mRF_peaks$TSS
p <- ggplot(mRF_peaks, aes(x = DistanceFromTSS, y = abs(Peak_GLMNET))) +
	geom_point(size = 1, shape = 1) +
	labs(x = "Distance From TSS", y = "|Regression coefficient|") +
	# xlim(c(-150000, 150000)) +
	theme_bw() +
	theme(text=element_text(size=14, color="black"))
ggsave(p, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1D TSS Distance v Regression Coefficient.pdf", width=5, height=2)

# Write out plot data
# dat <- mRF_peaks[ ,c("DistanceFromTSS", "Peak_GLMNET")]
# dat$Peak_GLMNET <- abs(dat$Peak_GLMNET)
# write.table(dat, file=paste0(dir.path, "ifnar_epigenome/Figures Data/Source/Figure 1D.csv"), sep=",", quote=F, col.names=T, row.names=F)

#################################################
### Fig. 1E: Rank genes by variance explained ###
#################################################

# Does ATAC explain ISG.RS gene expression?
# Plot ranked genes for KIRP

mRFobj <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/ISG.RS_mRFobj_KIRP.rds")
Rsq <- sapply(mRFobj, function(x) x$Rsq_GLMNET)
df.Rsq <- data.frame(
	Gene = gs.ISG.RS, 
	Rsq = as.numeric(Rsq[match(gs.ISG.RS, names(Rsq))]))
df.Rsq <- df.Rsq[match(gs.ISG.RS, df.Rsq$Gene), ]
# df.Rsq[is.na(df.Rsq)] <- 0
mat.plot <- na.omit(df.Rsq)

mat.plot <- mat.plot[sort.int(mat.plot$Rsq, decreasing=T, index.return=T)$ix, ]
mat.plot$rank <- 1:nrow(mat.plot)
mat.plot$label <- ifelse(mat.plot$Gene %in% c("OAS1", "BST2", "IFI44", "MX1", "LGALS3BP", "OASL", "STAT1", "ISG15"), "firebrick", "black")
mat.plot$size <- ifelse(mat.plot$Gene %in% c("OAS1", "BST2", "IFI44", "MX1", "LGALS3BP", "OASL", "STAT1", "ISG15"), 3, 1)

fig <- ggplot(mat.plot, aes(x=rank, Rsq)) +
	geom_point(size = mat.plot$size, color = mat.plot$label) +
    geom_label_repel(data=mat.plot[mat.plot$label == "firebrick", ], aes(label=Gene), point.padding=0.2, segment.size=0.5, segment.color="firebrick", col="firebrick", fontface="bold", size=3, box.padding=0.75) +
    ylim(c(0,1)) +
    xlab("Rank") +
    ylab("Variance explained (R2)") +
    theme_bw() +
    theme(axis.ticks = element_line(color = "black"),
          axis.text = element_text(color = "black", size = 12),
          axis.title.x = element_text(size=14),
          axis.title.y = element_text(size=14),
          plot.title = element_text(size = 10),
          aspect.ratio = 1)
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1E ISG.RS variance explained ranked.pdf", width=4, height=3)

# Write out plot data
# dat <- mat.plot[ ,c("Gene", "Rsq", "rank")]
# write.table(dat, file=paste0(dir.path, "ifnar_epigenome/Figures Data/Source/Figure 1E.csv"), sep=",", quote=F, col.names=T, row.names=F)

################################################################################
### Fig. 1F: Pan-cancer ISG RNA/ATAC effect on CD8 activity (paired samples) ###
################################################################################

# Get ATAC ISG metafeatures (cancer-specific)
df.isg_peaks <- lapply(1:length(cancer_types.paired), function(x) {
	cancer_type <- cancer_types.paired[x]
	metadat <- metadat.int[metadat.int$cancer_type == cancer_types.paired[x], ]
	atac.dat_ISG.RS <- atac.dat.int[match(ISG.RS_peaks_list[[cancer_types.paired[x]]], rownames(atac.dat.int)), match(metadat$TCGA_barcode_short, colnames(atac.dat.int))]
	# atac.dat_ISG.RS_zmad <- t(apply(atac.dat_ISG.RS, 1, calculate_zmad_score))

	atac.dat_IFNG.GS <- atac.dat.int[match(Immune_CD8_peaks_list[[cancer_types.paired[x]]]$Peak_Name, rownames(atac.dat.int)), match(metadat$TCGA_barcode_short, colnames(atac.dat.int))]
	# atac.dat_IFNG.GS_zmad <- t(apply(atac.dat_IFNG.GS, 1, calculate_zmad_score))
	df <- data.frame(
		TCGA_barcode_short = metadat$TCGA_barcode_short,
		caseUUID = metadat$caseUUID, 
		ATAC_ISG.RS = colMeans(atac.dat_ISG.RS), 
		ATAC_IFNG.GS = colMeans(atac.dat_IFNG.GS),
		cancer_type = cancer_types.paired[x])
	return(df)
})
df.isg_peaks <- do.call(rbind, df.isg_peaks)

# Add RNA and ATAC ISG metadat

rna.dat_ISG.RS <- t(apply(rna.dat.int[match(gs.ISG.RS, rownames(rna.dat.int)), ], 1, calculate_zmad_score))
rna.dat_IFNG.GS <- t(apply(rna.dat.int[na.omit(match(gs.IFNG, rownames(rna.dat.int))), ], 1, calculate_zmad_score))
metadat.int$RNA_ISG.RS <- colMeans(rna.dat_ISG.RS)
metadat.int$RNA_IFNG.GS <- colMeans(rna.dat_IFNG.GS)
metadat.int$ATAC_ISG.RS <- df.isg_peaks$ATAC_ISG.RS[match(metadat.int$TCGA_barcode_short, df.isg_peaks$TCGA_barcode_short)]
metadat.int$ATAC_IFNG.GS <- df.isg_peaks$ATAC_IFNG.GS[match(metadat.int$TCGA_barcode_short, df.isg_peaks$TCGA_barcode_short)]

### PAN-CANCER PLOTS ###

use.metric <- "cytolytic_score"
plot.cancers <- cancer_types.paired

set.seed(1)
idx.subsample <- unlist(sapply(1:length(plot.cancers), function(x) {
	idx <- which(metadat.int$cancer_type == plot.cancers[x])
	if(length(idx) > 15) idx <- idx[sample(1:length(idx), 15)]
}))
metadat.subsample <- metadat.int[idx.subsample, ]
print(table(metadat.subsample$cancer_type))

assay <- "RNA"
predictors <- c(paste0(assay, "_IFNG.GS"), paste0(assay, "_ISG.RS"))

fit.stats_combined <- fit_lm_regression(
	mat.plot = metadat.subsample, 
	use.metric = use.metric, 
	predictors = predictors,
	name = "ALL")
fit.stats_combined$gs <- factor(fit.stats_combined$gs, levels=rev(predictors))

p1 <- ggplot(data=fit.stats_combined, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
	geom_pointrange(aes(color = gs), fatten = 1, size = 1) +
	geom_hline(yintercept=0, lty = 2) +
	scale_color_manual(values = rev(c("navyblue", "firebrick"))) +
	coord_flip() +
	labs(x="Cancer type", y="Standardized regression coefficient (95% CI)", title=paste0("Pan-cancer ", assay, " (paired tumors)")) +
	# ggtitle(paste0("Pan-cancer ", assay, " (n=", nrow(metadat.subsample), ")")) +
	# ylim(min(fit.stats_combined$lower) - 0.05, max(fit.stats_combined$upper) + 0.01) +
	ylim(-1.1, 1.5) +
	theme_bw() +
	theme(axis.ticks = element_line(color = "black"),
		axis.text = element_text(color = "black", size=12),
		axis.title.x=element_text(size=13),
		axis.title.y=element_blank(),
		plot.title = element_text(size = 13),
		legend.position = "none",
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())

assay <- "ATAC"
predictors <- c(paste0(assay, "_IFNG.GS"), paste0(assay, "_ISG.RS"))

fit.stats_combined <- fit_lm_regression(
	mat.plot = metadat.subsample, 
	use.metric = use.metric, 
	predictors = predictors,
	name = "ALL")
fit.stats_combined$gs <- factor(fit.stats_combined$gs, levels=rev(predictors))
p2 <- ggplot(data=fit.stats_combined, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
	geom_pointrange(aes(color = gs), fatten = 1, size = 1) +
	geom_hline(yintercept=0, lty = 2) +
	scale_color_manual(values = rev(c("navyblue", "firebrick"))) +
	coord_flip() +
	labs(x="Cancer type", y="Standardized regression coefficient (95% CI)", title=paste0("Pan-cancer ", assay, " (paired tumors)")) +
	ylim(-1.1, 1.5) +
	theme_bw() +
	theme(axis.ticks = element_line(color = "black"),
		axis.text = element_text(color = "black", size=12),
		axis.title.x=element_text(size=13),
		axis.title.y=element_blank(),
		plot.title = element_text(size = 13),
		legend.position = "none",
		panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank())
fig <- p1 / p2
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1F Pancancer forest plot RNA ATAC.pdf", width=6, height=4)

# # Write out plot data
# dat <- fit.stats_combined[ ,c("est", "lower", "upper", "gs")]
# write.table(dat, file=paste0(dir.path, "ifnar_epigenome/Figures Data/Source/Figure 1F.csv"), sep=",", quote=F, col.names=T, row.names=F)

##########################################
### FIG. 1G: MOTIF ENRICHMENT ANALYSIS ###
##########################################

mat.combined <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Motif Enrichment/HOMER.enriched_cancer_v_immune_peaks.txt", sep="\t", stringsAsFactors=F, header=T)
sig <- 1e-10
mat.combined <- mat.combined[mat.combined$pval < sig, ]
mat.combined$ID <- factor(mat.combined$ID, levels=unique(mat.combined$ID))
mat.combined$Label <- factor(mat.combined$Label, levels=rev(unique(mat.combined$Label)))

fig <- ggplot(mat.combined, aes(x=Label, y=ID, color=-log10(pval), size=-log10(pval))) +
    geom_point(alpha=0.8) +
    scale_color_gradientn(colors=brewer.pal(9, "Greens")[3:9]) +
    ylab("Motif Archetype") +
    xlab("") +
    labs(size="") +
    guides(size = FALSE) +
    scale_size(range = c(2,8)) +
    coord_flip() +
    theme_bw(base_size=12) +
    theme(axis.text.y = element_text(size=10, color="black"),
        # axis.text.x = element_text(size=10, color="black"),
        axis.text.x = element_text(size=10, angle=30, hjust = 1, color="black"),
        axis.ticks = element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave(fig, file="~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 1/Fig 1G Cancer Immune ISG Motif Enrichment.pdf", width=8, height=3)

# # Write out plot data
# dat <- mat.combined[ ,c("Label", "ID", "pval")]
# write.table(dat, file=paste0(dir.path, "ifnar_epigenome/Figures Data/Source/Figure 1G.csv"), sep=",", quote=F, col.names=T, row.names=F)

######################################################
######################################################

### FUNCTIONS ###

# Fit lm model using ISG gene sets as predictors and CD8 activity as dependent variable
fit_lm_regression <- function(
	mat.plot,
	use.metric = "cytolytic_score",
	predictors = c("RNA_ISG.RS", "RNA_IFNG.GS"),
	name=NA) {
	
	# Y variable is CD8 infiltration (CIBERSORTx T.cells.CD8 or cytolytic_score)
	mat.plot$Score <- mat.plot[,use.metric]
	# if(use.metric == "cytolytic_score") mat.plot$Score <- (mat.plot$Score - min(mat.plot$Score)) / (max(mat.plot$Score) - min(mat.plot$Score))

	# Fit lm model
	dat <- mat.plot[ ,c(predictors, "Score")]
	dat <- data.frame(scale(dat)) # Get standardized regression coefficients
	fit <- lm(formula = paste0("Score ~ ", paste(predictors, collapse=" + ")), data=dat)

	stats_summary <- lapply(1:length(predictors), function(x) get_lm_stats(fit, name=predictors[x]))
	df.stats <- do.call(rbind, stats_summary)
	df.stats$cancer_type <- name
	df.stats$gs <- predictors
	return(df.stats)
}

# Extract lm stats from fit
get_lm_stats <- function(fit, name) {
	summary <- summary(fit)
	coefs <- summary$coefficients
	confint <- confint(fit, level=0.95)

	model_summary <- glance(fit)
	rsq <- model_summary$adj.r.squared
	pval <- as.numeric(model_summary$p.value)
	df.stats <- data.frame(coefs, confint, Rsq = rsq, model_pval = pval)
	colnames(df.stats) <- c("est", "se", "t.value", "beta_pval", "lower", "upper", "Rsq", "model_pval")
	return(df.stats[name,])
}

# Plot regression coefficients for lm regression in forest plot
visualize_forest_plot <- function(
	fit.stats, 
	metadat, 
	predictors = c("RNA_IFNG.GS", "RNA_ISG.RS"),
	colors = c("navyblue", "firebrick"),
	names = c("Immune ISGs", "Cancer ISGs"),
	plot.cancers) {

	metadat <- metadat[metadat$cancer_type %in% plot.cancers, ]

	# Sample size for each cancer
	sample_sizes <- table(metadat$cancer_type)[match(plot.cancers, names(table(metadat$cancer_type)))]

	p_list <- lapply(1:length(predictors), function(x) {
		df <- fit.stats[which(as.character(fit.stats$gs) == predictors[x]), ]
		df <- df[match(plot.cancers, df$cancer_type), ]
		df$name <- paste0(df$cancer_type, " (n=", sample_sizes, ")")
		df.sorted <- df[sort.int(df$est, decreasing = TRUE, index.return = TRUE)$ix, ]
		df.sorted$name <- factor(df.sorted$name, levels=rev(as.character(df.sorted$name)))

		p <- ggplot(data=df.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
			geom_pointrange(color = colors[x], fatten = 1, size = 1) +
			ggtitle(names[x]) +
			geom_hline(yintercept=0, lty=2) +
			coord_flip() +
			xlab("Cancer type") +
			ylab("Beta coefficient (95% CI)") +
			ylim(min(fit.stats$lower) - 0.05, max(fit.stats$upper) + 0.01) +
			theme_bw() +
			theme(axis.ticks = element_line(color = "black"),
				axis.text = element_text(color = "black", size = 7.5),
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				plot.title = element_text(size = 10),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank())
		return(p)
	})
	names(p_list) <- predictors

	# All samples combined
	set.seed(42)
	idx.subsample <- unlist(sapply(1:length(plot.cancers), function(x) {
		idx <- which(metadat$cancer_type == plot.cancers[x])
		if(length(idx) > min(table(metadat$cancer_type))) idx <- idx[sample(1:length(idx), min(table(metadat$cancer_type)))]
	}))
	metadat.subsample <- metadat[idx.subsample, ]
	print(table(metadat.subsample$cancer_type))

	fit.stats_combined <- fit_lm_regression(
		mat.plot = metadat.subsample, 
		use.metric = use.metric, 
		predictors = predictors,
		name = "ALL")
	fit.stats_combined$gs <- factor(fit.stats_combined$gs, levels=rev(predictors))

	p_coeff.comb <- ggplot(data=fit.stats_combined, aes(x=gs, y=est, ymin=lower, ymax=upper)) +
		geom_pointrange(aes(color = gs), fatten = 1, size = 1) +
		geom_hline(yintercept=0, lty = 2) +
		scale_color_manual(values = rev(colors)) +
		coord_flip() +
		labs(x="Cancer type", y="Standardized regression coefficient (95% CI)", title=paste0("Pan-cancer ", assay)) +
		# ylim(min(fit.stats_combined$lower) - 0.05, max(fit.stats_combined$upper) + 0.01) +
		ylim(min(fit.stats$lower) - 0.05, max(fit.stats$upper) + 0.01) +
		theme_bw() +
		theme(axis.ticks = element_line(color = "black"),
	        axis.text = element_text(color = "black", size = 7.5),
	        axis.title.x=element_text(size=10),
	        axis.title.y=element_blank(),
	        plot.title = element_text(size = 10),
	        legend.position = "none",
	        panel.grid.major = element_blank(), 
	        panel.grid.minor = element_blank())
	if(length(predictors) == 1) {
		fig <- p_list[[1]] / p_coeff.comb + plot_layout(heights=c(1,0.25))
	} else if(length(predictors) == 2) {
		fig <- p_list[[predictors[1]]] / p_list[[predictors[2]]] / p_coeff.comb + plot_layout(heights=c(1,1,0.25))
	} else if(length(predictors) == 3) {
		fig <- p_list[[predictors[1]]] / p_list[[predictors[2]]] / p_list[[predictors[3]]] / p_coeff.comb + plot_layout(heights=c(1,1,1,0.25))
	}
	return(fig)
}

# Calculate mean absolute deviation modified Z-score (for RNA/ATAC ISG scores)
calculate_zmad_score <- function(x) {
	mad <- mad(x)
	Zmad <- (x - median(x)) / (1.4826 * mad)
	return(Zmad)
}


