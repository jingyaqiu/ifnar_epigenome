### Link cis-regulatory elements to genes with mRF ###

library(caret)
library(GenomicRanges)

make_HOMER_bed <- function(gr, file=NULL) {
	bed <- data.frame(gr)[,1:3]
	bed$ID <- paste(bed[,1], bed[,2], bed[,3], sep="_")
	bed$col5 <- NA
	bed$strand <- as.character(strand(gr))
	if(!is.null(file)) write.table(bed, file=file, sep="\t", row.names=F, col.names=F, quote=F)
}

write_fasta <- function(gr, file=NULL) {
    fa <- getSeq(Hsapiens, gr)
    names(fa) <- paste0(as.character(gr@seqnames), ":", start(gr), "-", end(gr))
    writeXStringSet(fa, filepath=file, append=FALSE)
}

###############################
### Load in necessary files ###
###############################

### Paired RNA/ATAC ###

rna.dat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_RNA_Paired_Samples_Merged_Data.rds")
atac.dat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_ATAC_Paired_Samples_Merged_Data.rds")
metadat.int <- readRDS("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA_Paired_Samples_Metadata.rds")

### Cancer types ###

cancer_types <- c("SKCM", "ACCx", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "ESCA", "GBMx", "HNSC", "KIRC", "KIRP", "LGGx", "LIHC", "LUAD", "LUSC", "MESO", "PCPG", "PRAD", "STAD", "TGCT", "THCA", "UCEC")

# Remove cancer types with fewer than <8 samples, remove brain cancers with low immune infiltration
cancer_types.retain <- cancer_types[!cancer_types %in% c("CESC", "CHOL", "MESO", "GBMx", "LGGx")] 

# Keep cancer types with >15 paired samples
cancer_types.paired <- c("BRCA", "COAD", "KIRP", "PRAD", "LUAD",
						 "LIHC", "STAD", "LUSC", "KIRC", "ESCA")

### CANCER-SPECIFIC PEAK CALLS ###

# Download from NCI GDC: https://gdc.cancer.gov/about-data/publications/ATACseq-AWG

gr.atac_list <- lapply(1:length(cancer_types), function(x) {
	df <- read.csv(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA-ATAC_Cancer_Type-specific_PeakCalls/", gsub("x", "", cancer_types[x]), "_peakCalls.txt"), stringsAsFactors=F, header=T, sep="\t")
	gr <- makeGRangesFromDataFrame(df, keep.extra.columns = T)
	return(gr)
})
names(gr.atac_list) <- cancer_types

### PAN-CANCER PEAK SET ###

df.atac <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/TCGA-ATAC_PanCancer_PeakSet.txt", stringsAsFactors=F, header=T, sep="\t")
gr.atac <- makeGRangesFromDataFrame(df.atac, keep.extra.columns = T)
gr.atac <- gr.atac[match(rownames(atac.dat.int), gr.atac$name)]
print(all(rownames(atac.dat.int) == gr.atac$name))

### Gene annotations ###

bm <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Resource Files/human_biomart_annotations_2-20-20.txt", sep="\t", header=T, stringsAsFactors=F)
remove_idx <- grep(bm$chromosome_name, pattern="MT|GL|JH|X|Y")
bm <- bm[-remove_idx, ]
bm$strand <- ifelse(bm$strand == 1, "+", "-")
gr.bm <- makeGRangesFromDataFrame(bm, keep.extra.columns=TRUE, start.field="start_position", end.field="end_position")
# strand(gr.bm) <- "*"
gr.tss <- makeGRangesFromDataFrame(data.frame(
	chr = as.character(seqnames(gr.bm)), 
	start = gr.bm$tss, 
	end = gr.bm$tss, 
	Gene = gr.bm$hgnc_symbol), keep.extra.columns = T)

### ISG.RS ###

gs.ISG.RS <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Resource Files/ISG.RS.txt", stringsAsFactors=F)
gs.ISG.RS <- gs.ISG.RS$V1

### Immune peak-to-gene links ###

immune.links <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Corces2018_peak-to-gene_links_immune_response.txt", sep="\t", stringsAsFactors=F, header=T)
gr.immune.links <- makeGRangesFromDataFrame(immune.links, keep.extra.columns=T)
gr.immune.links_CD8 <- gr.immune.links[which(gr.immune.links$enriched_cell_type == "CD8")]

######################################################
### RUN mRF, IDENTIFY CANCER-SPECIFIC ISG.RS PEAKS ###
######################################################

print(all(colnames(rna.dat.int) == colnames(atac.dat.int))) # Check

W <- 92000 # 50000, 250000
gs.list <- c(gs.ISG.RS) # gs.ISG.RS, Random, c(gs.ISG.RS, gs.IFNG)
name <- "ISG.RS" # ISG.RS, Random, ALL_ISGs

for(i in 1:length(cancer_types.paired)) {
	cancer_type <- cancer_types.paired[i]
	print(cancer_type)
	if(cancer_type %in% c("CESC", "CHOL", "GBMx", "MESO", "LGGx")) next
	if(file.exists(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRFobj_", cancer_type, ".rds"))) next
	
	# Cancer-specific peak calls
	gr <- gr.atac_list[[cancer_type]]
	ol <- findOverlaps(gr, gr.atac)
	gr <- gr.atac[unique(to(ol))]

	# Paired RNA/ATAC data
	ids <- metadat.int$TCGA_barcode_short[metadat.int$cancer_type == cancer_type]
	rna.dat <- rna.dat.int[ ,match(ids, colnames(rna.dat.int))]
	atac.dat <- atac.dat.int[ ,match(ids, colnames(atac.dat.int))]
	print(all(colnames(rna.dat) == colnames(atac.dat)))

	# Run mRF
	mRFobj_list <- run_mRF(
		gs.list = gs.list,
		W = W,
		gr.bm = gr.bm,
		rna.dat = rna.dat,
		atac.dat = atac.dat,
		atac.gr = gr,
		seed = 42)
	mRFobj_list <- mRFobj_list[!is.na(mRFobj_list)]
	names(mRFobj_list) <- sapply(mRFobj_list, function(x) x$Gene)
	saveRDS(mRFobj_list, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRFobj_", cancer_type, ".rds"))

	# Write out ALL peaks
	mRFobj_list <- readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRFobj_", cancer_type, ".rds"))
	mRF_peaks <- extract_mRF_peaks(
		mRFobj_list = mRFobj_list, 
		use.metric = "glmnet", 
		th.glmnet = -100, 
		th.Rsq = 0)
	mRF_peaks <- mRF_peaks[mRF_peaks$Peak_Name != "(Intercept)", ]
	write.table(mRF_peaks, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRF_peaks_ALL_", cancer_type, ".txt"), sep="\t", col.names=T, row.names=F)
}

### IDENTIFY mRF LINKED PEAKS ###

name <- "ISG.RS" # ISG.RS, IFNG.GS
max_n <- 500 # 125, 500

for(i in 1:length(cancer_types.paired)) {
	cancer_type <- cancer_types.paired[i]
	print(cancer_type)

	# mRF_peaks <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/TCGA/DAR/", name, "/", name, "_mRFobj_", cancer_type, "2.4c.txt"), sep="\t", header=T)
	mRF_peaks <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRFobj_", cancer_type, ".txt"), sep="\t", header=T)
	mRF_peaks <- mRF_peaks[(mRF_peaks$Peak_GLMNET > 0.5 | mRF_peaks$Peak_Cor > 0.6) & mRF_peaks$Rsq_GLMNET > 0.5, ]
	gr <- gr.atac[gr.atac$name %in% unique(mRF_peaks$Peak_Name), ]
	print(table(!gr$name %in% gr.immune.links$Peak_Name))
	# if(name == "ISG.RS") gr <- gr[!gr$name %in% Immune_CD8_peaks_list[[cancer_type]]$Peak_Name] # Remove known immune-associated loci
	if(name == "ISG.RS") gr <- gr[!gr$name %in% gr.immune.links$Peak_Name]
	print(paste0("Number of ", name, " peaks: ", length(gr)))
	
	# Write out mRF peaks
	make_HOMER_bed(gr, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRF_peaks_", cancer_type, ".bed"))
	saveRDS(gr, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRF_peaks_", cancer_type, ".rds"))
}

name <- "ISG.RS"
ISG.RS_peaks_list <- lapply(1:length(cancer_types.paired), function(x) {
	cancer_type <- cancer_types.paired[x]
	gr <- readRDS(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/mRF/", name, "_mRF_peaks_Significant_", cancer_type, ".rds"))
	return(gr$name)
})
names(ISG.RS_peaks_list) <- cancer_types.paired
sapply(ISG.RS_peaks_list, function(x) length(x))

gr.isg.rs <- gr.atac[match(unique(unlist(ISG.RS_peaks_list)), gr.atac$name)]
make_HOMER_bed(gr.isg.rs, 
	file = "~/Dropbox/Minn/ifnar_epigenome/Processed Data/TCGA/Motif Enrichment/cancer_ISG.RS_ALL_peaks.bed")

#################
### Functions ###
#################

# Run mRF to identify putative gene-to-peak links
run_mRF <- function(
	gs.list,
	W = 92000,
	gr.bm,
	rna.dat,
	atac.dat,
	atac.gr,
	seed=42) {

	mRFobj_list <- lapply(1:length(gs.list), function(x) {
		print(gs.list[x])

		# Gene window
		gs <- gs.list[x]
		gr.bm_gs <- unique(gr.bm[gr.bm$hgnc_symbol %in% gs])
		start(gr.bm_gs) <- start(gr.bm_gs) - W
		end(gr.bm_gs) <- end(gr.bm_gs) + W

		# Peaks in gene window
		ol <- findOverlaps(atac.gr, gr.bm_gs)
		peak.anno <- atac.gr[unique(from(ol))]

		if(length(peak.anno) == 0) {
			print("No peaks")
			return(NA)
		}

		# Paired ATAC/RNA data
		# rna <- rna.dat[gs, match(caseuuids, colnames(rna.dat))]
		# atac <- t(atac.dat[match(peak.anno$name, rownames(atac.dat)), match(caseuuids, colnames(atac.dat))])
		rna <- rna.dat[gs, ]
		atac <- t(atac.dat[match(peak.anno$name, rownames(atac.dat)), ])
		if(length(peak.anno) == 1) {
			atac <- t(atac)
			colnames(atac) <- peak.anno$name
		}
		if(all(rownames(atac) == names(rna)) == FALSE) break

		# Leave-one-out CV is better than K-fold CV for small datasets?
		train_control <- trainControl(
			method = "repeatedcv", 
			number = round(nrow(atac) / 3), 
			repeats = 10)

		# Fit models
		if(length(peak.anno) == 1 | var(rna) == 0) {
			print("Only one peak")
			peak.rf <- var_exp_rf <- rmse_rf <- peak.glmnet <- var_exp_glmnet <- rmse_glmnet <- peak.lm <- var_exp_lm <- rmse_lm <- summary_lm <- NA
		} else {

			# RF
			# set.seed(seed)
			# model_rf <- train(atac, rna, 
			# 	method="rf", 
			# 	trControl=train_control, 
			# 	tuneGrid=expand.grid(mtry=c(2, floor(sqrt(ncol(atac))), ncol(atac)*0.5, ncol(atac))), 
			# 	importance=T)
			# mtry <- model_rf$bestTune$mtry
			# var_exp_rf <- model_rf$results$Rsquared[which(model_rf$results$mtry == mtry)]
			# rmse_rf <- model_rf$results$RMSE[which(model_rf$results$mtry == mtry)]
			# vimps_rf <- varImp(model_rf, scale=F)[[1]]
			# peak.rf <- vimps_rf[,1]
			# names(peak.rf) <- rownames(vimps_rf)

			# Glmnet
			# lasso <- glmnet(x = atac, y = rna)
			set.seed(seed)
			model_glmnet <- train(atac, rna, 
				method="glmnet", 
				trControl=train_control,
				# tuneGrid = expand.grid(alpha = c(0.1,0.5,1), lambda = lasso$lambda),
				tuneGrid = expand.grid(alpha = c(0.1,0.5,1), lambda = seq(0.001, 0.1, by = 0.001)),
				importance=T)
			lambda <- model_glmnet$bestTune$lambda
			alpha <- model_glmnet$bestTune$alpha
			var_exp_glmnet <- model_glmnet$results$Rsquared[model_glmnet$results$alpha == alpha & model_glmnet$results$lambda == lambda]
			rmse_glmnet <- model_glmnet$results$RMSE[model_glmnet$results$alpha == alpha & model_glmnet$results$lambda == lambda]
			# vimps_glmnet <- varImp(model_glmnet, scale=F)[[1]]
			vimps_glmnet <- coef(model_glmnet$finalModel, model_glmnet$finalModel$lambdaOpt)
			peak.glmnet <- vimps_glmnet[,1]
			# names(peak.glmnet) <- rownames(vimps_glmnet)

			# # lm
			# set.seed(seed)
			# model_lm <- train(atac, rna, 
			# 	method="lm", 
			# 	trControl=train_control,
			# 	importance=T)
			# var_exp_lm <- model_lm$results$Rsquared[1]
			# rmse_lm <- model_lm$results$RMSE[1]
			# # vimps_lm <- varImp(model_lm, scale=F)[[1]]
			# # peak.lm <- vimps_lm[,1]
			# # names(peak.lm) <- rownames(vimps_lm)
			# peak.lm <- summary(model_lm$finalModel)$coefficients
			# summary_lm <- summary(model_lm$finalModel)
		}

		# Pearson correlation
		peak.cors <- apply(atac, 2, cor, rna)

		# mRF object
		mRFobj <- list(Gene = gs,
			Window = paste0(as.character(seqnames(gr.bm_gs)), ":", start(gr.bm_gs), "-", end(gr.bm_gs)),
			# VIMPS_RF = peak.rf,
			# Rsq_RF = var_exp_rf,
			# Error_RF = rmse_rf,
			VIMPS_GLMNET = peak.glmnet,
			Rsq_GLMNET = var_exp_glmnet,
			Error_GLMNET = rmse_glmnet,
			# COEFFICIENTS_LM = peak.lm,
			# Rsq_LM = var_exp_lm,
			# Error_LM = rmse_lm,
			# SUMMARY_LM = summary_lm,
			CORRELATIONS = peak.cors)
		return(mRFobj)
	})
	names(mRFobj_list) <- gs.list
	return(mRFobj_list)
}

# Extract peak-to-gene links from mRF object
extract_mRF_peaks <- function(
	mRFobj_list,
	use.metric = "rf", # rf, cor, lm, glmnet
	th.rf = 5,
	th.glmnet = 0.2,
	th.lm = 0.05,
	th.cor = 0.5,
	th.Rsq = 0.3) {

	if(use.metric == "rf") {
		Rsq <- sapply(mRFobj_list, function(x) x$Rsq_RF)
		keep <- na.omit(names(Rsq)[Rsq > th.Rsq])
	} else if(use.metric == "glmnet") {
		Rsq <- sapply(mRFobj_list, function(x) x$Rsq_GLMNET)
		keep <- na.omit(names(Rsq)[Rsq > th.Rsq])
	} else if(use.metric == "lm") {
		Rsq <- sapply(mRFobj_list, function(x) x$Rsq_LM)
		keep <- na.omit(names(Rsq)[Rsq > th.Rsq])
	} else {
		keep <- names(mRFobj_list)
	}
	mRFobj_list <- mRFobj_list[names(mRFobj_list) %in% keep]
	peak2gene_links <- lapply(1:length(mRFobj_list), function(x) {
		mRFobj <- mRFobj_list[[x]]
		# if(use.metric == "rf" & is.na(mRFobj$VIMPS_RF[1])) return(NA)
		# peak.rf <- mRFobj$VIMPS_RF
		peak.glmnet <- mRFobj$VIMPS_GLMNET
		peak.cor <- mRFobj$CORRELATIONS
		# peak.lm <- mRFobj$COEFFICIENTS_LM
		# peak.lm <- peak.lm[2:nrow(peak.lm), ]
		# colnames(peak.lm) <- paste0("LM_", colnames(peak.lm))
		if(use.metric == "rf") {
			peak.select <- na.omit(names(peak.rf)[as.numeric(peak.rf) > th.rf])
		} else if(use.metric == "glmnet") {
			peak.select <- na.omit(names(peak.glmnet)[as.numeric(peak.glmnet) > th.glmnet])
		} else if(use.metric == "lm") {
			peak.select <- na.omit(rownames(peak.lm)[peak.lm[,"LM_Pr(>|t|)"] < th.lm])
			if(length(peak.select) == 1) {
				peak.lm <- data.frame(t(peak.lm[match(peak.select, rownames(peak.lm)), ]))
			} else if(length(peak.select) > 1) {
				peak.lm <- data.frame(peak.lm[match(peak.select, rownames(peak.lm)), ])
			}
		} else if(use.metric == "cor") {
			peak.select <- names(peak.cor)[abs(as.numeric(peak.cor)) > th.cor]
		} else if(use.metric == "both") {
			peak.select <- intersect(na.omit(names(peak.rf)[as.numeric(peak.rf) > th.rf]), na.omit(names(peak.glmnet)[as.numeric(peak.glmnet) > th.glmnet]))
		}

		if(length(peak.select) > 0) {
			peak2gene <- data.frame(
				Peak_Name = peak.select, 
				# Peak_RF = peak.rf[match(peak.select, names(peak.rf))], 
				Peak_GLMNET = peak.glmnet[match(peak.select, names(peak.glmnet))],
				# peak.lm[match(peak.select, rownames(peak.lm)), ],
				Peak_Cor = peak.cor[match(peak.select, names(peak.cor))],
				Gene = rep(mRFobj$Gene, length(peak.select)),
				TSS = rep(unique(gr.bm[which(gr.bm$hgnc_symbol == mRFobj$Gene)])$tss, length(peak.select)),
				# Rsq_RF = rep(mRFobj$Rsq_RF, length(peak.select)),
				Rsq_GLMNET = rep(mRFobj$Rsq_GLMNET, length(peak.select)))
				# Rsq_LM = rep(mRFobj$Rsq_LM, length(peak.select)))
			# colnames(peak2gene) <- c("Peak_Name", "Peak_RF", "Peak_GLMNET", paste0("LM_", colnames(peak.lm)), "Peak_Cor", "Gene", "TSS", "Rsq_RF", "Rsq_GLMNET", "Rsq_LM")
			colnames(peak2gene) <- c("Peak_Name", "Peak_GLMNET", "Peak_Cor", "Gene", "TSS", "Rsq_GLMNET")
		} else {
			peak2gene <- NA
		}
		return(peak2gene)
	})
	peak2gene_links <- do.call(rbind, na.omit(peak2gene_links))
	peak2gene_links <- peak2gene_links[!is.na(peak2gene_links$Peak_Name), ]
	return(peak2gene_links)
}
