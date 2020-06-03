
# loocv_rfSRC <- function(x, y, mtrys=c(2,5,9,15,25), nodes=c(2,5,10)) {
# 	mat <- data.frame(cbind(x, y))

# 	for(i in 1:length(mtrys)) {
# 		for(j in 1:length(nodes)) {
# 			errors <- sapply(1:nrow(mat), function(k) {
# 				train <- mat[-k, ]
# 				test <- mat[k, ]
# 				rf.form <- as.formula(paste0("Multivar(", paste(colnames(y), collapse=", "), ") ~ ."))
# 				train.rfsrc <- rfsrc(rf.form, data=train, importance=TRUE)
# 				pred <- predict(train.rfsrc, newdata=test)
				
# 				errors <- get.mv.error(pred)
# 				err <- errors[which(names(errors) == gene_id)]
# 				return(err)
# 			})
# 			cv.error <- mean(errors)
# 		}
# 	}
# }

mRFAR <- function(gs, W=50000, gr_rna, rna_anno, rna_dat, gr_atac, atac_dat, scaled=F, seed=T, version="v1") {

	mRFARobj <- list()
	for(i in 1:length(gs)) {
		print(i)
		gene_id <- gs[i]
		gene_symb <- as.character(rna_anno$mgi_symbol[which(rna_anno$ensembl_gene_id == gene_id)])
		gene_idx <- which(rna_anno$ensembl_gene_id == gene_id)
		tss <- rna_anno$tss[which(rna_anno$ensembl_gene_id == gene_id)]

		# Define window around TSS
		window.gs <- gr_rna[gene_idx]
		start(window.gs) <- start(window.gs) - W
		end(window.gs) <- end(window.gs) + W

		window <- c(start(window.gs), end(window.gs))

		# Find ATAC peaks in window
		o.atac <- findOverlaps(gr_atac, window.gs)
		o.atac.indices <- from(o.atac)

		# Find genes in window
		o.rna <- findOverlaps(gr_rna, window.gs)
		o.rna.indices <- from(o.rna)

		if(length(o.atac.indices) > 1) {
			if(version == "v1") {
				x <- t(as.matrix(atac_dat[o.atac.indices, ]))
				y <- as.numeric(rna_dat[gene_idx, ])

				if(seed == T) {
					set.seed(42)
				}

				train_control <- trainControl(method="LOOCV", number=3)
				model_rf <- train(x=x, y=y, trControl=train_control, method="rf", importance=TRUE)
				model_glmnet <- train(x=x, y=y, trControl=train_control, method="glmnet")

				# RF stats
				mtry <- model_rf$bestTune$mtry
				var_exp_rf <- model_rf$results$Rsquared[which(model_rf$results$mtry == mtry)]
				rmse_rf <- model_rf$results$RMSE[which(model_rf$results$mtry == mtry)]

				# # glmnet stats
				lambda <- model_glmnet$bestTune$lambda
				alpha <- model_glmnet$bestTune$alpha
				var_exp_glmnet <- model_glmnet$results$Rsquared[intersect(which(model_glmnet$results$alpha == alpha), which(model_glmnet$results$lambda == lambda))]
				rmse_glmnet <- model_glmnet$results$RMSE[intersect(which(model_glmnet$results$alpha == alpha), which(model_glmnet$results$lambda == lambda))]

				# Model-specific VIMPs seem to be better than model-independent loess vimps, based on looking at univariate relationships for i=4 (gs=499res.up)
				# glmnet varImp() is the same as coef(model_glmnet$finalModel, s=lambda)
				if(scaled == F) {
					vimps_rf <- varImp(model_rf, scale=F)[[1]]
					vimps_glmnet <- varImp(model_glmnet, scale=F)[[1]]
				} else {
					vimps_rf <- varImp(model_rf, scale=T)[[1]]
					vimps_glmnet <- varImp(model_glmnet, scale=T)[[1]]
				}
				correlation <- NA
			} else if(version == "v2") {
				x <- t(as.matrix(atac_dat[o.atac.indices, ]))
				y_multi <- t(as.matrix(rna_dat[o.rna.indices, ]))
				mat <- data.frame(cbind(x, y_multi))

				if(seed == T) {
					set.seed(42)
				}

				# Run randomForestSRC
	    		rf.form <- as.formula(paste0("Multivar(", paste(colnames(y_multi), collapse=", "), ") ~ ."))
				model_rfsrc <- rfsrc(rf.form, data=mat, importance=TRUE)
				# model_rfsrc <- rfsrc(rf.form, data=mat, ntree=1000, nodesize=2, importance=TRUE)
				output <- model_rfsrc$regrOutput[[which(names(model_rfsrc$regrOutput) == gene_id)]]
				# plot.variable(model_rfsrc, m.target=gene_id, partial=TRUE, nvar=12)

				errors <- get.mv.error(model_rfsrc)
				vimps <- get.mv.vimp(model_rfsrc)
				rmse_rf <- errors[which(names(errors) == gene_id)]
				vimps_rf <- vimps[,which(colnames(vimps) == gene_id)]
				
				# Extract variance explained
				if(model_rfsrc$family == "regr+") {
					v <- var(model_rfsrc$yvar[,which(colnames(model_rfsrc$yvar) == gene_id)], na.rm=TRUE)
					var_exp_rf <- 100 * (1 - rmse_rf/v)
				} else if(model_rfsrc$family == "regr") {
					v <- var(model_rfsrc$yvar, na.rm=TRUE)
					var_exp_rf <- 100 * (1 - rmse_rf/v)
				}

				var_exp_glmnet <- rmse_glmnet <- vimps_glmnet <- correlation <- NA
			}

		} else if(length(o.atac.indices) == 1) {
			x <- as.matrix(atac_dat[o.atac.indices, ])
			y <- as.matrix(rna_dat[gene_idx, ])
			# x_normalized <- scale(x, center=TRUE, scale=TRUE)

			# Calculate Pearson correlation between ATAC peak counts (normalized?) and gene expression
			correlation <- as.numeric(cor(t(x), t(y)))

			vimps_rf <- vimps_glmnet <- rmse_rf <- rmse_glmnet <- var_exp_rf <- var_exp_glmnet <- NA
		} else {
			print("No ATAC peaks in window")
			vimps_rf <- vimps_glmnet <- rmse_rf <- rmse_glmnet <- var_exp_rf <- var_exp_glmnet <- correlation <- NA
		}
		out <- list(geneid = gene_id, geneSymb=gene_symb, geneIndex=gene_idx, window=window, tss=tss, rnaIndices=o.rna.indices, atacIndices=o.atac.indices, vimps_rf=vimps_rf, vimps_glmnet=vimps_glmnet, error_rf=rmse_rf, error_glmnet=rmse_glmnet, Rsq_rf=var_exp_rf, Rsq_glmnet=var_exp_glmnet, cor=correlation)
		mRFARobj[[i]] <- out
	}
	return(mRFARobj)
}

extractPeaks <- function(mRFARobj, gr_atac) {
	gene_IDs <- sapply(mRFARobj, function(x) x$geneid)
	gene_symbols <- sapply(mRFARobj, function(x) as.character(x$geneSymb))

	vimps_rf <- sapply(mRFARobj, function(x) x$vimps_rf)
	vimps_glmnet <- sapply(mRFARobj, function(x) x$vimps_glmnet)
	atac_idx <- lapply(mRFARobj, function(x) x$atacIndices)
	window <- sapply(mRFARobj, function(x) x$window, simplify=FALSE)
	tss <- sapply(mRFARobj, function(x) x$tss)
	Rsq_rf <- sapply(mRFARobj, function(x) x$Rsq_rf)
	Rsq_glmnet <- sapply(mRFARobj, function(x) x$Rsq_glmnet)

	peaks_out_all <- data.frame(peak_ID=character(), VIMP_rf=numeric(), VIMP_glmnet=numeric(), chr=character(), start=integer(), end=integer(), gene=character(), TSS=integer(), Rsq_rf=numeric(), Rsq_glmnet=numeric())
	for(i in 1:length(mRFARobj)) {
		if(length(atac_idx[[i]]) > 1) {
			temp_df <- data.frame(matrix(nrow=length(vimps_rf[[i]]), ncol=10))

			temp_df[,1] <- as.character(atac_idx[[i]])
			temp_df[,2] <- as.numeric(vimps_rf[[i]])
			temp_df[,3] <- as.numeric(vimps_glmnet[[i]])
			temp_df[,4] <- as.character(seqnames(gr_atac[atac_idx[[i]]]))
			temp_df[,5] <- start(gr_atac[atac_idx[[i]]]) - 1
			temp_df[,6] <- end(gr_atac[atac_idx[[i]]])
			temp_df[,7] <- gene_symbols[i]
			temp_df[,8] <- tss[i]
			temp_df[,9] <- Rsq_rf[i]
			temp_df[,10] <- Rsq_glmnet[i]

			peaks_out_all <- rbind(peaks_out_all, temp_df)
		}	
	}
	colnames(peaks_out_all) <- c("peak_ID", "VIMP_rf", "VIMP_glmnet", "chr", "start", "end", "gene", "TSS", "Rsq_rf", "Rsq_glmnet")
	return(peaks_out_all)
}

# Calculate distance of each peak from TSS 
# Input data frame of peaks generated by extractPeaks
calcTSS_dist <- function(peaks) {
	if(as.numeric(peaks[11]) == 1) {
		dist <- mean(as.numeric(peaks[5]), as.numeric(peaks[6])) - as.numeric(peaks[8])
	} else {
		dist <- -(mean(as.numeric(peaks[5]), as.numeric(peaks[6])) - as.numeric(peaks[8]))
	}
	return(dist)
}

# Plot importance of peak vs distance from TSS
plot_distanceFromTSS_v_VIMP <- function(peaks) {
	peaks_out_all_distanceFromTSS <- apply(peaks, 1, calcTSS_dist)

	df_rf <- data.frame(distance=peaks_out_all_distanceFromTSS, vimp=peaks$VIMP_rf)
	df_glmnet <- data.frame(distance=peaks_out_all_distanceFromTSS, vimp=peaks$VIMP_glmnet)

	fig_rf <- ggplot(df_rf, aes(x=distance, y=vimp)) +
		geom_point(aes(colour=vimp), size=1, na.rm=TRUE) +
		xlab("Distance from TSS") +
		ylab("VIMP score") +
		xlim(-190000, 190000) +
		# ylim(-10, 100)
		# ylim(-10, 40) +
		theme_classic()

	fig_glmnet <- ggplot(df_glmnet, aes(x=distance, y=vimp)) +
		geom_point(aes(colour=vimp), size=1, na.rm=TRUE) +
		xlab("Distance from TSS") +
		ylab("VIMP score") +
		xlim(-190000, 190000) +
		# ylim(-10, 100)
		# ylim(-1, 7) +
		theme_classic()

	return(list(fig_rf, fig_glmnet))
}

# Convert GRanges object to fasta sequences
GRanges_toFasta <- function(gr.peaks) {
	print("Convert to mouse fasta sequence")
	gr.peaks_dna <- getSeq(Mmusculus, gr.peaks)
	names(gr.peaks_dna) <- paste("chromosome:GRCm38:", as.character(gr.peaks@seqnames), ":", start(gr.peaks), ":", end(gr.peaks), sep="")
	return(gr.peaks_dna)
}

# Trim DNAStringSet
# Output trimmed genomic ranges in GRanges object
trimPeakset <- function(gr.peaks, N) {
	gr.peaks_trimmed <- gr.peaks
	gr.peaks_mean <- (start(gr.peaks) + end(gr.peaks))/2
	start(gr.peaks_trimmed) <- gr.peaks_mean - N
	end(gr.peaks_trimmed) <- gr.peaks_mean + N
	return(gr.peaks_trimmed)
}

extract_mRFARobj_peaks <- function(mRFARobj, gr.atac, dir_path, peak_type, trimN=500, rna.anno, version) {

	# Load in mRFAR object
	print(paste(peak_type, "_", length(mRFARobj), sep=""))

	# Create directory if it doesn't exist
	if(!dir.exists(dir_path)) {
		dir.create(dir_path)
		dir.create(paste(dir_path, "/fasta", sep=""))
		dir.create(paste(dir_path, "/peaks", sep=""))
		dir.create(paste(dir_path, "/figures", sep=""))
		dir.create(paste(dir_path, "/bed", sep=""))
	}

	# Extract all peaks
	peaks_out_all <- extractPeaks(mRFARobj, gr.atac)
	peaks_out_all$gene_strand <- rna.anno$strand[match(peaks_out_all$gene, rna.anno$mgi_symbol)]
	write.table(peaks_out_all, file=paste(dir_path, "/peaks/peaks_", peak_type, "_ALL.txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T)
    
	# Convert to GRanges of unique peaks only
	gr.peakset <- gr.atac[as.numeric(unique(peaks_out_all$peak_ID))]
	start(gr.peakset) <- start(gr.peakset) - 1

	# Plot importance of peak vs distance from TSS
	fig_TSS_distance <- plot_distanceFromTSS_v_VIMP(peaks_out_all)
	ggsave(filename=paste("distanceFromTSS_v_VIMP_RF_", peak_type, ".pdf", sep=""), plot=fig_TSS_distance[[1]], device="pdf", path=paste(dir_path, "/figures", sep=""), width=5, height=5)
	ggsave(filename=paste("distanceFromTSS_v_VIMP_glmnet_", peak_type, ".pdf", sep=""), plot=fig_TSS_distance[[2]], device="pdf", path=paste(dir_path, "/figures", sep=""), width=5, height=5)

	# Extract FASTA sequence of peaks
	gr.peakset_dna <- GRanges_toFasta(gr.peakset)
	writeXStringSet(gr.peakset_dna, filepath=paste(dir_path, "/fasta/peaks_", peak_type, "_ALL.fasta", sep=""), append=FALSE)

	if(version == "v1") {
		# Extract only peaks with VIMP score > 0 (impt)
		type <- c("rf", "glmnet")

		cutoff <- 0
		peaks_out_impt_rf <- peaks_out_all[which(peaks_out_all$VIMP_rf > cutoff), ]
		peaks_out_impt_glmnet <- peaks_out_all[which(peaks_out_all$VIMP_glmnet > cutoff), ]
		peaks_out_impt_list <- list(peaks_out_impt_rf, peaks_out_impt_glmnet)
		lapply(1:length(peaks_out_impt_list), function(x) write.table(peaks_out_impt_list[[x]], file=paste(dir_path, "/peaks/peaks_", peak_type, "_impt_", type[x], ".txt", sep=""), sep="\t", quote=F, row.names=F, col.names=T))

		gr.peakset_impt_rf <- gr.atac[unique(as.numeric(peaks_out_impt_rf$peak_ID))]
		gr.peakset_impt_glmnet <- gr.atac[unique(as.numeric(peaks_out_impt_glmnet$peak_ID))]

		gr.peakset_impt_dna_list <- list(GRanges_toFasta(gr.peakset_impt_rf), GRanges_toFasta(gr.peakset_impt_glmnet))
		lapply(1:length(gr.peakset_impt_dna_list), function(x) writeXStringSet(gr.peakset_impt_dna_list[[x]], filepath=paste(dir_path, "/fasta/peaks_", peak_type, "_impt_", type[x], ".fasta", sep=""), append=FALSE))

		# Write out bed files
		gr.peakset_df <- data.frame(gr.peakset)[,1:3]
		gr.peakset_impt_df_rf <- data.frame(gr.peakset_impt_rf)[,1:3]
		gr.peakset_impt_df_glmnet <- data.frame(gr.peakset_impt_glmnet)[,1:3]
		gr.peakset_impt_df_list <- list(gr.peakset_impt_df_rf, gr.peakset_impt_df_glmnet)

		write.table(gr.peakset_df, file=paste(dir_path, "/bed/peaks_", peak_type, "_ALL.bed", sep=""), quote=F, row.names=F, col.names=F, sep="\t")
		lapply(1:length(gr.peakset_impt_df_list), function(x) write.table(gr.peakset_impt_df_list[[x]], file=paste(dir_path, "/bed/peaks_", peak_type, "_impt_", type[x], ".bed", sep=""), quote=F, row.names=F, col.names=F, sep="\t"))

		# Trim datasets
		gr.peakset_trimmed <- trimPeakset(gr.peakset, trimN)
		gr.peakset_impt_rf_trimmed <- trimPeakset(gr.peakset_impt_rf, trimN)
		gr.peakset_impt_glmnet_trimmed <- trimPeakset(gr.peakset_impt_glmnet, trimN)

		gr.peakset_trimmed_dna <- GRanges_toFasta(gr.peakset_trimmed)
		gr.peakset_impt_rf_trimmed_dna <- GRanges_toFasta(gr.peakset_impt_rf_trimmed)
		gr.peakset_impt_glmnet_trimmed_dna <- GRanges_toFasta(gr.peakset_impt_glmnet_trimmed)

		writeXStringSet(gr.peakset_trimmed_dna, filepath=paste(dir_path, "/fasta/peaks_", peak_type, "_ALL_trimmed", trimN, ".fasta", sep=""), append=FALSE)
		writeXStringSet(gr.peakset_impt_rf_trimmed_dna, filepath=paste(dir_path, "/fasta/peaks_", peak_type, "_impt_rf_trimmed", trimN, ".fasta", sep=""), append=FALSE)
		writeXStringSet(gr.peakset_impt_glmnet_trimmed_dna, filepath=paste(dir_path, "/fasta/peaks_", peak_type, "_impt_glmnet_trimmed", trimN, ".fasta", sep=""), append=FALSE)
	}
}	




