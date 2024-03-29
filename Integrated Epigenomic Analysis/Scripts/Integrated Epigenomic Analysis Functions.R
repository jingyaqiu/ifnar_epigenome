#################
### Functions ###
#################

# Convert peak IDs (chr_start_end) to GRanges
peakids2GRanges <- function(peakids, delim="_") {
  df <- data.frame(
    chr = sapply(strsplit(peakids, split = delim), function(x) x[1]), 
    start = as.numeric(sapply(strsplit(peakids, split = delim), function(x) x[2])), 
    end = as.numeric(sapply(strsplit(peakids, split = delim), function(x) x[3])))
  gr <- makeGRangesFromDataFrame(df, starts.in.df.are.0based = FALSE)
  gr$ID <- peakids
  gr$Coordinates <- paste0(as.character(seqnames(gr)), ":", start(gr), "-", end(gr))
  return(gr)
}

make_HOMER_bed <- function(gr, file=NULL) {
    bed <- data.frame(gr)[,1:3]
    bed$ID <- paste(bed[,1], bed[,2], bed[,3], sep="_")
    bed$col5 <- NA
    bed$strand <- as.character(strand(gr))
    if(!is.null(file)) write.table(bed, file=file, sep="\t", row.names=F, col.names=F, quote=F)
}

# Plot PCA, color by annotation
plot_PCA <- function(dat, n_var_features=2500, annotation, color.pal, title) {

	# Use most variable features for PCA
	rv <- rowVars(as.matrix(dat))
	select <- order(rv, decreasing=TRUE)[seq_len(min(n_var_features, length(rv)))]
	pca <- prcomp(t(dat[select,]), scale=T, center=T)
	df.pca <- data.frame(pca$x, sample=colnames(dat), anno=annotation)

	# Extract variance explained
	percentVar <- pca$sdev^2 / sum( pca$sdev^2 ) 
	pca <- list(df.pca, percentVar*100)
	names(pca) <- c("x", "percentVar")

	# Plot
	fig.pca <- ggplot(pca[[1]], aes_string("PC1", "PC2", color="anno", label="anno")) +
        geom_point(size=10, alpha=1, shape=1, stroke=3) +
        scale_color_manual(values=color.pal) +
	    theme_bw() + 
	    geom_vline(xintercept=0, linetype="dashed") +
	    geom_hline(yintercept=0, linetype="dashed") +
        xlab(paste0("PC1 (", round(pca[[2]][1], 2), "%)")) +
        ylab(paste0("PC2 (", round(pca[[2]][2], 2), "%)")) +
        labs(col = "Condition") +
	    ggtitle(title) +
	    coord_fixed() +
	    # geom_text_repel(aes(label=rownames(dat)), size=3) +
	    theme(axis.line=element_line(colour="black"),
	          panel.border=element_rect(colour="black", fill=NA, size=1),
	          legend.title=element_text(size=22),
	          legend.text=element_text(size=22),
	          text = element_text(size=22),
	          aspect.ratio=1,
	          plot.title = element_text("Helvetica"))
	return(fig.pca)
}

# Extract normalized data for peak so interest (gr.peaks)
extract_signal_matrix <- function(
	gr.peaks, 
	mat.assay, 
	gr.assay,
	minoverlap = 0) {
	
	ol <- findOverlaps(gr.peaks, gr.assay, minoverlap=minoverlap)
	gr.ol <- gr.assay[to(ol)]
	gr.ol$name <- paste(data.frame(gr.ol)[,1], data.frame(gr.ol)[,2], data.frame(gr.ol)[,3], sep="_")
	mat <- mat.assay[match(unique(gr.ol$name), rownames(mat.assay)), ]
	mat <- t(scale(t(mat)))
	return(mat)
}

# Calculate average signal at each region of interest
calculate_signal_density <- function(
	gr.peaks, 
	assay.plot=c("H3K27Ac", "H3K4me1", "ATAC"), 
	mat_list) {

	density_signals_list <- lapply(assay.plot, function(assay) {
		print(assay)

		# Extract signal matrix at peaks overlapping regions of interest
		mat.peaks <- extract_signal_matrix(
			gr.peaks=gr.peaks, 
			mat.assay=mat_list[[assay]], 
			gr.assay=peakids2GRanges(rownames(mat_list[[assay]]), delim="_"))

		# Average signal at each peak within each condition
		mat.peaks <- data.frame(t(mat.peaks))
		mat.peaks$Condition <- sapply(strsplit(rownames(mat.peaks), split="_|\\."), function(x) paste(x[1:2], collapse="_"))
		mat.means <- mat.peaks %>%
			group_by(Condition) %>%
			summarize_all(mean) %>%
			pivot_longer(cols = !Condition, names_to = "Peak", values_to = "avg_signal")
		mat.means$Assay <- assay
		return(mat.means)
	})
	mat.density <- do.call(rbind, density_signals_list)
	return(mat.density)
}

# Plot signal at regions of interest #
plot_signals_density <- function(
	gr.peaks, 
	assay.plot=c("H3K27Ac", "H3K4me1", "ATAC"), 
	plot.conditions,  
	mat_list, 
	colors) {

	mat.density <- calculate_signal_density(
		gr.peaks=gr.peaks, 
		assay.plot=assay.plot, 
		# gr_list=gr_list, 
		mat_list=mat_list)

	# Density plot
	mat.plot <- mat.density[mat.density$Condition %in% plot.conditions, ]
	mat.plot$Condition <- droplevels(factor(mat.plot$Condition, levels=plot.conditions))
	mat.plot$Assay <- droplevels(factor(mat.plot$Assay, levels=assay.plot))
	mat.summary <- mat.plot %>%
		group_by(Assay, Condition) %>%
		summarize(mean=mean(avg_signal))

	fig <- ggplot(mat.plot, aes(x=avg_signal, color=Condition, fill=Condition)) +
		geom_density(alpha=0.2) +
		facet_wrap(~Assay, nrow=1) +
		scale_color_manual(values=colors) +
		scale_fill_manual(values=colors) +
		geom_vline(data=mat.summary, aes(xintercept=mean, color=Condition), linetype="dashed", size=0.6) +
		labs(title="", x = "Average signal", y = "Density") +
		guides(color=guide_legend(nrow=2, byrow=TRUE)) +
		theme_bw() +
		theme(panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			axis.text.y = element_text(size=14, color="black"),
            axis.text.x = element_text(size=14, color="black"),
			aspect.ratio=0.9,
			text=element_text(size=18),
			strip.background=element_rect(color="black", fill="white", size=0.8),
            strip.text.x=element_text(size=18),
            legend.position="bottom",
            legend.title=element_blank())
	return(fig)
}

# Identify promoter and enhancer ATAC peaks for each gene
identify_cis_REs_ATAC <- function(
	gr.peaks,
	gs, 
	W=92000, 
	gr.bm, 
	gr.tss,
	adat, 
	edat, 
	cor_cutoff=0.2) {

	# Gene cis-regulatory windows
	gr.rna_gs <- gr.bm[na.omit(match(gs, gr.bm$mgi_symbol))]
	gr.rna_gs_window <- gr.rna_gs
	start(gr.rna_gs_window) <- start(gr.rna_gs) - W
	end(gr.rna_gs_window) <- end(gr.rna_gs) + W

	# ATAC peak coordinates
	gr.adat <- peakids2GRanges(rownames(adat), delim="_")

	# Identify promoter (TSS) and enhancers (r > cor_cutoff) for each gene of interest
	gs_linked_REs <- lapply(1:length(gs), function(x) {
		ol <- findOverlaps(gr.peaks, gr.rna_gs_window[x], minoverlap=100)
		gr.ol <- gr.peaks[unique(from(ol))]
		if(length(gr.ol) < 1) {
			idx <- list(NA, NA)
		} else {
			gr.tss_gs <- gr.tss[match(gr.rna_gs_window[x]$mgi_symbol, gr.tss$gene)]
			ol.tss <- findOverlaps(gr.ol, gr.tss_gs)
			if(length(ol.tss) < 1 & length(gr.ol[gr.ol$RE == "enhancer"]) < 1) {
				idx <- list(NA, NA)
			} else if(length(gr.ol[gr.ol$RE == "enhancer"]) < 1) {
				gr.p <- gr.ol[unique(from(ol.tss))]
				idx <- list(gr.p$idx, NA)
			} else {
				# Annotate promoters
				gr.p <- gr.ol[unique(from(ol.tss))]

				# Retain only "enhancers" that correlate with gene expression
				gr.e <- gr.ol[gr.ol$RE == "enhancer"]
				ol <- findOverlaps(gr.e, gr.adat)
				mat.e <- adat[to(ol), ]
				mat.e$ID <- from(ol)
				mat.e_mean <- mat.e %>% 
					group_by(ID) %>%
					summarize_all(mean)
				rna <- as.numeric(edat[match(gr.rna_gs_window[x]$ensembl_gene_id, rownames(edat)), ])
				atac <- mat.e_mean[,colnames(mat.e_mean) %in% colnames(edat)]
				cors <- sapply(1:nrow(atac), function(y) cor(rna, as.numeric(atac[y,])))
				gr.e$cors <- NA
				gr.e$cors[1:length(gr.e) %in% from(ol)] <- cors
				gr.e <- gr.e[!is.na(gr.e$cors)] # Remove closed REs - but maybe don't have to do this!!
				gr.e <- gr.e[abs(gr.e$cors) > cor_cutoff]
				idx <- list(gr.p$idx, gr.e$idx)
				# idx <- unique(c(gr.p$idx, gr.e$idx))
			}
		}
		return(idx)
	})
	names(gs_linked_REs) <- gs
	return(gs_linked_REs)
}

# Summarize signals in cis-regulatory windows
calculate_cis_genes_average_signal <- function(
	gs_linked_REs,
	gr.peaks,
	mat_list,
	assay.plot,
	plot.conditions,
	use.pal,
	summarize_genes = TRUE,
	method = "gsva",
	return_values = TRUE,
	plot.boxplot = FALSE,
	point_size = 1.5,
    stroke = 1,
    font_size = 14,
    font_size_strip = 14,
    title = "",
    legend_title = "Condition",
    aspect.ratio = 1,
    ylim = NULL) {

	gs.ensembl <- gr.bm$ensembl_gene_id[match(names(gs_linked_REs), gr.bm$mgi_symbol)]
	gr.assay <- gr.peaks[unique(unlist(gs_linked_REs))]
	out_list <- lapply(1:length(assay.plot), function(x) {
		assay <- assay.plot[x]
		print(assay)

		# Extract signal matrix at peaks overlapping regions of interest
		if(assay == "RNA") {
			# mat.peaks <- mat_list[["RNA"]][match(gs.ensembl, rownames(mat_list[["RNA"]])), ]
			# mat.peaks <- t(scale(t(mat.peaks)))
			ids <- gs.ensembl
		} else {
			ol <- findOverlaps(gr.assay, peakids2GRanges(rownames(mat_list[[assay]])))
			ids <- rownames(mat_list[[assay]])[unique(to(ol))]
		}
			# if(summarize_genes == TRUE) {
			# 	mat.peaks <- lapply(1:length(gs.ensembl), function(y) {
			# 		gr <- gr.peaks[gs_linked_REs[[y]]]
			# 		mat <- extract_signal_matrix(
			# 			gr.peaks = gr, 
			# 			mat.assay = mat_list[[assay]], 
			# 			gr.assay = peakids2GRanges(rownames(mat_list[[assay]]), delim="_"))
			# 		return(as.numeric(colMeans(mat)))
			# 	})
			# 	mat.peaks <- do.call(rbind, mat.peaks)
			# 	rownames(mat.peaks) <- gs.ensembl
			# 	colnames(mat.peaks) <- colnames(mat_list[[assay]])
			# 	mat.peaks <- mat.peaks[!is.na(mat.peaks[,1]), ]
			# 	mat.peaks <- t(scale(t(mat.peaks)))
			# } else {
			# 	# mat.peaks <- extract_signal_matrix(
			# 	# 	gr.peaks=gr.peaks[unique(unlist(gs_linked_REs))], 
			# 	# 	mat.assay=mat_list[[assay]], 
			# 	# 	gr.assay=peakids2GRanges(rownames(mat_list[[assay]]), delim="_"))
				
			# }
		# df <- data.frame(
		# 	Sample=colnames(mat.peaks), 
		# 	Signal=as.numeric(colMeans(mat.peaks)))
		# df$Condition <- sapply(strsplit(as.character(df$Sample), split="_|\\."), function(x) paste(x[1:2], collapse="_"))
		# df <- df[df$Condition %in% plot.conditions, ]
		# df$Condition <- factor(df$Condition, levels=plot.conditions)
		# df$Assay <- assay

		# Plot average signal
		gsva <- gsva(
	        expr = as.matrix(mat_list[[assay]]),
	        gset.idx.list = list(ids),
	        method = ifelse(method == "both", ifelse(assay == "RNA", "gsva", "zscore"), method),
	        kcdf = "Gaussian",
	        min.sz = 10,
	        max.sz = 1000000)
	    gsva <- t(scale(t(gsva)))

	    df <- data.frame(
            Sample = colnames(gsva),
            Signal = as.numeric(gsva),
            Assay = assay)
	    df$Condition <- sapply(strsplit(as.character(df$Sample), split="_|\\."), function(x) paste(x[1:2], collapse="_"))
		df <- df[df$Condition %in% plot.conditions, ]
		df$Condition <- factor(df$Condition, levels=plot.conditions)

		if(return_values == TRUE) {
			return(df)
		} else {
			p <- ggplot(df, aes(x=Condition, y=Signal, color=Condition)) +
				facet_wrap(~Assay, nrow=1, drop=TRUE)
			if(plot.boxplot == TRUE) {
		        p <- p +
		            geom_boxplot(width = 0.4, outlier.shape = NULL, fatten=NULL)
		    }
		    p <- p +
		    	geom_jitter(shape=1, size = point_size, stroke = stroke, position=position_jitter(0.1)) +
		        # facet_wrap(~Feature, nrow = round(length(geneset_list) / 5), labeller = labeller(Feature = label_wrap_gen(20))) +
		        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
		        scale_colour_manual(values=use.pal) +
		        labs(x = "", y="Signal", title = title, color = legend_title) +
		        guides(color = guide_legend(override.aes = list(shape = 16), ncol = 1)) +
		        theme_bw(base_size=font_size) +
		        theme(
		            axis.text.y = element_text(size = font_size, color = "black"),
		            axis.title.y=element_text(size=font_size),
		            axis.title.x=element_text(size=font_size),
		            aspect.ratio = aspect.ratio,
		            strip.background = element_rect(fill="white"),
		            strip.text.x = element_text(margin = ggplot2::margin(5, 0, 5, 0, "pt"), size = font_size_strip),
		            legend.position="none",
		            axis.text.x = element_text(size=font_size, angle = 45, hjust=1, color = "black"))

		    if(x == 1) {
				p <- p + labs(y="Mean ISG signal", x="")
			} else {
				p <- p + labs(y="", x="")
			}

			if(!is.null(ylim)) {
				p <- p +
					ylim(ylim)
			}
			return(p)
		}
	})
	return(out_list)
}

plot_summary <- function(
    mat.plot, 
    geneset_list, 
    geneset_names, 
    metadat,
    groupBy = "Condition",
    colorBy = "Condition",
    facetBy = NULL,
    normalize = NULL,
    use.pal = pal_uchicago()(6),
    method = "gsva",
    plot.boxplot = TRUE,
    plot.legend = TRUE,
    ylim = NULL,
    show.xaxis.labels = TRUE,
    font_size = 14,
    font_size_strip = 14,
    point_size = 1.5,
    stroke = 1,
    aspect.ratio = 1,
    nrow = 1,
    ncol = 1,
    title = "",
    legend_title = "Condition",
    return_values = FALSE) {

    print(paste0("Count data and sample metadata match: ", all(metadat$Sample == colnames(mat.plot))))
    
    # Summary
    gsva <- gsva(
        expr = as.matrix(mat.plot),
        gset.idx.list = geneset_list,
        method = method,
        kcdf = "Gaussian",
        min.sz = 10,
        max.sz = 1000000)
    gsva <- t(scale(t(gsva)))

    if(length(geneset_list) > 1) {
        mat.plot_Summary <- data.frame(
            Sample = colnames(gsva),
            Group = metadat[match(colnames(gsva), metadat$Sample), groupBy],
            Color_Group = metadat[match(colnames(gsva), metadat$Sample), colorBy],
            t(gsva))
        if(!is.null(facetBy)) mat.plot_Summary$Facet_Group <- metadat[match(colnames(gsva), metadat$Sample), facetBy]

        mat.plot_Summary <- mat.plot_Summary %>% 
            pivot_longer(
                cols = make.names(names(geneset_list)), 
                names_to = "Feature", 
                values_to = "Signal")
        mat.plot_Summary$Feature <- factor(mat.plot_Summary$Feature, levels = make.names(names(geneset_list)))
        levels(mat.plot_Summary$Feature) <- geneset_names
    } else {
        mat.plot_Summary <- data.frame(
            Sample = colnames(gsva),
            Group = metadat[match(colnames(gsva), metadat$Sample), groupBy],
            Color_Group = metadat[match(colnames(gsva), metadat$Sample), colorBy],
            Signal = as.numeric(gsva),
            Feature = geneset_names)
        if(!is.null(facetBy)) mat.plot_Summary$Facet_Group <- metadat[match(colnames(gsva), metadat$Sample), facetBy]
    }

    if(!is.null(normalize)) {
        # mat.plot_Summary <- mat.plot_Summary %>%
        #     group_by(get(normalize)) %>%
        #     mutate(Signal_Normalized = (Signal - min(Signal)) / (max(Signal) - min(Signal)))
        ref_group <- levels(mat.plot_Summary$Facet_Group)[1]
        mat.plot_Summary <- mat.plot_Summary %>%
            group_by(get(normalize)) %>%
            mutate(Signal_Normalized = Signal - min(Signal) + 1) %>%
            mutate(Signal_Normalized = Signal_Normalized / min(Signal_Normalized[Facet_Group == ref_group]))
    }
    
    # Visualize
    fig <- ggplot(mat.plot_Summary, aes_string(x="Group", y=ifelse(!is.null(normalize), "Signal_Normalized", "Signal"), color = "Color_Group"))
    if(plot.boxplot == TRUE) {
        fig <- fig +
            geom_boxplot(width = 0.4, outlier.shape = NULL, fatten=NULL)
    }
    if(!is.null(facetBy) & length(geneset_list) == 1) {
        fig <- fig +
            facet_wrap(~Facet_Group, nrow = nrow)
    } else if(!is.null(facetBy) & length(geneset_list) > 1) {
        fig <- fig +
            facet_wrap(~Facet_Group + Feature, nrow = nrow, ncol = ncol)
    } else {
        fig <- fig +
            facet_wrap(~Feature, nrow = nrow)
    }
    fig <- fig +
        geom_jitter(shape=1, size = point_size, stroke = stroke, position=position_jitter(0.1)) +
        # facet_wrap(~Feature, nrow = round(length(geneset_list) / 5), labeller = labeller(Feature = label_wrap_gen(20))) +
        stat_summary(fun.y = "mean", fun.ymin = "mean", fun.ymax= "mean", width=0.5, size=0.5, geom = "crossbar") +
        scale_colour_manual(values=use.pal) +
        labs(x = "", y="Signal", title = title, color = legend_title) +
        guides(color = guide_legend(override.aes = list(shape = 16), ncol = 1)) +
        theme_bw(base_size=font_size) +
        theme(
            axis.text.y = element_text(size = font_size, color = "black"),
            axis.title.y=element_text(size=font_size),
            axis.title.x=element_text(size=font_size),
            aspect.ratio = aspect.ratio,
            strip.background = element_rect(fill="white"),
            strip.text.x = element_text(margin = ggplot2::margin(5, 0, 5, 0, "pt"), size = font_size_strip))
    if(!is.null(ylim)) {
        fig <- fig +
            ylim(ylim)
    }
    if(show.xaxis.labels == TRUE) {
        fig <- fig +
            theme(axis.text.x = element_text(size=font_size, angle = 45, hjust=1, color = "black"))
    } else {
        fig <- fig + 
            theme(axis.text.x = element_blank())
    }
    if(plot.legend == FALSE) {
        fig <- fig +
            theme(legend.position = "none")
    }

    if(return_values == TRUE) {
        out <- mat.plot_Summary
    } else {
        out <- fig
    }
    return(out)
}
