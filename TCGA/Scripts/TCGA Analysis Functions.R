### FUNCTIONS ###

# Fit lm model using ISG gene sets as predictors and CD8 activity as dependent variable
fit_lm_regression <- function(
	df,
	use.metric = "cytolytic_score",
	predictors = c("RNA_ISG.RS", "RNA_IFNG.GS"),
	cancer_type = NA) {
	
	# Fit lm model
	df$Score <- df[ ,use.metric] # Outcome is CD8 cytolytic_score
	dat <- df[ ,c(predictors, "Score")]
	dat <- data.frame(scale(dat)) # Get standardized regression coefficients
	fit <- lm(formula = paste0("Score ~ ", paste(predictors, collapse=" + ")), data = dat)

	# Extract model stats
	stats_summary <- lapply(1:length(predictors), function(x) get_lm_stats(fit, predictor = predictors[x]))
	df.stats <- do.call(rbind, stats_summary)
	df.stats$cancer_type <- cancer_type
	df.stats$gs <- predictors
	return(df.stats)
}

# Extract lm stats from fit
get_lm_stats <- function(fit, predictor) {
	summary <- summary(fit)
	coefs <- summary$coefficients
	confint <- confint(fit, level=0.95)
	model_summary <- glance(fit)
	adjRsq <- model_summary$adj.r.squared
	model_pval <- as.numeric(model_summary$p.value)
	df.stats <- data.frame(coefs, confint, adjRsq = adjRsq, model_pval = model_pval)
	colnames(df.stats) <- c("est", "se", "t.value", "coef_pval", "lower", "upper", "adjRsq", "model_pval")
	return(df.stats[predictor,])
}

# Plot regression coefficients for lm regression in forest plot
visualize_forest_plot <- function(
	fit.stats, 
	metadat, 
	predictors = c("RNA_IFNG.GS", "RNA_ISG.RS"),
	colors = c("navyblue", "firebrick"),
	labels = c("Immune ISGs", "Cancer ISGs"),
	plot.cancers) {

	metadat <- metadat[metadat$cancer_type %in% plot.cancers, ]

	# Sample size for each cancer
	sample_sizes <- table(metadat$cancer_type)[match(plot.cancers, names(table(metadat$cancer_type)))]

	p_list <- lapply(1:length(predictors), function(x) {
		df <- fit.stats[which(as.character(fit.stats$gs) == predictors[x]), ]
		df <- df[match(plot.cancers, df$cancer_type), ]
		df$name <- paste0(df$cancer_type, " (n=", as.numeric(sample_sizes), ")")
		df.sorted <- df[sort.int(df$est, decreasing = TRUE, index.return = TRUE)$ix, ]
		df.sorted$name <- factor(df.sorted$name, levels=rev(as.character(df.sorted$name)))

		p <- ggplot(data=df.sorted, aes(x=name, y=est, ymin=lower, ymax=upper)) +
			geom_pointrange(color = colors[x], fatten = 1, size = 1) +
			geom_hline(yintercept=0, lty=2) +
			coord_flip() +
			labs(x = "Cancer type", y = "Beta coefficient (95% CI)", title = labels[x]) +
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
	# print(table(metadat.subsample$cancer_type))

	fit.stats_combined <- fit_lm_regression(
		df = metadat.subsample, 
		use.metric = use.metric, 
		predictors = predictors,
		cancer_type = "ALL")
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