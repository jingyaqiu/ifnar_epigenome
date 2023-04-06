

color_2 <- brewer.pal(11, "Spectral")[c(10,2,4)]

p_theme2 <- theme_bw(base_size = 12) +
    theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        legend.text = element_text(size = 8))

p_theme3 <- theme_bw(base_size = 12) +
    theme(legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#############################################################
### Figure 4C: GSEA on type I ISGs in Res 499 and Res 237 ###
#############################################################

library(DESeq2)
library(fgsea)

plotEnrichment2 <- function (pathway, stats, leadingEdge = NULL, showLE = TRUE, 
                             le.text.size = 4, gseaParam = 1, nudge.y = 0.05) {
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, 
                          returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) +
    geom_point(color = "green", size = 0.1) + 
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(color = "green") + theme_bw() + 
    geom_segment(data = data.frame(x = pathway),
                 mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = 0.2) + 
    theme(panel.grid.minor = element_blank()) + 
    labs(x = "rank", y = "enrichment score")
  if (showLE) {
    LEdat <- data.frame(x = pathway, y = 0.1, gene = names(gseaRes$tops))
    if (is.null(leadingEdge)) leadingEdge <- LEdat$gene
    g <- g + geom_text_repel(data = subset(LEdat, LEdat$gene %in% leadingEdge),
                             aes(x = x, y = diff/2, label = gene),
                             direction = "x", angle = 90, size = le.text.size,
                             segment.size = 0.2, segment.alpha = 0.5,
                             nudge_y = nudge.y)
  }
  return(g)
}

### In vitro RNA-seq in Res499 ###

load(file.path("~/Dropbox/Xu Qiu IFN paper/V2/18. In vitro Oas1 RNA seq/anno.ens.mm10.RData"))
anno <- anno.ens
load(file.path("~/Dropbox/Xu Qiu IFN paper/V2/18. In vitro Oas1 RNA seq/dds.vitro.Rdata"))

# Getting DEGs, specifying after "condition", "expt", "ctrl"
res.R499 <- results(dds.vitro, contrast = c("condition", "R499_OasKO", "R499"))
res.R499.all <- as.data.frame(res.R499)
z <- match(anno$GeneID, rownames(res.R499.all))
z1 <- anno[!is.na(z),]
res.R499.all.anno.filt <- res.R499.all[-which(is.na(match(rownames(res.R499.all), z1$GeneID)) == TRUE),]
res.R499.all.anno.filt$GeneID <- rownames(res.R499.all.anno.filt)
res.R499.all.anno <- merge(res.R499.all.anno.filt, z1, by = "GeneID")

# Map mouse gene ID to human symbol
load(file.path("~/Dropbox/Xu Qiu IFN paper/V2/18. In vitro Oas1 RNA seq/mouse2human.Rdata"))

# Join human symbol to res.all.anno
res.tibble <- res.R499.all.anno %>%
    as_tibble()
res.mouse2human<- inner_join(res.tibble, bm, by=c("GeneID" = "ensembl_gene_id"))
res.mouse2human.filt <- res.mouse2human %>% 
    dplyr::select(hsapiens_homolog_associated_gene_name, stat) %>% 
    na.omit() %>%
    distinct() %>% 
    group_by(hsapiens_homolog_associated_gene_name) %>% 
    dplyr::summarize(stat=mean(stat))

# Main fGSEA call
ranks <- deframe(res.mouse2human.filt)
pathways.hallmark <- gmtPathways("~/Dropbox/Xu Qiu IFN paper/V2/18. In vitro Oas1 RNA seq/h.all.v6.2.symbols.gmt")
genesets <- vector("list", 1) # if only want to check one geneset
genesets [[1]] <- pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE
names(genesets) <- "IFN.I"
fgseaRes <- fgsea(pathways=genesets, stats=ranks, nperm=1000)

# pval: 0.001980198; NES: 1.657026

p_R499.GSA <- plotEnrichment2(
    pathways.hallmark$HALLMARK_INTERFERON_ALPHA_RESPONSE, 
    showLE = FALSE,
    stats = ranks, 
    le.text.size = 4,
    nudge.y = 0.1) + 
    ylim(-0.5, 0.1) + 
    xlab("Gene Rank") + ylab("Enrichment Score") +
    theme(axis.title = element_text(size = 12)) +
    p_theme3
p_R499.GSA

### In vitro R237 GSA ###

res.R237 <- results(dds.vitro, contrast = c("condition", "R237_OasKO", "R237"))
res.R237.all <- as.data.frame(res.R237)
z <- match(anno$GeneID, rownames(res.R237.all))
z1 <- anno[!is.na(z),]
res.R237.all.anno.filt <- res.R237.all[-which(is.na(match(rownames(res.R237.all), z1$GeneID)) == TRUE),]
res.R237.all.anno.filt$GeneID <- rownames(res.R237.all.anno.filt)
res.R237.all.anno <- merge(res.R237.all.anno.filt, z1, by = "GeneID")

# Join human symbol to res.all.anno
res.tibble <- res.R237.all.anno %>%
    as_tibble()
res.mouse2human <- inner_join(res.tibble, bm, by=c("GeneID" = "ensembl_gene_id"))
res.mouse2human.filt <- res.mouse2human %>% 
    dplyr::select(hsapiens_homolog_associated_gene_name, stat) %>% 
    na.omit() %>%
    distinct() %>% 
    group_by(hsapiens_homolog_associated_gene_name) %>% 
    dplyr::summarize(stat=mean(stat))

# Main fGSEA call
ranks <- deframe(res.mouse2human.filt)
fgseaRes <- fgsea(pathways=genesets, stats=ranks, nperm=1000)

# pval: 0.001703578; NES: 2.69285

p_R237.GSA <- plotEnrichment2(
    pathways.hallmark[["HALLMARK_INTERFERON_ALPHA_RESPONSE"]], 
    showLE = FALSE,
    stats = ranks, 
    le.text.size = 4,
    nudge.y = 0.1) + 
    ylim(-0.85, 0.1) + 
    xlab("Gene Rank") + 
    ylab("Enrichment Score") +
    theme(axis.title = element_text(size = 12)) +
    p_theme3
p_R237.GSA

# Replicates: Res499 (n=2), Res499 OAS1 KO (n=2), Res 237 (n=2), Res 237 OAS1 KO (n=2)

######################################################################################################
### Figure 4D: PolyIC treatment induces higher IFNa/b, ISGs in Res 499, blunted by OAS1 KO (Bihui) ###
######################################################################################################

## Res499 produces more IFNa and IFNb compared to B16 after pIC; Oas1 blunts this

dat <- read.delim("~/Dropbox/Xu Qiu IFN paper/V2/1. pIC B16 R499 Oas1 IFN/B16_R499_Oas1.txt")
dat$Cell <- factor(dat$Cell, levels = c("B16", "Res499", "Oas1KO"))
levels(dat$Cell) <- c("B16", "Res499", "Oas1KO")
Step1 <- dat %>% filter(IFN %in% "IFNa")
p_R499.IFN <- ggplot(data = Step1, aes(x = Cell, y = Conc, color = Cell)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(w = 0.2), size = 1.5) +
    facet_wrap(~IFN, scales = "free") +
    scale_color_manual(values = color_2) +
    ylab("Concentration (pg/ml)") +
    p_theme2 +
    theme(legend.position = "none")
p_R499.IFN

with(subset(dat, subset = IFN %in% "IFNa" & Cell %in% c("B16", "Res499")),
     t.test(Conc ~ Cell, alternative = "two.sided"))   
with(subset(dat, subset = IFN %in% "IFNa" & Cell %in% c("Oas1KO", "Res499")),
     t.test(Conc ~ Cell, alternative = "two.sided"))
with(subset(dat, subset = IFN %in% "IFNb" & Cell %in% c("B16", "Res499")),
     t.test(Conc ~ Cell, alternative = "two.sided"))  
with(subset(dat, subset = IFN %in% "IFNb" & Cell %in% c("Oas1KO", "Res499")),
     t.test(Conc ~ Cell, alternative = "two.sided"))  

# p-value: 1.294e-07, 6.671e-07,  0.04117, 0.009958
# Replicates: B16 (n=11), Res499 (n=23), Oas1KO (n=17)

## Res499 induces higher ISGs (Isg15, Mx1) compared to B16 after pIC; Oas1 blunts it

dat <- read.delim("~/Dropbox/Xu Qiu IFN paper/V2/2. pIC B16 R499 ISGs/B16 R499 ISG induction.txt")
dat$Cell <- factor(dat$Cell, levels = c("B16", "R499", "Oas1KO"))
p_ISG <- ggplot(data = dat, aes(x = Cell, y = Expression, color = Cell)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(w = 0.2), size = 1.5) +
    facet_wrap(~ISG, scales = "free") +
    scale_color_manual(values = color_2) +
    ylab("Relative Expression") +
    p_theme2 +
    theme(legend.position = "none")
p_ISG

with(subset(dat, subset = ISG %in% "Isg15" & Cell %in% c("B16", "R499")),
     t.test(Expression ~ Cell, alternative = "two.sided"))   
with(subset(dat, subset = ISG %in% "Isg15" & Cell %in% c("Oas1KO", "R499")),
     t.test(Expression ~ Cell, alternative = "two.sided"))
with(subset(dat, subset = ISG %in% "Mx1" & Cell %in% c("B16", "R499")),
     t.test(Expression ~ Cell, alternative = "two.sided"))  
with(subset(dat, subset = ISG %in% "Mx1" & Cell %in% c("Oas1KO", "R499")),
     t.test(Expression ~ Cell, alternative = "two.sided"))  

# p-val: 0.0003026, 0.004987, 0.001725, 0.002594
# Replicates: B16 (n=8), Res499 (n=22), Oas1KO (n=14)

##################################################################
### Figure 4F: STAT1/IRF3 KO reverses memory at ISG.RS cis-REs ###
##################################################################

print("Load normalized count matrices")

mat_list <- lapply(1:length(assays), function(x) {
    if(assays[x] == "ATAC") {
        mat <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_original_consensus_mat_tn5_insertion_counts_IDR_rlog.txt", sep="\t", header=T, stringsAsFactors=F)
        # mat <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/consensus_mat_tn5_insertion_counts_IDR_vst.txt", sep="\t", header=T, stringsAsFactors=F)
        colnames(mat) <- gsub("cas", "WT", colnames(mat))
    } else {
        mat <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/", assays[x], "_consensus_mat_insertion_counts_WT_pooled_vst.txt"), sep="\t", header=T, stringsAsFactors=F)
    }
    return(mat)
})
names(mat_list) <- assays
mat_list[["RNA"]] <- mat.RNA_original

# Revisions ATAC data
mat_list1 <- mat_list
mat_list1[["ATAC"]] <- read.table(paste0("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/ATAC_revisions_consensus_mat_tn5_insertion_counts_IDR_vst.txt"), sep="\t", stringsAsFactors=F, header=T)
colnames(mat_list1[["ATAC"]]) <- gsub("\\.", "_", colnames(mat_list1[["ATAC"]]))

# Combine annotated promoter and enhancer REs
gr.promoter <- RE_repertoire_list[[1]]
gr.promoter$RE <- "promoter"
gr.enhancer <- RE_repertoire_list[[2]]
gr.enhancer$RE <- "enhancer"
ol <- findOverlaps(RE_list[["activated_enhancer_ATAC_inclusive"]], gr.enhancer)
gr.enhancer <- gr.enhancer[unique(to(ol))]
RE.combined <- c(gr.promoter, gr.enhancer)
RE.combined$idx <- 1:length(RE.combined)

gs.name <- "ISG.RS" # ISG.RS, IFN.I
gs <- gs_list[[gs.name]]
gs <- gs[gs %in% rownames(mat.RNA_original)]
gs <- gr.bm$mgi_symbol[match(gs, gr.bm$ensembl_gene_id)]

W <- 92000

# Link REs to target genes
# Promoter - overlap with TSS
# Enhancer - ATAC correlation with RNA with r > 0.1
gs_linked_REs <- identify_cis_REs_ATAC(
    gr.peaks = RE.combined,
    gs = gs, 
    W = W, 
    gr.bm = gr.bm, 
    gr.tss = gr.tss_window,
    adat = mat.ATAC_original, 
    edat = mat.RNA_original,
    cor_cutoff = 0.2)
gs_linked_REs_filt <- sapply(gs_linked_REs, function(x) na.omit(unique(c(x[[1]], x[[2]]))))

# Remove genes with no associated REs
remove_idx <- unlist(sapply(1:length(gs_linked_REs_filt), function(x) {
    if(all(is.na(gs_linked_REs_filt[[x]])) | length(gs_linked_REs_filt[[x]])==0) return(x)
}))
if(length(remove_idx) > 0) gs_linked_REs_filt <- gs_linked_REs_filt[-remove_idx]

# ISG.RS DKO
p <- calculate_cis_genes_average_signal(
    gs_linked_REs = gs_linked_REs_filt, 
    gr.peaks = RE.combined,
    mat_list = mat_list1,
    assay.plot = "ATAC",
    plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO"),
    use.pal = pal_nejm()(6),
    method = "gsva",
    summarize_genes = TRUE,
    return_values = FALSE,
    plot.boxplot = FALSE,
    point_size = 5,
    stroke = 1.2,
    title = "ISG.RS cis-REs",
    ylim = c(-2,2))
# ggsave(p[[1]], file=paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", gs.name, " signals summary DKO.pdf"), width=5, height=4)

df <- calculate_cis_genes_average_signal(
    gs_linked_REs = gs_linked_REs_filt, 
    gr.peaks = RE.combined,
    mat_list = mat_list1,
    assay.plot = "ATAC",
    plot.conditions = c("B16_WT", "B16_SKO", "R499_WT", "R499_SKO", "R499_IRF3KO", "R499_DKO"),
    use.pal = pal_nejm()(6),
    method = "gsva",
    summarize_genes = TRUE,
    return_values = TRUE,
    plot.boxplot = FALSE,
    point_size = 5,
    stroke = 1.2,
    title = "ISG.RS cis-REs",
    ylim = c(-2,2))
df <- df[[1]]
t.test(df$Signal[df$Condition == "B16_WT"], df$Signal[df$Condition == "R499_WT"])
t.test(df$Signal[df$Condition == "B16_WT"], df$Signal[df$Condition == "B16_SKO"])
t.test(df$Signal[df$Condition == "B16_SKO"], df$Signal[df$Condition == "R499_SKO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_SKO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_IRF3KO"])
t.test(df$Signal[df$Condition == "R499_WT"], df$Signal[df$Condition == "R499_DKO"])

# # Write out plot data
# dat <- df[ ,c("Sample", "Condition", "Signal")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 4F.csv"), sep=",", quote=F, col.names=T, row.names=F)

############################################################
### Figure 4G: STAT1/IRF3 KO reverses memory at IFN-IMDs ###
############################################################

set.seed(12)
gr.atac_revisions <- peakids2GRanges(rownames(mat_list1[["ATAC"]]), delim="_")
gr.atac_revisions$ID <- rownames(mat_list1[["ATAC"]])

metadat.ATAC_revisions <- read.csv("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Feature Matrices/metadata_ATAC_revisions.csv", header=T)
metadat.ATAC_revisions$Label <- factor(metadat.ATAC_revisions$Label, levels = c("B16 WT", "B16 STAT1 KO", "Res499 WT", "Res499 STAT1 KO", "Res499 IRF3 KO", "Res499 DKO"))

# IFN-IMD summary
ol <- findOverlaps(gr.memory, gr.atac_revisions)
gs <- gr.atac_revisions$ID[unique(to(ol))]
p <- plot_summary(
    mat.plot = mat_list1[["ATAC"]],
    geneset_list = list(gs),
    geneset_names = "ATAC",
    metadat = metadat.ATAC_revisions,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(6),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")
ggsave(p, file = paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", name, " Domains DKO Dependence Summary.pdf"), width = 4.5, height = 4.5)

df <- plot_summary(
    mat.plot = mat_list1[["ATAC"]],
    geneset_list = list(gs),
    geneset_names = "ATAC",
    metadat = metadat.ATAC_revisions,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(6),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "Res499 WT"])
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "B16 STAT1 KO"])
t.test(df$Signal[df$Group == "B16 STAT1 KO"], df$Signal[df$Group == "Res499 STAT1 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 STAT1 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 IRF3 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 DKO"])

# # Write out plot data
# dat <- df[ ,c("Sample", "Group", "Signal")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 4G.1.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Heatmap
use.pal_revisions <- pal_nejm()(6)
names(use.pal_revisions) <- levels(metadat.ATAC_revisions$Label)
ha_col_revisions <- HeatmapAnnotation(
    Condition=factor(metadat.ATAC_revisions$Label, levels=names(use.pal_revisions)), 
    col=list(Condition=use.pal_revisions), 
    simple_anno_size=unit(0.7, "cm"),
    border = TRUE)

set.seed(12)
ol <- findOverlaps(gr.memory_filtered, gr.atac_revisions)
gs <- gr.atac_revisions$ID[unique(to(ol))]
mat.plot <- mat_list1[["ATAC"]][match(gs, rownames(mat_list1[["ATAC"]])), ]
mat.plot <- t(scale(t(mat.plot)))
hm <- Heatmap(
    mat.plot,
    name = "ATAC",
    use_raster = TRUE, 
    raster_quality = 1, 
    cluster_columns = F, 
    cluster_rows = T, 
    show_column_names = F, 
    show_row_names = F, 
    row_names_gp = gpar(fontsize = 7),
    col = colorRamp2(seq(-2,2,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)), 
    column_split = metadat.ATAC_revisions$Label[match(colnames(mat.plot), metadat.ATAC_revisions$Sample)],
    bottom_annotation = ha_col_revisions,
    column_title = NULL,
    row_title = NULL,
    border = TRUE)
ht <- draw(hm)

pdf(paste0("~/Dropbox/Minn/ifnar_epigenome_old/results/Reviewer Response Figures/", name, " Domains DKO Dependence Heatmap Activated REs.pdf"), width = 4.5, height = 3.25)
ht
dev.off()

# # Write out plot data
# dat <- mat.plot
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 4G.2.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Resolved Domains summary
ol <- findOverlaps(gr.resolved, gr.atac_revisions)
gs <- gr.atac_revisions$ID[unique(to(ol))]
p <- plot_summary(
    mat.plot = mat_list1[["ATAC"]],
    geneset_list = list(gs),
    geneset_names = "ATAC",
    metadat = metadat.ATAC_revisions,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(6),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "")

df <- plot_summary(
    mat.plot = mat_list1[["ATAC"]],
    geneset_list = list(gs),
    geneset_names = "ATAC",
    metadat = metadat.ATAC_revisions,
    groupBy = "Label",
    colorBy = "Label",
    use.pal = pal_nejm()(6),
    method = "ssgsea",
    plot.boxplot = FALSE,
    show.xaxis.labels = TRUE,
    plot.legend = FALSE,
    font_size = 14,
    font_size_strip = 15,
    point_size = 5,
    stroke = 1.2,
    aspect.ratio = 0.9,
    title = "",
    legend_title = "",
    return_values = TRUE)
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "Res499 WT"])
t.test(df$Signal[df$Group == "B16 WT"], df$Signal[df$Group == "B16 STAT1 KO"])
t.test(df$Signal[df$Group == "B16 STAT1 KO"], df$Signal[df$Group == "Res499 STAT1 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 STAT1 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 IRF3 KO"])
t.test(df$Signal[df$Group == "Res499 WT"], df$Signal[df$Group == "Res499 DKO"])

# # Write out plot data
# dat <- df[ ,c("Sample", "Group", "Signal")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 4G.3.csv"), sep=",", quote=F, col.names=T, row.names=F)

########################################
### Fig 4H: Oas1 KO Survival (Bihui) ###
########################################

library(survival)
library(survminer)
library(pzfx)

# read.prizm: read prizm exported data
read.prizm <- function(file){
    require(reshape2)
    dat <- read.table(file, header = T, sep="\t")
    rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    dat <- melt(dat, id.vars="Day", na.rm=T)
    colnames(dat) <- c("Time","Predictor","Event")
    return(dat)
}

### Oas1KO sensitizes Res499 survival

path <- "7. R499 Oas1 in vivo"
fileName <- file.path("~/Dropbox/Xu Qiu IFN paper/V2/7. R499 Oas1 in vivo/499 Oas repeat survival_pool.txt")
survdat <- read.prizm(fileName)
levels(survdat$Predictor)
Step1 <- survdat
levels(Step1$Predictor) <- c("WT NT", "WT CP", "Oas1KO NT", "Oas1KO CP")
Step2 <- Step1 %>% filter(Predictor %in% c("WT NT", "WT CP", "Oas1KO NT",
                                           "Oas1KO CP"))
Step3 <- separate(Step2, Predictor, c("Cell", "Treatment"), " ")
Step3$Cell <- factor(Step3$Cell, levels = c("WT", "Oas1KO"))
Step3$Treatment <- factor(Step3$Treatment, levels = c("NT", "CP"))

res <- pairwise_survdiff(Surv(Time, Event) ~ Predictor, 
    data = Step2,
    p.adjust.method = "none")
res$p.value # pvals

# pval: CP cas vs CP Oas: 0.00570
sfit <- survfit(Surv(Time, Event) ~ Treatment + Cell, data = Step3)
p_R499Oas1KO.surv <- ggsurv2(sfit, facetby = "Treatment", facet.order = c("NT", "CP"), size.est = 1,
                             surv.col = color_2, main = "",
                             reflevel = "WT") + p_theme3
p_R499Oas1KO.surv

######################################################
### Figure 4I: RIG-I KO Survival in B16 or Res 499 ###
######################################################

### Survival ###

library(survival)
library(survminer)

source("~/Dropbox/Minn/resources/useful_protocols/Custom ggsurv function.R")

input_file <- "~/Dropbox/Minn/ifnar_epigenome_old/Darwin Revision Expts/RIG KO Survival.csv"
output_file <- "~/Dropbox/Minn/ifnar_epigenome_old/Darwin Revision Expts/Figures/RIG KO Survival.pdf"

dat <- read.csv(input_file, header=T)
dat$Treatment <- factor(gsub("nt", "NT", dat$Treatment), levels=c("NT", "dICB"))
dat$Genotype <- factor(dat$Genotype, levels=c("b16", "b16rigi", "res499", "res499rigi"))
levels(dat$Genotype) <- c("B16 WT", "B16 RIG-I KO", "Res499 WT", "Res499 RIG-I KO")

# # Write out plot data
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 4I.csv"), sep=",", quote=F, col.names=T, row.names=F)

# Compute survival curve for censored data
fit <- survfit( Surv(Days, Survival) ~ Treatment + Genotype, data = dat)
res <- pairwise_survdiff(Surv(Days, Survival) ~ Treatment + Genotype, data = dat)
res$p.value # pvals

# Visualize

use.pal <- brewer.pal(11, "Spectral")[c(10,2,4)]

dat.B16 <- dat[grep(dat$Genotype, pattern="B16"), ]
dat.B16$Genotype <- droplevels(dat.B16$Genotype)
dat.R499 <- dat[grep(dat$Genotype, pattern="Res499"), ]
dat.R499$Genotype <- droplevels(dat.R499$Genotype)

fit.B16 <- survfit( Surv(Days, Survival) ~ Treatment + Genotype, data = dat.B16)
fit.R499 <- survfit( Surv(Days, Survival) ~ Treatment + Genotype, data = dat.R499)

p.survival_B16 <- ggsurv2(
    s = fit.B16, 
    facetby = "Treatment", 
    facet.order = c("NT", "dICB"),
    surv.col = use.pal,
    size.est = 1,
    reflevel = "B16 WT") + 
    theme_bw(base_size = 16) +
    theme(legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(color = "black"))

p.survival_R499 <- ggsurv2(
    s = fit.R499, 
    facetby = "Treatment", 
    facet.order = c("NT", "dICB"),
    surv.col = use.pal,
    size.est = 1,
    reflevel = "Res499 WT") + 
    theme_bw(base_size = 16) +
    theme(legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        axis.text = element_text(color = "black"))

p.survival <- ggarrange(p.survival_B16, p.survival_R499, nrow=2)
p.survival

ggsave(p.survival, file=output_file, width=6.5, height=5.5)

###################################################################################
### Extended Data Figure 5I: STAT1+IRF3 DKO reverses memory at Res499 enhancers ###
###################################################################################

# Heatmap
use.pal_revisions <- pal_nejm()(6)
names(use.pal_revisions) <- levels(metadat.ATAC_revisions$Label)
ha_col_revisions <- HeatmapAnnotation(
    Condition=factor(metadat.ATAC_revisions$Label, levels=names(use.pal_revisions)), 
    col=list(Condition=use.pal_revisions), 
    simple_anno_size=unit(0.7, "cm"),
    border = TRUE)

set.seed(12)
ol <- findOverlaps(gr.memory_filtered, gr.atac_revisions)
gs <- gr.atac_revisions$ID[unique(to(ol))]
mat.plot <- mat_list1[["ATAC"]][match(gs, rownames(mat_list1[["ATAC"]])), ]
mat.plot <- t(scale(t(mat.plot)))
hm <- Heatmap(
    mat.plot,
    name = "ATAC",
    use_raster = TRUE, 
    raster_quality = 1, 
    cluster_columns = F, 
    cluster_rows = T, 
    show_column_names = F, 
    show_row_names = F, 
    row_names_gp = gpar(fontsize = 7),
    col = colorRamp2(seq(-2,2,length=12), colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(12)), 
    column_split = metadat.ATAC_revisions$Label[match(colnames(mat.plot), metadat.ATAC_revisions$Sample)],
    bottom_annotation = ha_col_revisions,
    column_title = NULL,
    row_title = NULL,
    border = TRUE)
ht <- draw(hm)

# Write out plot data
dat <- mat.plot
write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Figure 5I.csv"), sep=",", quote=F, col.names=T, row.names=T)
