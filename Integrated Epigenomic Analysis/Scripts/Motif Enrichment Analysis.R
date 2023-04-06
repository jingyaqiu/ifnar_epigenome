
### Fig 2G: Motif enrichment in annotated REs ###

# Run HOMER, find enriched motifs
~/Dropbox/Minn/ifnar_epigenome/scripts/Motif enrichment analysis/motif_analysis.sh

###########################################
### Plot HOMER motif enrichment results ###
###########################################

rm(list=ls())

library(tidyverse)
library(RColorBrewer)

RE.names <- c("gained_TSS_promoters", "lost_TSS_promoters", "activated_enhancer_ATAC_inclusive", "deactivated_enhancer_ATAC_inclusive", "unchanged_REs")

denovo_motifs <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Motif Enrichment/denovo_motifs_combined.txt", sep="\t", stringsAsFactors=F, header=T)
known_motifs <- read.table("~/Dropbox/Minn/ifnar_epigenome/Processed Data/Integrated Epigenomic Analysis/Motif Enrichment/known_motifs_combined.txt", sep="\t", stringsAsFactors=F, header=T)
denovo_motifs <- denovo_motifs[denovo_motifs$name %in% RE.names, ]
known_motifs <- known_motifs[known_motifs$name %in% RE.names, ]

mat.filt.denovo <- filter_motifs(
    mat.motifs=denovo_motifs, 
    pval_threshold=1e-10, 
    max_n=3)
mat.filt.known <- filter_motifs(
    mat.motifs=known_motifs, 
    pval_threshold=1e-6, 
    max_n=3)

# Combined
mat.combined <- rbind(mat.filt.denovo[,c("Label", "ID", "pval", "OE", "Perc_target")], mat.filt.known[,c("Label", "ID", "pval", "OE", "Perc_target")])
mat.combined <- mat.combined %>%
    group_by(Label, ID) %>%
    slice_min(pval, n=1) %>%
    arrange(Label, pval)
mat.combined$Label <- factor(mat.combined$Label, levels=rev(unique(mat.combined$Label)))
mat.combined$ID <- factor(mat.combined$ID, levels=unique(mat.combined$ID))

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
    theme(axis.text.y = element_text(size=14, color="black"),
        axis.text.x = element_text(size=13.5, angle=30, hjust = 1, color="black"),
        axis.ticks = element_blank(),
        plot.margin=unit(c(1,1,1,1),"cm"),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave(fig, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Final Figures/Figure 2/Fig 2G Annotated REs Motif Enrichment.pdf"), width=11, height=3.5)

# dat <- data.frame(fig$data)[ ,c("Label", "ID", "pval")]
# write.table(dat, file=paste0("~/Dropbox/Minn/ifnar_epigenome/Figures Data/Source/Figure 2G.csv"), sep=",", quote=F, col.names=T, row.names=F)

#################
### Functions ###
#################

# Filter enriched motifs on p-threshold, rename for plotting
filter_motifs <- function(mat.motifs, pval_threshold=1e-12, max_n=5) {
    mat.motifs <- mat.motifs[sort.int(mat.motifs$pval, decreasing=F, index.return=T)$ix, ]

    # Motifs with no 
    motifs <- mat.motifs$Archetype_Motif
    names(motifs) <- mat.motifs$Similar_motif
    print("No archetype match: ")
    print(names(motifs[is.na(motifs)]))

    mat <- mat.motifs %>%
        mutate(name = dplyr::recode(mat.motifs$name, 
        "gained_TSS_promoters"="Activated promoters",
        "lost_TSS_promoters"="Deactivated promoters",
        "gained_other_promoters"="Activated promoters (non-protein-coding)",
        "lost_other_promoters"="Deactivated promoters (non-protein-coding)",
        "denovo_enhancer"="Strongly activated enhancers",
        "activated_enhancer"="Activated enhancers",
        "activated_enhancer_ATAC"="Activated enhancers",
        "activated_enhancer_ATAC_inclusive"="Activated enhancers",
        "deactivated_enhancer"="Deactivated enhancers",
        "deactivated_enhancer_ATAC_inclusive"="Deactivated enhancers",
        "unchanged_REs"="Constitutive REs",
        "random_peaks"="Random peaks (n=2000)")) %>%
    dplyr::filter(!is.na(Archetype_Motif)) %>%
    dplyr::filter(pval < pval_threshold)

    # Select most significantly enriched motif within each motif archetype in a RE set
    mat.filt <- mat %>%
        group_by(name, Archetype_Motif) %>%
        slice_min(pval, n=1)
    mat.filt$Label <- str_wrap(mat.filt$name, width=50) # 6,25
    mat.filt$Label <- factor(mat.filt$Label, levels=c("Activated enhancers", "Deactivated enhancers", "Activated promoters", "Deactivated promoters", "Constitutive REs", "Random peaks (n=2000)"))
    mat.filt <- mat.filt %>%
        arrange(Label, pval) %>%
        ungroup() %>%
        group_by(Label) %>%
        slice_min(order_by=pval, n=max_n)
    mat.filt$ID <- factor(mat.filt$Archetype_Motif, levels=rev(unique(mat.filt$Archetype_Motif)))
    mat.filt$OE <- mat.filt$Perc_target / mat.filt$Perc_background
    print(mat.filt %>% arrange(pval))
    return(mat.filt)
}

