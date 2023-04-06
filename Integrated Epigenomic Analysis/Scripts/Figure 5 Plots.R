

color_2 <- brewer.pal(11, "Spectral")[c(10,2,4)]

p_theme2 <- theme_bw(base_size = 12) +
    theme(axis.ticks = element_line(color = "black"),
        axis.text = element_text(color = "black", size = 10),
        legend.text = element_text(size = 8))

p_theme3 <- theme_bw(base_size = 12) +
    theme(legend.title=element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())

##############################################################################################
### Figure 5A: IFNAR sensitizes to PD1 survival through CD8, B2m, prf dependent mechanisms ###
##############################################################################################

library(survival)
library(survminer)
library(pzfx)

# read.prizm: read prizm exported data
read.prizm <- function(file){
    require(reshape2)
    dat <- read.table(file, header = T, sep="\t")
    # rownames(dat) <- dat[,1]
    dat <- dat[,-1]
    dat <- melt(dat, id.vars="Day", na.rm=T)
    colnames(dat) <- c("Time","Predictor","Event")
    return(dat)
}

use.pal <- c("royalblue", "orange", "firebrick", "grey50", "gold3", "skyblue3")

# ~/Dropbox/Xu Qiu IFN paper/V2/depletion master/depletion master.txt
fileName <- "~/Dropbox/Xu Qiu IFN paper/V2/depletion master/depletion master.txt"
survdat <- read.prizm(fileName)
Step1 <- survdat
# levels(Step1$Predictor) <- c("WT no_dep", "WTPD1 no_dep", "IFNARKO no_dep", "IFNARKOPD1 no_dep",
#                              "IFNARKO aNK", "IFNARKOPD1 aNK", "IFNARKO aCD8", "IFNARKOPD1 aCD8",
#                              "IFNARB2MKO no_dep", "IFNARB2MKOPD1 no_dep", "IFNARKO Prf", "IFNARKOPD1 Prf")
# Step2 <- separate(Step1, Predictor, c("Cell", "Treatment"), " ")
# Step2$Treatment <- factor(Step2$Treatment, levels = c("no_dep", "aNK", "aCD8", "Prf"))
# Step2$Cell <- factor(Step2$Cell, levels = c("WT", "WTPD1", "IFNARKO", "IFNARKOPD1", "IFNARB2MKO", "IFNARB2MKOPD1"))

res <- pairwise_survdiff(Surv(Time, Event) ~ Predictor, data = Step1,
                         p.adjust.method = "none")
# NT: p= 0.00462, PD1: p=0.00038
sfit <- survfit(Surv(Time, Event) ~ Cell + Treatment, data = Step2)
p_Ifnar.PD1.surv <- ggsurv2(sfit, facetby = "Treatment", facet.order = c("no_dep", "aNK", "aCD8", "Prf"), size.est = 1,
                            surv.col = use.pal, main = "",
                            strata.order = c("WT", "WTPD1", "IFNARKO", "IFNARKOPD1", "IFNARB2MKO", "IFNARB2MKOPD1")) + p_theme3
p_Ifnar.PD1.surv

################################################
### Fig 5B: Ifnar KO increases CD8 frequency ###
################################################

DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/cd8 percentage_pool.txt"), header=T, sep="\t")
colnames(DataSet) <- c("WT NT", "WT PD1", "IFNARKO NT", "IFNARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("WT NT", "IFNARKO NT", "WT PD1", "IFNARKO PD1"))
Step2 <- separate(Step1, Group, c("Cell", "Treatment"), " ")
Step2$Cell <- factor(Step2$Cell, levels = c("WT", "IFNARKO"))
Step2$Treatment <- factor(Step2$Treatment, levels = c("NT", "PD1"))
Step2 <- Step2[complete.cases(Step2),]

# Boxplots
p_CD8 <- ggplot(data = Step2, aes(x = Cell, y = Percentage, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("CD8 T cells/Total CD45") +
  facet_wrap(~ Treatment, nrow = 1) + 
  #labs(title = "IFNAR CD8") +
  p_theme2 +
  theme(legend.position = "none")
p_CD8

with(Step2[Step2$Treatment == "NT",], t.test(Percentage ~ Cell, alternative = "two.sided"))
with(Step2[Step2$Treatment == "PD1",], t.test(Percentage ~ Cell, alternative = "two.sided"))

# NT: p=0.0004197, p=8.879e-05
# Replicates: WT NT (n=22), IFNAR KO NT (n=25), WT PD1 (n=32), IFNAR KO PD1 (n=29)

#######################################################
### Figure 5C: IFNAR KO has enhanced MHC I and PDL1 ###
#######################################################

color_7 <- brewer.pal(11, "Spectral")[c(2,3,4,6,8,10,11)]
DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/24. ARKO tumor MHC/R499 ARKO MHC.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
colnames(DataSet) <- c("CTRL IgG", "CTRL", "IFNARKO IgG", "IFNARKO", "IFNGRKO IgG", "IFNGRKO")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("CTRL", "CTRL IgG", "IFNARKO", "IFNARKO IgG", "IFNGRKO", "IFNGRKO IgG"))
levels(Step1$Group) <- c("WT", "WT IgG", "IFNARKO", "IFNARKO IgG", "IFNGRKO", "IFNGRKO IgG")
Step2 <-  Step1
Step2 <- Step2[!is.na(Step2$Percentage),]

# Boxplots
p_MHC.R499 <- ggplot(data = Step2, aes(x = Group, y = Percentage, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_7) +
  ylab("MHC-I MFI") +
  p_theme2 +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p_MHC.R499

with(Step2[Step2$Group %in% c("WT", "IFNARKO"),], t.test(Percentage ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("WT", "IFNGRKO"),], t.test(Percentage ~ Group, alternative = "two.sided"))

# p= 0.001936, p = 0.2232
# Replicates: WT (n=4), WT IgG (n=5), IFNAR KO (n=5), IFNAR KO IgG (n=5), IFNGR KO (n=5), IFNGR KO IgG (n=5)

################################################################
### Figure 5D: Mixing IFNAR and R499 sensitizes tumor growth ###
################################################################

DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/27. IFNAR mixing/Mixing Day 15 Total volume.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
#colnames(DataSet) <- c("CTRL NT", "ARKO NT", "CTRL PD1", "ARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Volume")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("IFNAR.0_R499.100", "IFNAR.25_R499.75", "IFNAR.50_R499.50",
                                              "IFNAR.75_R499.25", "IFNAR.100_R499.0"))
Step2 <- Step1[complete.cases(Step1),]

bluePalette <- brewer.pal(9, "Blues")[c(9, 7, 6, 5, 3)]

# Boxplots
p_mix <- ggplot(data = Step2, aes(x = Group, y = Volume, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = bluePalette) +
  ylab("Total volume (cm^3)") +
  p_theme2 +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")
p_mix

with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.25_R499.75"),], t.test(Volume ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.50_R499.50"),], t.test(Volume ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.75_R499.25"),], t.test(Volume ~ Group, alternative = "two.sided"))
# pval: 0.005848; 0.0028; 0.0014

### Mixing IFNGR and R499 does not sensitize tumor growth

DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/27. IFNAR mixing/Mixing Day 15 Total volume_IFNGRKO.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
#colnames(DataSet) <- c("CTRL NT", "ARKO NT", "CTRL PD1", "ARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Volume")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("IFNGR.0_R499.100", "IFNGR.25_R499.75", "IFNGR.50_R499.50",
                                              "IFNGR.75_R499.25", "IFNGR.100_R499.0"))
Step2 <- Step1[complete.cases(Step1),]

bluePalette <- brewer.pal(9, "Blues")[c(9, 7, 6, 5, 3)]

# Boxplots
p_mix <- ggplot(data = Step2, aes(x = Group, y = Volume, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = bluePalette) +
  ylab("Total volume (cm^3)") +
  p_theme2 +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")
p_mix

with(Step2[Step2$Group %in% c("IFNGR.0_R499.100", "IFNGR.25_R499.75"),], t.test(Volume ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNGR.0_R499.100", "IFNGR.50_R499.50"),], t.test(Volume ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNGR.0_R499.100", "IFNGR.75_R499.25"),], t.test(Volume ~ Group, alternative = "two.sided"))
# pval: 0.8013; 0.3041; 0.0436

# Replicates: n=5 per group

###########################################################
### Figure 5E: CD8+ flow with different IFNAR KO mixing ###
###########################################################

### Flow shoing WT/KO mixing leads to increased CD8 infiltrate by mixing
DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/Ifnar mixing CD8 percent.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
#colnames(DataSet) <- c("CTRL NT", "ARKO NT", "CTRL PD1", "ARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percent")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
levels(Step1$Group) <- c("IFNAR.100_R499.0", "IFNAR.75_R499.25", "IFNAR.0_R499.100", "IFNAR.50_R499.50", "IFNAR.25_R499.75")
Step1$Group <- factor(Step1$Group, levels = c("IFNAR.0_R499.100", "IFNAR.25_R499.75", "IFNAR.50_R499.50",
                                              "IFNAR.75_R499.25", "IFNAR.100_R499.0"))
Step2 <- Step1[complete.cases(Step1),]

bluePalette <- brewer.pal(9, "Blues")[c(9, 7, 6, 5, 3)]

# Boxplots
p_mix.flow <- ggplot(data = Step2, aes(x = Group, y = Percent, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = bluePalette) +
  ylab("Percent CD8 T Cells") +
  p_theme2 +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  theme(legend.position = "none")
p_mix.flow

with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.25_R499.75"),], t.test(Percent ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.50_R499.50"),], t.test(Percent ~ Group, alternative = "two.sided"))
with(Step2[Step2$Group %in% c("IFNAR.0_R499.100", "IFNAR.75_R499.25"),], t.test(Percent ~ Group, alternative = "two.sided"))
# pval: 0.002975; 0.002571; 0.01012

# Replicates: n=5 per group

###############################################################
### Extended Figure 5B: Oas1 baseline expression B16/Res499 ###
###############################################################

### Res499 has higher Oas1 RNA expression base line

dat <- read.delim(file.path("~/Dropbox/Xu Qiu IFN paper/V2/4. B16 R499 baseline Oas1/B16 R499 baseline Oas1.txt"))
dat1 <- dat %>% filter(ISG %in% "Oas1g")
dat1$Cell <- factor(dat1$Cell, levels = c("B16", "R499"))

p_Oas1.R499base <- ggplot(data = dat1, aes(x = Cell, y = Expression, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("Relative Expression") +
  #labs(title = "R499 Oas1 expression") +
  p_theme2 +
  theme(legend.position = "none")

p_Oas1.R499base

with(dat[dat$ISG=="Oas1g", ],
     t.test(Expression ~ Cell, alternative = "two.sided"))
# p = 0.03013

write.table(dat1, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5B.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: n=5 per group

###############################################################
### Extended Figure 5D: Oas1 baseline expression TSA/Res237 ###
###############################################################

### Res237 have higher Oas1 in vitro compared to TSA
dat <- read.delim(file.path("~/Dropbox/Xu Qiu IFN paper/V2/5. TSA R237 baseline Oas1/TSA_R237 baseline Oas1.txt"))
dat1 <- dat %>% filter(ISG %in% "Oas1g")
dat1$Cell <- factor(dat1$Cell, levels = c("TSA", "R237"))

p_Oas1.R237base <- ggplot(data = dat1, aes(x = Cell, y = Expression, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("Relative Expression") +
  #labs(title = "R237 Oas1 expression") +
  p_theme2 +
  theme(legend.position = "none")
p_Oas1.R237base

with(dat[dat$ISG=="Oas1g", ],
     t.test(Expression ~ Cell, alternative = "two.sided"))
# pval: 5.394e-05

write.table(dat1, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5D.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: n=6 per group

#######################################################################################
### Extended Figure 5G: Res237 produces more IFNa/b after polyIC, Oas1 KO blunts it ###
#######################################################################################

# Res237 produces more IFNa and IFNb compared to TSA after pIC; Oas1 KO blunts it
dat <- read.delim(file.path("~/Dropbox/Xu Qiu IFN paper/V2/3. pIC TSA R237 Oas1 IFN/TSA_R237_Oas1KO scatter.txt"))
dat$Cell <- factor(dat$Cell, levels = c("TSA", "R237", "cl2", "cl4"))
# Add detection limit line
hline.data <- data.frame(dtl = c(2.38, 15.6), IFN = c("IFNa", "IFNb"))

p_R237.IFN <- ggplot(data = dat, aes(x = Cell, y = Conc, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  geom_hline(data = hline.data, aes(yintercept = dtl), linetype="dashed") + # add detection limit line
  facet_wrap(~IFN, scales = "free") +
  scale_color_manual(values = color_4) +
  ylab("Concentration (pg/ml)") +
  #labs(title = "R237 IFN concentration") +
  p_theme2 +
  theme(legend.position = "none")
p_R237.IFN 

with(subset(dat, subset = IFN %in% "IFNa" & Cell %in% c("TSA", "R237")),
     t.test(Conc ~ Cell, alternative = "two.sided"))   
with(subset(dat, subset = IFN %in% "IFNa" & Cell %in% c("R237", "cl2")),
     t.test(Conc ~ Cell, alternative = "two.sided"))
with(subset(dat, subset = IFN %in% "IFNa" & Cell %in% c("R237", "cl4")),
     t.test(Conc ~ Cell, alternative = "two.sided"))
with(subset(dat, subset = IFN %in% "IFNb" & Cell %in% c("TSA", "R237")),
     t.test(Conc ~ Cell, alternative = "two.sided"))   
with(subset(dat, subset = IFN %in% "IFNb" & Cell %in% c("R237", "cl2")),
     t.test(Conc ~ Cell, alternative = "two.sided"))
with(subset(dat, subset = IFN %in% "IFNb" & Cell %in% c("R237", "cl4")),
     t.test(Conc ~ Cell, alternative = "two.sided"))
# pval: 0.0006783; 0.0007384; 0.001683; 4.456e-12; 1.565e-06; 1.225e-05
# Replicates IFNa: TSA(n=8), Res237 (n=8), cl2 (n=7), cl4 (n=4)
# Replicates IFNa: TSA(n=6), Res237 (n=6), cl2 (n=5), cl4 (n=5)

write.table(dat, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5G.csv", sep=",", col.names=T, row.names=F, quote=F)

#####################################################
### Extended Data Figure 5H: STAT1 + IRF3 KO qPCR ###
#####################################################

### STAT1 downregulates Oas1a expression in Res499
path <- "siIRF"
dat <- read.delim(file.path("~/Dropbox/Xu Qiu IFN paper/V2/siIRF/master STAT1KO.txt"))
Step1 <- dat %>% filter(ISG %in% "OAS1A")

p_SKO <- ggplot(data = Step1, aes(x = KO, y = Expression, color = KO)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  #facet_wrap(~ISG, scales = "free") +
  scale_color_manual(values = color_2) +
  ylab("Relative Oas1 Expression") +
  p_theme2 +
  theme(legend.position = "none")
p_SKO 

with(subset(dat, subset = ISG %in% "OAS1A"),
     t.test(Expression ~ KO, alternative = "two.sided"))
with(subset(dat, subset = ISG %in% "OAS1G"),
     t.test(Expression ~ KO, alternative = "two.sided"))
# Oas1A pval: 6.107e-07; Oas1G pval: 5.357e-06
# Replicates: n=10 per group

write.table(dat, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5H.1.csv", sep=",", col.names=T, row.names=F, quote=F)

### siIRF3 downregulates Oas1a expression in STAT1 KO

dat <- read.delim(file.path("~/Dropbox/Xu Qiu IFN paper/V2/siIRF/master siIRF.txt"))
Step1 <- dat %>% filter(ISG %in% "Oas1A" & !(KD %in% "siIRF9"))

p_siIRF3 <- ggplot(data = Step1, aes(x = KD, y = Expression, color = KD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  #facet_wrap(~ISG, scales = "free") +
  scale_color_manual(values = color_4) +
  ylab("Relative OAS1 Expression") +
  p_theme2 +
  theme(legend.position = "none")
p_siIRF3

with(subset(dat, subset = KD %in% c("siC", "siIRF3") & ISG %in% "Oas1A"),
     t.test(Expression ~ KD, alternative = "two.sided"))
with(subset(dat, subset = KD %in% c("siC", "siIRF3") & ISG %in% "Oas1G"),
     t.test(Expression ~ KD, alternative = "two.sided"))
# Oas1A pval: 0.0001873; Oas1G pval: 0.003086

write.table(dat[dat$ISG == "Oas1A", ], file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5H.2.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: siControl (n=10), siIRF1 (n=7), siIRF3 (n=10), siIRF7 (n=7), siIRF9 (n=7)

###########################################################
### Extended Data Figure 6A: Oas1 KO sensitizes Res 237 ###
###########################################################

# Oas1KO sensitizes Res237 survival
fileName <- file.path("~/Dropbox/Xu Qiu IFN paper/V2/8. R237 Oas1 in vivo/Res237 OasKO survival_pool.txt")
survdat <- read.prizm(fileName)
levels(survdat$Predictor)
Step1 <- survdat
levels(Step1$Predictor) <- c("Cas9 NT", "Cas9 aCTLA4", "cl2 NT",
                             "cl2 aCTLA4", "cl4 NT", "cl4 aCTLA4")
Step2 <- separate(Step1, Predictor, c("Cell", "Treatment"), " ")
Step3.1 <- Step2 %>% filter(Step2$Treatment %in% "NT")
Step3.2 <- Step2 %>% filter(Step2$Treatment %in% "aCTLA4")

res.1 <- pairwise_survdiff(Surv(Time, Event) ~ Cell, data = Step3.1,
                           p.adjust.method = "none")
res.2 <- pairwise_survdiff(Surv(Time, Event) ~ Cell, data = Step3.2,
                           p.adjust.method = "none")
# pval: NT: cl2: 0.0009, cl4:0.0020; C4: cl2: 0.032, cl4: 0.034

sfit <- survfit(Surv(Time, Event) ~ Treatment + Cell, data = Step2)
p_R237Oas1KO.surv <- ggsurv2(sfit, 
    facetby = "Treatment", 
    facet.order = c("NT", "aCTLA4"), 
    ize.est = 1,
    surv.col = color_2, 
    main = "",
    reflevel = "Cas9") +
    p_theme3
p_R237Oas1KO.surv

write.table(Step2, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 6A.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: n=10 per group for NT, n=20 per group for aCTLA4

##############################################################
### Extended Data Figure 6B: Oas1 KO increases CD8 T-cells ###
##############################################################

### R237 Oas1KO increases CD8 T cells

DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/10. R237 Oas1 flow/Res237_Oas1 immune flow.txt"), header=T, sep="\t")
Step1 <- DataSet %>% filter(Cell %in% c("CD8", "Ki67_GzmB_PD1_Eomes"))
Step1 <- droplevels(Step1)
Step1$Genotype<- factor(Step1$Genotype, levels = c("Res237", "Oas1KO"))
Step2.1 <- Step1 %>% filter(Cell %in% "CD8")

# Boxplots
p_R237.CD8 <- ggplot(data = Step2.1, aes(x = Genotype, y = Percentage_100x, color = Genotype)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("CD8 T cells/Total CD45+") +
  #facet_wrap(~ Treatment, nrow = 1) + 
  #labs(title = "R237 Oas1 CD8") +
  p_theme2 +
  theme(legend.position = "none")
p_R237.CD8

with(Step2.1, t.test(Percentage_100x ~ Genotype, alternative = "two.sided"))   
# p-value = 1.411e-05

write.table(Step2.1, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 6B.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: Res237 (n=10), Oas1 KO (n=9)

#################################################################
### Extended Data Fig 7B: Ifng in CD8 increased with IFNAR KO ###
#################################################################

### ARKO leads to increased IFN gamma in CD8
DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/cd8ifng.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
colnames(DataSet) <- c("WT NT", "KO NT", "WT PD1", "KO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("WT NT", "KO NT", "WT PD1", "KO PD1"))
Step2 <- separate(Step1, Group, c("Cell", "Treatment"), " ")
Step2$Treatment <- as.factor(Step2$Treatment)
Step2$Cell <- factor(Step2$Cell, levels = c("WT", "KO"))
Step3 <- Step2[!is.na(Step2$Percentage),]

# Boxplots
p_CD8.gamma <- ggplot(data = Step3, aes(x = Cell, y = Percentage, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  facet_wrap(~ Treatment, nrow = 1) + 
  ylab("Percent of CD8 T cells") +
  # labs(title = "IFNg positive") +
  p_theme2
p_CD8.gamma

with(Step2[Step2$Treatment %in% "NT",], t.test(Percentage ~ Cell, alternative = "two.sided"))
# p= 0.01949
with(Step2[Step2$Treatment %in% "PD1",], t.test(Percentage ~ Cell, alternative = "two.sided"))
# p= 0.001479
with(Step2[Step2$Cell %in% "KO",], t.test(Percentage ~ Treatment, alternative = "two.sided"))
# p= 0.007295

write.table(Step2, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 7B.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: n=13 per group

#############################################################
### Extended Data Fig 5A: Intratumoral immune cells Ifnb1 ###
#############################################################

### Violin plot of IFNg, IFNg score and IFNb

umap.dat <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/20. ARKO scRNA/CD45_umap.txt"), header=T, sep="\t")
umap.dat$Cell <- factor(umap.dat$Cell, levels = c("CTRL", "KO"))
umap.dat$Treatment <- factor(umap.dat$Treatment, levels = c("NT", "PD1"))
umap.dat$condition <- factor(umap.dat$condition, levels = c("R499", "ARKO", "R499P", "ARKOP"))

# Remove outlier for violin plots (IFNg,  IFNG.GS, IFNb)
remove.outlier <- T
nout <- 2.5
violin.dat <- umap.dat %>% dplyr::select(condition, ifng.score, ifng, ifnb)
violin.dat <- gather(violin.dat, key = geneset, value = exp, 2:ncol(violin.dat))
#violin.dat$geneset <- as.factor(violin.dat$geneset)

if (remove.outlier) {
  violin.dat <- violin.dat %>% group_by(geneset, condition) %>%
    filter((exp < median(exp) + nout * IQR(exp)) & (exp > median(exp) - nout * IQR(exp))) %>%
    ungroup()
}

ifnb.violin.dat <- violin.dat %>% filter(geneset %in% "ifnb")
p_ifnb.violin <- ggplot(ifnb.violin.dat, aes(x = condition, y = exp, fill = condition)) +
  geom_violin(scale = "width", size = 0.3) +
  #geom_boxplot(width = 0.2, outlier.shape = NA, size = 0.5) +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color="black", size = 0.1) +
  scale_fill_manual(values = color_4) +
  #scale_y_continuous(trans='log2', labels = function(x) format(x, digits = 2, scientific = TRUE)) +
  ylab("ifnb expression") +
  p_theme2 +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(ifnb.violin.dat$condition))) +
  theme(legend.position="none")

with(umap.dat[which(umap.dat$condition == "R499" | umap.dat$condition == "ARKO"), ],
     wilcox.test(ifnb ~ condition))
with(umap.dat[which(umap.dat$condition == "R499P" | umap.dat$condition == "ARKOP"), ],
     wilcox.test(ifnb ~ condition))
# NT: p < 2.2e-16, p < 2.2e-16

write.table(data.frame(ifnb.violin.dat), file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 5A.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: R499 (n=7692 cells), R499P (n=8076 cells), ARKO (n=6758 cells), ARKOP (n=7334 cells); n=2 biological replicates per group
# n=29,584 cells total

#######################################################################
### Extended Data Fig 6E: ARKO increases percent of Perforin in CD8 ###
#######################################################################

DataSet <- read.table("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/4-20-17 IFNAR PD1 CD8 perforin.txt", header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
colnames(DataSet) <- c("WT NT", "WT PD1", "IFNARKO NT", "IFNARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("WT NT", "IFNARKO NT", "WT PD1", "IFNARKO PD1"))
Step2 <- separate(Step1, Group, c("Cell", "Treatment"), " ")
Step2$Cell <- factor(Step2$Cell, levels = c("WT", "IFNARKO"))
Step2$Treatment <- factor(Step2$Treatment, levels = c("NT", "PD1"))

# Boxplots
p_perf <- ggplot(data = Step2, aes(x = Cell, y = Percentage, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("Percent of Perforin+ CD8 T Cells") +
  facet_wrap(~ Treatment, nrow = 1) + 
  p_theme2 +
  theme(legend.position = "none")

p_perf

with(Step2[Step2$Treatment == "NT",], t.test(Percentage ~ Cell, alternative = "two.sided"))
with(Step2[Step2$Treatment == "PD1",], t.test(Percentage ~ Cell, alternative = "two.sided"))
# NT: p= 0.08347, p=0.01871

##################################################
### Extended Data Fig 6F: B16 IFNAR KO control ###
##################################################

survdat <- read.csv("~/Dropbox/Xu Qiu IFN paper/V2/B16 IFNAR KO control/data.csv", header = T)
Step1 <- survdat
Step1$Treatment <- as.factor(Step1$Treatment)
levels(Step1$Treatment)
Step1$Genotype <- as.factor(Step1$Genotype)
Step1 <- Step1 %>% dplyr::filter(Genotype %in% c("B16", "B16_IFNAR")) %>% droplevels()

levels(Step1$Genotype)
levels(Step1$Genotype) <- c("WT", "KO")

Step2.1 <- Step1 %>% dplyr::filter(Treatment %in% "PD1")
res.1 <- pairwise_survdiff(Surv(Days, Survival) ~ Genotype, data = Step2.1,
                           p.adjust.method = "none")

Step2.2 <- StExtep1 %>% dplyr::filter(Treatment %in% "NT")
res.2 <- pairwise_survdiff(Surv(Days, Survival) ~ Genotype, data = Step2.2,
                           p.adjust.method = "none")

sfit <- survfit(Surv(Days, Survival) ~ Genotype + Treatment, data = Step1)
p_B16.surv <- ggsurv2(sfit, facetby = "Treatment", facet.order = c("NT", "PD1"), size.est = 1,
                      surv.col = color_2, main = "",
                      strata.order = c("WT", "KO")) + p_theme3
p_B16.surv

# Replicates: PD1 n=10; NT n=5

######################################################
### Extended Data Fig 6G: CT26 IFNAR KO sensitizes ###
######################################################

survdat <- read.csv("~/Dropbox/Xu Qiu IFN paper/V2/CT26 IFNAR KO/data.csv", header = T)
Step1 <- survdat
Step1$Treatment <- as.factor(Step1$Treatment)
levels(Step1$Treatment)
Step1$Genotype <- as.factor(Step1$Genotype)
levels(Step1$Genotype)
levels(Step1$Genotype) <- c("WT", "KO")

Step2.1 <- Step1 %>% dplyr::filter(Treatment %in% "NT")
Step2.2 <- Step1 %>% dplyr::filter(Treatment %in% "PD1")

res.1 <- pairwise_survdiff(Surv(Days, Survival) ~ Genotype, data = Step2.1,
                           p.adjust.method = "none")
res.2 <- pairwise_survdiff(Surv(Days, Survival) ~ Genotype, data = Step2.2,
                           p.adjust.method = "none")
# NT: p= 0.023, PD1: p=0.03

sfit <- survfit(Surv(Days, Survival) ~ Genotype + Treatment, data = Step1)
p_CT26KO.surv <- ggsurv2(sfit, facetby = "Treatment", facet.order = c("NT", "PD1"), size.est = 1,
                         surv.col = color_2, main = "",
                         strata.order = c("WT", "KO")) + p_theme3 + ggtitle("CT26")
p_CT26KO.surv

# Replicates: PD1 n=10; NT n=5

############################################################
### Extended Data Fig 7A: Intratumoral immune cells Ifng ###
############################################################

ifng.violin.dat <- violin.dat %>% filter(geneset %in% "ifng")
p_ifng.violin <- ggplot(ifng.violin.dat, aes(x = condition, y = exp, fill = condition)) +
  geom_violin(scale = "width", size = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 0.5) +
  # stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  # geom="pointrange", color="black") +
  scale_fill_manual(values = color_4) +
  #scale_y_continuous(trans='log2', labels = function(x) format(x, digits = 2, scientific = TRUE)) +
  ylab("ifng expression") +
  p_theme2 +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(ifng.violin.dat$condition))) +
  theme(legend.position="none")

gscore.violin.dat <- violin.dat %>% filter(geneset %in% "ifng.score")
p_ifng.score.violin <- ggplot(gscore.violin.dat, aes(x = condition, y = exp, fill = condition)) +
  geom_violin(scale = "width", size = 0.3) +
  geom_boxplot(width = 0.2, outlier.shape = NA, size = 0.5) +
  #stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #            geom="pointrange", color="black", size = 1) +
  scale_fill_manual(values = color_4) +
  ylab("ifng.score") +
  # scale_y_continuous(trans='log2') +
  p_theme2 +
  coord_flip()+
  scale_x_discrete(limits = rev(levels(gscore.violin.dat$condition))) +
  theme(legend.position="none")

# Wilcoxon test
with(umap.dat[which(umap.dat$condition == "R499" | umap.dat$condition == "ARKO"), ],
     wilcox.test(ifng.score ~ condition))
with(umap.dat[which(umap.dat$condition == "R499P" | umap.dat$condition == "ARKOP"), ],
     wilcox.test(ifng.score ~ condition))
# NT: p < 2.2e-16; treated: p < 2.2e-16

with(umap.dat[which(umap.dat$condition == "R499" | umap.dat$condition == "ARKO"), ],
     wilcox.test(ifng ~ condition))
with(umap.dat[which(umap.dat$condition == "R499P" | umap.dat$condition == "ARKOP"), ],
     wilcox.test(ifng ~ condition))
# NT: p = 0.0001108, p < 2.2e-16

# Replicates: R499 (n=7692 cells), R499P (n=8076 cells), ARKO (n=6758 cells), ARKOP (n=7334 cells); n=2 biological replicates per group
# n=29,584 cells total

write.table(data.frame(ifng.violin.dat), file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 7A Ifng.csv", sep=",", col.names=T, row.names=F, quote=F)
write.table(data.frame(gscore.violin.dat), file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 7A IFNG.GS.csv", sep=",", col.names=T, row.names=F, quote=F)

#############################################
### Extended Data Fig 7F: CD11c+ DCs flow ###
#############################################

# Ifnar immune flow CD86 on cDC1 DCs

DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/CD11c CD86 MFI pool.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
colnames(DataSet) <- c("CTRL NT", "ARKO NT", "CTRL PD1", "ARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("CTRL NT", "ARKO NT", "CTRL PD1", "ARKO PD1"))
Step2 <- separate(Step1, Group, c("Cell", "Treatment"), " ")
Step2$Cell <- factor(Step2$Cell, levels = c("CTRL", "ARKO"))
Step2$Treatment <- factor(Step2$Treatment, levels = c("NT", "PD1"))
Step2 <- Step2[complete.cases(Step2),]
Step3 <- Step2 %>% filter(Treatment == "NT")

# Boxplots
p_DC.act <- ggplot(data = Step3, aes(x = Cell, y = Percentage, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("CD86 MFI") +
  facet_wrap(~ Treatment, nrow = 1) + 
  #labs(title = "IFNAR CD8") +
  p_theme2 +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p_DC.act

with(Step2[Step2$Treatment == "NT",], t.test(Percentage ~ Cell, alternative = "two.sided"))
with(Step2[Step2$Treatment == "PD1",], t.test(Percentage ~ Cell, alternative = "two.sided"))
# NT: p= 0.003179, PD1: p=NS

write.table(Step2, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 7F.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: n=9 per group

################################################
### Extended Data Fig 7G: Batf3 -/- survival ###
################################################

# Ifnar Batf3 dependent survival
fileName <- file.path("~/Dropbox/Xu Qiu IFN paper/V2/26. IFNAR in vivo Batf3/Ifnar in Batf3.txt")
survdat <- read.prizm(fileName)
Step1 <- survdat
levels(Step1$Predictor)
levels(Step1$Predictor) <- c("WT NT", "WT PD1", "Batf3-/- NT", "Batf3-/- PD1")
Step2 <- separate(Step1, Predictor, c("Genotype", "Treatment"), " ")
Step2$Treatment <- factor(Step2$Treatment, levels= c("NT", "PD1"))
Step2$Genotype <- factor(Step2$Genotype, levels = c("WT", "Batf3-/-"))
Step3 <- Step2 %>% filter(Treatment %in% c("NT", "PD1"))

res <- pairwise_survdiff(Surv(Time, Event) ~ Genotype + Treatment, data = Step3,
                         p.adjust.method = "none")
# NT: p= 0.1029, PD1: p=  0.0065
sfit <- survfit(Surv(Time, Event) ~ Predictor, data = Step1)
p_Batf3.surv <- ggsurv2(sfit, size.est = 1,
                        surv.col = use.pal, main = "",
                        strata.order = c("WT NT", "WT PD1", "Batf3-/- NT", "Batf3-/- PD1")) + p_theme3
p_Batf3.surv

write.table(Step2, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 7G.csv", sep=",", col.names=T, row.names=F, quote=F)

############################################################
### Extended Data Fig 8D: Effector-like CD8 T-cells flow ###
############################################################

### Flow showing increased Teff-like
DataSet <- read.table(file.path("~/Dropbox/Xu Qiu IFN paper/V2/16. IFNAR immune flow/Cx3cr1+CD62L- Teff.txt"), header=T, sep="\t",
                      fileEncoding="UTF-8-BOM")
colnames(DataSet) <- c("CTRL NT", "ARKO NT", "CTRL PD1","ARKO PD1")
Step1 <- gather(DataSet, key= "Group", value= "Percentage")
Step1$Group <- as.factor(Step1$Group)
levels(Step1$Group)
Step1$Group <- factor(Step1$Group, levels = c("ARKO NT", "ARKO PD1", "CTRL NT", "CTRL PD1"))
Step2 <- separate(Step1, Group, c("Cell", "Treatment"), " ")
Step2$Cell <- factor(Step2$Cell, levels = c("CTRL", "ARKO"))
Step2$Treatment <- factor(Step2$Treatment, levels = c("NT", "PD1"))
Step2 <- Step2[complete.cases(Step2),]

# Boxplots
p_Teff <- ggplot(data = Step2, aes(x = Cell, y = Percentage, color = Cell)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(w = 0.2), size = 1.5) +
  scale_color_manual(values = color_2) +
  ylab("Percent of Tim3- CD8 T cells") +
  facet_wrap(~ Treatment, nrow = 1) + 
  #labs(title = "IFNAR CD8") +
  p_theme2  +
  theme(legend.position="none")

p_Teff

with(Step2[Step2$Treatment == "NT",], t.test(Percentage ~ Cell, alternative = "two.sided"))
with(Step2[Step2$Treatment == "PD1",], t.test(Percentage ~ Cell, alternative = "two.sided"))
# NT: 0.001485, p=  0.001733

write.table(Step2, file = "~/Dropbox/Minn/ifnar_epigenome/Supplementary Data Files/Source/Extended Data Fig 8D.csv", sep=",", col.names=T, row.names=F, quote=F)

# Replicates: WT NT (n=8), IFNAR KO NT (n=10), WT PD1 (n=10), IFNAR KO PD1 (n=10)

