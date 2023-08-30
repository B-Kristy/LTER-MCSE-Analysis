rm(list=ls())
# MAC OS
knitr::opts_knit$set(root.dir = "~/Dropbox/LTER_MCSE/data/seq_data")

# Load Required Libraries
library(bbmle)
library(cowplot)
library(data.table)
library(DESeq2)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(hillR)
library(indicspecies)
library(knitr)
library(lazyeval)
library(metacoder)
library(multcompView)
library(pander)
library(permute)
library(phyloseq)
library(plyr)
library(lattice)
library(magrittr)
library(microbiome)
library(MetBrewer)
library(nlme)
library(parallel)
library(reshape2)
library(tidyverse)
library(vegan)
load("phylo_filtered.RData")
load("phylo_rarefied.RData")

## Figure 1: Alpha-Diversity ----------------------------

# Subset only the in-situ soil core samples
bact_mcse <- subset_samples(bact.rfy, incubation == "in-situ")

# Calculate Shannon Diversity
otu_table <- t(as(bact_mcse@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

shan <- hill_taxa(otu_table, q = 1)
shan <- merge(bact_mcse@sam_data, shan, by = 0)

colnames(shan)[colnames(shan) == "y"] ="shan"

# Seasonal Variation - Statistics
mod.aov <- aov(shan ~ treatment + sample_date, data = shan)
mod <- lm(shan ~ sample_date + 0 , data = shan)

# Summer Plot
shan_summer <- shan %>%
               subset(sample_date == "summer")

# Add Tukey's Post-Hoc Comparisons
mod.aov <- aov(shan ~ treatment, data = shan_summer)
mod <- lm(shan ~ treatment, data = shan_summer)
mod <- lm(shan ~ treatment + 0, data = shan_summer)
mod <- lm(shan ~ system, data = shan_summer)
mod <- lm(shan ~ system + 0, data = shan_summer)

tHSD <- TukeyHSD(mod.aov, ordered = FALSE, conf.level = 0.95)

letters.df <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = shan_summer))$treatment[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$treatment <- rownames(letters.df) 
letters.df$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)


# Reorder Treatments
shan_summer$treatment <- factor(shan_summer$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_summer$system <- factor(shan_summer$system, levels = c("annual", "perennial", "forest"))

shan_summer_p <- shan_summer %>%
                  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
                  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
                  geom_jitter(width = 0.1) +
                  geom_text(data = letters.df, aes(x=treatment, y=placement, label=Letter), 
                  size =5.5, color="black", fontface="bold") +
                  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
                  ylab("Shannon Diversity") +
                  xlab("") +
                  ggtitle("Summer Samples") +
                  theme_bw() +
                  theme(axis.title.x = element_text(size = 14), 
                        title = element_text(size = 14),
                        axis.title.y = element_text(size = 14), 
                        axis.text.x = element_text(size = 12, color = "black"), 
                        axis.text.y = element_text(size = 12, color = "black"), 
                        panel.border = element_rect(colour = "black", fill=NA, size=2),
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size = 14), 
                        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Fall plot
shan_fall <- shan %>%
  subset(sample_date == "fall")

# Add Tukey's Post-Hoc Comparisons
mod.aov <- aov(shan ~ treatment, data = shan_fall)
mod <- lm(shan ~ treatment + 0, data = shan_fall)
mod <- lm(shan ~ treatment, data = shan_fall)

tHSD <- TukeyHSD(mod.aov, ordered = FALSE, conf.level = 0.95)

letters.df <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = shan_fall))$treatment[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$treatment <- rownames(letters.df) 
letters.df$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)


# Reorder Treatments
shan_fall$treatment <- factor(shan_fall$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_fall$system <- factor(shan_fall$system, levels = c("annual", "perennial", "forest"))

shan_fall_p <- shan_fall %>%
  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_text(data = letters.df, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Shannon Diversity") +
  xlab("") +
  ggtitle("Fall Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

Fig_1 <- ggarrange(
        shan_summer_p, 
        shan_fall_p,
        labels = c("A", "B"),
        ncol = 1,
        nrow = 2,
        common.legend = T,
        legend = "right"
      )

ggsave("Fig_1.png", plot = Fig_1, dpi = 400, w = 10, h = 7)

mod <- lm(n.fix ~ shan, data = shan)

## Figure 2: Beta-diversity ordination ----------------------------------

# Ordination
set.seed(37920)
bact.mcse.rel <- transform_sample_counts(bact_mcse, function(x) x / sum(x) )
bray.NMDS.mcse <- ordinate(bact.mcse.rel, method = "NMDS", distance = "bray") # Stress = 0.1731933

# Beta Diversity Plot
bray.points <- as.data.frame(bray.NMDS.mcse$points)
bray.points <- rownames_to_column(bray.points, "id")

sam_data <- as(sample_data(bact_mcse), "data.frame")
sam_data <- rownames_to_column(sam_data, "id")

bray.points.m <- merge(bray.points, sam_data, by = "id")

bray.points.mean <- bray.points.m %>%
  group_by(treatment, sample_date) %>% 
  summarise(Mean_MDS1 = mean(MDS1), 
            SD_MDS1 = sd(MDS1), 
            N = n(),
            SE_MDS1 = SD_MDS1/sqrt(N), 
            Mean_MDS2 = mean(MDS2), 
            SD_MDS2 = sd(MDS2), 
            N = n(),
            SE_MDS2 = SD_MDS2/sqrt(N), 
  ) 

bray.points.mean$treatment <- factor(bray.points.mean$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))

# PLOT ordination
div_plot <- bray.points.mean %>%
  ggplot(aes(x=Mean_MDS1,y=Mean_MDS2)) +
  geom_point(aes(fill = treatment, shape = sample_date), size = 6.75) +
  geom_errorbar(aes(ymax = Mean_MDS2 + SE_MDS2, ymin = Mean_MDS2 - SE_MDS2, width = .001))+ 
  geom_errorbarh(aes(xmax = Mean_MDS1 + SE_MDS1, xmin = Mean_MDS1 - SE_MDS1, height = .001)) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = 0.12, y = 0.07, label = "stress = 0.17") +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  scale_shape_manual(values =c(21, 22), name = "Sampling Date", 
                     labels = c("Fall", "Summer")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21, 22)))) +
  theme_bw() +
  labs(x = "NMDS1", y = "NMDS2") +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 10, color = "black"))  + 
  theme(axis.text.y = element_text(size = 10, color = "black"))  + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(legend.text = element_text(size = 10)) + 
  theme(legend.position="right", legend.title=element_text(size=10))+
  theme(legend.position="right")

# permANOVA
set.seed(37920)
bray.NMDS.mcse <- phyloseq::distance(bact.mcse.rel,"bray")
set.seed(37920)
permanova <- adonis(bray.NMDS.mcse~sample_data(bact.mcse.rel)$treatment *
                      sample_data(bact.mcse.rel)$sample_date, permutations = 9999, method = "bray")
permanova["aov.tab"]

# GGsave
ggsave("Fig_2.png", plot = div_plot, dpi = 400, w = 10, h = 7)


## Figure 3: Differential Abundance -------------------------------------------------------

# Summer plot 
bact_mcse <- subset_samples(bact.p_filtered, incubation == "in-situ")
bact_mcse <- subset_samples(bact.p_filtered, sample_date == "summer")

## DIFFERENTIAL ABUNDANCE: FORESTED SITES VS. MCSE SITES
sample_data(bact_mcse)$agro <- ifelse(bact_mcse@sam_data$treatment == "SF" |
                                        bact_mcse@sam_data$treatment == "DF" |
                                        bact_mcse@sam_data$treatment == "CF" , "Forest", "Agro")

diagdds <- phyloseq_to_deseq2(bact_mcse, ~ agro) 
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")

results <- results(diagdds, alpha = 0.01, contrast = c("agro", "Forest","Agro"))

sigtab_results <- results[which(results$padj < 0.01), ]

sigtab_results_TAX <- cbind(as(sigtab_results, "data.frame"), as(tax_table(bact_mcse)[row.names(sigtab_results), ], "matrix"))

sigtab_results_TAX_genra <- subset(sigtab_results_TAX, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Class <- factor(as.character(sigtab_results_TAX_genra$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Genus <- factor(as.character(sigtab_results_TAX_genra$Genus), levels=names(x))


summer_p <- sigtab_results_TAX_genra %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  theme_classic() +
  geom_point(size = 4.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_text(x=37.5, y=1.5, label="Forest-Enriched", size = 5) +
  geom_text(x=37.5, y=-1.5, label="MCSE-Enriched", size = 5) +
  ggtitle("Summer Samples") +
  scale_fill_manual(values=met.brewer("Egypt", 14)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =14), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15)
  ) 
summary(sigtab_results)

# Fall Plot 
bact_mcse <- subset_samples(bact.p_filtered, incubation == "in-situ")
bact_mcse <- subset_samples(bact.p_filtered, sample_date == "fall")

## DIFFERENTIAL ABUNDANCE: FORESTED SITES VS. MCSE SITES
sample_data(bact_mcse)$agro <- ifelse(bact_mcse@sam_data$treatment == "SF" |
                                        bact_mcse@sam_data$treatment == "DF" |
                                        bact_mcse@sam_data$treatment == "CF" , "Forest", "Agro")

diagdds <- phyloseq_to_deseq2(bact_mcse, ~ agro) 
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")

results <- results(diagdds, alpha = 0.01, contrast = c("agro", "Forest","Agro"))

sigtab_results <- results[which(results$padj < 0.01), ]

sigtab_results_TAX <- cbind(as(sigtab_results, "data.frame"), as(tax_table(bact_mcse)[row.names(sigtab_results), ], "matrix"))

sigtab_results_TAX_genra <- subset(sigtab_results_TAX, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Class <- factor(as.character(sigtab_results_TAX_genra$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Genus <- factor(as.character(sigtab_results_TAX_genra$Genus), levels=names(x))


fall_p <- sigtab_results_TAX_genra %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  theme_classic() +
  geom_point(size = 4.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_text(x=5.0, y=-0.25, label="MCSE-Enriched", size = 5) +
  ggtitle("Fall Samples") +
  scale_fill_manual(values=met.brewer("Egypt", 14)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =14), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15)
  ) 
summary(sigtab_results)

# Arrange summer and fall plots together
Fig_3 <- ggarrange(
  summer_p, 
  fall_p,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)

ggsave("Fig_3.png", plot = Fig_3, dpi = 400, w = 10, h = 12)
## Figure 4: Fumigation on Diazotroph Alpha-diversity and N-fixation ----------------------
# Subset only fumigated samples
bact_fum <- subset_samples(bact.rfy, incubation != "in-situ")

# Calculate Shannon Diversity
otu_table <- t(as(bact_fum@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

shan <- hill_taxa(otu_table, q = 1)
shan <- merge(bact_fum@sam_data, shan, by = 0)

colnames(shan)[colnames(shan) == "y"] ="shan"

# shannon diversity plot
shan$incubation <- as.numeric(shan$incubation)
shan$treatment <- factor(shan$treatment, levels = c("T1", "T3", "T4", "T5",
                                                            "T7", "SF", "DF", "CF"))
shan$system <- factor(shan$system, levels = c("annual", "perennial", "forest"))

shan_plot <- shan %>%
            ggplot(aes(x=incubation, y = shan, color = treatment)) +
            geom_point(size = 5.5) +
            facet_grid(~system) +
            geom_smooth(method = 'lm') +
            scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                                        '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                               labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                          "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                          "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
            ylab("Shannon Diversity") +
            xlab("Fumigation Time (Hr)") +
            theme_bw() +
            theme(axis.title.x = element_text(size = 14), 
                  axis.title.y = element_text(size = 14), 
                  strip.text.x = element_text(size = 14),
                  axis.text.x = element_text(size = 12, color = "black"), 
                  axis.text.y = element_text(size = 12, color = "black"), 
                  panel.border = element_rect(colour = "black", fill=NA, size=2),
                  legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 14), 
                  plot.margin=unit(c(.5,1,.5,.5),"cm"))

mod <- aov(shan ~ incubation*system, data = shan)
mod <- aov(shan ~ incubation*treatment, data = shan)
mod <- lm(shan ~ incubation*system, data = shan)

# Partial correlation test
library(ppcor)
shan_pcor <- shan %>%
  filter(doc != 'NA')

shan_pcor <- shan_pcor %>%
  filter(don != 'NA')

shan_pcor <- shan_pcor %>%
  filter(n.fix != 'NA')
# Partial correlation between shannon diversity and fumigation given "doc" and "don"
pcor.test(shan_pcor$shan, shan_pcor$incubation, shan_pcor[,c("doc", "don")])
pcor.test(shan_pcor$shan, shan_pcor$n.fix, shan_pcor[,c("doc", "don")])



mod <- aov(n.fix ~ shan * system, data =shan)


# Nitrogen fixation potential
shan <- shan %>%
  filter(n.fix != 'NA')

n.fix_plot <- shan %>%
  ggplot(aes(x=incubation, y = n.fix, color = treatment)) +
  geom_point(size = 5.5) +
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Potential (\u00b5g/g soil)") +
  xlab("Fumigation Time (Hr)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

mod <- aov(n.fix ~ incubation*treatment, data = shan)
mod <- lm(n.fix ~ incubation, data = shan)

Fig_4 <- ggarrange(
  shan_plot, 
  n.fix_plot,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)


ggsave("Fig_4.png", plot = Fig_4, dpi = 400, w = 16, h = 14)

## Figure 5: Fumigation Dataset Ordination ---------------------------
set.seed(37920)
bact.fum.rel <- transform_sample_counts(bact_fum, function(x) x / sum(x) )
bray.NMDS.fum <- ordinate(bact.fum.rel, method = "NMDS", distance = "bray")

env <- as.data.frame(sample_data(bact_fum))
set.seed(37920)
ord.fit.sig <- envfit(bray.NMDS.fum ~ env$treatment + env$incubation + 
                        env$sample_date + env$n.fix + env$mbc + env$doc + env$don , 
                      perm = 9999, na.rm = TRUE)

data.scores <- as.data.frame(scores(bray.NMDS.fum)$sites)
data.scores$treatment <- env$treatment
data.scores$incubation <- env$incubation
data.scores$sample_date <- env$sample_date
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
#data.scores$mbn <- env$mbn (Not significant)
data.scores$doc <- env$doc
data.scores$don <- env$don

# permANOVA
set.seed(37920)
bray.NMDS.fum <- phyloseq::distance(bact.fum.rel,"bray")
set.seed(37920)
permanova <- adonis(bray.NMDS.fum~sample_data(bact.fum.rel)$treatment * 
                                   sample_data(bact.fum.rel)$incubation, 
                    permutations = 9999, method = "bray")
permanova["aov.tab"]

en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors")) * ordiArrowMul(ord.fit.sig)
en_coord_cont$name <- c("N-Fixation", "MBC", "DOC", "DON")

data.scores$treatment <- factor(data.scores$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
Fig_5 <- data.scores %>%
          ggplot(aes(x=NMDS1,y=NMDS2)) +
          geom_point(aes(fill = treatment, shape = incubation), size = 5.5) +
          geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
          geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
          geom_text(x = 1.2, y = 1.3, label = "stress = 0.19") +
          scale_fill_manual(values=c('#EEAC4E', '#D7785D','#9C7C50', '#577753','#1A4C76', 
                                     '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                            labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",  
                                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
          scale_shape_manual(values =c(21, 22, 23, 24), name = "Fumigation Time (Hr)", 
                             labels = c("0 Hours", "24 Hours", "72 Hours", "Reference")) +
          geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
                       data = en_coord_cont, size =1.25, alpha = 0.7, colour = "black") +
          geom_label(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "black", 
                     fontface = "bold", size = 5.5, label = en_coord_cont$name, fill = "white") + 
          guides(fill = guide_legend(override.aes = list(shape=c(21, 22, 23,24)))) +
          theme_bw() +
          theme(axis.title.x = element_text(size = 15)) + 
          theme(axis.title.y = element_text(size = 15)) + 
          theme(axis.text.x = element_text(size = 12, color = "black"))  + 
          theme(axis.text.y = element_text(size = 12, color = "black"))  + 
          theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
          theme(legend.text = element_text(size = 12)) +
          theme(legend.position="right", legend.title=element_text(size=14))+
          theme(legend.position="right")

ggsave("Fig_5.png", plot = Fig_5, dpi = 400, w = 10, h = 7)

## Figure 6: Biodiv vs. Ecosystem Function ----------------------------------
# Calculate Shannon Diversity
otu_table <- t(as(bact_fum@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

shan <- hill_taxa(otu_table, q = 1)
shan <- merge(bact_fum@sam_data, shan, by = 0)

colnames(shan)[colnames(shan) == "y"] ="shan"

# Remove na nfixation
shan <- shan %>%
  filter(n.fix != 'NA')

shan$treatment <- factor(shan$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan$system <- factor(shan$system, levels = c("annual", "perennial", "forest"))

Fig_6 <- shan %>%
  ggplot(aes(x=shan, y = n.fix, col = treatment)) +
  geom_point(aes(shape = incubation), size = 4.5) +  
  facet_grid(~system, scales = 'free') +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                     labels = c(
                       "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  xlab("Shannon Diversity") + 
  ylab("Nitrogen Fixation Potential (\u00b5g/g)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

ggsave("Fig_6.png", plot = Fig_6, dpi = 400, w = 12, h = 7)

# Effects parameterization for biodiversity-ecosystem function relationships
mod <- lm(n.fix ~ shan * system, data = shan)
emtrends(mod, ~ system, var = "shan")
mod <- lm(n.fix ~ shan*treatment, data = shan)
emtrends(mod, ~ treatment, var = "shan")








## Supplementary Figure 1: Rarefaction Plot --------------------------
# Load  phyloseq object
load("phylo_filtered.RData")

#Obtain OTU Table
bact.tab<-t(as(otu_table(bact.p_filtered), "matrix"))

# Perform Rarefaction with rarecurve() function from vegan
vegan_rarefaction <- rarecurve(bact.tab, step = 100, sample = raremax, col = "blue", cex = 0.6, tidy = TRUE)

# merge with metadata 
sam_data <- as(sample_data(bact.p_filtered), "data.frame")
sam_data <- rownames_to_column(sam_data, "Site")

vegan_rarefaction_m <- merge(vegan_rarefaction, sam_data, by = "Site")

# Reorder Treatment Factors
vegan_rarefaction_m$treatment <- factor(vegan_rarefaction_m$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))


# Plot Rarefaction Curve
Fig_S1 <- vegan_rarefaction_m %>%
  ggplot(aes(x=Sample, y = Species, col = treatment)) +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592','#577B96','#B0381D'), 
                     name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat", "T4: Bio-Based Corn/Soy/Wheat", 
                                "T5: Perennial Poplar", "T7: Early Successional Grassland", 
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest" )) +
  geom_point() +
  xlab("Sequencing Depth") +
  ylab("Diazotroph Species Richness") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Save High Resolution Plot for Publication
ggsave("Fig_S1.png", plot = Fig_S1, dpi = 400, w = 10, h =5)



## Supplementary Figure 2: Nitrogen Fixation Potentials --------------------------------
# Remove NAs
shan <- shan %>%
  subset(n.fix != 'NA')

# Seasonal Variation - Statistics
mod.aov <- aov(n.fix ~ treatment + sample_date, data = shan)
mod <- lm(n.fix ~ sample_date + 0 , data = shan)

cor.test(shan$shan, shan$n.fix, method = "spearman")

# Summer plot
shan_summer <- shan %>%
  subset(sample_date == "summer")

# Add Tukey's Post-Hoc Comparisons
mod.aov <- aov(n.fix ~ treatment, data = shan_summer)

tHSD <- TukeyHSD(mod.aov, ordered = FALSE, conf.level = 0.95)

letters.df <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = shan_summer))$treatment[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$treatment <- rownames(letters.df) 
letters.df$placement <- c(11, 11, 11, 11, 11, 11, 11, 11)


# Reorder Treatments
shan_summer$treatment <- factor(shan_summer$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_summer$system <- factor(shan_summer$system, levels = c("annual", "perennial", "forest"))

shan_summer_p <- shan_summer %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_text(data = letters.df, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Potential (\u00b5g/g soil)") +
  xlab("") +
  ggtitle("Summer Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Fall plot
shan_fall <- shan %>%
  subset(sample_date == "fall")

# Add Tukey's Post-Hoc Comparisons
mod.aov <- aov(n.fix ~ treatment, data = shan_fall)

tHSD <- TukeyHSD(mod.aov, ordered = FALSE, conf.level = 0.95)

letters.df <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = shan_fall))$treatment[,4])$Letters)
colnames(letters.df)[1] <- "Letter" #Reassign column name
letters.df$treatment <- rownames(letters.df) 
letters.df$placement <- c(18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5, 18.5)


# Reorder Treatments
shan_fall$treatment <- factor(shan_fall$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_fall$system <- factor(shan_fall$system, levels = c("annual", "perennial", "forest"))

shan_fall_p <- shan_fall %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  geom_text(data = letters.df, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Potential (\u00b5g/g soil)") +
  xlab("") +
  ggtitle("Fall Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Arrange figures
Fig_S2 <- ggarrange(
  shan_summer_p, 
  shan_fall_p,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)

ggsave("Fig_S2.png", plot = Fig_S2, dpi = 400, w = 10, h = 7)

## Supplementary Figure 3: Differential Abundance ------------------------
## DIFFERENTIAL ABUNDANCE: Summer vs. Fall
bact_mcse <- subset_samples(bact.p_filtered, incubation == "in-situ")
diagdds <- phyloseq_to_deseq2(bact_mcse, ~ sample_date) 
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")


results <- results(diagdds, alpha = 0.01, contrast = c("sample_date", "summer","fall"))

sigtab_results <- results[which(results$padj < 0.01), ]

sigtab_results_TAX <- cbind(as(sigtab_results, "data.frame"), as(tax_table(bact_mcse)[row.names(sigtab_results), ], "matrix"))

sigtab_results_TAX_genra <- subset(sigtab_results_TAX, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Class <- factor(as.character(sigtab_results_TAX_genra$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Genus <- factor(as.character(sigtab_results_TAX_genra$Genus), levels=names(x))


Fig_S3 <- sigtab_results_TAX_genra %>%
          ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
          theme_classic() +
          geom_point(size = 5.75, pch = 21) +
          geom_hline(yintercept=0, linetype="dashed", color = "red") +
          geom_text(x=2.3, y=1.5, label="Summer-Enriched", size = 5) +
          scale_fill_manual(values=met.brewer("Egypt", 12)) +
          theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =20), 
                axis.text.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size =15),
                legend.text = element_text(size = 15),
                legend.title = element_text(size =15))
summary(sigtab_results)

##GGsave
ggsave("Fig_S3.png", plot = Fig_S3, dpi = 400, w = 10, h = 7)


## Supplementary Figure 4: DON/DOC --------------------------------
bact_fum <- subset_samples(bact.rfy, incubation != "in-situ")
don_data <- bact_fum@sam_data
don_data <- as.matrix(don_data)
don_data <- as.data.frame(don_data)

# DON
don_data <- don_data %>%
  filter(don != 'NA')
don_data$incubation <- as.numeric(don_data$incubation)
don_data$don <- as.numeric(don_data$don)
don_data$treatment <- factor(don_data$treatment, levels = c("T1", "T3", "T4", "T5",
                             "T7", "SF", "DF", "CF"))
don_data$system <- factor(don_data$system, levels = c("annual", "perennial", "forest"))

don_fig <- don_data %>%
          ggplot(aes(x=incubation, y = don, color = treatment)) +
          geom_point(size = 5.5) +
          facet_grid(~system) +
          geom_smooth(method = 'lm') +
          ggtitle("Dissolvable Organic Nitrogen") +
          scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                                      '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                             labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                        "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                        "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
          ylab("DON (\u00b5g/g soil)") +
          xlab("Fumigation Time (Hr)") +
          theme_bw() +
          theme(axis.title.x = element_text(size = 14), 
                axis.title.y = element_text(size = 14), 
                strip.text.x = element_text(size = 14),
                axis.text.x = element_text(size = 12, color = "black"), 
                axis.text.y = element_text(size = 12, color = "black"), 
                panel.border = element_rect(colour = "black", fill=NA, size=2),
                legend.text = element_text(size = 12), 
                legend.title = element_text(size = 14), 
                plot.margin=unit(c(.5,1,.5,.5),"cm"))

# DOC
doc_data <- bact_fum@sam_data
doc_data <- as.matrix(doc_data)
doc_data <- as.data.frame(doc_data)
doc_data$incubation <- as.numeric(doc_data$incubation)
doc_data$doc <- as.numeric(doc_data$doc)

doc_data <- doc_data %>%
  filter(doc != 'NA')

doc_data$treatment <- factor(doc_data$treatment, levels = c("T1", "T3", "T4", "T5",
                                                            "T7", "SF", "DF", "CF"))
doc_data$system <- factor(doc_data$system, levels = c("annual", "perennial", "forest"))

doc_fig <- doc_data %>%
  ggplot(aes(x=incubation, y = doc, color = treatment)) +
  geom_point(size = 5.5) +
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  ggtitle("Dissolvable Organic Carbon") +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("DOC (\u00b5g/g soil)") +
  xlab("Fumigation Time (Hr)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Arrange Figures
Fig_S4 <- ggarrange(
  don_fig, 
  doc_fig,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)


ggsave("Fig_S4.png", plot = Fig_S4, dpi = 400, w = 16, h = 14)



mod <- lm(don ~ incubation, data = don_data)
mod <- lm(doc ~ incubation, data = doc_data)
# Other soil biogeochemical tests
mbc_data <- bact_fum@sam_data
mbc_data <- as.matrix(mbc_data)
mbc_data <- as.data.frame(mbc_data)
mbc_data$incubation <- as.numeric(mbc_data$incubation)
mbc_data$mbc <- as.numeric(mbc_data$mbc)

mbc_data <- mbc_data %>%
  filter(mbc != 'NA')
mod <- lm(mbc ~ incubation, data = mbc_data)

mbn_data <- bact_fum@sam_data
mbn_data <- as.matrix(mbn_data)
mbn_data <- as.data.frame(mbn_data)
mbn_data$incubation <- as.numeric(mbn_data$incubation)
mbn_data$mbn <- as.numeric(mbn_data$mbn)

mbn_data <- mbn_data %>%
  filter(mbn != 'NA')
mod <- lm(mbn ~ incubation, data = mbn_data)


## Supplementary Figure 5: Differential Abundance Fumigation---------------------------------
bact_fum_ext <- subset_samples(bact.p_filtered, incubation == "0" | incubation == "24")

## DIFFERENTIAL ABUNDANCE: 0 Hr. vs 24Hr 
diagdds <- phyloseq_to_deseq2(bact_fum_ext, ~ incubation) 
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")


results <- results(diagdds, alpha = 0.01, contrast = c("incubation", "0","24"))


sigtab_results <- results[which(results$padj < 0.01), ]
sigtab_results_0_vs_24 <- results[which(results$padj < 0.01), ]


sigtab_results_TAX <- cbind(as(sigtab_results, "data.frame"), as(tax_table(bact_fum_ext)[row.names(sigtab_results), ], "matrix"))

sigtab_results_TAX_genra <- subset(sigtab_results_TAX, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Class <- factor(as.character(sigtab_results_TAX_genra$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Genus <- factor(as.character(sigtab_results_TAX_genra$Genus), levels=names(x))


p1 <- sigtab_results_TAX_genra %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  theme_classic() +
  geom_point(size = 5.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_text(x=24.5, y=1.5, label="0 Hr. Enriched", size = 5) +
  geom_text(x=24.5, y=-1.5, label="24 Hr. Enriched", size = 5) +
  scale_fill_manual(values=met.brewer("Egypt", 13)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, 
                                   vjust=0.5, size =20), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15)) 

## DIFFERENTIAL ABUNDANCE: 0 Hr. vs 72Hr
bact_fum_ext <- subset_samples(bact.p_filtered, incubation == "0" | incubation == "72")
diagdds <- phyloseq_to_deseq2(bact_fum_ext, ~ incubation) 
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans = apply(counts(diagdds), 1, gm_mean)
diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds = DESeq(diagdds, fitType="parametric")


results <- results(diagdds, alpha = 0.01, contrast = c("incubation", "0","72"))

sigtab_results <- results[which(results$padj < 0.01), ]
sigtab_results_0_vs_72 <- results[which(results$padj < 0.01), ]

sigtab_results_TAX <- cbind(as(sigtab_results, "data.frame"), as(tax_table(bact_fum_ext)[row.names(sigtab_results), ], "matrix"))

sigtab_results_TAX_genra <- subset(sigtab_results_TAX, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Class <- factor(as.character(sigtab_results_TAX_genra$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra$log2FoldChange, sigtab_results_TAX_genra$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra$Genus <- factor(as.character(sigtab_results_TAX_genra$Genus), levels=names(x))


p2 <- sigtab_results_TAX_genra %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  theme_classic() +
  geom_point(size = 5.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "red") +
  geom_text(x=16.5, y=1.5, label="0 Hr. Enriched", size = 5) +
  geom_text(x=16.5, y=-1.5, label="72 Hr. Enriched", size = 5) +
  scale_fill_manual(values=met.brewer("Egypt", 12)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, 
                                   vjust=0.5, size =20), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15)
  ) 

Fig_S5 <- ggarrange(
  p1, 
  p2,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)

ggsave("Fig_S5.png", plot = Fig_S5, dpi = 400, w = 10, h = 16)
summary(sigtab_results_0_vs_24)
summary(sigtab_results_0_vs_72)

