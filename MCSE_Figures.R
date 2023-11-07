### Code and data analysis for: [Insert Manuscript Title Here]
### DOI:
### Please contact Brandon Kristy (kristybr@msu.edu) for any questions
### All analyses are organized as presented in the manuscript.

rm(list=ls()) # Remove all objects from Environment
# MAC OS
# Note -- Change this to the GitHub Repository Directory once downloading onto your PC
setwd("~/Dropbox/Projects/LTER_MCSE/data/seq_data/")
# Required Libraries
library(bbmle)
library(caret)
library(cowplot)
library(data.table)
library(DESeq2)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(grid)
library(hillR)
library(knitr)
library(lazyeval)
library(metacoder)
library(multcompView)
library(pander)
library(permute)
library(phyloseq)
library(plyr)
library(ppcor)
library(lattice)
library(magrittr)
library(microbiome)
library(MetBrewer)
library(nlme)
library(parallel)
library(reshape2)
library(tidyverse)
library(vegan)


# Load phyloseq objects - filtered and rarefied data
load("phylo_filtered.RData") # Filtered, unrarefied data -- We will use this object for DESEQ2, which performs its own internal normalization
load("phylo_rarefied.RData") # Filtered, rarefied data -- We will use this object for all other analyses

## Figure 1: Alpha-Diversity ----------------------------

# Subset only the in-situ soil core samples from rarefied data
bact_mcse <- subset_samples(bact.rfy, incubation == "in-situ")

# Obtain OTU table from phyloseq object, transpose, and covert to data frame
otu_table <- t(as(bact_mcse@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

# Calculate hill number estimates for each sample at q=1 -- proxy to shannon diversity
shan.mcse <- hill_taxa(otu_table, q = 1)
# Add metadata from the phyloseq object to the dataframe with our hill number estimates
shan.mcse <- merge(bact_mcse@sam_data, shan.mcse, by = 0)
# Change the name of the Hill number column from 'y' to 'shan'
colnames(shan.mcse)[colnames(shan.mcse) == "y"] ="shan"

# Statistical Analyses -- Effects of factors (Land management, Sample Date) on hill number estimate (q=1)

# Total model: How does sampling date and land management treatment effect alpha diversity?
#ANOVA Model
mod.aov <- aov(shan ~ treatment * sample_date, data = shan.mcse)
summary(mod.aov)

# What is the effect of sample date (Fall vs. Summer) on alpha diversity?
mod <- lm(shan ~ sample_date + 0, data = shan.mcse)
summary(mod)
confint(mod) # obtain confidence intervals for publication

# Subset in-situ samples by sampling date - Summer or Fall 
shan_summer <- shan.mcse %>%
  subset(sample_date == "summer")

shan_fall <- shan.mcse %>%
  subset(sample_date == "fall")

# We will now evaluate the individual effects of land management treatment, or more broadly system, within each sampling date:
# What is the effect of treatment -within- each sampling date on alpha diversity?
# Summer
mod.aov <- aov(shan ~ treatment, data = shan_summer)
summary(mod.aov)

# Obtain individual effects of each treatment
mod <- lm(shan ~ treatment + 0, data = shan_summer)
summary(mod)

# What is the effect by system (Annual, Perennial, Forest) on alpha diversity? 
mod.aov <- aov(shan ~ system, data = shan_summer)
summary(mod.aov)
confint(mod.aov) # obtain confidence intervals


# Obtain mean alpha diversity values for each system (annual, perennial, treatment) in the summer sampling date 
mod <- lm(shan ~ system + 0, data = shan_summer)
summary(mod)

# Add model estimates into dataframe for GGPLOT
est <- c(244.85, 168.82, 87.54)
system <- c("annual", "perennial", "forest")
summer_systems_means <- data.frame(est, system)

# Fall
mod.aov <- aov(shan ~ treatment, data = shan_fall)
summary(mod.aov)

# Obtain individual effects of each treatment
mod <- lm(shan ~ treatment + 0, data = shan_fall)
summary(mod)
confint(mod)

# What is the effect by system (Annual, Perennial, Forest) on alpha diversity? 
mod.aov <- aov(shan ~ system, data = shan_fall)
summary(mod.aov)

# Obtain mean alpha diversity values for each system (annual, perennial, treatment) in the fall sampling date 
mod <- lm(shan ~ system + 0, data = shan_fall)
summary(mod)
# Add model information into dataframe
est <- c(281.24, 212.20, 193.25)
system <- c("annual", "perennial", "forest")
fall_systems_means <- data.frame(est, system)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters.df.summer <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = shan_summer))$treatment[,4])$Letters)
colnames(letters.df.summer)[1] <- "Letter" #Reassign column name
letters.df.summer$treatment <- rownames(letters.df.summer) 
letters.df.summer$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
letters.df.summer$system <- value_map[letters.df.summer$treatment]
letters.df.summer$system <- factor(letters.df.summer$system, levels = c("annual", "perennial", "forest"))

# Add Tukey's Post-Hoc Comparisons fir each MCSSE Treatment: Fall sampless
letters.df.fall <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = shan_fall))$treatment[,4])$Letters)
colnames(letters.df.fall)[1] <- "Letter" #Reassign column name
letters.df.fall$treatment <- rownames(letters.df.fall) 
letters.df.fall$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
letters.df.fall$system <- value_map[letters.df.fall$treatment]
letters.df.fall$system <- factor(letters.df.fall$system, levels = c("annual", "perennial", "forest"))


# Summer Alpha-Diversity Plot - Figure 1
# Reorder Treatments
shan_summer$treatment <- factor(shan_summer$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_summer$system <- factor(shan_summer$system, levels = c("annual", "perennial", "forest"))
summer_systems_means$system <- factor(summer_systems_means$system, levels = c("annual", "perennial", "forest"))

# Plot alpha diversity across all treatments in the Summer -- with Post-hoc comparisons across treatment and system-model mean estimates
shan_summer_p <- shan_summer %>%
                  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
                  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
                  geom_text(data = letters.df.summer, aes(x=treatment, y=placement, label=Letter), 
                            size =5.5, color="black", fontface="bold") +
                  facet_grid(~system, scales = 'free') +
                  geom_hline(aes(yintercept=est), linetype = 'dashed', col = 'black', size = 1.25,
                             data=summer_systems_means) +
                  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
                  ylab("Hill Number (q=1)") +
                  xlab("") +
                  ggtitle("Summer Samples") +
                  theme_bw() +
                  theme(axis.title.x = element_text(size = 14), 
                        title = element_text(size = 14),
                        axis.title.y = element_text(size = 14), 
                        strip.text.x = element_text(size = 14),
                        axis.text.x = element_text(size = 12, color = "black"), 
                        axis.text.y = element_text(size = 12, color = "black"), 
                        legend.text = element_text(size = 12), 
                        legend.title = element_text(size = 14), 
                        plot.margin=unit(c(.5,1,.5,.5),"cm"))



# Fall Alpha-Diversity Plot
# Reorder Treatments
shan_fall$treatment <- factor(shan_fall$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_fall$system <- factor(shan_fall$system, levels = c("annual", "perennial", "forest"))
fall_systems_means$system <- factor(fall_systems_means$system, levels = c("annual", "perennial", "forest"))



# Plot alpha diversity across all treatments in the Fall -- with Post-hoc comparisons across treatment and system-model mean estimates
shan_fall_p <- shan_fall %>%
  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(~system, scales = 'free') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', col = 'black', size = 1.25,
             data=fall_systems_means) +
  geom_text(data = letters.df.fall, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Hill Number (q=1)") +
  xlab("") +
  ggtitle("Fall Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        strip.text.x = element_text(size = 14),
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Combine Summer and Fall plots into Figure 1 (Panel A and B)
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


## Figure 2: Beta-diversity ordinations ----------------------------------
## Total Dataset Ordination -- This is NOT in the publication -- as we separated the analyses into our in-situ soil samples and our biodiversity reduction samples -- our reduction treatment has a strong affect on community composition displayed in the first ordniation

# Create envfit model to visualize effects of MBC/DOC/DON as vectors on the NMDS ordination
set.seed(37920)
# Convert to relative abundance 
bact.rel <- transform_sample_counts(bact.rfy, function(x) x / sum(x) )
# Perform NMDS ordination with Bray curtis distance using the ordinate() function
bray.NMDS.fum.ord <- ordinate(bact.rel, method = "NMDS", distance = "bray") # stress = 0.218

# Exctract sample metadata from the rarefied phyloseq object 
env <- as.data.frame(sample_data(bact.rfy))
set.seed(37920)
# Perform the envfit() function to quantify effects of continuous variables on community composition
ord.fit.sig <- envfit(bray.NMDS.fum.ord ~ env$n.fix + env$mbc + env$mbn + env$doc + env$don , 
                      perm = 9999, na.rm = TRUE)

# Obtain NMDS point coordinates for each sample, and append relevant metadata from the env dataframe 
data.scores <- as.data.frame(scores(bray.NMDS.fum.ord)$sites)
data.scores$treatment <- env$treatment
data.scores$incubation <- env$incubation
data.scores$sample_date <- env$sample_date
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
data.scores$mbn <- env$mbn
data.scores$doc <- env$doc
data.scores$don <- env$don

# Extract vector coordinates of metadata variables 
en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors"))
# Extract p-values of each vector
en_coord_cont$pval <- ord.fit.sig$vectors$pvals
# Subset vectors that are statistically significant (p <0.05)
en_coord_cont_sig <- en_coord_cont %>% filter(pval < 0.05)
# Add names of relevant vectors to the dataframe -- for the ggplot
en_coord_cont_sig$name <- c("N.fix", "MBC", "MBN", "DOC", "DON")

# Reorganize treatment factor 
data.scores$treatment <- factor(data.scores$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))

# Display total ordination of all data -- this is not in publication 
data.scores %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(fill = treatment, shape = incubation), size = 5.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = 1.2, y = 1.3, label = "stress = 0.197") +
  scale_fill_manual(values=c('#EEAC4E', '#D7785D','#9C7C50', '#577753','#1A4C76', 
                             '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",  
                               "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22, 23, 24), name = "Biodiversity Reduction Time") +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont_sig, size =1.25, alpha = 0.7, colour = "black") +
  geom_label(data = en_coord_cont_sig, aes(x = NMDS1, y = NMDS2), colour = "black", 
             fontface = "bold", size = 5.5, label = en_coord_cont_sig$name, fill = "white") + 
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

## In-situ Analysis -- What are the effects of management, seasonality, and soil edaphic variables on diazotroph composition 
set.seed(37920)
# Convert to relative abundance
bact.mcse.rel <- transform_sample_counts(bact_mcse, function(x) x / sum(x) )
# Perform NMDS ordination with Bray curtis distance using the ordinate() function
bray.NMDS.mcse.ord <- ordinate(bact.mcse.rel, method = "NMDS", distance = "bray") # Stress = 0.185

# Calculate Bray-curtis distance using phyloseq::distance to identify significance of our factors: treatment and season/sample date
set.seed(37920)
bray.NMDS.mcse <- phyloseq::distance(bact.mcse.rel,"bray")

# What are the interactive effects of sampling date and treatment on diazotroph Beta-diversity? -- perform permanova
set.seed(37920)
permanova <- adonis(bray.NMDS.mcse~sample_data(bact.mcse.rel)$treatment *
                      sample_data(bact.mcse.rel)$sample_date, permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication

# Create envfit model to visualize effects of MBC/DOC/DON as vectors on the NMDS ordination
# Exctract sample metadata from the rarefied phyloseq object 
env <- as.data.frame(sample_data(bact_mcse))
set.seed(37920)
# Perform the envfit() function to quantify effects of continuous variables on community composition
ord.fit.sig <- envfit(bray.NMDS.mcse.ord ~ env$soil_temp + 
                        env$soil_moisture +
                        env$n.fix + env$mbc + env$mbn + 
                        env$doc + env$don, 
                      perm = 9999, na.rm = TRUE)

# Obtain NMDS point coordinates for each sample, and append relevant metadata from the env dataframe 
data.scores <- as.data.frame(scores(bray.NMDS.mcse.ord)$sites)
data.scores$treatment <- env$treatment
data.scores$sample_date <- env$sample_date
data.scores$soil_temp <- env$soil_temp
data.scores$soil_moisture <- env$soil_moisture
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
data.scores$mbn <- env$mbn
data.scores$doc <- env$doc
data.scores$don <- env$don

# Extract vector coordinates of metadata variables 
en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors"))
# Extract p-values of each vector
en_coord_cont$pval <- ord.fit.sig$vectors$pvals
# Subset vectors that are statistically significant (p <0.05)
en_coord_cont_sig <- en_coord_cont %>% filter(pval < 0.05)
# Add names of relevant vectors to the dataframe -- for the ggplot
en_coord_cont_sig$name <- c("Soil.Moisture", "N.fix", "DOC", "DON")

# Reorganize treatment factor 
data.scores$treatment <- factor(data.scores$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))

# Display NMDS ordination of our in-situ dataset
Fig_2a <- data.scores %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(fill = treatment, shape = sample_date), size = 5.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = 0.2, y = 0.15, label = "stress = 0.185") +
  scale_fill_manual(values=c('#EEAC4E', '#D7785D','#9C7C50', '#577753','#1A4C76', 
                             '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",  
                               "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22), name = "Sample Date", 
                     labels = c("Fall", "Summer")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont_sig, size =1.25, alpha = 0.7, colour = "black") +
  geom_label(data = en_coord_cont_sig, aes(x = NMDS1, y = NMDS2), colour = "black", 
             fontface = "bold", size = 5.5, label = en_coord_cont_sig$name, fill = "white") + 
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  ggtitle("In-situ Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 12, color = "black"))  + 
  theme(axis.text.y = element_text(size = 12, color = "black"))  + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(legend.text = element_text(size = 12)) +
  theme(legend.position="right", legend.title=element_text(size=14))+
  theme(legend.position="right")



## In-vitro Analysis -- What are the effects of management, biodiversity reduction, and soil edaphic variables on diazotroph community composition 
# Subset fumigated samples from rarefied phyloseq object
bact_fum <- subset_samples(bact.rfy, incubation != "in-situ")
# Conver the abundance values to relative abundance 
bact.fum.rel <- transform_sample_counts(bact_fum, function(x) x / sum(x) )
# Perform NMDS ordination with Bray curtis distance using the ordinate() function
set.seed(37920)
bray.NMDS.fum.ord <- ordinate(bact.fum.rel, method = "NMDS", distance = "bray")

# Calculate Bray-curtis distance using phyloseq::distance to identify significance of our factors: treatment and biodiversity reduction
bray.NMDS.fum <- phyloseq::distance(bact.fum.rel,"bray")
set.seed(37920)
# What are the interactive effects of biodiversity reduction and treatment on diazotroph Beta-diversity? -- perform permanova
permanova <- adonis(bray.NMDS.fum~sample_data(bact.fum.rel)$treatment * 
                      sample_data(bact.fum.rel)$incubation, 
                    permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication 


# Create envfit model to visualize effects of MBC/DOC/DON as vectors on the NMDS ordination
# Extract metadata from biodiversity reduction phyloseq object
env <- as.data.frame(sample_data(bact_fum))
set.seed(37920)
# Perform envfit() function to quantify effects of continuous variables on diazotroph community composition
ord.fit.sig <- envfit(bray.NMDS.fum.ord ~ env$n.fix + env$mbc + env$mbn + env$doc + env$don , 
                      perm = 9999, na.rm = TRUE)

# Obtain NMDS coordinate locations for every sample in the dataframe -- and add relevant metadata
data.scores <- as.data.frame(scores(bray.NMDS.fum.ord)$sites)
data.scores$treatment <- env$treatment
data.scores$incubation <- env$incubation
data.scores$sample_date <- env$sample_date
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
data.scores$mbn <- env$mbn
data.scores$doc <- env$doc
data.scores$don <- env$don

# Extract coordinates of our continuous vectors 
en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors"))
# Extract p-values for each vector
en_coord_cont$pval <- ord.fit.sig$vectors$pvals
# Subset and only keep significant variables  (p < 0.05)
en_coord_cont_sig <- en_coord_cont %>% filter(pval < 0.05)
# Add the names of each vector onto the dataframe -- for ggplot
en_coord_cont_sig$name <- c("N.fix", "MBC", "MBN", "DOC", "DON")

# Reorganize treatment factor 
data.scores$treatment <- factor(data.scores$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))

Fig_2b <- data.scores %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(fill = treatment, shape = incubation), size = 5.5) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = 1.2, y = 1.3, label = "stress = 0.197") +
  scale_fill_manual(values=c('#EEAC4E', '#D7785D','#9C7C50', '#577753','#1A4C76', 
                             '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",  
                               "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22, 23), name = "Biodiversity Reduction Time", 
                     labels = c("0 Hours", "24 Hours", "72 hours")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont_sig, size =1.25, alpha = 0.7, colour = "black") +
  geom_label(data = en_coord_cont_sig, aes(x = NMDS1, y = NMDS2), colour = "black", 
             fontface = "bold", size = 5.5, label = en_coord_cont_sig$name, fill = "white") + 
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  ggtitle("In-vitro Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 12, color = "black"))  + 
  theme(axis.text.y = element_text(size = 12, color = "black"))  + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=2))+
  theme(legend.text = element_text(size = 12)) +
  theme(legend.position="right", legend.title=element_text(size=14))+
  theme(legend.position="right")

## Arrange both ordinations together
Fig_2 <- ggarrange(
  Fig_2a, 
  Fig_2b,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)



# GGsave
ggsave("Fig_2.png", plot = Fig_2, dpi = 400, w = 11, h = 14)


## Figure 3: FLNF Rates across in-situ Dataset --------------------------------------------------------------
# Remove NAs
shan.mcse <- shan.mcse %>%
  subset(n.fix != 'NA')

# Statistical Analyses
# Total model: How does sampling date and treatment effect FLNF rates?
#ANOVA Model
mod.aov <- aov(n.fix ~ treatment * sample_date, data = shan.mcse)
summary(mod.aov)

# Total model: How does sampling date and system (grouped MCSE treatments) effect FLNF rates
mod.aov <- aov(n.fix ~ system * sample_date, data = shan.mcse)
summary(mod.aov)

# Obtain effect size of Sampling data /w effects parameterization 
mod <- lm(n.fix ~ sample_date + 0, data = shan.mcse)
summary(mod)
confint(mod)


# Obtain effect size of System /w effects paranmeterization
mod <- lm(n.fix ~ system + 0, data = shan.mcse)
summary(mod)

# Subset in-situ samples by sampling date - Summer or Fall 
shan_summer <- shan.mcse %>%
  subset(sample_date == "summer")

shan_fall <- shan.mcse %>%
  subset(sample_date == "fall")

# What is the effect of treatment -within- each sampling date on mitrogen fixation?
# Summer
mod.aov <- aov(n.fix ~ treatment, data = shan_summer)
summary(mod.aov)

# Obtain individual effects of each treatment
mod <- lm(n.fix ~ treatment + 0, data = shan_summer)
summary(mod)



# What is the effect by system (Annual, Perennial, Forest) on nitrogen fixation? 
mod.aov <- aov(n.fix ~ system, data = shan_summer)
summary(mod.aov)

# Obtain individual effects of each system 
mod <- lm(n.fix ~ system + 0, data = shan_summer)
summary(mod)

est <- c(3.5283, 2.6012, 1.4306)
system <- c("annual", "perennial", "forest")
summer_systems_means <- data.frame(est, system)


# Fall
mod.aov <- aov(n.fix ~ treatment, data = shan_fall)
summary(mod.aov)

# Obtain individual effects of each treatment
mod <- lm(n.fix ~ treatment + 0, data = shan_fall)
summary(mod)

# What is the effect by system (Annual, Perennial, Forest) on alpha diversity? 
mod.aov <- aov(n.fix ~ system, data = shan_fall)
summary(mod.aov)

# Obtain individual effects of each system 
mod <- lm(n.fix ~ system + 0, data = shan_fall)
summary(mod)

est <- c(5.004, 9.550, 2.940)
system <- c("annual", "perennial", "forest")
fall_systems_means <- data.frame(est, system)

# Add Tukey's Post-Hoc Comparisons: Summer
letters.df.summer <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = shan_summer))$treatment[,4])$Letters)
colnames(letters.df.summer)[1] <- "Letter" #Reassign column name
letters.df.summer$treatment <- rownames(letters.df.summer) 
letters.df.summer$placement <- c(11, 11, 11, 11, 11, 11, 11, 11)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
letters.df.summer$system <- value_map[letters.df.summer$treatment]
letters.df.summer$system <- factor(letters.df.summer$system, levels = c("annual", "perennial", "forest"))

# Add Tukey's Post-Hoc Comparisons: Fall
letters.df.fall <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = shan_fall))$treatment[,4])$Letters)
colnames(letters.df.fall)[1] <- "Letter" #Reassign column name
letters.df.fall$treatment <- rownames(letters.df.fall) 
letters.df.fall$placement <- c(15, 15, 15, 15, 15, 15, 15, 15)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
letters.df.fall$system <- value_map[letters.df.fall$treatment]
letters.df.fall$system <- factor(letters.df.fall$system, levels = c("annual", "perennial", "forest"))

# Summer Plot
# Reorder Treatments
shan_summer$treatment <- factor(shan_summer$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_summer$system <- factor(shan_summer$system, levels = c("annual", "perennial", "forest"))
summer_systems_means$system <- factor(summer_systems_means$system, levels = c("annual", "perennial", "forest"))

Fig_3A <- shan_summer %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(~system, scales = 'free') +
  geom_text(data = letters.df.summer, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', col = 'black', size = 1.25,
             data=summer_systems_means) +
  scale_y_continuous(limits = c(0, 15)) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  xlab("") +
  ggtitle("Summer Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))


# Fall Plot
# Reorder Treatments
shan_fall$treatment <- factor(shan_fall$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
shan_fall$system <- factor(shan_fall$system, levels = c("annual", "perennial", "forest"))
fall_systems_means$system <- factor(fall_systems_means$system, levels = c("annual", "perennial", "forest"))

Fig_3B <- shan_fall %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(~system, scales = 'free') +
  geom_text(data = letters.df.fall, aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', col = 'black', size = 1.25,
             data=fall_systems_means) +
  scale_y_continuous(limits = c(0, 15)) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  xlab("") +
  ggtitle("Fall Samples") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.position="none",
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

# Interactive plot 
shan.mcse$sample_date <- factor(shan.mcse$sample_date, levels = c("summer", "fall"))
shan.mcse$system <- factor(shan.mcse$system, levels = c("annual", "perennial", "forest"))

library(interactions)
mod <- lm(n.fix ~ system*sample_date, data = shan.mcse)
Fig_3C <- cat_plot(mod, pred = sample_date, 
                  modx = system, 
                  geom = "line", 
                  point.shape = TRUE) +
                  scale_y_continuous(limits = c(0, 15)) +
                  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
                  xlab("Sample Date") +
                  theme_bw() +
                  theme(axis.title.x = element_text(size = 14),
                  title = element_text(size = 14),
                  axis.title.y = element_text(size = 14), 
                  axis.text.x = element_text(size = 12, color = "black"), 
                  axis.text.y = element_text(size = 12, color = "black"), 
                  legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 14), 
                  plot.margin=unit(c(.5,1,.5,.5),"cm"))


Fig_3AB <- ggarrange(
  Fig_3A + rremove("ylab") + rremove("legend"),
  Fig_3B + rremove("ylab"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T, 
  legend = "right"
)

Fig_3 <- ggarrange(
  Fig_3AB,
  Fig_3C + rremove("ylab"), 
  labels = c("", "C"),
  ncol = 1, 
  nrow = 2,
  common.legend = F,
  legend = "right"
)

Fig_3 <- annotate_figure(Fig_3, left = textGrob("Nitrogen Fixation Rate (\u00b5g/g soil)", 
                                                  rot = 90, gp = gpar(cex = 1.5), hjust = 0.25))

Fig_3

ggsave("Fig_3.png", plot = Fig_3, dpi = 600, w = 10, h = 14)

## Figure 4: Fumigation on Diazotroph Alpha-diversity, N-fixation, and MBC (community size) ----------------------
# Subset only fumigated samples
bact_fum <- subset_samples(bact.rfy, incubation != "in-situ")

# Calculate alpha-diveristy (Hill number q=1 -- 'Shannon Diversity')
otu_table <- t(as(bact_fum@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

shan.fum <- hill_taxa(otu_table, q = 1)
shan.fum <- merge(bact_fum@sam_data, shan.fum, by = 0)

colnames(shan.fum)[colnames(shan.fum) == "y"] ="shan"

# Set incubation as a numeric variable for statistical tests
shan.fum$incubation <- as.numeric(shan.fum$incubation)

# Statistical Analyses

# What is the effect of biodiversity reduction on Diazotroph-Alpha Diversity?
# Perform partial correlation between (div~fum), given MBC/DOC/DON - controlling for community size-effects AND changes in soil biogeochemistry
# Remove 'NA' values for doc/don/mbc 
shan.fum_pcor <- shan.fum %>%
  filter(doc != 'NA')

shan.fum_pcor <- shan.fum_pcor %>%
  filter(don != 'NA')

shan.fum_pcor <- shan.fum_pcor %>%
  filter(mbc != 'NA')

shan.fum_pcor <- shan.fum_pcor %>%
  filter(n.fix != 'NA')

# Perform partial correlation analysis - identify relationship when controlling for covariates 
pcor.div.fum <- pcor.test(shan.fum_pcor$shan, shan.fum_pcor$incubation, shan.fum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.div.fum # No significant correlation when adjusting for soil biogeochem + community size

# How does fumigation time and treatment interact to influence diazotroph alpha diversity 
mod.aov <- aov(shan ~ incubation*treatment, data = shan.fum)
summary(mod.aov)

# How does fumigation time and system interact to influence diazotroph alpha diversity
mod.aov <- aov(shan ~ incubation*system, data = shan.fum)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(shan ~ incubation, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(shan ~ incubation, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(shan ~ incubation, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -1.6", "Est = -0.3", "Est = -0.8" )
R <- c("Adj-R.sq = 0.35", "Adj-R.sq = -0.03", "Adj-R.sq = 0.1")
p <- c("p = <0.001", "p = 0.6", "p=0.03")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# shannon diversity plot - with system model information
shan.fum$treatment <- factor(shan.fum$treatment, levels = c("T1", "T3", "T4", "T5",
                                                            "T7", "SF", "DF", "CF"))
shan.fum$system <- factor(shan.fum$system, levels = c("annual", "perennial", "forest"))

shan_plot <- shan.fum %>%
            ggplot(aes(x=incubation, y = shan)) +
            geom_point(aes(color=treatment), size = 5.5) +
            geom_text(x=50, y = 480, aes(label = est), size = 5, data = model.df) +
            geom_text(x=50, y = 460, aes(label = R), size = 5, data = model.df) +
            geom_text(x=50, y = 440, aes(label = p), size = 5, data = model.df) +
            facet_grid(~system) +
            geom_smooth(method = 'lm') +
            scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                                        '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                               labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                          "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                          "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
            ylab("Hill Number (q=1)") +
            xlab("Fumigation Time (Hr)") +
            theme_bw() +
            theme(axis.title.x = element_text(size = 14), 
                  axis.title.y = element_text(size = 14), 
                  strip.text.x = element_text(size = 14),
                  axis.text.x = element_text(size = 12, color = "black"), 
                  axis.text.y = element_text(size = 12, color = "black"), 
                  legend.text = element_text(size = 12), 
                  legend.title = element_text(size = 14), 
                  plot.margin=unit(c(.5,1,.5,.5),"cm"))



## Nitrogen Fixation Rate analyses/plot
shan.fum <- shan.fum %>%
  filter(n.fix != 'NA')

# Statistical Analyses

# What is the effect of fumigation time on Nitrogen Fixation Rate
# Perform partial correlation between (div~fum), given MBC/DOC/DON - controlling for community size-effects AND changes in soil biogeochemistry
# Remove 'NA' values for doc/don/mbc 
shan.fum_pcor <- shan.fum_pcor %>%
  filter(n.fix != 'NA')

shan.fum_pcor <- shan.fum_pcor %>%
  filter(n.fix.mbc != 'NA')

# Perform partial correlation analysis
pcor.nfix.fum <- pcor.test(shan.fum_pcor$n.fix, shan.fum_pcor$incubation, shan.fum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.fum # significant correlation when adjusting for soil biogeochem + community size (Rho = -0.34; p = 0.002)

# How does fumigation time and treatment interact to influence Nitrogen Fixation Rate
mod.aov <- aov(n.fix ~ incubation*treatment, data = shan.fum)
summary(mod.aov)

# How does biodiversity reduction  and system interact to influence diazotroph alpha diversity
mod.aov <- aov(n.fix ~ incubation*system, data = shan.fum)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(n.fix ~ incubation, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix ~ incubation, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(n.fix ~ incubation, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -0.09", "Est = -0.05", "Est = -0.06" )
R <- c("Adj-R.sq = 0.42", "Adj-R.sq = 0.38", "Adj-R.sq = 0.27")
p <- c("p = <0.001", "p = <0.001", "p=0.003")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)

model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# Plot Nitrogen Fixation Rate - with system model information
n.fix_plot <- shan.fum %>%
  ggplot(aes(x=incubation, y = n.fix)) +
  geom_point(aes(color=treatment), size = 5.5) +
  geom_text(x=50, y = 14, aes(label = est), size =5, data = model.df) +
  geom_text(x=50, y = 13, aes(label = R), size =5, data = model.df) +
  geom_text(x=50, y = 12, aes(label = p), size = 5, data = model.df) +
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  xlab("Fumigation Time (Hr)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

## Microbial biomass C analyses/plot
shan.fum <- shan.fum %>%
  filter(mbc != 'NA')

# Statistical Analyses

# What is the effect of fumigation time on microbial biomass C
# Remove 'NA' values for doc/don/mbc 


# Perform partial correlation analysis
pcor.mbc.fum <- pcor.test(shan.fum_pcor$mbc, shan.fum_pcor$incubation, shan.fum_pcor[,c("doc", "don")], method = "spearman")
pcor.mbc.fum # significant correlation when adjusting for soil biogeochem + community size (Rho = -0.14. p = 0.2)

# How does fumigation time and treatment interact to influence MBC
mod.aov <- aov(mbc ~ incubation*treatment, data = shan.fum)
summary(mod.aov)

# How does biodiversity reduction  and system interact to influence MBC
mod.aov <- aov(mbc ~ incubation*system, data = shan.fum)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(mbc ~ incubation, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(mbc ~ incubation, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(mbc ~ incubation, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.22", "Est = -3.54", "Est = -0.85" )
R <- c("Adj-R.sq = -0.02", "Adj-R.sq = 0.34", "Adj-R.sq = -0.02")
p <- c("p = 0.85", "p = 0.002", "p=0.48")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)

model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# Plot Nitrogen Fixation Rate - with system model information
mbc_plot <- shan.fum %>%
  ggplot(aes(x=incubation, y = mbc)) +
  geom_point(aes(color=treatment), size = 5.5) +
  geom_text(x=50, y = 770, aes(label = est), size = 5, data = model.df) +
  geom_text(x=50, y = 720, aes(label = R), size = 5, data = model.df) +
  geom_text(x=50, y = 670, aes(label = p), size = 5, data = model.df) +
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("Microbial Biomass Carbon (\u00b5g/g soil)") +
  xlab("Fumigation Time (Hr)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))




# Arrange both Shannon Diversity, N-fixation, and MBC plots together
Fig_4 <- ggarrange(
  shan_plot + rremove("xlab"), 
  mbc_plot + rremove("xlab"),
  n.fix_plot + rremove ("xlab"),
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3,
  common.legend = T,
  legend = "right"
)
require(grid)   # for the textGrob() function
Fig_4<- annotate_figure(Fig_4, bottom = textGrob("Biodiversity Reduction Time (Hr)", 
                                         vjust = -0.75, hjust = 0.85, gp = gpar(cex = 1.5)))

ggsave("Fig_4.png", plot = Fig_4, dpi = 600, w = 16, h = 21)

## Figure 5: Biodiv vs. Ecosystem Function ----------------------------------

# Remove na nfixation
shan.fum <- shan.fum %>%
  filter(n.fix != 'NA')

# Statistical Analyses

# Perform partial correlation between diazotroph alpha-diversity~n.fixation when controlling for (doc.don.mbc)
pcor.nfix.div.fum <- pcor.test(shan.fum_pcor$shan, shan.fum_pcor$n.fix, shan.fum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.fum # Significant positive association when controlling for covariance in soil biogechemistry and community size

# How does this correlation compare to only unfumigated samples (0 hour control)
shan.unfum <- shan.fum %>%
  subset(incubation == "0")

shan.unfum_pcor <- shan.unfum %>%
  filter(doc != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(don != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(mbc != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(n.fix != 'NA')


pcor.nfix.div.unfum <- pcor.test(shan.unfum_pcor$shan, shan.unfum_pcor$n.fix, shan.unfum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.unfum # No significant association
# How does the biodiversity-ecosystem function relationship vary by treatment and system - generate models


# How does fumigation time and treatment interact to influence diazotroph alpha diversity 
mod.aov <- aov(n.fix ~ shan*treatment, data = shan.fum)
summary(mod.aov)

# How does fumigation time and system interact to influence diazotroph alpha diversity
mod.aov <- aov(n.fix ~ shan*system, data = shan.fum)
summary(mod.aov)

# What is the range in diazotroph alpha diversity estimates between in-situ and in-vitro datasets?
min(shan.mcse$shan) #33.59
max(shan.mcse$shan) #442.29

min(shan.fum$shan) #37.89
max(shan.fum$shan) # 501.27
mean(shan.fum$shan) # 209.2404


min(shan.unfum$shan) # 37.89
max(shan.unfum$shan) # 501.27
mean(shan.unfum$shan) # 263.2013

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(n.fix ~ shan, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix ~ shan, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(n.fix ~ shan, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.01", "Est = 0.007", "Est = 0.04" )
R <- c("Adj-R.sq = 0.087", "Adj-R.sq = 0.027", "Adj-R.sq = 0.49")
p <- c("p = 0.045", "p = 0.21", "p<0.001")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))
# Create Plot with system-model estimates
shan.fum$treatment <- factor(shan.fum$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.fum$system <- factor(shan.fum$system, levels = c("annual", "perennial", "forest"))
shan.fum$incubation <- as.factor(shan.fum$incubation)


Fig_5a <- shan.fum %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(shape = incubation, fill = treatment), size = 4.5) +  
  facet_grid(~system) +
  geom_text(x=420, y = 12, aes(label = est), size = 5, data = model.df) +
  geom_text(x=420, y = 11, aes(label = R), size = 5, data = model.df) +
  geom_text(x=420, y = 10, aes(label = p), size = 5, data = model.df) +
  geom_smooth(method = 'lm') +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                     labels = c(
                       "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22, 23), name = "Biodiversity Reduction Time", 
                     labels = c("0 Hours", "24 Hours", "72 hours")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Hill Number (q=1)") + 
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  ggtitle("In-vitro Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))


## In-situ Div-Fun analysis
shan.mcse <- shan.mcse %>%
  filter(n.fix != 'NA')

shan.mcse <- shan.mcse %>%
  filter(doc != 'NA')

# Perform partial correlation between diazotroph alpha-diversity~n.fixation when controlling for (doc.don.mbc)
pcor.nfix.div.mcse <- pcor.test(shan.mcse$shan, shan.mcse$n.fix, shan.mcse[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.mcse # Significant positive association when controlling for covariance in soil biogechemistry and community size


# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.mcse %>%
  subset(system =="annual")
lm.annual <- lm(n.fix ~ shan, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.mcse %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix ~ shan, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.mcse %>%
  subset(system =="forest")
lm.forest <- lm(n.fix ~ shan, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -0.003", "Est = 0.02", "Est = 0.001" )
R <- c("Adj-R.sq = -0.046", "Adj-R.sq = 0.11", "Adj-R.sq = -0.066")
p <- c("p = 0.78", "p = 0.12", "p=0.91")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))



shan.mcse$treatment <- factor(shan.mcse$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.mcse$system <- factor(shan.mcse$system, levels = c("annual", "perennial", "forest"))

Fig_5b <- shan.mcse %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(fill=treatment), shape =21, size = 4.5) +  
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  ggtitle("In-situ Experiment") +
  geom_text(x=340, y = 18, aes(label = est), size = 5, data = model.df) +
  geom_text(x=340, y = 17, aes(label = R), size = 5, data = model.df) +
  geom_text(x=340, y = 16, aes(label = p), size = 5, data = model.df) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                     labels = c(
                       "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  xlab("Hill Number (q=1)") + 
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))



Fig_5 <- ggarrange(
  Fig_5a + rremove("xlab") + rremove("ylab"), 
  Fig_5b+ rremove("xlab") + rremove("ylab"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)
require(grid)   # for the textGrob() function
Fig_5 <- annotate_figure(Fig_5, bottom = textGrob("Hill Number (q=1)", 
                                                 vjust = -0.75, hjust = 1.25, gp = gpar(cex = 1.5)))
Fig_5 <- annotate_figure(Fig_5, left = textGrob("Nitrogen Fixation Rate (\u00b5g/g soil)", 
                                                   rot = 90, gp = gpar(cex = 1.5)))



ggsave("Fig_5.png", plot = Fig_5, dpi = 400, w = 16, h = 14)

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


## Supplementary Figure 2: Differential Abundance -------------------------------------------------------

# Because alternative normalization is used in DESEQ2 -- we use unrarefied phyloseq objects for these analyses 
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
Fig_S2 <- ggarrange(
  summer_p, 
  fall_p,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)

ggsave("Fig_S2.png", plot = Fig_S2, dpi = 400, w = 10, h = 12)


#Differential Abundance ------------------------ Summer vs fall
## DIFFERENTIAL ABUNDANCE: Summer vs. Fall
bact_mcse <- subset_samples(bact.p_filtered, incubation == "in-situ")
diagdds <- phyloseq_to_deseq2(bact_mcse, ~ sample_date) 
# calculate geometric means prior to estimate size factors
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

## Supplementary Figure 3: Stepwise Regression Analysis --------------------------------
## Identify alternative factors that influence in-situ nitrogen fixation with stepwise regression modeling 
## Evaluating predictors of nitrogen fixation across the in-situ land management gradient
# Subset data for in-situ samples
situ <- subset_samples(bact.rfy, incubation == "in-situ")

# Obtain alpha-diversity measurement via Hill Numbers
# Convert otu table to data frame
otu_table <- t(as(situ@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

# Obtain Hill Number estimate
shan_mcse <- hill_taxa(otu_table, q = 1)

#Merge with metadata
n_fix_data <- merge(situ@sam_data, shan_mcse, by = 0)
colnames(n_fix_data)[colnames(n_fix_data) == "y"] ="shan"
n_fix_data$s.soil_moisture <- as.numeric(scale(n_fix_data$soil_moisture))
n_fix_data$s.soil_temp <- as.numeric(scale(n_fix_data$soil_temp))

n_fix_data_fig <- n_fix_data %>%
  select(n.fix, sample_date, system, s.soil_moisture, soil_moisture, shan, s.soil_temp, soil_temp, mbc, mbn, doc, don ) %>%
  na.omit()
# Perform a Step-wise regression analysis to identify which edaphic/biological factors contribute to nitrogen fixation 
# Subset data 
n_fix_data <- n_fix_data %>%
  select(n.fix, shan, s.soil_temp, s.soil_moisture, mbc, mbn, doc, don ) %>%
  na.omit()

n_fix_data$n.fix <- as.numeric(n_fix_data$n.fix)
n_fix_data$shan <- as.numeric(n_fix_data$shan)
n_fix_data$s.soil_temp <- as.numeric(n_fix_data$s.soil_temp)
n_fix_data$s.soil_moisture <- as.numeric(n_fix_data$s.soil_moisture)
n_fix_data$mbc <- as.numeric(n_fix_data$mbc)
n_fix_data$mbn <- as.numeric(n_fix_data$mbn)
n_fix_data$doc <- as.numeric(n_fix_data$doc)
n_fix_data$don <- as.numeric(n_fix_data$don)

# Fit the full model 
full.model <- lm(n.fix ~., data = n_fix_data)
# Set seed for reproducibility
set.seed(37920)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(n.fix ~., data = n_fix_data,
                    method = "leapBackward", 
                    tuneGrid = data.frame(nvmax = 1:7),
                    trControl = train.control
)
step.model$results

step.model$bestTune
summary(step.model$finalModel)
coef(step.model$finalModel, 2) # Soil moisture and soil temperature predict nitrogen fixation rates


mod <- lm(n.fix ~ s.soil_temp + s.soil_moisture, data = n_fix_data)
summary(mod)
confint(mod)

# Add model information into dataframe
est <- c("Est = -1.571")
R <- c("Adj-R.sq = 0.2487")
p <- c("p = 0.00361")
model.df <- data.frame(est, R, p)

# look at soil temperature co-variance with sample date
mod <- aov(s.soil_temp ~ sample_date, data = n_fix_data_fig)
summary(mod)

# Visualize significant model -- Soil Temperature
Fig_S3A<- n_fix_data_fig %>%
  ggplot(aes(x=soil_temp, y = n.fix)) +
  geom_point(aes(fill=sample_date), shape = 21, size = 5) +
  geom_text(x=27, y = 18, aes(label = est), size = 5, data = model.df) +
  geom_text(x=27, y = 17, aes(label = R), size = 5, data = model.df) +
  geom_text(x=27, y = 16, aes(label = p), size = 5, data = model.df) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 20)) +
  xlab("Soil Temperature (C)") +
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  geom_smooth(method='lm') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))


# Look at soil moisture co-variance with system
mod <- aov(s.soil_moisture ~ system, data = n_fix_data_fig)
summary(mod)

# Perform effects parameterization 
mod <- lm(soil_moisture ~ system + 0, data = n_fix_data_fig)
summary(mod)
confint(mod)









# Add model parameters to soil moisture figure 
est <- c("Est = 1.287")
R <- c("Adj-R.sq = 0.2487")
p <- c("p = 0.01731")
model.df <- data.frame(est, R, p)
# Visualize significant model -- Soil Moisture 
Fig_S3B<- n_fix_data_fig %>%
  ggplot(aes(x=soil_moisture, y = n.fix)) +
  geom_point(aes(fill=system), shape = 21, size = 5) +
  geom_text(x=0.165, y = 18, aes(label = est), size = 5, data = model.df) +
  geom_text(x=0.165, y = 17, aes(label = R), size = 5, data = model.df) +
  geom_text(x=0.165, y = 16, aes(label = p), size = 5, data = model.df) +
  scale_y_continuous(limits = c(0, 20)) +
  xlab("Soil Moisture") +
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  geom_smooth(method='lm') +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))





# Arrange figures

Fig_S3 <- ggarrange(
  Fig_S3A + rremove("ylab"),
  Fig_S3B + rremove("ylab"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)
Fig_S3 <- annotate_figure(Fig_S3, left = textGrob("Nitrogen Fixation Rate (\u00b5g/g soil)", 
                                                rot = 90, gp = gpar(cex = 1.5)))

ggsave("Fig_S3.png", plot = Fig_S3, dpi = 600, w = 10, h = 12)

## Supplementary Figure 4: DON/DOC/MBC/MBN --------------------------------
bact_fum <- subset_samples(bact.rfy, incubation != "in-situ")
don_data <- bact_fum@sam_data
don_data <- as.matrix(don_data)
don_data <- as.data.frame(don_data)

# DON
don_data <- don_data %>%
  filter(don != 'NA')
don_data$incubation <- as.numeric(don_data$incubation)
don_data$don <- as.numeric(don_data$don)

# Statistical Analyses
# How does fumigation time and treatment interact to influence dissolvable organic nitrogen? 
mod.aov <- aov(don ~ incubation*treatment, data = don_data)
summary(mod.aov)

# How does fumigation time and system interact to influence dissolvable organic nitrogen? 
mod.aov <- aov(don ~ incubation*system, data = don_data)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
don_annual <- don_data %>%
  subset(system =="annual")
lm.annual <- lm(don ~ incubation, data = don_annual)
summary(lm.annual)

# Perennial
don_perennial <- don_data %>%
  subset(system =="perennial")
lm.perennial <- lm(don ~ incubation, data = don_perennial)
summary(lm.perennial)

# Forest
don_forest <- don_data %>%
  subset(system =="forest")
lm.forest <- lm(don ~ incubation, data = don_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.055", "Est = 1.02", "Est = 0.26" )
R <- c("Adj-R.sq = 0.27", "Adj-R.sq = 0.27", "Adj-R.sq = 0.021")
p <- c("p = <0.001", "p = 0.005", "p=0.23")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)

model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))


#DON Plot
don_data$treatment <- factor(don_data$treatment, levels = c("T1", "T3", "T4", "T5",
                             "T7", "SF", "DF", "CF"))
don_data$system <- factor(don_data$system, levels = c("annual", "perennial", "forest"))


don_fig <- don_data %>%
          ggplot(aes(x=incubation, y = don)) +
          geom_point(aes(color=treatment),size = 5.5) +
          geom_text(x=50, y = 165, aes(label = est), size = 5, data = model.df) +
          geom_text(x=50, y = 155, aes(label = R), size = 5, data = model.df) +
          geom_text(x=50, y = 145, aes(label = p), size = 5, data = model.df) +
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

# Statistical Analyses
# How does fumigation time and treatment interact to influence dissolvable organic nitrogen? 
mod.aov <- aov(don ~ incubation*treatment, data = don_data)
summary(mod.aov)

# How does fumigation time and system interact to influence dissolvable organic nitrogen? 
mod.aov <- aov(don ~ incubation*system, data = don_data)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
doc_annual <- doc_data %>%
  subset(system =="annual")
lm.annual <- lm(doc ~ incubation, data = doc_annual)
summary(lm.annual)

# Perennial
doc_perennial <- doc_data %>%
  subset(system =="perennial")
lm.perennial <- lm(doc ~ incubation, data = doc_perennial)
summary(lm.perennial)

# Forest
doc_forest <- doc_data %>%
  subset(system =="forest")
lm.forest <- lm(doc ~ incubation, data = doc_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.73", "Est = 1.91", "Est = 1.33" )
R <- c("Adj-R.sq = 0.21", "Adj-R.sq = 0.21", "Adj-R.sq = 0.055")
p <- c("p = 0.004", "p = 0.01", "p=0.13")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)

model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# DOC plot
doc_data$treatment <- factor(doc_data$treatment, levels = c("T1", "T3", "T4", "T5",
                                                            "T7", "SF", "DF", "CF"))
doc_data$system <- factor(doc_data$system, levels = c("annual", "perennial", "forest"))

doc_fig <- doc_data %>%
  ggplot(aes(x=incubation, y = doc)) +
  geom_point(aes(color=treatment), size = 5.5) +
  geom_text(x=50, y = 770, aes(label = est), size = 5, data = model.df) +
  geom_text(x=50, y = 720, aes(label = R), size = 5, data = model.df) +
  geom_text(x=50, y = 670, aes(label = p), size = 5, data = model.df) +
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
## MBN
mbn_data <- bact_fum@sam_data
mbn_data <- as.matrix(mbn_data)
mbn_data <- as.data.frame(mbn_data)
mbn_data$incubation <- as.numeric(mbn_data$incubation)
mbn_data$mbn <- as.numeric(mbn_data$mbn)

mbn_data <- mbn_data %>%
  filter(mbn != 'NA')

# Statistical Analyses
# How does fumigation time and treatment interact to influence microbial biomass N?
mod.aov <- aov(mbn ~ incubation*treatment, data = mbn_data)
summary(mod.aov)

# How does fumigation time and system interact to influence microbial biomass N?
mod.aov <- aov(mbn ~ incubation*system, data = mbn_data)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
mbn_annual <- mbn_data %>%
  subset(system =="annual")
lm.annual <- lm(mbn ~ incubation, data = mbn_annual)
summary(lm.annual)

# Perennial
mbn_perennial <- mbn_data %>%
  subset(system =="perennial")
lm.perennial <- lm(mbn ~ incubation, data = mbn_perennial)
summary(lm.perennial)

# Forest
mbn_forest <- mbn_data %>%
  subset(system =="forest")
lm.forest <- lm(mbn ~ incubation, data = mbn_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.25", "Est = -0.36", "Est = -0.003" )
R <- c("Adj-R.sq = 0.10", "Adj-R.sq = 0.11", "Adj-R.sq = -0.04")
p <- c("p = 0.03", "p = 0.07", "p=0.99")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)

model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# MBN Plot
mbn_data$treatment <- factor(mbn_data$treatment, levels = c("T1", "T3", "T4", "T5",
                                                            "T7", "SF", "DF", "CF"))
mbn_data$system <- factor(mbn_data$system, levels = c("annual", "perennial", "forest"))

mbn_fig <- mbn_data %>%
  ggplot(aes(x=incubation, y = mbn)) +
  geom_point(aes(color = treatment), size = 5.5) +
  geom_text(x=50, y = 100, aes(label = est), size = 5, data = model.df) +
  geom_text(x=50, y = 90, aes(label = R), size = 5, data = model.df) +
  geom_text(x=50, y = 80, aes(label = p), size = 5, data = model.df) +
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  ggtitle("Microbial Biomass Nitrogen") +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76','#86A592','#577B96','#B0381D'), name = "MCSE Treatment",
                     labels = c("T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                                "T4: Bio-Based Corn/Soy/Wheat", "T5: Poplar", "T7: Early Successional",
                                "SF: Successional Forest", "DF: Deciduous Forest", "CF: Coniferous Forest")) +
  ylab("MBN (\u00b5g/g soil)") +
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


# Arrange all figures
Fig_S4 <- ggarrange(
  don_fig + rremove("xlab"),
  doc_fig + rremove("xlab"),
  mbn_fig + rremove("xlab"),
  labels = c("A", "B", "C"),
  ncol = 1,
  nrow = 3,
  common.legend = T,
  legend = "right"
)

require(grid)   # for the textGrob() function
Fig_S4<- annotate_figure(Fig_S4, bottom = textGrob("Biodiversity Reduction Time (Hr)", 
                                                 vjust = -0.75, hjust = 0.85, gp = gpar(cex = 1.5)))

ggsave("Fig_S4.png", plot = Fig_S4, dpi = 600, w = 16, h = 21)

## In-situ analysis of MBC/MBN/DOC/DON
# MBC
shan.mcse_mbc <- shan.mcse %>%
  filter(mbc != 'NA')
mod.aov <- aov(mbc ~ treatment*sample_date, data = shan.mcse_mbc)
summary(mod.aov)

# MBN
shan.mcse_mbn <- shan.mcse %>%
  filter(mbn != 'NA')
mod.aov <- aov(mbn ~ treatment*sample_date, data = shan.mcse_mbn)
summary(mod.aov)

# DOC
shan.mcse_doc <- shan.mcse %>%
  filter(doc != 'NA')
mod.aov <- aov(doc ~ treatment*sample_date, data = shan.mcse_doc)
summary(mod.aov)

# DON
shan.mcse_don <- shan.mcse %>%
  filter(don != 'NA')
mod.aov <- aov(don ~ treatment*sample_date, data = shan.mcse_don)
summary(mod.aov)

## Supplementary Figure 5: Differential Abundance Fumigation---------------------------------
bact_fum_ext <- subset_samples(bact.p_filtered, incubation == "0" | incubation == "24")

## DIFFERENTIAL ABUNDANCE: 0 Hr. vs 24Hr 
diagdds <- phyloseq_to_deseq2(bact_fum_ext, ~ incubation) 
# calculate geometric means prior to estimate size factors
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



## Supplementary Figure 6: Biodiv vs. Function on Normalized MBC data ---------------------
# Remove na nfixation
shan.fum <- shan.fum %>%
  filter(n.fix.mbc != 'NA')

# Statistical Analyses
shan.fum_pcor <- shan.fum_pcor %>%
  filter(n.fix.mbc != 'NA')

# Perform partial correlation between diazotroph alpha-diversity~n.fixation when controlling for (doc.don.mbc)
pcor.nfix.div.fum <- pcor.test(shan.fum_pcor$shan, shan.fum_pcor$n.fix.mbc, shan.fum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.fum # Significant positive association when controlling for covariance in soil biogechemistry

# How does this correlation compare to only unfumigated samples (0 hour control)
shan.unfum <- shan.fum %>%
  subset(incubation == "0")

shan.unfum_pcor <- shan.unfum %>%
  filter(doc != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(don != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(mbc != 'NA')

shan.unfum_pcor <- shan.unfum_pcor %>%
  filter(n.fix != 'NA')


pcor.nfix.div.unfum <- pcor.test(shan.unfum_pcor$shan, shan.unfum_pcor$n.fix.mbc, shan.unfum_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.unfum # No significant association
# How does the biodiversity-ecosystem function relationship vary by treatment and system - generate models
# How does fumigation time and treatment interact to influence diazotroph alpha diversity 
mod.aov <- aov(n.fix.mbc ~ shan*treatment, data = shan.fum)
summary(mod.aov)

# How does fumigation time and system interact to influence diazotroph alpha diversity
mod.aov <- aov(n.fix.mbc ~ shan*system, data = shan.fum)
summary(mod.aov)

# What is the range in diazotroph alpha diversity estimates between in-situ and in-vitro datasets?
min(shan.mcse$shan) #33.59
max(shan.mcse$shan) #442.29

min(shan.fum$shan) #37.89
max(shan.fum$shan) # 501.27
mean(shan.fum$shan) # 209.2404


min(shan.unfum$shan) # 37.89
max(shan.unfum$shan) # 501.27
mean(shan.unfum$shan) # 263.2013

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(n.fix.mbc ~ shan, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix.mbc ~ shan, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(n.fix.mbc ~ shan, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = < 0.0001", "Est = <0.0001", "Est = 0.00015" )
R <- c("Adj-R.sq = 0.0033", "Adj-R.sq = -0.022", "Adj-R.sq = 0.45")
p <- c("p = 0.3", "p = 0.5", "p<0.001")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))
# Create Plot with system-model estimates
shan.fum$treatment <- factor(shan.fum$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.fum$system <- factor(shan.fum$system, levels = c("annual", "perennial", "forest"))
shan.fum$incubation <- as.factor(shan.fum$incubation)


Fig_S6a <- shan.fum %>%
  ggplot(aes(x=shan, y = n.fix.mbc)) +
  geom_point(aes(shape = incubation, fill = treatment), size = 4.5) +  
  facet_grid(~system) +
  geom_text(x=350, y = 0.090, aes(label = est), size = 5, data = model.df) +
  geom_text(x=350, y = 0.085, aes(label = R), size = 5, data = model.df) +
  geom_text(x=350, y = 0.080, aes(label = p), size = 5, data = model.df) +
  geom_smooth(method = 'lm') +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22, 23), name = "Biodiversity Reduction Time", 
                     labels = c("0 Hours", "24 Hours", "72 hours")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Hill Number (q=1)") + 
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  ggtitle("In-vitro Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))




## In-situ Div-Fun analysis
shan.mcse <- shan.mcse %>%
  filter(n.fix.mbc != 'NA')

shan.mcse_pcor <- shan.mcse %>%
  filter(doc != 'NA')

shan.mcse_pcor <- shan.mcse_pcor %>%
  filter(don != 'NA')

shan.mcse_pcor <- shan.mcse_pcor %>%
  filter(mbc != 'NA')

shan.mcse_pcor <- shan.mcse_pcor %>%
  filter(n.fix != 'NA')


pcor.nfix.div.mcse <- pcor.test(shan.mcse_pcor$shan, shan.mcse_pcor$n.fix.mbc, shan.mcse_pcor[,c("doc", "don", "mbc")], method = "spearman")
pcor.nfix.div.mcse
# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.mcse %>%
  subset(system =="annual")
lm.annual <- lm(n.fix.mbc ~ shan, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.mcse %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix.mbc ~ shan, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.mcse %>%
  subset(system =="forest")
lm.forest <- lm(n.fix.mbc ~ shan, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = <0.001", "Est = <0.001", "Est = <0.001" )
R <- c("Adj-R.sq = 0.037", "Adj-R.sq = 0.001", "Adj-R.sq = 0.079")
p <- c("p = 0.19", "p = 0.33", "p=0.14")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))



shan.mcse$treatment <- factor(shan.mcse$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.mcse$system <- factor(shan.mcse$system, levels = c("annual", "perennial", "forest"))

Fig_S6b <- shan.mcse %>%
  ggplot(aes(x=shan, y = n.fix.mbc)) +
  geom_point(aes(fill=treatment), shape =21, size = 4.5) +  
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  ggtitle("In-situ Experiment") +
  geom_text(x=330, y = 0.10, aes(label = est), size = 5, data = model.df) +
  geom_text(x=330, y = 0.095, aes(label = R), size = 5, data = model.df) +
  geom_text(x=330, y = 0.090, aes(label = p), size = 5, data = model.df) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  xlab("Hill Number (q=1)") + 
  ylab("Normalized Nitrogen Fixation Rate (\u00b5g/g soil/MBC)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

Fig_S6 <- ggarrange(
  Fig_S6a + rremove("xlab") + rremove("ylab"), 
  Fig_S6b+ rremove("xlab") + rremove("ylab"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)

require(grid)   # for the textGrob() function
Fig_S6<- annotate_figure(Fig_S6, bottom = textGrob("Hill Number (q=1)", 
                                                  vjust = -0.75, hjust = 1.25, gp = gpar(cex = 1.5)))
Fig_S6 <- annotate_figure(Fig_S6, left = textGrob("Normalized Nitrogen Fixation Rate (\u00b5g/g soil/MBC)", 
                                                rot = 90, gp = gpar(cex = 1.5)))




ggsave("Fig_S6.png", plot = Fig_S6, dpi = 600, w = 16, h = 14)


## Supplementary Figure 7: MBC vs. Function ---------------------

# Remove na nfixation
shan.fum <- shan.fum %>%
  filter(mbc != 'NA')

# Statistical Analyses


# Perform partial correlation between diazotroph alpha-diversity~mbc when controlling for (doc.don.mbc)
pcor.mbc.div.fum <- pcor.test(shan.fum_pcor$mbc, shan.fum_pcor$n.fix, shan.fum_pcor[,c("doc", "don","shan")], method = "spearman")
pcor.mbc.div.fum # No sig correlation

# How does fumigation time and treatment interact to influence diazotroph alpha diversity 
mod.aov <- aov(n.fix ~ mbc*treatment, data = shan.fum)
summary(mod.aov)

# How does fumigation time and system interact to influence diazotroph alpha diversity
mod.aov <- aov(n.fix ~ mbc*system, data = shan.fum)
summary(mod.aov)

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.fum %>%
  subset(system =="annual")
lm.annual <- lm(n.fix ~ mbc, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.fum %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix ~ mbc, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.fum %>%
  subset(system =="forest")
lm.forest <- lm(n.fix ~ mbc, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -0.008", "Est = 0.009", "Est = < 0.0001")
R <- c("Adj-R.sq = 0.16", "Adj-R.sq = 0.42", "Adj-R.sq = -0.040")
p <- c("p = 0.01", "p = 0.0004", "p = 0.83")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))
# Create Plot with system-model estimates
shan.fum$treatment <- factor(shan.fum$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.fum$system <- factor(shan.fum$system, levels = c("annual", "perennial", "forest"))
shan.fum$incubation <- as.factor(shan.fum$incubation)


Fig_S7a <- shan.fum %>%
  ggplot(aes(x=mbc, y = n.fix)) +
  geom_point(aes(shape = incubation, fill = treatment), size = 4.5) +  
  facet_grid(~system) +
  geom_text(x=420, y = 12, aes(label = est), size = 5, data = model.df) +
  geom_text(x=420, y = 11, aes(label = R), size = 5, data = model.df) +
  geom_text(x=420, y = 10, aes(label = p), size = 5, data = model.df) +
  geom_smooth(method = 'lm') +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(21, 22, 23), name = "Biodiversity Reduction Time", 
                     labels = c("0 Hours", "24 Hours", "72 hours")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Microbial Biomass Carbon (\u00b5g/g soil)") + 
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  ggtitle("In-vitro Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))


## In-situ MBC-Fun analysis
shan.mcse <- shan.mcse %>%
  filter(mbc != 'NA')

# How does fumigation time and treatment interact to influence diazotroph alpha diversity 
mod.aov <- aov(n.fix ~ mbc*treatment, data = shan.mcse)
summary(mod.aov)

# How does fumigation time and system interact to influence diazotroph alpha diversity
mod.aov <- aov(n.fix ~ mbc*system, data = shan.mcse)
summary(mod.aov)


# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
shan_annual <- shan.mcse %>%
  subset(system =="annual")
lm.annual <- lm(n.fix ~ mbc, data = shan_annual)
summary(lm.annual)

# Perennial
shan_perennial <- shan.mcse %>%
  subset(system =="perennial")
lm.perennial <- lm(n.fix ~ mbc, data = shan_perennial)
summary(lm.perennial)

# Forest
shan_forest <- shan.mcse %>%
  subset(system =="forest")
lm.forest <- lm(n.fix ~ mbc, data = shan_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -0.004", "Est = -0.019", "Est = -0.0026" )
R <- c("Adj-R.sq = -0.019", "Adj-R.sq = 0.087", "Adj-R.sq = -0.064")
p <- c("p = 0.45", "p = 0.14", "p=0.84")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))



shan.mcse$treatment <- factor(shan.mcse$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
shan.mcse$system <- factor(shan.mcse$system, levels = c("annual", "perennial", "forest"))

Fig_S7b <- shan.mcse %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(fill=treatment), shape =21, size = 4.5) +  
  facet_grid(~system) +
  geom_smooth(method = 'lm') +
  ggtitle("In-situ Experiment") +
  geom_text(x=340, y = 18, aes(label = est), size = 5, data = model.df) +
  geom_text(x=340, y = 17, aes(label = R), size = 5, data = model.df) +
  geom_text(x=340, y = 16, aes(label = p), size = 5, data = model.df) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  xlab("Microbial Biomass Carbon (\u00b5g/g soil)") + 
  ylab("Nitrogen Fixation Rate (\u00b5g/g soil)") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12),
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))



Fig_S7 <- ggarrange(
  Fig_S7a + rremove("xlab") + rremove("ylab"), 
  Fig_S7b+ rremove("xlab") + rremove("ylab"),
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)
require(grid)   # for the textGrob() function
Fig_S7 <- annotate_figure(Fig_S7, bottom = textGrob("Microbial Biomass Carbon (\u00b5g/g soil)", 
                                                  vjust = -0.75, hjust = 0.75, gp = gpar(cex = 1.5)))
Fig_S7 <- annotate_figure(Fig_S7, left = textGrob("Nitrogen Fixation Rate (\u00b5g/g soil)", 
                                                rot = 90, gp = gpar(cex = 1.5)))



ggsave("Fig_S7.png", plot = Fig_S7, dpi = 400, w = 16, h = 14)

