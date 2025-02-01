### Code and data analysis for: [Insert Manuscript Title Here]
### DOI:
### Please contact Branetn Kristy (kristybr@msu.edu) for any questions
### All analyses are organized as presented in the manuscript.

rm(list=ls()) # Remove all objects from Environment

# MAC OS
# Note -- Change this to the GitHub Repository Directory once downloading onto your PC
setwd("~/Dropbox/Projects/LTER_MCSE/data/seq_data/") # Mac OS directory location 
#setwd("C:/Users/Brand/Dropbox/Projects/LTER_MCSE/data/seq_data/") # Windows OS directory location



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
library(ggResidpanel)
library(hillR)
library(knitr)
library(lattice)
library(lazyeval)
library(leaps)
library(lme4)
library(magrittr)
library(metacoder)
library(MetBrewer)
library(microbiome) 
library(multcompView)
library(nlme)
library(pander)
library(parallel)
library(pairwiseAdonis)
library(permute)
library(phyloseq) 
library(plyr)
library(ppcor)
library(purrr)
library(reshape2)
library(tidyverse)
library(vegan)

# Load phyloseq objects - filtered and rarefied data
load("phylo_filtered.RData") # Filtered, unrarefied data -- We will use this object for DESEQ2, which performs its own internal normalization
load("phylo_rarefied.RData") # Filtered, rarefied data -- We will use this object for all other analyses

# Subset data 
# Field Experiment
bact_mcse <- subset_samples(bact.rfy, incubation == "in-situ")
# Biodiversity Reduction Experiment
bact_microcosm <- subset_samples(bact.rfy, incubation != "in-situ")





# Field Experiment Statistics####
# Calculate alpha-diversity:
# Obtain OTU table from phyloseq object, transpose, and covert to data frame
otu_table <- t(as(bact_mcse@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

# Calculate hill number estimates for each sample at q=1 -- proxy to shannon diversity
shan_div <- hill_taxa(otu_table, q = 1)

# Add metadata from the phyloseq object to the dataframe with our hill number estimates
mcse_data <- merge(bact_mcse@sam_data, shan_div, by = 0)
# Change the name of the Hill number column from 'y' to 'shan'
colnames(mcse_data)[colnames(mcse_data) == "y"] ="shan"

# Test for variance homogeneity and normal distribution with general linear mixed effect modesl for each (Div, FLNF, MBC, MBN, EOC, ETN)
# Alpha diversity
div_glm <- lmer(shan ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(div_glm)

# FLNF
flnf_glm <- lmer(n.fix ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(flnf_glm)

# MBC
mbc_glm <- lmer(mbc ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(mbc_glm)

# MBN
mbn_glm <- lmer(mbn ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(mbn_glm)

# EOC
eoc_glm <- lmer(eoc ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(eoc_glm)

# ETN
etn_glm <- lmer(etn ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(etn_glm)

# Moisture
moisture_glm <- lmer(soil_moisture ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(moisture_glm)

# Temperature
temp_glm <- lmer(soil_temp ~ treatment * sample_date + (1 | Rep), data = mcse_data)
resid_panel(temp_glm)



### 2-way ANOVA: Management, Season, Management*Season ---------
# Alpha Diversity 
two_way_anova_diversity <- aov(shan ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_diversity)

# Obtain effect sizes: 
diversity_season_effects <- lm(shan~sample_date + 0, data = mcse_data)
summary(diversity_season_effects)
confint(diversity_season_effects)

diversity_treatment_effects <- lm(shan~treatment + 0, data = mcse_data)
summary(diversity_treatment_effects)
confint(diversity_treatment_effects)


# Post-Hoc Test
tukey_diversity <- TukeyHSD(two_way_anova_diversity, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_diversity

# FLNF
two_way_anova_flnf <- aov(n.fix ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_flnf)

flnf_season_effects <- lm(n.fix~sample_date + 0, data = mcse_data)
summary(flnf_season_effects)
confint(flnf_season_effects)

# MBC
two_way_anova_mbc <- aov(mbc ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_mbc)

# Subset summer and fall samples:
mcse_data_summer <- mcse_data %>% subset(sample_date == "summer")
mcse_data_fall <- mcse_data %>% subset(sample_date == "fall")

#Post-Hoc Test: Both Summer and Fall due to significant interaction
one_way_anova_mbc_summer <- aov(mbc ~ treatment, data = mcse_data_summer)
tukey_mbc_summer <- TukeyHSD(one_way_anova_mbc_summer, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbc_summer

one_way_anova_mbc_fall <- aov(mbc ~ treatment, data = mcse_data_fall)
tukey_mbc_fall <- TukeyHSD(one_way_anova_mbc_fall, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbc_fall

# MBN
two_way_anova_mbn <- aov(mbn ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_mbn)

#Post-Hoc Test: Both Summer and Fall due to significant interaction
one_way_anova_mbn_summer <- aov(mbn ~ treatment, data = mcse_data_summer)
tukey_mbn_summer <- TukeyHSD(one_way_anova_mbn_summer, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbn_summer

one_way_anova_mbn_fall <- aov(mbn ~ treatment, data = mcse_data_fall)
tukey_mbn_fall <- TukeyHSD(one_way_anova_mbn_fall, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbn_fall

# EOC
two_way_anova_eoc <- aov(eoc ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_eoc)

#Post-Hoc Test: Both Summer and Fall due to significant interaction
one_way_anova_eoc_summer <- aov(eoc ~ treatment, data = mcse_data_summer)
tukey_eoc_summer <- TukeyHSD(one_way_anova_eoc_summer, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_eoc_summer

one_way_anova_eoc_fall <- aov(eoc ~ treatment, data = mcse_data_fall)
tukey_eoc_fall <- TukeyHSD(one_way_anova_eoc_fall, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_eoc_fall

eoc_effects <- lm(eoc~sample_date + 0, mcse_data)
# ETN
two_way_anova_etn <- aov(etn ~ treatment * sample_date, data = mcse_data)
summary(two_way_anova_etn)

etn_effects <- lm(etn~sample_date + 0, mcse_data)

#Post-Hoc Test: Both Summer and Fall due to significant interaction
one_way_anova_etn_summer <- aov(etn ~ treatment, data = mcse_data_summer)
tukey_etn_summer <- TukeyHSD(one_way_anova_etn_summer, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_etn_summer

one_way_anova_etn_fall <- aov(etn ~ treatment, data = mcse_data_fall)
tukey_etn_fall <- TukeyHSD(one_way_anova_etn_fall, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_etn_fall

### 2-way ANOVA: Category, Season, Management*Season ---------
# Alpha Diversity 
two_way_anova_diversity_category <- aov(shan ~ system * sample_date, data = mcse_data)
summary(two_way_anova_diversity_category)

# Post-Hoc Test
tukey_diversity_category <- TukeyHSD(two_way_anova_diversity_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_diversity_category

# FLNF
two_way_anova_flnf_category <- aov(n.fix ~ system * sample_date, data = mcse_data)
summary(two_way_anova_flnf_category)

# Post-Hoc Test between Summer and Fall (due to significant interaction)
one_way_anova_flnf_category_summer <- aov(n.fix ~ system, data = mcse_data_summer)
tukey_flnf_category_summer <- TukeyHSD(one_way_anova_flnf_category_summer, "system", ordered = FALSE, conf.level = 0.95)
tukey_flnf_category_summer

one_way_anova_flnf_category_fall <- aov(n.fix ~ system, data = mcse_data_fall)
tukey_flnf_category_fall <- TukeyHSD(one_way_anova_flnf_category_fall, "system", ordered = FALSE, conf.level = 0.95)
tukey_flnf_category_fall

# MBC
two_way_anova_mbc <- aov(mbc ~ system * sample_date, data = mcse_data)
summary(two_way_anova_mbc)

tukey_mbc_category <- TukeyHSD(two_way_anova_mbc, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbc_category
  
  
# MBN
two_way_anova_mbn <- aov(mbn ~ system * sample_date, data = mcse_data)
summary(two_way_anova_mbn)

# Post-Hoc Test between Summer and Fall (due to significant interaction)
one_way_anova_mbn_category_summer <- aov(mbn ~ system, data = mcse_data_summer)
tukey_mbn_category_summer <- TukeyHSD(one_way_anova_mbn_category_summer, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbn_category_summer

one_way_anova_mbn_category_fall <- aov(mbn ~ system, data = mcse_data_fall)
tukey_mbn_category_fall <- TukeyHSD(one_way_anova_mbn_category_fall, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbn_category_fall


# EOC
two_way_anova_eoc <- aov(eoc ~ system * sample_date, data = mcse_data)
summary(two_way_anova_eoc)

# Post-Hoc Test between Summer and Fall (due to significant interaction)
one_way_anova_eoc_category_summer <- aov(eoc ~ system, data = mcse_data_summer)
tukey_eoc_category_summer <- TukeyHSD(one_way_anova_eoc_category_summer, "system", ordered = FALSE, conf.level = 0.95)
tukey_eoc_category_summer

one_way_anova_eoc_category_fall <- aov(eoc ~ system, data = mcse_data_fall)
tukey_eoc_category_fall <- TukeyHSD(one_way_anova_eoc_category_fall, "system", ordered = FALSE, conf.level = 0.95)
tukey_eoc_category_fall

# ETN
two_way_anova_etn <- aov(etn ~ system * sample_date, data = mcse_data)
summary(two_way_anova_etn)

# Post-hoc test (system and sample date)
etn_season_effects <- lm(etn ~ sample_date + 0, data = mcse_data)
summary(etn_season_effects)

etn_category_effects <- lm(etn ~ system + 0, data = mcse_data)
summary(etn_category_effects)


tukey_etn_category <- TukeyHSD(two_way_anova_etn, "system", ordered = FALSE, conf.level = 0.95)
tukey_etn_category


### FLNF: Stepwise Regression---------
mcse_data$s.soil_moisture <- as.numeric(scale(mcse_data$soil_moisture))
mcse_data$s.soil_temp <- as.numeric(scale(mcse_data$soil_temp))
mcse_data$shan <- as.numeric(mcse_data$shan)
mcse_data$n.fix <- as.numeric(mcse_data$n.fix)
mcse_data$mbc <- as.numeric(mcse_data$mbc)
mcse_data$mbn <- as.numeric(mcse_data$mbn)
mcse_data$eoc <- as.numeric(mcse_data$eoc)
mcse_data$etn <- as.numeric(mcse_data$etn)

mcse_data_stepwise_regression <- mcse_data %>%
  subset(select = c("s.soil_moisture", "s.soil_temp", "shan", "n.fix",
                    "mbc", "mbn", "eoc", "etn")) %>%
  na.omit()

full.model <- lm(n.fix ~., data = mcse_data_stepwise_regression)

# Set seed for reproducibility
set.seed(37920)
# Set up repeated k-fold cross-validation
train.control <- trainControl(method = "cv", number = 10)
# Train the model
step.model <- train(n.fix ~., data = mcse_data_stepwise_regression,
                    method = "leapBackward", 
                    tuneGrid = data.frame(nvmax = 1:7),
                    trControl = train.control
)
step.model$results

step.model$bestTune
summary(step.model$finalModel)

coef(step.model$finalModel, 2) # Soil moisture and soil temperature predict nitrogen fixation rates

mod <- lm(n.fix ~ s.soil_temp + s.soil_moisture, data = mcse_data_stepwise_regression)
summary(mod)
confint(mod)

### permANOVA: Management, Season ---------------------------------------------------------
set.seed(37920)
# Convert to relative abundance
bact_mcse_rel <- transform_sample_counts(bact_mcse, function(x) x / sum(x) )
# Perform NMDS ordination with Bray curtis distance using the ordinate() function
bray_NMDS_mcse_ord <- ordinate(bact_mcse_rel, method = "NMDS", distance = "bray") # Stress = 0.185

# Calculate Bray-curtis distance using phyloseq::distance to identify significance of our factors: treatment and season/sample date
set.seed(37920)
bray_NMDS_mcse <- phyloseq::distance(bact_mcse_rel,"bray")

# What are the interactive effects of sampling date and treatment on diazotroph Beta-diversity? -- perform permanova
set.seed(37920)
permanova <- adonis(bray_NMDS_mcse~sample_data(bact_mcse_rel)$treatment *
                      sample_data(bact_mcse_rel)$sample_date, permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication

# Pairwise adonis2 between each season 
# Summer
bact_mcse_rel_summer <- subset_samples(bact_mcse_rel, sample_date == "summer")
set.seed(37920)
bray_NMDS_mcse_summer <- phyloseq::distance(bact_mcse_rel_summer,"bray")

pairwise.adonis(bray_NMDS_mcse_summer, sample_data(bact_mcse_rel_summer)$treatment)

# Fall
bact_mcse_rel_fall <- subset_samples(bact_mcse_rel, sample_date == "fall")
set.seed(37920)
bray_NMDS_mcse_fall <- phyloseq::distance(bact_mcse_rel_fall,"bray")
pairwise.adonis(bray_NMDS_mcse_fall, sample_data(bact_mcse_rel_fall)$treatment)

# Create envfit model to visualize effects of MBC/eoc/etn as vectors on the NMDS ordination
# Exctract sample metadata from the rarefied phyloseq object 
env <- as.data.frame(sample_data(bact_mcse))
set.seed(37920)
# Perform the envfit() function to quantify effects of continuous variables on community composition
ord.fit.sig <- envfit(bray_NMDS_mcse_ord ~ env$soil_temp + 
                        env$soil_moisture +
                        env$n.fix + env$mbc + env$mbn + 
                        env$eoc + env$etn, 
                      perm = 9999, na.rm = TRUE)

### permANOVA: Category, Season ---------------------------------------------------------
# What are the interactive effects of sampling date and treatment on diazotroph Beta-diversity? -- perform permanova
set.seed(37920)
permanova <- adonis(bray_NMDS_mcse~sample_data(bact_mcse_rel)$system *
                      sample_data(bact_mcse_rel)$sample_date, permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication

# Pairwise adonis2 between each season 
# Summer
bact_mcse_rel_summer <- subset_samples(bact_mcse_rel, sample_date == "summer")
set.seed(37920)
bray_NMDS_mcse_summer <- phyloseq::distance(bact_mcse_rel_summer,"bray")

pairwise.adonis(bray_NMDS_mcse_summer, sample_data(bact_mcse_rel_summer)$system)
# Fall
bact_mcse_rel_fall <- subset_samples(bact_mcse_rel, sample_date == "fall")
set.seed(37920)
bray_NMDS_mcse_fall <- phyloseq::distance(bact_mcse_rel_fall,"bray")
pairwise.adonis(bray_NMDS_mcse_fall, sample_data(bact_mcse_rel_fall)$system)

#### Partial Correlation Analysis: FLNF vs. Diversity ------------------------------
# Perform partial correlation between diazotroph alpha-diversity~n.fixation when controlling for (eoc.etn.mbc)
mcse_data_nona <- mcse_data %>% drop_na(eoc)

pcor_fixation_diversity_mcse <- pcor.test(mcse_data_nona$shan, mcse_data_nona$n.fix, mcse_data_nona[,c("eoc", "etn", "mbc", "mbn")], method = "spearman")

pcor_fixation_diversity_mcse
#### DESEQ --------------------------------------------------------------------
##### Summer Vs. Fall ----------------------------------------------------------
bact_mcse_unrarefied <- subset_samples(bact.p_filtered, incubation == "in-situ")

diagdds_summer_v_fall <- phyloseq_to_deseq2(bact_mcse_unrarefied, ~ sample_date) 

# calculate geometric means prior to estimate size factors
# gm_mean function
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

geoMeans <-  apply(counts(diagdds_summer_v_fall), 1, gm_mean)
diagdds_summer_v_fall <- estimateSizeFactors(diagdds_summer_v_fall, geoMeans = geoMeans)
diagdds_summer_v_fall <- DESeq(diagdds_summer_v_fall, fitType="parametric")

results_summer_v_fall <- results(diagdds_summer_v_fall, alpha = 0.01, contrast = c("sample_date", "summer","fall"))

sigtab_results_summer_v_fall <- results_summer_v_fall[which(results_summer_v_fall$padj < 0.01), ]

sigtab_results_TAX_summer_v_fall <- cbind(as(sigtab_results_summer_v_fall, "data.frame"), as(tax_table(bact_mcse_unrarefied)[row.names(sigtab_results_summer_v_fall), ], "matrix"))
summary(sigtab_results_summer_v_fall)


sigtab_results_TAX_genra_summer_v_fall <- subset(sigtab_results_TAX_summer_v_fall, !is.na(Genus))
# Order Genus 
x <- tapply(sigtab_results_TAX_genra_summer_v_fall$log2FoldChange, sigtab_results_TAX_genra_summer_v_fall$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_summer_v_fall$Class <- factor(as.character(sigtab_results_TAX_genra_summer_v_fall$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra_summer_v_fall$log2FoldChange, sigtab_results_TAX_genra_summer_v_fall$Genus, function(x) max(x))
x <- sort(x, TRUE)

sigtab_results_TAX_genra_summer_v_fall$Genus <- factor(as.character(sigtab_results_TAX_genra_summer_v_fall$Genus), levels=names(x))

#####Annual + Perennial vs. Forest ---------------------------------------
# Summer
bact_mcse_unrarefied_summer <- subset_samples(bact_mcse_unrarefied, sample_date == "summer")
## DIFFERENTIAL ABUNDANCE: FORESTED SITES VS. MCSE SITES
sample_data(bact_mcse_unrarefied_summer)$agro <- ifelse(bact_mcse_unrarefied_summer@sam_data$treatment == "SF" |
                                        bact_mcse_unrarefied_summer@sam_data$treatment == "DF" |
                                        bact_mcse_unrarefied_summer@sam_data$treatment == "CF" , "Forest", "Agro")

diagdds_mcse_vs_forest_summer <- phyloseq_to_deseq2(bact_mcse_unrarefied_summer, ~ agro) 

# calculate geometric means prior to estimate size factors
geoMeans <- apply(counts(diagdds_mcse_vs_forest_summer), 1, gm_mean)
diagdds_mcse_vs_forest_summer <- estimateSizeFactors(diagdds_mcse_vs_forest_summer, geoMeans = geoMeans)
diagdds_mcse_vs_forest_summer <- DESeq(diagdds_mcse_vs_forest_summer, fitType="parametric")

results_mcse_vs_forest_summer <- results(diagdds_mcse_vs_forest_summer, alpha = 0.01, contrast = c("agro", "Forest","Agro"))

sigtab_results_mcse_vs_forest_summer <- results_mcse_vs_forest_summer[which(results_mcse_vs_forest_summer$padj < 0.01), ]
summary(sigtab_results_mcse_vs_forest_summer)


sigtab_results_TAX_mcse_vs_forest_summer <- cbind(as(sigtab_results_mcse_vs_forest_summer, "data.frame"), as(tax_table(bact_mcse_unrarefied_summer)[row.names(sigtab_results_mcse_vs_forest_summer), ], "matrix"))

sigtab_results_TAX_genra_mcse_vs_forest_summer <- subset(sigtab_results_TAX_mcse_vs_forest_summer, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra_mcse_vs_forest_summer$log2FoldChange, sigtab_results_TAX_genra_mcse_vs_forest_summer$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_mcse_vs_forest_summer$Class <- factor(as.character(sigtab_results_TAX_genra_mcse_vs_forest_summer$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra_mcse_vs_forest_summer$log2FoldChange, sigtab_results_TAX_genra_mcse_vs_forest_summer$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_mcse_vs_forest_summer$Genus <- factor(as.character(sigtab_results_TAX_genra_mcse_vs_forest_summer$Genus), levels=names(x))
sigtab_results_TAX_genra_mcse_vs_forest_summer$sample_date <- "summer"


# Fall
bact_mcse_unrarefied_fall <- subset_samples(bact_mcse_unrarefied, sample_date == "fall")

sample_data(bact_mcse_unrarefied_fall)$agro <- ifelse(bact_mcse_unrarefied_fall@sam_data$treatment == "SF" |
                                                          bact_mcse_unrarefied_fall@sam_data$treatment == "DF" |
                                                          bact_mcse_unrarefied_fall@sam_data$treatment == "CF" , "Forest", "Agro")

diagdds_mcse_vs_forest_fall <- phyloseq_to_deseq2(bact_mcse_unrarefied_fall, ~ agro) 

# calculate geometric means prior to estimate size factors
geoMeans <- apply(counts(diagdds_mcse_vs_forest_fall), 1, gm_mean)
diagdds_mcse_vs_forest_fall <- estimateSizeFactors(diagdds_mcse_vs_forest_fall, geoMeans = geoMeans)
diagdds_mcse_vs_forest_fall <- DESeq(diagdds_mcse_vs_forest_fall, fitType="parametric")

results_mcse_vs_forest_fall <- results(diagdds_mcse_vs_forest_fall, alpha = 0.01, contrast = c("agro", "Forest","Agro"))

sigtab_results_mcse_vs_forest_fall <- results_mcse_vs_forest_fall[which(results_mcse_vs_forest_fall$padj < 0.01), ]
summary(sigtab_results_mcse_vs_forest_fall)

sigtab_results_TAX_mcse_vs_forest_fall <- cbind(as(sigtab_results_mcse_vs_forest_fall, "data.frame"), as(tax_table(bact_mcse_unrarefied_fall)[row.names(sigtab_results_mcse_vs_forest_fall), ], "matrix"))

sigtab_results_TAX_genra_mcse_vs_forest_fall <- subset(sigtab_results_TAX_mcse_vs_forest_fall, !is.na(Genus))


# Order Genus 
x <- tapply(sigtab_results_TAX_genra_mcse_vs_forest_fall$log2FoldChange, sigtab_results_TAX_genra_mcse_vs_forest_fall$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_mcse_vs_forest_fall$Class <- factor(as.character(sigtab_results_TAX_genra_mcse_vs_forest_fall$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra_mcse_vs_forest_fall$log2FoldChange, sigtab_results_TAX_genra_mcse_vs_forest_fall$Genus, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_mcse_vs_forest_fall$Genus <- factor(as.character(sigtab_results_TAX_genra_mcse_vs_forest_fall$Genus), levels=names(x))
sigtab_results_TAX_genra_mcse_vs_forest_fall$sample_date <- "fall"

# Combine both DESEQ objects into one dataframe for the plot 
sigtab_results_TAX_genra_mcse_vs_forest_combined <- rbind(sigtab_results_TAX_genra_mcse_vs_forest_summer,
                                                          sigtab_results_TAX_genra_mcse_vs_forest_fall)




# Field Experiment Figures ####
## Diversity ####
# Obtain mean alpha diversity values for each system (annual, perennial, forest) in the summer sampling date 
lm_diversity_category_summer <- lm(shan ~ system + 0, data = mcse_data_summer)
summary(lm_diversity_category_summer)

lm_diversity_category_fall <- lm(shan ~ system + 0, data = mcse_data_fall)
summary(lm_diversity_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(244.85, 168.82, 87.54)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_diversity_effects_summer <- data.frame(est, system, sample_date)

est <- c(281.24, 212.20, 193.25)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_diversity_effects_fall <- data.frame(est, system, sample_date)

merge_diversity_effects <- list(category_diversity_effects_summer,
                                category_diversity_effects_fall)

merge_diversity_effects <- merge_diversity_effects %>% 
  purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_diversity_summer <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_diversity_summer)[1] <- "Letter" #Reassign column name
letters_df_diversity_summer$treatment <- rownames(letters_df_diversity_summer) 
letters_df_diversity_summer$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_diversity_summer$system <- value_map[letters_df_diversity_summer$treatment]
letters_df_diversity_summer$system <- factor(letters_df_diversity_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_diversity_summer$sample_date <- value_map_summer[letters_df_diversity_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_diversity_fall <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_diversity_fall)[1] <- "Letter" #Reassign column name
letters_df_diversity_fall$treatment <- rownames(letters_df_diversity_fall) 
letters_df_diversity_fall$placement <- c(475, 475, 475, 475, 475, 475, 475, 475)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                      "T3" = "fall",
                      "T4" = "fall",
                      "T5" = "fall",
                      "T7" = "fall",
                      "DF" = "fall",
                      "SF" = "fall",
                      "CF" = "fall")
letters_df_diversity_fall$system <- value_map[letters_df_diversity_fall$treatment]
letters_df_diversity_fall$system <- factor(letters_df_diversity_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_diversity_fall$sample_date <- value_map_fall[letters_df_diversity_fall$treatment]

# Merge all letters together
merge_tukey_diversity_mcse <- list(letters_df_diversity_summer,
                                   letters_df_diversity_fall)
merge_tukey_diversity_mcse <- merge_tukey_diversity_mcse  %>% purrr::reduce(full_join)
  
# Reorder Treatments
mcse_data$treatment <- factor(mcse_data$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
mcse_data$system <- factor(mcse_data$system, levels = c("annual", "perennial", "forest"))
merge_diversity_effects$system <- factor(merge_diversity_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_diversity_mcse$system <- factor(merge_tukey_diversity_mcse$system, levels = c("annual", "perennial", "forest"))

mcse_data$sample_date <- factor(mcse_data$sample_date, levels = c("summer", "fall"))
merge_diversity_effects$sample_date <- factor(merge_diversity_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_diversity_mcse$sample_date <- factor(merge_tukey_diversity_mcse$sample_date, levels = c("summer", "fall"))


# Plot alpha diversity across all treatments in the Summer -- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_diversity_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_diversity_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_diversity_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab("Hill Number (q=1)") +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_diversity_plot


## FLNF ####
# Obtain mean flnf values for each system (annual, perennial, forest)
lm_flnf_category_summer <- lm(n.fix ~ system + 0, data = mcse_data_summer)
summary(lm_flnf_category_summer)

lm_flnf_category_fall <- lm(n.fix ~ system + 0, data = mcse_data_fall)
summary(lm_flnf_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(3.5283, 2.6012, 1.4306)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_flnf_effects_summer <- data.frame(est, system, sample_date)

est <- c(5.004, 9.550, 2.940)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_flnf_effects_fall <- data.frame(est, system, sample_date)

merge_flnf_effects <- list(category_flnf_effects_summer,
                                category_flnf_effects_fall)
merge_flnf_effects <- merge_flnf_effects %>% purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_flnf_summer <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_flnf_summer)[1] <- "Letter" #Reassign column name
letters_df_flnf_summer$treatment <- rownames(letters_df_flnf_summer) 
letters_df_flnf_summer$placement <- c(15, 15, 15, 15, 15, 15, 15, 15)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_flnf_summer$system <- value_map[letters_df_flnf_summer$treatment]
letters_df_flnf_summer$system <- factor(letters_df_flnf_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_flnf_summer$sample_date <- value_map_summer[letters_df_flnf_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_flnf_fall <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_flnf_fall)[1] <- "Letter" #Reassign column name
letters_df_flnf_fall$treatment <- rownames(letters_df_flnf_fall) 
letters_df_flnf_fall$placement <- c(15, 15, 15, 15, 15, 15, 15, 15)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                    "T3" = "fall",
                    "T4" = "fall",
                    "T5" = "fall",
                    "T7" = "fall",
                    "DF" = "fall",
                    "SF" = "fall",
                    "CF" = "fall")
letters_df_flnf_fall$system <- value_map[letters_df_flnf_fall$treatment]
letters_df_flnf_fall$system <- factor(letters_df_flnf_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_flnf_fall$sample_date <- value_map_fall[letters_df_flnf_fall$treatment]

# Merge all letters together
merge_tukey_flnf_mcse <- list(letters_df_flnf_summer,
                                   letters_df_flnf_fall)
merge_tukey_flnf_mcse <- merge_tukey_flnf_mcse  %>% purrr::reduce(full_join)

# Reorder Treatments
merge_flnf_effects$system <- factor(merge_flnf_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_flnf_mcse$system <- factor(merge_tukey_flnf_mcse$system, levels = c("annual", "perennial", "forest"))

merge_flnf_effects$sample_date <- factor(merge_flnf_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_flnf_mcse$sample_date <- factor(merge_tukey_flnf_mcse$sample_date, levels = c("summer", "fall"))


# Plot flnf across all treatments-- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_flnf_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_flnf_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_flnf_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_flnf_plot

## MBC ####
# Obtain mean flnf values for each system (annual, perennial, forest)
lm_mbc_category_summer <- lm(mbc ~ system + 0, data = mcse_data_summer)
summary(lm_mbc_category_summer)

lm_mbc_category_fall <- lm(mbc ~ system + 0, data = mcse_data_fall)
summary(lm_mbc_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(240.00, 471.57, 262.63)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_mbc_effects_summer <- data.frame(est, system, sample_date)

est <- c(299.80, 437.88, 315.14)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_mbc_effects_fall <- data.frame(est, system, sample_date)

merge_mbc_effects <- list(category_mbc_effects_summer,
                           category_mbc_effects_fall)
merge_mbc_effects <- merge_mbc_effects %>% purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_mbc_summer <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_mbc_summer)[1] <- "Letter" #Reassign column name
letters_df_mbc_summer$treatment <- rownames(letters_df_mbc_summer) 
letters_df_mbc_summer$placement <- c(700, 700, 700, 700, 700, 700, 700, 700)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_mbc_summer$system <- value_map[letters_df_mbc_summer$treatment]
letters_df_mbc_summer$system <- factor(letters_df_mbc_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_mbc_summer$sample_date <- value_map_summer[letters_df_mbc_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_mbc_fall <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_mbc_fall)[1] <- "Letter" #Reassign column name
letters_df_mbc_fall$treatment <- rownames(letters_df_mbc_fall) 
letters_df_mbc_fall$placement <- c(615, 615, 615, 615, 615, 615, 615, 615)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                    "T3" = "fall",
                    "T4" = "fall",
                    "T5" = "fall",
                    "T7" = "fall",
                    "DF" = "fall",
                    "SF" = "fall",
                    "CF" = "fall")
letters_df_mbc_fall$system <- value_map[letters_df_mbc_fall$treatment]
letters_df_mbc_fall$system <- factor(letters_df_mbc_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_mbc_fall$sample_date <- value_map_fall[letters_df_mbc_fall$treatment]

# Merge all letters together
merge_tukey_mbc_mcse <- list(letters_df_mbc_summer,
                              letters_df_mbc_fall)
merge_tukey_mbc_mcse <- merge_tukey_mbc_mcse  %>% purrr::reduce(full_join)

# Reorder Treatments
merge_mbc_effects$system <- factor(merge_mbc_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_mbc_mcse$system <- factor(merge_tukey_mbc_mcse$system, levels = c("annual", "perennial", "forest"))

merge_mbc_effects$sample_date <- factor(merge_mbc_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_mbc_mcse$sample_date <- factor(merge_tukey_mbc_mcse$sample_date, levels = c("summer", "fall"))


# Plot mbc across all treatments-- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_mbc_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = mbc, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_mbc_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_mbc_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(MBC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_mbc_plot


## MBN ####
# Obtain mean flnf values for each system (annual, perennial, forest)
lm_mbn_category_summer <- lm(mbn ~ system + 0, data = mcse_data_summer)
summary(lm_mbn_category_summer)

lm_mbn_category_fall <- lm(mbn ~ system + 0, data = mcse_data_fall)
summary(lm_mbn_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(32.100, 81.211, 37.941)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_mbn_effects_summer <- data.frame(est, system, sample_date)

est <- c(41.823, 70.536, 36.169)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_mbn_effects_fall <- data.frame(est, system, sample_date)

merge_mbn_effects <- list(category_mbn_effects_summer,
                          category_mbn_effects_fall)
merge_mbn_effects <- merge_mbn_effects %>% purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_mbn_summer <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_mbn_summer)[1] <- "Letter" #Reassign column name
letters_df_mbn_summer$treatment <- rownames(letters_df_mbn_summer) 
letters_df_mbn_summer$placement <- c(120, 120, 120, 120, 120, 120, 120, 120)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_mbn_summer$system <- value_map[letters_df_mbn_summer$treatment]
letters_df_mbn_summer$system <- factor(letters_df_mbn_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_mbn_summer$sample_date <- value_map_summer[letters_df_mbn_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_mbn_fall <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_mbn_fall)[1] <- "Letter" #Reassign column name
letters_df_mbn_fall$treatment <- rownames(letters_df_mbn_fall) 
letters_df_mbn_fall$placement <- c(120, 120, 120, 120, 120, 120, 120, 120)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                    "T3" = "fall",
                    "T4" = "fall",
                    "T5" = "fall",
                    "T7" = "fall",
                    "DF" = "fall",
                    "SF" = "fall",
                    "CF" = "fall")
letters_df_mbn_fall$system <- value_map[letters_df_mbn_fall$treatment]
letters_df_mbn_fall$system <- factor(letters_df_mbn_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_mbn_fall$sample_date <- value_map_fall[letters_df_mbn_fall$treatment]

# Merge all letters together
merge_tukey_mbn_mcse <- list(letters_df_mbn_summer,
                             letters_df_mbn_fall)
merge_tukey_mbn_mcse <- merge_tukey_mbn_mcse  %>% purrr::reduce(full_join)

# Reorder Treatments
merge_mbn_effects$system <- factor(merge_mbn_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_mbn_mcse$system <- factor(merge_tukey_mbn_mcse$system, levels = c("annual", "perennial", "forest"))

merge_mbn_effects$sample_date <- factor(merge_mbn_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_mbn_mcse$sample_date <- factor(merge_tukey_mbn_mcse$sample_date, levels = c("summer", "fall"))


# Plot mbn across all treatments-- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_mbn_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = mbn, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_mbn_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_mbn_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(MBN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_mbn_plot

## EOC ####
# Obtain mean eoc values for each system (annual, perennial, forest)
# Use removed outliers dataframe for EOC due to high T1R1 estimates 
lm_eoc_category_summer <- lm(eoc ~ system + 0, data = mcse_data_summer)
summary(lm_eoc_category_summer)

lm_eoc_category_fall <- lm(eoc ~ system + 0, data = mcse_data_fall)
summary(lm_eoc_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(94.517, 99.599, 116.055)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_eoc_effects_summer <- data.frame(est, system, sample_date)

est <- c(47.474, 52.099, 46.124)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_eoc_effects_fall <- data.frame(est, system, sample_date)

merge_eoc_effects <- list(category_eoc_effects_summer,
                          category_eoc_effects_fall)
merge_eoc_effects <- merge_eoc_effects %>% purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_eoc_summer <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_eoc_summer)[1] <- "Letter" #Reassign column name
letters_df_eoc_summer$treatment <- rownames(letters_df_eoc_summer) 
letters_df_eoc_summer$placement <- c(220, 220, 220, 220, 220, 220, 220, 220)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_eoc_summer$system <- value_map[letters_df_eoc_summer$treatment]
letters_df_eoc_summer$system <- factor(letters_df_eoc_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_eoc_summer$sample_date <- value_map_summer[letters_df_eoc_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_eoc_fall <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_eoc_fall)[1] <- "Letter" #Reassign column name
letters_df_eoc_fall$treatment <- rownames(letters_df_eoc_fall) 
letters_df_eoc_fall$placement <- c(90, 90, 90, 90, 90, 90, 90, 90)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                    "T3" = "fall",
                    "T4" = "fall",
                    "T5" = "fall",
                    "T7" = "fall",
                    "DF" = "fall",
                    "SF" = "fall",
                    "CF" = "fall")
letters_df_eoc_fall$system <- value_map[letters_df_eoc_fall$treatment]
letters_df_eoc_fall$system <- factor(letters_df_eoc_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_eoc_fall$sample_date <- value_map_fall[letters_df_eoc_fall$treatment]

# Merge all letters together
merge_tukey_eoc_mcse <- list(letters_df_eoc_summer,
                             letters_df_eoc_fall)
merge_tukey_eoc_mcse <- merge_tukey_eoc_mcse  %>% purrr::reduce(full_join)

# Reorder Treatments
merge_eoc_effects$system <- factor(merge_eoc_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_eoc_mcse$system <- factor(merge_tukey_eoc_mcse$system, levels = c("annual", "perennial", "forest"))

merge_eoc_effects$sample_date <- factor(merge_eoc_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_eoc_mcse$sample_date <- factor(merge_tukey_eoc_mcse$sample_date, levels = c("summer", "fall"))


# Plot eoc across all treatments-- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_eoc_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = eoc, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_eoc_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_eoc_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(EOC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_eoc_plot

## ETN ####
# Obtain mean etn values for each system (annual, perennial, forest)
lm_etn_category_summer <- lm(etn ~ system + 0, data = mcse_data_summer)
summary(lm_etn_category_summer)

lm_etn_category_fall <- lm(etn ~ system + 0, data = mcse_data_fall)
summary(lm_etn_category_fall)

# Add model estimates into dataframe for GGPLOT
est <- c(10.945, 6.489, 15.216)
system <- c("annual", "perennial", "forest")
sample_date <- c("summer", "summer", "summer")
category_etn_effects_summer <- data.frame(est, system, sample_date)

est <- c(16.2137, 14.7978, 15.0729)
system <- c("annual", "perennial", "forest")
sample_date <- c("fall", "fall", "fall")
category_etn_effects_fall <- data.frame(est, system, sample_date)

merge_etn_effects <- list(category_etn_effects_summer,
                          category_etn_effects_fall)
merge_etn_effects <- merge_etn_effects %>% purrr::reduce(full_join)


# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Summer Samples 
letters_df_etn_summer <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ treatment, data = mcse_data_summer))$treatment[,4])$Letters)

colnames(letters_df_etn_summer)[1] <- "Letter" #Reassign column name
letters_df_etn_summer$treatment <- rownames(letters_df_etn_summer) 
letters_df_etn_summer$placement <- c(30, 30, 30, 30, 30, 30, 30, 30)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_summer <- c("T1" = "summer",
                      "T3" = "summer",
                      "T4" = "summer",
                      "T5" = "summer",
                      "T7" = "summer",
                      "DF" = "summer",
                      "SF" = "summer",
                      "CF" = "summer")
letters_df_etn_summer$system <- value_map[letters_df_etn_summer$treatment]
letters_df_etn_summer$system <- factor(letters_df_etn_summer$system, levels = c("annual", "perennial", "forest"))
letters_df_etn_summer$sample_date <- value_map_summer[letters_df_etn_summer$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Fall Samples 
letters_df_etn_fall <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ treatment, data = mcse_data_fall))$treatment[,4])$Letters)

colnames(letters_df_etn_fall)[1] <- "Letter" #Reassign column name
letters_df_etn_fall$treatment <- rownames(letters_df_etn_fall) 
letters_df_etn_fall$placement <- c(30, 30, 30, 30, 30, 30, 30, 30)
value_map <- c("T1" = "annual",
               "T3" = "annual",
               "T4" = "annual",
               "T5" = "perennial",
               "T7" = "perennial",
               "DF" = "forest",
               "SF" = "forest",
               "CF" = "forest")
value_map_fall <- c("T1" = "fall",
                    "T3" = "fall",
                    "T4" = "fall",
                    "T5" = "fall",
                    "T7" = "fall",
                    "DF" = "fall",
                    "SF" = "fall",
                    "CF" = "fall")
letters_df_etn_fall$system <- value_map[letters_df_etn_fall$treatment]
letters_df_etn_fall$system <- factor(letters_df_etn_fall$system, levels = c("annual", "perennial", "forest"))
letters_df_etn_fall$sample_date <- value_map_fall[letters_df_etn_fall$treatment]

# Merge all letters together
merge_tukey_etn_mcse <- list(letters_df_etn_summer,
                             letters_df_etn_fall)
merge_tukey_etn_mcse <- merge_tukey_etn_mcse  %>% purrr::reduce(full_join)

# Reorder Treatments
merge_etn_effects$system <- factor(merge_etn_effects$system, levels = c("annual", "perennial", "forest"))
merge_tukey_etn_mcse$system <- factor(merge_tukey_etn_mcse$system, levels = c("annual", "perennial", "forest"))

merge_etn_effects$sample_date <- factor(merge_etn_effects$sample_date, levels = c("summer", "fall"))
merge_tukey_etn_mcse$sample_date <- factor(merge_tukey_etn_mcse$sample_date, levels = c("summer", "fall"))


# Plot etn across all treatments-- with Post-hoc comparisons across treatment and system-model mean estimates
mcse_etn_plot <- mcse_data %>%
  ggplot(aes(x=treatment, y = etn, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  geom_text(data = merge_tukey_etn_mcse, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  facet_grid(sample_date ~ system, scales = 'free_x') +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_etn_effects) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(ETN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_etn_plot

# Ordination ####
# Obtain NMDS point coordinates for each sample, and append relevant metadata from the env dataframe 
data.scores <- as.data.frame(scores(bray_NMDS_mcse_ord)$sites)
data.scores$treatment <- env$treatment
data.scores$system <- env$system
data.scores$treatment <- env$treatment
data.scores$sample_date <- env$sample_date
data.scores$soil_temp <- env$soil_temp
data.scores$soil_moisture <- env$soil_moisture
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
data.scores$mbn <- env$mbn
data.scores$eoc <- env$eoc
data.scores$etn <- env$etn

# Extract vector coordinates of metadata variables 
en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors"))
# Extract p-values of each vector
en_coord_cont$pval <- ord.fit.sig$vectors$pvals
# Subset vectors that are statistically significant (p <0.05)
en_coord_cont_sig <- en_coord_cont %>% filter(pval < 0.05)
# Add names of relevant vectors to the dataframe -- for the ggplot
en_coord_cont_sig$name <- c("Soil.Moisture", "FLNF", "EOC")

# Add ellipse for system level variable
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the system factor
data.scores$system <- as.factor(data.scores$system)
df_ell.system <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores$system)){
  df_ell.system <- rbind(df_ell.system, cbind(as.data.frame(with(data.scores [data.scores$system==g,],
                                                                                   veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                                          wt=rep(1/length(NMDS1),
                                                                                                                 length(NMDS1)))$cov,
                                                                                                   center=c(mean(NMDS1),
                                                                                                            mean(NMDS2))))),system=g))
}

# data for labelling the ellipse
NMDS.mean.system<-aggregate(data.scores[ ,c("NMDS1", "NMDS2")], 
                         list(group = data.scores$system), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(data.scores[,c("NMDS1", "NMDS2")], 
                    list(group = data.scores$system), mean)






# Reorganize treatment factor 
data.scores$treatment <- factor(data.scores$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
data.scores$system <- factor(data.scores$system, levels = c("annual", "perennial", "forest"))


df_ell.system$system <- factor(df_ell.system$system, levels = c("annual", "perennial", "forest"))

# Display NMDS ordination of our in-situ dataset
mcse_nmds_ordination <- data.scores %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(color = system, shape = sample_date), size = 5.5, alpha = 1.0) +
  geom_path(data = df_ell.system, aes(x = NMDS1, y = NMDS2, group = system, color = system), linetype = 1, size = 1, alpha = 1.0) +
  geom_polygon(data = df_ell.system, aes(x = NMDS1, y = NMDS2, fill = system), alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = -0.35, y = 0.16, label = "stress = 0.185") +
  scale_color_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                     labels = c("Annual", "Perennial", "Forest")) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                     labels = c("Annual", "Perennial", "Forest")) +
  scale_shape_manual(values =c(16, 17), name = "Sample Date", 
                     labels = c("Fall", "Summer")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont_sig, size =1.25, alpha = 0.7, colour = "black") +
  geom_label(data = en_coord_cont_sig, aes(x = NMDS1, y = NMDS2), colour = "black", 
             fontface = "bold", size = 5.5, label = en_coord_cont_sig$name, fill = "white") + 
  guides(fill = guide_legend(override.aes = list(shape=c(16)))) +
  ggtitle("Field Experiment") +
  xlim(-0.4, 0.4) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 12, color = "black"))  + 
  theme(axis.text.y = element_text(size = 12, color = "black"))  + 
  theme(legend.text = element_text(size = 12)) +
  theme(plot.margin=unit(c(.5,1,.5,.5),"cm")) +
  theme(legend.position="right", legend.title=element_text(size=14))+
  theme(legend.position="right")

mcse_nmds_ordination

# Differential Abundance ####
### Summer Vs. Fall ####
mcse_summer_vs_fall_plot <- sigtab_results_TAX_genra_summer_v_fall %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  theme_classic() +
  geom_point(size = 5.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth = 0.95) +
  geom_text(x=2.3, y=1.5, label="Summer-Enriched", size = 5) +
  scale_fill_manual(values=met.brewer("Egypt", 2)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =20), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15))
mcse_summer_vs_fall_plot

### Annual + Perennial Vs. Forest Summer ####
mcse_vs_forest_plot <- sigtab_results_TAX_genra_mcse_vs_forest_combined %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  facet_grid(~sample_date, scales = "free", space = "free") +
  theme_bw() +
  geom_point(size = 4.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth = 0.95) +
  scale_fill_manual(values=met.brewer("Egypt", 5)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =14), 
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        plot.margin=unit(c(.5,1,.5,.5),"cm"),
        legend.title = element_text(size =15))
mcse_vs_forest_plot
# FLNF vs. Diversity Regression ####
mcse_flnf_div <- mcse_data %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(color = treatment), size = 4.5) +  
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Hill Number (q=1)") + 
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  ggtitle("Field Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
mcse_flnf_div

# Biodiversity Reduction Experiment Statistics ####
# Calculate alpha-diversity:
# Obtain OTU table from phyloseq object, transpose, and covert to data frame
otu_table <- t(as(bact_microcosm@otu_table, "matrix"))
otu_table <- as.data.frame(otu_table)

# Calculate hill number estimates for each sample at q=1 -- proxy to shannon diversity
shan_div <- hill_taxa(otu_table, q = 1)
# Add metadata from the phyloseq object to the dataframe with our hill number estimates
microcosm_data <- merge(bact_microcosm@sam_data, shan_div, by = 0)
# Change the name of the Hill number column from 'y' to 'shan'
colnames(microcosm_data)[colnames(microcosm_data) == "y"] ="shan"

# Change the name of the incubation factors 
microcosm_data$incubation <- as.factor(microcosm_data$incubation)
levels(microcosm_data$incubation) <- c("control", "moderate", "high")

# Test for variance homogeneity and normal distribution
# Alpha diversity
div_glm <- lmer(shan ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(div_glm)

# FLNF
flnf_glm <- lmer(n.fix ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(flnf_glm)

# MBC
mbc_glm <- lmer(mbc ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(mbc_glm)

# MBN
mbn_glm <- lmer(mbn ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(mbn_glm)

# EOC
eoc_glm <- lmer(eoc ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(eoc_glm)

# ETN
etn_glm <- lmer(etn ~ treatment * incubation + (1 | Rep), data = microcosm_data)
resid_panel(etn_glm)

#### 2-way ANOVA: Management, Reduction, Management*Reduction ####
# Alpha Diversity
two_way_anova_diversity_microcosm <- aov(shan ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_diversity_microcosm)





# Post-Hoc testing: Both treatment and incubation
tukey_diversity_microcosm_treatment <- TukeyHSD(two_way_anova_diversity_microcosm, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_diversity_microcosm_treatment

tukey_diversity_microcosm_incubation <- TukeyHSD(two_way_anova_diversity_microcosm, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_diversity_microcosm_incubation

# FLNF 
two_way_anova_flnf_microcosm <- aov(n.fix ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_flnf_microcosm)

# Post-Hoc testing: Both treatment and incubation
tukey_flnf_microcosm_treatment <- TukeyHSD(two_way_anova_flnf_microcosm, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_flnf_microcosm_treatment

tukey_flnf_microcosm_incubation <- TukeyHSD(two_way_anova_flnf_microcosm, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_flnf_microcosm_incubation

# MBC
two_way_anova_mbc_microcosm <- aov(mbc ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_mbc_microcosm)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
microcosm_data_control <- microcosm_data %>% subset(incubation == "control")
microcosm_data_moderate <- microcosm_data %>% subset(incubation == "moderate")
microcosm_data_high <- microcosm_data %>% subset(incubation == "high")

one_way_anova_mbc_microcosm_control <- aov(mbc ~ treatment, data = microcosm_data_control)
tukey_mbc_microcosm_control_treatment <- TukeyHSD(one_way_anova_mbc_microcosm_control, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_control_treatment

one_way_anova_mbc_microcosm_moderate <- aov(mbc ~ treatment, data = microcosm_data_moderate)
tukey_mbc_microcosm_moderate_treatment <- TukeyHSD(one_way_anova_mbc_microcosm_moderate, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_moderate_treatment

one_way_anova_mbc_microcosm_high <- aov(mbc ~ treatment, data = microcosm_data_high)
tukey_mbc_microcosm_high_treatment <- TukeyHSD(one_way_anova_mbc_microcosm_high, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_high_treatment

# MBN
two_way_anova_mbn_microcosm <- aov(mbn ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_mbn_microcosm)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
microcosm_data_control <- microcosm_data %>% subset(incubation == "control")
microcosm_data_moderate <- microcosm_data %>% subset(incubation == "moderate")
microcosm_data_high <- microcosm_data %>% subset(incubation == "high")

one_way_anova_mbn_microcosm_control <- aov(mbn ~ treatment, data = microcosm_data_control)
tukey_mbn_microcosm_control_treatment <- TukeyHSD(one_way_anova_mbn_microcosm_control, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_control_treatment
  
one_way_anova_mbn_microcosm_moderate <- aov(mbn ~ treatment, data = microcosm_data_moderate)
tukey_mbn_microcosm_moderate_treatment <- TukeyHSD(one_way_anova_mbn_microcosm_moderate, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_moderate_treatment
  
one_way_anova_mbn_microcosm_high <- aov(mbn ~ treatment, data = microcosm_data_high)
tukey_mbn_microcosm_high_treatment <- TukeyHSD(one_way_anova_mbn_microcosm_high, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_high_treatment
  
# EOC
two_way_anova_eoc_microcosm <- aov(eoc ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_eoc_microcosm)

# Post-Hoc testing: Both treatment and incubation
tukey_eoc_microcosm_treatment <- TukeyHSD(two_way_anova_eoc_microcosm, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_eoc_microcosm_treatment

tukey_eoc_microcosm_incubation <- TukeyHSD(two_way_anova_eoc_microcosm, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_eoc_microcosm_incubation

# ETN
two_way_anova_etn_microcosm <- aov(etn ~ treatment * incubation, data = microcosm_data)
summary(two_way_anova_etn_microcosm)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
one_way_anova_etn_microcosm_control <- aov(etn ~ treatment, data = microcosm_data_control)
tukey_etn_microcosm_control_treatment <- TukeyHSD(one_way_anova_etn_microcosm_control, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_control_treatment

one_way_anova_etn_microcosm_moderate <- aov(etn ~ treatment, data = microcosm_data_moderate)
tukey_etn_microcosm_moderate_treatment <- TukeyHSD(one_way_anova_etn_microcosm_moderate, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_moderate_treatment

one_way_anova_etn_microcosm_high <- aov(etn ~ treatment, data = microcosm_data_high)
tukey_etn_microcosm_high_treatment <- TukeyHSD(one_way_anova_etn_microcosm_high, "treatment", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_high_treatment

#### 2-way ANOVA: Category, Reduction, Category*Reduction ######
# Alpha Diversity
two_way_anova_diversity_microcosm_category <- aov(shan ~ system * incubation, data = microcosm_data)
summary(two_way_anova_diversity_microcosm_category)

# Post-Hoc testing: Both category and incubation
tukey_diversity_microcosm_category <- TukeyHSD(two_way_anova_diversity_microcosm_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_diversity_microcosm_category

tukey_diversity_microcosm_incubation_category <- TukeyHSD(two_way_anova_diversity_microcosm_category, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_diversity_microcosm_incubation_category

# FLNF
two_way_anova_flnf_microcosm_category <- aov(n.fix ~ system * incubation, data = microcosm_data)
summary(two_way_anova_flnf_microcosm_category)

# Post-Hoc testing: Both category and incubation
tukey_flnf_microcosm_category <- TukeyHSD(two_way_anova_flnf_microcosm_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_flnf_microcosm_category

tukey_flnf_microcosm_incubation_category <- TukeyHSD(two_way_anova_flnf_microcosm_category, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_flnf_microcosm_incubation_category

# MBC
two_way_anova_mbc_microcosm_category <- aov(mbc ~ system * incubation, data = microcosm_data)
summary(two_way_anova_mbc_microcosm_category)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
one_way_anova_mbc_microcosm_control_category <- aov(mbc ~ system, data = microcosm_data_control)
tukey_mbc_microcosm_control_category <- TukeyHSD(one_way_anova_mbc_microcosm_control_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_control_category

one_way_anova_mbc_microcosm_moderate_category <- aov(mbc ~ system, data = microcosm_data_moderate)
tukey_mbc_microcosm_moderate_category <- TukeyHSD(one_way_anova_mbc_microcosm_moderate_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_moderate_category

one_way_anova_mbc_microcosm_high_category <- aov(mbc ~ system, data = microcosm_data_high)
tukey_mbc_microcosm_high_category <- TukeyHSD(one_way_anova_mbc_microcosm_high_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbc_microcosm_high_category




# MBN
two_way_anova_mbn_microcosm_category <- aov(mbn ~ system * incubation, data = microcosm_data)
summary(two_way_anova_mbn_microcosm_category)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
one_way_anova_mbn_microcosm_control_category <- aov(mbc ~ incubation, data = microcosm_data)
tukey_mbn_microcosm_control_category <- TukeyHSD(one_way_anova_mbn_microcosm_control_category, "incubation", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_control_category

one_way_anova_mbn_microcosm_moderate_category <- aov(mbn ~ system, data = microcosm_data_moderate)
tukey_mbn_microcosm_moderate_category <- TukeyHSD(one_way_anova_mbn_microcosm_moderate_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_moderate_category

one_way_anova_mbn_microcosm_high_category <- aov(mbn ~ system, data = microcosm_data_high)
tukey_mbn_microcosm_high_category <- TukeyHSD(one_way_anova_mbn_microcosm_high_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_mbn_microcosm_high_category

# EOC
two_way_anova_eoc_microcosm_category <- aov(eoc ~ system * incubation, data = microcosm_data)
summary(two_way_anova_eoc_microcosm_category)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
one_way_anova_eoc_microcosm_control_category <- aov(eoc ~ system, data = microcosm_data_control)
tukey_eoc_microcosm_control_category <- TukeyHSD(one_way_anova_eoc_microcosm_control_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_eoc_microcosm_control_category

one_way_anova_eoc_microcosm_moderate_category <- aov(eoc ~ system, data = microcosm_data_moderate)
tukey_eoc_microcosm_moderate_category <- TukeyHSD(one_way_anova_eoc_microcosm_moderate_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_eoc_microcosm_moderate_category

one_way_anova_eoc_microcosm_high_category <- aov(eoc ~ system, data = microcosm_data_high)
tukey_eoc_microcosm_high_category <- TukeyHSD(one_way_anova_eoc_microcosm_high_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_eoc_microcosm_high_category




# ETN
two_way_anova_etn_microcosm_category <- aov(etn ~ system * incubation, data = microcosm_data)
summary(two_way_anova_etn_microcosm_category)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
one_way_anova_etn_microcosm_control_category <- aov(etn ~ system, data = microcosm_data_control)
tukey_etn_microcosm_control_category <- TukeyHSD(one_way_anova_etn_microcosm_control_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_control_category

one_way_anova_etn_microcosm_moderate_category <- aov(etn ~ system, data = microcosm_data_moderate)
tukey_etn_microcosm_moderate_category <- TukeyHSD(one_way_anova_etn_microcosm_moderate_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_moderate_category

one_way_anova_etn_microcosm_high_category <- aov(etn ~ system, data = microcosm_data_high)
tukey_etn_microcosm_high_category <- TukeyHSD(one_way_anova_etn_microcosm_high_category, "system", ordered = FALSE, conf.level = 0.95)
tukey_etn_microcosm_high_category

### permANOVA: Management, Biodiv reduction, Management*Biodiv reduction ---------------------------------------------------------
set.seed(37920)
# Convert to relative abundance
bact_microcosm_rel <- transform_sample_counts(bact_microcosm, function(x) x / sum(x) )

# Change the name of the incubation factors 
sample_data(bact_microcosm_rel)$incubation <- as.factor(sample_data(bact_microcosm_rel)$incubation)
levels(sample_data(bact_microcosm_rel)$incubation) <- c("control", "moderate", "high")

# Perform NMDS ordination with Bray curtis distance using the ordinate() function
bray_NMDS_microcosm_ord <- ordinate(bact_microcosm_rel, method = "NMDS", distance = "bray") # Stress = 0.194

# Calculate Bray-curtis distance using phyloseq::distance to identify significance of our factors: treatment and season/sample date
set.seed(37920)
bray_NMDS_microcosm <- phyloseq::distance(bact_microcosm_rel,"bray")

# What are the interactive effects of management and biodiversity reduction on diazotroph Beta-diversity? -- perform permanova
set.seed(37920)
permanova <- adonis(bray_NMDS_microcosm~sample_data(bact_microcosm_rel)$treatment *
                      sample_data(bact_microcosm_rel)$incubation, permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication

# Pairwise adonis2 between each biodiversity reduction level 
# Control
bact_microcosm_rel_control <- subset_samples(bact_microcosm_rel, incubation == "control")
set.seed(37920)
bray_NMDS_microcosm_control <- phyloseq::distance(bact_microcosm_rel_control,"bray")

pairwise.adonis(bray_NMDS_microcosm_control, sample_data(bact_microcosm_rel_control)$treatment)

# Moderate
bact_microcosm_rel_moderate <- subset_samples(bact_microcosm_rel, incubation == "moderate")
set.seed(37920)
bray_NMDS_microcosm_moderate <- phyloseq::distance(bact_microcosm_rel_moderate,"bray")

pairwise.adonis(bray_NMDS_microcosm_moderate, sample_data(bact_microcosm_rel_moderate)$treatment)

# High
bact_microcosm_rel_high <- subset_samples(bact_microcosm_rel, incubation == "high")
set.seed(37920)
bray_NMDS_microcosm_high <- phyloseq::distance(bact_microcosm_rel_high,"bray")

pairwise.adonis(bray_NMDS_microcosm_high, sample_data(bact_microcosm_rel_high)$treatment)
# Create envfit model to visualize effects of MBC/eoc/etn as vectors on the NMDS ordination
# Exctract sample metadata from the rarefied phyloseq object 
env <- as.data.frame(sample_data(bact_microcosm))
set.seed(37920)
# Perform the envfit() function to quantify effects of continuous variables on community composition
ord.fit.sig <- envfit(bray_NMDS_microcosm_ord ~ 
                        env$n.fix + env$mbc + env$mbn + 
                        env$eoc + env$etn, 
                        perm = 9999, na.rm = TRUE)
### permANOVA: Category, Biodiv reduction, Category*Biodiv reduction ---------------------------------------------------------
# What are the interactive effects of land-use category and biodiversity reduction on diazotroph Beta-diversity? -- perform permanova
set.seed(37920)
permanova <- adonis(bray_NMDS_microcosm~sample_data(bact_microcosm_rel)$system *
                      sample_data(bact_microcosm_rel)$incubation, permutations = 9999, method = "bray")
permanova["aov.tab"] # Obtain for publication

# Pairwise adonis2 between each biodiversity reduction level 
# Control
pairwise.adonis(bray_NMDS_microcosm_control, sample_data(bact_microcosm_rel_control)$system)

# Moderate
pairwise.adonis(bray_NMDS_microcosm_moderate, sample_data(bact_microcosm_rel_moderate)$system)

# High
pairwise.adonis(bray_NMDS_microcosm_high, sample_data(bact_microcosm_rel_high)$system)

# Pairwise adonis2 for biodiversity reduction 
pairwise.adonis(bray_NMDS_microcosm, sample_data(bact_microcosm_rel)$incubation)
### Partial Correlation Analysis: FLNF vs. Diversity  -------------
# Perform partial correlation between diazotroph alpha-diversity~n.fixation when controlling for (eoc.etn.mbc)
microcosm_data_nona <- microcosm_data %>% drop_na(eoc)

pcor_fixation_diversity_microcosm <- pcor.test(microcosm_data_nona$shan, microcosm_data_nona$n.fix, microcosm_data_nona[,c("etn", "mbc", "mbn")], method = "spearman")

pcor_fixation_diversity_microcosm

### Partial Correlation Analysis: MBC normalized FLNF vs. Diversity ---------------------------------------
microcosm_data_nona <- microcosm_data %>% drop_na(eoc, etn, mbc, mbn, n.fix.mbc)

pcor_normfixation_diversity_microcosm <- pcor.test(microcosm_data_nona$shan, microcosm_data_nona$n.fix.mbc, microcosm_data_nona[,c("eoc", "etn", "mbn")], method = "spearman")

pcor_normfixation_diversity_microcosm

### 2-way ANOVA: Diversity and land-use category on nitrogen fixation: Comparison between field and lab ---------------------------
FLNF_div_mcse_lm <- aov(n.fix ~ shan*system, data = mcse_data)
summary(FLNF_div_mcse_lm)

FLNF_div_microcosm_lm <- aov(n.fix ~ shan*system, data = microcosm_data)
summary(FLNF_div_microcosm_lm)

# Post-Hoc testing: Both treatment by incubation (due to significant interaction)
microcosm_data_annual <- microcosm_data %>% subset(system == "annual")
microcosm_data_perennial <- microcosm_data %>% subset(system == "perennial")
microcosm_data_forest <- microcosm_data %>% subset(system == "forest")

one_way_anova_FLNF_div_microcosm_annual <- aov(n.fix ~ shan, data = microcosm_data_annual)
summary(one_way_anova_FLNF_div_microcosm_annual)

one_way_anova_FLNF_div_microcosm_perennial <- aov(n.fix ~ shan, data = microcosm_data_perennial)
summary(one_way_anova_FLNF_div_microcosm_perennial)

one_way_anova_FLNF_div_microcosm_forest <- aov(n.fix ~ shan, data = microcosm_data_forest)
summary(one_way_anova_FLNF_div_microcosm_forest)

### 2-way ANOVA: Diversity and land-use category on MBC normalized nitrogen fixation: Comparison between field and lab ---------------------------
FLNF_div_mcse_lm <- aov(n.fix.mbc ~ shan*system, data = mcse_data)
summary(FLNF_div_mcse_lm)

FLNF_div_microcosm_lm <- aov(n.fix.mbc ~ shan*system, data = microcosm_data)
summary(FLNF_div_microcosm_lm)

### DESEQ -------------------------------------
##### Control vs. Moderate --------------------
bact_microcosm_unrarefied <- subset_samples(bact.p_filtered, incubation == "0" | incubation == "24")

# Change the name of the incubation factors 
sample_data(bact_microcosm_unrarefied)$incubation <- as.factor(sample_data(bact_microcosm_unrarefied)$incubation)
levels(sample_data(bact_microcosm_unrarefied)$incubation) <- c("control", "moderate")


diagdds_control_v_moderate <- phyloseq_to_deseq2(bact_microcosm_unrarefied, ~ incubation) 

# calculate geometric means prior to estimate size factors
geoMeans <-  apply(counts(diagdds_control_v_moderate), 1, gm_mean)
diagdds_control_v_moderate <- estimateSizeFactors(diagdds_control_v_moderate, geoMeans = geoMeans)
diagdds_control_v_moderate <- DESeq(diagdds_control_v_moderate, fitType="parametric")

results_control_v_moderate <- results(diagdds_control_v_moderate, alpha = 0.01, contrast = c("incubation", "control","moderate"))

sigtab_results_control_v_moderate <- results_control_v_moderate[which(results_control_v_moderate$padj < 0.01), ]

sigtab_results_TAX_control_v_moderate <- cbind(as(sigtab_results_control_v_moderate, "data.frame"), as(tax_table(bact_microcosm_unrarefied)[row.names(sigtab_results_control_v_moderate), ], "matrix"))
summary(sigtab_results_control_v_moderate)


sigtab_results_TAX_genra_control_v_moderate <- subset(sigtab_results_TAX_control_v_moderate, !is.na(Genus))

# Order Genus 
x <- tapply(sigtab_results_TAX_genra_control_v_moderate$log2FoldChange, sigtab_results_TAX_genra_control_v_moderate$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_control_v_moderate$Class <- factor(as.character(sigtab_results_TAX_genra_control_v_moderate$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra_control_v_moderate$log2FoldChange, sigtab_results_TAX_genra_control_v_moderate$Genus, function(x) max(x))
x <- sort(x, TRUE)

sigtab_results_TAX_genra_control_v_moderate$Genus <- factor(as.character(sigtab_results_TAX_genra_control_v_moderate$Genus), levels=names(x))
sigtab_results_TAX_genra_control_v_moderate$incubation <- "Control vs. Moderate"
sigtab_results_TAX_genra_control_v_moderate$Genus <- gsub("Candidatus_Azobacteroides", "Candidatus", sigtab_results_TAX_genra_control_v_moderate$Genus)



##### Control vs. High ------------------------
bact_microcosm_unrarefied <- subset_samples(bact.p_filtered, incubation == "0" | incubation == "72")

# Change the name of the incubation factors 
sample_data(bact_microcosm_unrarefied)$incubation <- as.factor(sample_data(bact_microcosm_unrarefied)$incubation)
levels(sample_data(bact_microcosm_unrarefied)$incubation) <- c("control", "high")


diagdds_control_v_high <- phyloseq_to_deseq2(bact_microcosm_unrarefied, ~ incubation) 

# calculate geometric means prior to estimate size factors
geoMeans <-  apply(counts(diagdds_control_v_high), 1, gm_mean)
diagdds_control_v_high <- estimateSizeFactors(diagdds_control_v_high, geoMeans = geoMeans)
diagdds_control_v_high <- DESeq(diagdds_control_v_high, fitType="parametric")

results_control_v_high <- results(diagdds_control_v_high, alpha = 0.01, contrast = c("incubation", "control","high"))

sigtab_results_control_v_high <- results_control_v_high[which(results_control_v_high$padj < 0.01), ]

sigtab_results_TAX_control_v_high <- cbind(as(sigtab_results_control_v_high, "data.frame"), as(tax_table(bact_microcosm_unrarefied)[row.names(sigtab_results_control_v_high), ], "matrix"))
summary(sigtab_results_control_v_high)


sigtab_results_TAX_genra_control_v_high <- subset(sigtab_results_TAX_control_v_high, !is.na(Genus))

# Order Genus 
x <- tapply(sigtab_results_TAX_genra_control_v_high$log2FoldChange, sigtab_results_TAX_genra_control_v_high$Class, function(x) max(x))
x <- sort(x, TRUE)
sigtab_results_TAX_genra_control_v_high$Class <- factor(as.character(sigtab_results_TAX_genra_control_v_high$Class), levels=names(x))
# Genus order
x <- tapply(sigtab_results_TAX_genra_control_v_high$log2FoldChange, sigtab_results_TAX_genra_control_v_high$Genus, function(x) max(x))
x <- sort(x, TRUE)

sigtab_results_TAX_genra_control_v_high$Genus <- factor(as.character(sigtab_results_TAX_genra_control_v_high$Genus), levels=names(x))
sigtab_results_TAX_genra_control_v_high$incubation <- "Control vs. High"


sigtab_results_TAX_genra_deseq_combined <- rbind(sigtab_results_TAX_genra_control_v_moderate,                   sigtab_results_TAX_genra_control_v_high)

sigtab_results_TAX_genra_deseq_combined$Genus <- gsub("Candidatus_Azobacteroides", "Candidatus", sigtab_results_TAX_genra_deseq_combined$Genus)

sigtab_results_TAX_genra_control_v_high$Genus <- gsub("Candidatus_Azobacteroides", "Candidatus", sigtab_results_TAX_genra_control_v_high$Genus)



##### Moderate vs. High -----------------------
bact_microcosm_unrarefied <- subset_samples(bact.p_filtered, incubation == "24" | incubation == "72")

# Change the name of the incubation factors 
sample_data(bact_microcosm_unrarefied)$incubation <- as.factor(sample_data(bact_microcosm_unrarefied)$incubation)
levels(sample_data(bact_microcosm_unrarefied)$incubation) <- c("moderate", "high")


diagdds_moderate_v_high <- phyloseq_to_deseq2(bact_microcosm_unrarefied, ~ incubation) 

# calculate geometric means prior to estimate size factors
geoMeans <-  apply(counts(diagdds_moderate_v_high), 1, gm_mean)
diagdds_moderate_v_high <- estimateSizeFactors(diagdds_moderate_v_high, geoMeans = geoMeans)
diagdds_moderate_v_high <- DESeq(diagdds_moderate_v_high, fitType="parametric")

results_moderate_v_high <- results(diagdds_moderate_v_high, alpha = 0.01, contrast = c("incubation", "moderate","high"))

sigtab_results_moderate_v_high <- results_moderate_v_high[which(results_moderate_v_high$padj < 0.01), ]
# Biodiversity Reduction Experiment Figures #######
#### Diversity #####
# Obtain mean alpha diversity values for each system (annual, perennial, forest) in the summer sampling date 
lm_diversity_category_microcosm_control <- lm(shan ~ system + 0, data = microcosm_data_control)
summary(lm_diversity_category_microcosm_control)

lm_diversity_category_microcosm_moderate <- lm(shan ~ system + 0, data = microcosm_data_moderate)
summary(lm_diversity_category_microcosm_moderate)

lm_diversity_category_microcosm_high <- lm(shan ~ system + 0, data = microcosm_data_high)
summary(lm_diversity_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(368.49, 248.77, 135.65)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_diversity_effects_microcosm_annual <- data.frame(est, incubation, system)

est <- c(258.34, 191.16, 81.15)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_diversity_effects_microcosm_perennial <- data.frame(est, incubation, system)

est <- c(236.29, 216.19, 71.82)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_diversity_effects_microcosm_forest <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_diversity_effects_microcosm_annual,                         category_diversity_effects_microcosm_perennial,                      category_diversity_effects_microcosm_forest)

merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: Annual groups
letters_df_diversity_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_diversity_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm_control$treatment <- rownames(letters_df_diversity_microcosm_control) 

letters_df_diversity_microcosm_control$placement <- c(515, 515, 515, 515, 515, 515, 515, 515)
value_map_system <- c("T1" = "annual",
                     "T3" = "annual",
                     "T4" = "annual",
                     "T5" = "perennial",
                     "T7" = "perennial",
                     "DF" = "forest",
                     "SF" = "forest",
                     "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_diversity_microcosm_control$system <- value_map_system[letters_df_diversity_microcosm_control$treatment]
letters_df_diversity_microcosm_control$system <- factor(letters_df_diversity_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_diversity_microcosm_control$incubation <- value_map_incubation[letters_df_diversity_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_diversity_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_diversity_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm_moderate$treatment <- rownames(letters_df_diversity_microcosm_moderate) 

letters_df_diversity_microcosm_moderate$placement <- c(515, 515, 515, 515, 515, 515, 515, 515)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_diversity_microcosm_moderate$system <- value_map_system[letters_df_diversity_microcosm_moderate$treatment]
letters_df_diversity_microcosm_moderate$system <- factor(letters_df_diversity_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_diversity_microcosm_moderate$incubation <- value_map_incubation[letters_df_diversity_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_diversity_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_diversity_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm_high$treatment <- rownames(letters_df_diversity_microcosm_high) 

letters_df_diversity_microcosm_high$placement <- c(515, 515, 515, 515, 515, 515, 515, 515)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_diversity_microcosm_high$system <- value_map_system[letters_df_diversity_microcosm_high$treatment]
letters_df_diversity_microcosm_high$system <- factor(letters_df_diversity_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_diversity_microcosm_high$incubation <- value_map_incubation[letters_df_diversity_microcosm_high$treatment]

# Merge all letters together
merge_tukey_diversity_microcosm <- list(letters_df_diversity_microcosm_control,
                                        letters_df_diversity_microcosm_moderate,
                                        letters_df_diversity_microcosm_high)
merge_tukey_diversity_microcosm <- merge_tukey_diversity_microcosm %>% purrr::reduce(full_join)

# Reorder Treatments
microcosm_data$treatment <- factor(microcosm_data$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "DF", "CF"))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_diversity_microcosm$incubation <- factor(merge_tukey_diversity_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_diversity_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = shan, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_diversity_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab("Hill Number (q=1)") +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_diversity_plot

## Diversity by Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(shan ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(450, 450, 450, 450, 450, 450, 450, 450, 450)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                           "perennial_control",
                                                                                           "forest_control",
                                                                                           "annual_moderate",
                                                                                           "perennial_moderate",
                                                                                           "forest_moderate",
                                                                                           "annual_high",
                                                                                           "perennial_high",
                                                                                           "forest_high"
                                                                                          ))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                           "perennial_control",
                                                                                           "forest_control",
                                                                                           "annual_moderate",
                                                                                           "perennial_moderate",
                                                                                           "forest_moderate",
                                                                                           "annual_high",
                                                                                           "perennial_high",
                                                                                           "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_div <- microcosm_data %>%
      ggplot(aes(x=system_incubation, y = shan)) +
      geom_boxplot(aes(fill = system), outlier.shape = NA) +
      geom_text(data = letters_df_diversity_microcosm, 
                aes(x = system_incubation, y=placement, label=Letter), 
                size =5.5, color="black", fontface="bold") +
      #geom_hline(aes(yintercept=est), linetype = 'dashed', 
                 #col = 'black', linewidth = 1.25,
                 #data=merge_effects_category_microcosm) +
      scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                        labels = c("Annual", "Perennial", "Forest")) +
      ylab("Hill Number (q=1)") +
      xlab("Fumigation Exposure Treatment") +
      scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
      geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
      theme_bw() +
      theme(axis.title.x = element_text(size = 14), 
            title = element_text(size = 14),
            axis.title.y = element_text(size = 14), 
            strip.text.x = element_text(size = 14),
            strip.text.y = element_text(size = 14),
            axis.text.x = element_text(size = 12, color = "black"),
            axis.text.y = element_text(size = 12, color = "black"), 
            legend.text = element_text(size = 12), 
            legend.title = element_text(size = 14), 
            plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_div



#### FLNF #####
# Obtain mean alpha diversity values for each system (annual, perennial, forest) in the summer sampling date 
lm_flnf_category_microcosm_control <- lm(n.fix ~ system + 0, data = microcosm_data_control)
summary(lm_flnf_category_microcosm_control)

lm_flnf_category_microcosm_moderate <- lm(n.fix ~ system + 0, data = microcosm_data_moderate)
summary(lm_flnf_category_microcosm_moderate)

lm_flnf_category_microcosm_high <- lm(n.fix ~ system + 0, data = microcosm_data_high)
summary(lm_flnf_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(8.1296, 4.7024, 5.9548)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_flnf_effects_microcosm_control <- data.frame(est, incubation, system)

est <- c(2.6519, 0.4735, 1.1257)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_flnf_effects_microcosm_moderate <- data.frame(est, incubation, system)

est <- c(1.0689, 0.4720, 0.9239)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_flnf_effects_microcosm_high <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_flnf_effects_microcosm_control,
                                         category_flnf_effects_microcosm_moderate,
                                         category_flnf_effects_microcosm_high)
merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: control hour
letters_df_flnf_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_flnf_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_flnf_microcosm_control$treatment <- rownames(letters_df_flnf_microcosm_control) 

letters_df_flnf_microcosm_control$placement <- c(15, 15, 15, 15, 15, 15, 15, 15)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_flnf_microcosm_control$system <- value_map_system[letters_df_flnf_microcosm_control$treatment]
letters_df_flnf_microcosm_control$system <- factor(letters_df_flnf_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_flnf_microcosm_control$incubation <- value_map_incubation[letters_df_flnf_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_flnf_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_flnf_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_flnf_microcosm_moderate$treatment <- rownames(letters_df_flnf_microcosm_moderate) 

letters_df_flnf_microcosm_moderate$placement <- c(10, 10, 10, 10, 10, 10, 10, 10)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_flnf_microcosm_moderate$system <- value_map_system[letters_df_flnf_microcosm_moderate$treatment]
letters_df_flnf_microcosm_moderate$system <- factor(letters_df_flnf_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_flnf_microcosm_moderate$incubation <- value_map_incubation[letters_df_flnf_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_flnf_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_flnf_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_flnf_microcosm_high$treatment <- rownames(letters_df_flnf_microcosm_high) 

letters_df_flnf_microcosm_high$placement <- c(10, 10, 10, 10, 10, 10, 10, 10)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_flnf_microcosm_high$system <- value_map_system[letters_df_flnf_microcosm_high$treatment]
letters_df_flnf_microcosm_high$system <- factor(letters_df_flnf_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_flnf_microcosm_high$incubation <- value_map_incubation[letters_df_flnf_microcosm_high$treatment]

# Merge all letters together
merge_tukey_flnf_microcosm <- list(letters_df_flnf_microcosm_control,
                                        letters_df_flnf_microcosm_moderate,
                                        letters_df_flnf_microcosm_high)
merge_tukey_flnf_microcosm <- merge_tukey_flnf_microcosm %>% purrr:::reduce(full_join)

# Reorder Treatments
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_flnf_microcosm$incubation <- factor(merge_tukey_flnf_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_flnf_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = n.fix, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_flnf_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_flnf_plot


## FLFNF By Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(n.fix ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(15, 15, 15, 15, 15, 15, 15, 15, 15)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                        "perennial_control",
                                                                                        "forest_control",
                                                                                        "annual_moderate",
                                                                                        "perennial_moderate",
                                                                                        "forest_moderate",
                                                                                        "annual_high",
                                                                                        "perennial_high",
                                                                                        "forest_high"
))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                                                        "perennial_control",
                                                                                                                        "forest_control",
                                                                                                                        "annual_moderate",
                                                                                                                        "perennial_moderate",
                                                                                                                        "forest_moderate",
                                                                                                                        "annual_high",
                                                                                                                        "perennial_high",
                                                                                                                        "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_flnf <- microcosm_data %>%
  ggplot(aes(x=system_incubation, y = n.fix)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  geom_text(data = letters_df_diversity_microcosm, 
            aes(x = system_incubation, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  #geom_hline(aes(yintercept=est), linetype = 'dashed', 
  #col = 'black', linewidth = 1.25,
  #data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("Fumigation Exposure Treatment") +
  scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_flnf



#### MBC #####
# Obtain mean mbc values for each system (annual, perennial, forest) in the summer sampling date 
lm_mbc_category_microcosm_control <- lm(mbc ~ system + 0, data = microcosm_data_control)
summary(lm_mbc_category_microcosm_control)

lm_mbc_category_microcosm_moderate <- lm(mbc ~ system + 0, data = microcosm_data_moderate)
summary(lm_mbc_category_microcosm_moderate)

lm_mbc_category_microcosm_high <- lm(mbc ~ system + 0, data = microcosm_data_high)
summary(lm_mbc_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(240.00, 471.57, 262.63)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_mbc_effects_microcosm_control <- data.frame(est, incubation, system)

est <- c(265.10, 309.11, 138.47)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_mbc_effects_microcosm_moderate <- data.frame(est, incubation, system)

est <- c(252.70, 244.87, 387.92)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_mbc_effects_microcosm_high <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_mbc_effects_microcosm_control,
                                         category_mbc_effects_microcosm_moderate,
                                         category_mbc_effects_microcosm_high)
merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: control hour
letters_df_mbc_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_mbc_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_mbc_microcosm_control$treatment <- rownames(letters_df_mbc_microcosm_control) 

letters_df_mbc_microcosm_control$placement <- c(650, 650, 650, 650, 650, 650, 650, 650)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_mbc_microcosm_control$system <- value_map_system[letters_df_mbc_microcosm_control$treatment]
letters_df_mbc_microcosm_control$system <- factor(letters_df_mbc_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_mbc_microcosm_control$incubation <- value_map_incubation[letters_df_mbc_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_mbc_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_mbc_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_mbc_microcosm_moderate$treatment <- rownames(letters_df_mbc_microcosm_moderate) 

letters_df_mbc_microcosm_moderate$placement <- c(650, 650, 650, 650, 650, 650, 650, 650)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_mbc_microcosm_moderate$system <- value_map_system[letters_df_mbc_microcosm_moderate$treatment]
letters_df_mbc_microcosm_moderate$system <- factor(letters_df_mbc_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_mbc_microcosm_moderate$incubation <- value_map_incubation[letters_df_mbc_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_mbc_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_mbc_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_mbc_microcosm_high$treatment <- rownames(letters_df_mbc_microcosm_high) 

letters_df_mbc_microcosm_high$placement <- c(650, 650, 650, 650, 650, 650, 650, 650)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_mbc_microcosm_high$system <- value_map_system[letters_df_mbc_microcosm_high$treatment]
letters_df_mbc_microcosm_high$system <- factor(letters_df_mbc_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_mbc_microcosm_high$incubation <- value_map_incubation[letters_df_mbc_microcosm_high$treatment]

# Merge all letters together
merge_tukey_mbc_microcosm <- list(letters_df_mbc_microcosm_control,
                                   letters_df_mbc_microcosm_moderate,
                                   letters_df_mbc_microcosm_high)
merge_tukey_mbc_microcosm <- merge_tukey_mbc_microcosm %>% purrr::reduce(full_join)

# Reorder Treatments
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_mbc_microcosm$incubation <- factor(merge_tukey_mbc_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_mbc_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = mbc, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_mbc_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(MBC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_mbc_plot

## MBC By Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(mbc ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(650, 650, 650, 650, 650, 650, 650, 650, 650)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                        "perennial_control",
                                                                                        "forest_control",
                                                                                        "annual_moderate",
                                                                                        "perennial_moderate",
                                                                                        "forest_moderate",
                                                                                        "annual_high",
                                                                                        "perennial_high",
                                                                                        "forest_high"
))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                                                        "perennial_control",
                                                                                                                        "forest_control",
                                                                                                                        "annual_moderate",
                                                                                                                        "perennial_moderate",
                                                                                                                        "forest_moderate",
                                                                                                                        "annual_high",
                                                                                                                        "perennial_high",
                                                                                                                        "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_mbc <- microcosm_data %>%
  ggplot(aes(x=system_incubation, y = mbc)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  geom_text(data = letters_df_diversity_microcosm, 
            aes(x = system_incubation, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  #geom_hline(aes(yintercept=est), linetype = 'dashed', 
  #col = 'black', linewidth = 1.25,
  #data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(MBC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("Fumigation Exposure Treatment") +
  scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 11), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_mbc




#### MBN #####
# Obtain mean mbn values for each system (annual, perennial, forest) in the summer sampling date 
lm_mbn_category_microcosm_control <- lm(mbn ~ system + 0, data = microcosm_data_control)
summary(lm_mbn_category_microcosm_control)

lm_mbn_category_microcosm_moderate <- lm(mbn ~ system + 0, data = microcosm_data_moderate)
summary(lm_mbn_category_microcosm_moderate)

lm_mbn_category_microcosm_high <- lm(mbn ~ system + 0, data = microcosm_data_high)
summary(lm_mbn_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(32.100, 81.211, 37.941)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_mbn_effects_microcosm_control <- data.frame(est, incubation, system)

est <- c(4.291, 8.850, 6.752)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_mbn_effects_microcosm_moderate <- data.frame(est, incubation, system)

est <- c(25.915, 29.265, 25.829)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_mbn_effects_microcosm_high <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_mbn_effects_microcosm_control,
                                         category_mbn_effects_microcosm_moderate,
                                         category_mbn_effects_microcosm_high)
merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: control hour
letters_df_mbn_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_mbn_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_mbn_microcosm_control$treatment <- rownames(letters_df_mbn_microcosm_control) 

letters_df_mbn_microcosm_control$placement <- c(100, 100, 100, 100, 100, 100, 100, 100)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_mbn_microcosm_control$system <- value_map_system[letters_df_mbn_microcosm_control$treatment]
letters_df_mbn_microcosm_control$system <- factor(letters_df_mbn_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_mbn_microcosm_control$incubation <- value_map_incubation[letters_df_mbn_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_mbn_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_mbn_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_mbn_microcosm_moderate$treatment <- rownames(letters_df_mbn_microcosm_moderate) 

letters_df_mbn_microcosm_moderate$placement <- c(100, 100, 100, 100, 100, 100, 100, 100)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_mbn_microcosm_moderate$system <- value_map_system[letters_df_mbn_microcosm_moderate$treatment]
letters_df_mbn_microcosm_moderate$system <- factor(letters_df_mbn_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_mbn_microcosm_moderate$incubation <- value_map_incubation[letters_df_mbn_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_mbn_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_mbn_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_mbn_microcosm_high$treatment <- rownames(letters_df_mbn_microcosm_high) 

letters_df_mbn_microcosm_high$placement <- c(100, 100, 100, 100, 100, 100, 100, 100)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_mbn_microcosm_high$system <- value_map_system[letters_df_mbn_microcosm_high$treatment]
letters_df_mbn_microcosm_high$system <- factor(letters_df_mbn_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_mbn_microcosm_high$incubation <- value_map_incubation[letters_df_mbn_microcosm_high$treatment]

# Merge all letters together
merge_tukey_mbn_microcosm <- list(letters_df_mbn_microcosm_control,
                                  letters_df_mbn_microcosm_moderate,
                                  letters_df_mbn_microcosm_high)
merge_tukey_mbn_microcosm <- merge_tukey_mbn_microcosm %>% purrr::reduce(full_join)

# Reorder Treatments
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_mbn_microcosm$incubation <- factor(merge_tukey_mbn_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_mbn_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = mbn, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_mbn_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(MBN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_mbn_plot

## MBN By Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(mbn ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(90, 90, 90, 90, 90, 90, 90, 90, 90)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                        "perennial_control",
                                                                                        "forest_control",
                                                                                        "annual_moderate",
                                                                                        "perennial_moderate",
                                                                                        "forest_moderate",
                                                                                        "annual_high",
                                                                                        "perennial_high",
                                                                                        "forest_high"
))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                                                        "perennial_control",
                                                                                                                        "forest_control",
                                                                                                                        "annual_moderate",
                                                                                                                        "perennial_moderate",
                                                                                                                        "forest_moderate",
                                                                                                                        "annual_high",
                                                                                                                        "perennial_high",
                                                                                                                        "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_mbn <- microcosm_data %>%
  ggplot(aes(x=system_incubation, y = mbn)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  geom_text(data = letters_df_diversity_microcosm, 
            aes(x = system_incubation, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  #geom_hline(aes(yintercept=est), linetype = 'dashed', 
  #col = 'black', linewidth = 1.25,
  #data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(MBN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("Fumigation Exposure Treatment") +
  scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 11), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_mbn





#### EOC #####
# Obtain mean eoc values for each system (annual, perennial, forest) in the summer sampling date 
lm_eoc_category_microcosm_control <- lm(eoc ~ system + 0, data = microcosm_data_control)
summary(lm_eoc_category_microcosm_control)

lm_eoc_category_microcosm_moderate <- lm(eoc ~ system + 0, data = microcosm_data_moderate)
summary(lm_eoc_category_microcosm_moderate)

lm_eoc_category_microcosm_high <- lm(eoc ~ system + 0, data = microcosm_data_high)
summary(lm_eoc_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(94.517, 99.60, 116.06)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_eoc_effects_microcosm_control <- data.frame(est, incubation, system)

est <- c(155.24, 271.82, 275.37)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_eoc_effects_microcosm_moderate <- data.frame(est, incubation, system)

est <- c(158.63, 247.91, 223.87)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_eoc_effects_microcosm_high <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_eoc_effects_microcosm_control,
                                         category_eoc_effects_microcosm_moderate,
                                         category_eoc_effects_microcosm_high)
merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: control hour
letters_df_eoc_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_eoc_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_eoc_microcosm_control$treatment <- rownames(letters_df_eoc_microcosm_control) 

letters_df_eoc_microcosm_control$placement <- c(400, 400, 400, 400, 400, 400, 400, 400)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_eoc_microcosm_control$system <- value_map_system[letters_df_eoc_microcosm_control$treatment]
letters_df_eoc_microcosm_control$system <- factor(letters_df_eoc_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_eoc_microcosm_control$incubation <- value_map_incubation[letters_df_eoc_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_eoc_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_eoc_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_eoc_microcosm_moderate$treatment <- rownames(letters_df_eoc_microcosm_moderate) 

letters_df_eoc_microcosm_moderate$placement <- c(650, 650, 650, 650, 650, 650, 650, 650)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_eoc_microcosm_moderate$system <- value_map_system[letters_df_eoc_microcosm_moderate$treatment]
letters_df_eoc_microcosm_moderate$system <- factor(letters_df_eoc_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_eoc_microcosm_moderate$incubation <- value_map_incubation[letters_df_eoc_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_eoc_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_eoc_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_eoc_microcosm_high$treatment <- rownames(letters_df_eoc_microcosm_high) 

letters_df_eoc_microcosm_high$placement <- c(650, 650, 650, 650, 650, 650, 650, 650)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_eoc_microcosm_high$system <- value_map_system[letters_df_eoc_microcosm_high$treatment]
letters_df_eoc_microcosm_high$system <- factor(letters_df_eoc_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_eoc_microcosm_high$incubation <- value_map_incubation[letters_df_eoc_microcosm_high$treatment]

# Merge all letters together
merge_tukey_eoc_microcosm <- list(letters_df_eoc_microcosm_control,
                                  letters_df_eoc_microcosm_moderate,
                                  letters_df_eoc_microcosm_high)
merge_tukey_eoc_microcosm <- merge_tukey_eoc_microcosm %>% purrr::reduce(full_join)

# Reorder Treatments
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_eoc_microcosm$incubation <- factor(merge_tukey_eoc_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_eoc_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = eoc, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_eoc_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(EOC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_eoc_plot


## EOC by Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(eoc ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(400, 400, 400, 400, 400, 400, 400, 400, 400)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                        "perennial_control",
                                                                                        "forest_control",
                                                                                        "annual_moderate",
                                                                                        "perennial_moderate",
                                                                                        "forest_moderate",
                                                                                        "annual_high",
                                                                                        "perennial_high",
                                                                                        "forest_high"
))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                                                        "perennial_control",
                                                                                                                        "forest_control",
                                                                                                                        "annual_moderate",
                                                                                                                        "perennial_moderate",
                                                                                                                        "forest_moderate",
                                                                                                                        "annual_high",
                                                                                                                        "perennial_high",
                                                                                                                        "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_eoc <- microcosm_data %>%
  ggplot(aes(x=system_incubation, y = eoc)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  geom_text(data = letters_df_diversity_microcosm, 
            aes(x = system_incubation, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  #geom_hline(aes(yintercept=est), linetype = 'dashed', 
  #col = 'black', linewidth = 1.25,
  #data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(EOC~('\u00b5g'~C~g^-1~dry~soil)) +
  xlab("Fumigation Exposure Treatment") +
  scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 11), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_eoc




#### ETN #####
# Obtain mean etn values for each system (annual, perennial, forest) in the summer sampling date 
lm_etn_category_microcosm_control <- lm(etn ~ system + 0, data = microcosm_data_control)
summary(lm_etn_category_microcosm_control)

lm_etn_category_microcosm_moderate <- lm(etn ~ system + 0, data = microcosm_data_moderate)
summary(lm_etn_category_microcosm_moderate)

lm_etn_category_microcosm_high <- lm(etn ~ system + 0, data = microcosm_data_high)
summary(lm_etn_category_microcosm_high)

# Add model estimates into dataframe for GGPLOT
est <- c(10.945, 6.489, 15.216)
system <- c("annual", "perennial", "forest")
incubation <- c("control", "control", "control")
category_etn_effects_microcosm_control <- data.frame(est, incubation, system)

est <- c(56.499, 108.709, 60.880)
system <- c("annual", "perennial", "forest")
incubation <- c("moderate", "moderate", "moderate")
category_etn_effects_microcosm_moderate <- data.frame(est, incubation, system)

est <- c(47.709, 76.727, 54.406)
system <- c("annual", "perennial", "forest")
incubation <- c("high", "high", "high")
category_etn_effects_microcosm_high <- data.frame(est, incubation, system)

merge_effects_category_microcosm <- list(category_etn_effects_microcosm_control,
                                         category_etn_effects_microcosm_moderate,
                                         category_etn_effects_microcosm_high)
merge_effects_category_microcosm <- merge_effects_category_microcosm %>% purrr::reduce(full_join)

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: control hour
letters_df_etn_microcosm_control <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ treatment, data = microcosm_data_control))$treatment[,4])$Letters)

colnames(letters_df_etn_microcosm_control)[1] <- "Letter" #Reassign column name
letters_df_etn_microcosm_control$treatment <- rownames(letters_df_etn_microcosm_control) 

letters_df_etn_microcosm_control$placement <- c(50, 50, 50, 50, 50, 50, 50, 50)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "control",
                           "T3" = "control",
                           "T4" = "control",
                           "T5" = "control",
                           "T7" = "control",
                           "DF" = "control",
                           "SF" = "control",
                           "CF" = "control")
letters_df_etn_microcosm_control$system <- value_map_system[letters_df_etn_microcosm_control$treatment]
letters_df_etn_microcosm_control$system <- factor(letters_df_etn_microcosm_control$system, levels = c("annual", "perennial", "forest"))
letters_df_etn_microcosm_control$incubation <- value_map_incubation[letters_df_etn_microcosm_control$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: moderate hour
letters_df_etn_microcosm_moderate <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ treatment, data = microcosm_data_moderate))$treatment[,4])$Letters)

colnames(letters_df_etn_microcosm_moderate)[1] <- "Letter" #Reassign column name
letters_df_etn_microcosm_moderate$treatment <- rownames(letters_df_etn_microcosm_moderate) 

letters_df_etn_microcosm_moderate$placement <- c(150, 150, 150, 150, 150, 150, 150, 150)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "moderate",
                           "T3" = "moderate",
                           "T4" = "moderate",
                           "T5" = "moderate",
                           "T7" = "moderate",
                           "DF" = "moderate",
                           "SF" = "moderate",
                           "CF" = "moderate")
letters_df_etn_microcosm_moderate$system <- value_map_system[letters_df_etn_microcosm_moderate$treatment]
letters_df_etn_microcosm_moderate$system <- factor(letters_df_etn_microcosm_moderate$system, levels = c("annual", "perennial", "forest"))
letters_df_etn_microcosm_moderate$incubation <- value_map_incubation[letters_df_etn_microcosm_moderate$treatment]

# Add Tukey's Post-Hoc Comparisons for each MCSE Treatment: high hour
letters_df_etn_microcosm_high <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ treatment, data = microcosm_data_high))$treatment[,4])$Letters)

colnames(letters_df_etn_microcosm_high)[1] <- "Letter" #Reassign column name
letters_df_etn_microcosm_high$treatment <- rownames(letters_df_etn_microcosm_high) 

letters_df_etn_microcosm_high$placement <- c(150, 150, 150, 150, 150, 150, 150, 150)
value_map_system <- c("T1" = "annual",
                      "T3" = "annual",
                      "T4" = "annual",
                      "T5" = "perennial",
                      "T7" = "perennial",
                      "DF" = "forest",
                      "SF" = "forest",
                      "CF" = "forest")
value_map_incubation <-  c("T1" = "high",
                           "T3" = "high",
                           "T4" = "high",
                           "T5" = "high",
                           "T7" = "high",
                           "DF" = "high",
                           "SF" = "high",
                           "CF" = "high")
letters_df_etn_microcosm_high$system <- value_map_system[letters_df_etn_microcosm_high$treatment]
letters_df_etn_microcosm_high$system <- factor(letters_df_etn_microcosm_high$system, levels = c("annual", "perennial", "forest"))
letters_df_etn_microcosm_high$incubation <- value_map_incubation[letters_df_etn_microcosm_high$treatment]

# Merge all letters together
merge_tukey_etn_microcosm <- list(letters_df_etn_microcosm_control,
                                  letters_df_etn_microcosm_moderate,
                                  letters_df_etn_microcosm_high)
merge_tukey_etn_microcosm <- merge_tukey_etn_microcosm %>% purrr::reduce(full_join)

# Reorder Treatments
merge_effects_category_microcosm$system <- factor(merge_effects_category_microcosm$system, levels = c("annual", "perennial", "forest"))
merge_effects_category_microcosm$incubation <- factor(merge_effects_category_microcosm$incubation, levels = c("control", "moderate", "high"))
merge_tukey_etn_microcosm$incubation <- factor(merge_tukey_etn_microcosm$incubation, levels = c("control", "moderate", "high"))

microcosm_etn_plot <- microcosm_data %>%
  ggplot(aes(x=treatment, y = etn, fill = treatment)) +
  geom_boxplot(alpha = 0.90, outlier.shape = NA) +
  facet_grid(incubation ~system, scales = 'free_x') +
  geom_text(data = merge_tukey_etn_microcosm, 
            aes(x=treatment, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  geom_hline(aes(yintercept=est), linetype = 'dashed', 
             col = 'black', linewidth = 1.25,
             data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76','#86A592','#577B96','#B0381D'), 
                    name = "MCSE Treatment",
                    labels = c("T1: Conventional Corn/Soy/Wheat", 
                               "T3: Reduced Input Corn/Soy/Wheat",
                               "T4: Bio-Based Corn/Soy/Wheat", 
                               "T5: Poplar", "T7: Early Successional",
                               "SF: Successional Forest", "DF: Deciduous Forest", 
                               "CF: Coniferous Forest")) +
  ylab(ETN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_etn_plot


## ETN By Land Use Category #######################################################
# Unite land-use category and biodiversity reduction into one column 
microcosm_data <- microcosm_data %>%
  unite(system_incubation, c(system, incubation), remove = FALSE)

letters_df_diversity_microcosm <- data.frame(multcompLetters(TukeyHSD(aov(etn ~ system_incubation, data = microcosm_data))$system_incubation[,4])$Letters)

colnames(letters_df_diversity_microcosm)[1] <- "Letter" #Reassign column name
letters_df_diversity_microcosm$system_incubation <- rownames(letters_df_diversity_microcosm) 
letters_df_diversity_microcosm$placement <- c(160, 160, 160, 160, 160, 160, 160, 160, 160)

microcosm_data$system_incubation <- factor(microcosm_data$system_incubation, levels = c("annual_control",
                                                                                        "perennial_control",
                                                                                        "forest_control",
                                                                                        "annual_moderate",
                                                                                        "perennial_moderate",
                                                                                        "forest_moderate",
                                                                                        "annual_high",
                                                                                        "perennial_high",
                                                                                        "forest_high"
))
letters_df_diversity_microcosm$system_incubation <- factor(letters_df_diversity_microcosm$system_incubation, levels = c("annual_control",
                                                                                                                        "perennial_control",
                                                                                                                        "forest_control",
                                                                                                                        "annual_moderate",
                                                                                                                        "perennial_moderate",
                                                                                                                        "forest_moderate",
                                                                                                                        "annual_high",
                                                                                                                        "perennial_high",
                                                                                                                        "forest_high"
))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))

microcosm_etn <- microcosm_data %>%
  ggplot(aes(x=system_incubation, y = etn)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  geom_text(data = letters_df_diversity_microcosm, 
            aes(x = system_incubation, y=placement, label=Letter), 
            size =5.5, color="black", fontface="bold") +
  #geom_hline(aes(yintercept=est), linetype = 'dashed', 
  #col = 'black', linewidth = 1.25,
  #data=merge_effects_category_microcosm) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(ETN~('\u00b5g'~N~g^-1~dry~soil)) +
  xlab("Fumigation Exposure Treatment") +
  scale_x_discrete(labels=c("", "Control", "", "", "Moderate (24h)", "", "", "High (72h)", "")) +
  geom_vline(xintercept = c(3.5, 6.5), linetype = "dotted", size = 0.75) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 11), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_etn






# Ordination ####
# Obtain NMDS point coordinates for each sample, and append relevant metadata from the env dataframe 
data.scores <- as.data.frame(scores(bray_NMDS_microcosm_ord)$sites)
data.scores$treatment <- env$treatment
data.scores$system <- env$system
data.scores$sample_date <- env$sample_date
data.scores$incubation <- env$incubation
data.scores$soil_temp <- env$soil_temp
data.scores$soil_moisture <- env$soil_moisture
data.scores$n.fix <- env$n.fix
data.scores$mbc <- env$mbc
data.scores$mbn <- env$mbn
data.scores$eoc <- env$eoc
data.scores$etn <- env$etn

# Extract vector coordinates of metadata variables 
en_coord_cont <- as.data.frame(scores(ord.fit.sig, "vectors"))
# Extract p-values of each vector
en_coord_cont$pval <- ord.fit.sig$vectors$pvals
# Subset vectors that are statistically significant (p <0.05)
en_coord_cont_sig <- en_coord_cont %>% filter(pval < 0.05)
# Add names of relevant vectors to the dataframe -- for the ggplot
en_coord_cont_sig$name <- c("FLNF", "MBN", "EOC", "ETN")

# Add ellipse for system level variable
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#data for ellipse, in this case using the system factor
data.scores$system <- as.factor(data.scores$system)
df_ell.system <- data.frame() #sets up a data frame before running the function.
for(g in levels(data.scores$system)){
  df_ell.system <- rbind(df_ell.system, cbind(as.data.frame(with(data.scores [data.scores$system==g,],
                                                                       veganCovEllipse(cov.wt(cbind(NMDS1,NMDS2),
                                                                                              wt=rep(1/length(NMDS1),
                                                                                                     length(NMDS1)))$cov,
                                                                                       center=c(mean(NMDS1),
                                                                                                mean(NMDS2))))),system=g))
}

# data for labelling the ellipse
NMDS.mean.system<-aggregate(data.scores[ ,c("NMDS1", "NMDS2")], 
                               list(group = data.scores$system), mean)

# data for labelling the ellipse
NMDS.mean=aggregate(data.scores[,c("NMDS1", "NMDS2")], 
                    list(group = data.scores$system), mean)






# Reorganize treatment factor 

data.scores$system<- factor(data.scores$system, levels = c("annual", "perennial", "forest"))

df_ell.system$system <- factor(df_ell.system$system, levels = c("annual", "perennial", "forest"))


# Change the name of the incubation factors 
data.scores$incubation <- as.factor(data.scores$incubation)
levels(data.scores$incubation) <- c("control", "moderate", "high") # CHECK AND SEE IF THIS MATCHES ORIGINAL NAMES

# Display NMDS ordination of our in-situ dataset
microcosm_nmds_ordination <- data.scores %>%
  ggplot(aes(x=NMDS1,y=NMDS2)) +
  geom_point(aes(color = system, shape = incubation), size = 5.5,  alpha = 1.0) +
  geom_path(data = df_ell.system, aes(x = NMDS1, y = NMDS2, group = system, color = system), linetype = 1, size = 1, alpha = 1.0) +
  geom_polygon(data = df_ell.system, aes(x = NMDS1, y = NMDS2, fill = system), alpha = 0.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "darkgrey", linewidth = 1) +
  geom_vline(xintercept=0, linetype = "dashed", color = "darkgrey", linewidth = 1) +
  geom_text(x = -1.8, y = 1.40, label = "stress = 0.197") +
  scale_color_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                     labels = c("Annual", "Perennial", "Forest")) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  scale_shape_manual(values =c(16, 17, 18), name = "Fumigation Exposure", 
                     labels = c("Control", "Moderate (24h)", "High (72h)")) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont_sig, size =1.25, alpha = 0.7, colour = "black") +
  geom_label(data = en_coord_cont_sig, aes(x = NMDS1, y = NMDS2), colour = "black", 
             fontface = "bold", size = 5.5, label = en_coord_cont_sig$name, fill = "white") + 
  guides(fill = guide_legend(override.aes = list(shape=c(16)))) +
  ggtitle("Biodiversity Manipulation Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15)) + 
  theme(axis.title.y = element_text(size = 15)) + 
  theme(axis.text.x = element_text(size = 12, color = "black"))  + 
  theme(axis.text.y = element_text(size = 12, color = "black"))  + 
  theme(legend.text = element_text(size = 12)) +
  theme(plot.margin=unit(c(.5,1,.5,.5),"cm")) +
  theme(legend.position="right", legend.title=element_text(size=14))+
  theme(legend.position="right")

microcosm_nmds_ordination

# Differential Abundance ####
##### Control v. High ---------------
microcosm_deseq <- sigtab_results_TAX_genra_control_v_high %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  facet_grid(rows = vars(incubation), scales = 'free_x') +
  theme_bw() +
  geom_point(size = 5.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth = 0.95) +
  #geom_text(x=24.5, y=1.5, label="Control Enriched", size = 5) +
  #geom_text(x=24.5, y=-1.5, label="Moderate Enriched", size = 5) +
  scale_fill_manual(values=met.brewer("Egypt", 11)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =13),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        plot.margin=unit(c(.5,1,.5,.5),"cm"),
        legend.title = element_text(size =15))
microcosm_deseq

##### Control v. Moderate ---------------
microcosm_deseq_control_v_moderate <- sigtab_results_TAX_genra_control_v_moderate %>%
  ggplot(aes(x=Genus, y=log2FoldChange, fill=Class)) + 
  facet_grid(rows = vars(incubation), scales = 'free_x') +
  theme_bw() +
  geom_point(size = 5.75, pch = 21) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", linewidth = 0.95) +
  #geom_text(x=24.5, y=1.5, label="Control Enriched", size = 5) +
  #geom_text(x=24.5, y=-1.5, label="Moderate Enriched", size = 5) +
  scale_fill_manual(values=met.brewer("Egypt", 13)) +
  theme(axis.text.x = element_text(angle = -75, hjust = 0, vjust=0.5, size =13), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size =15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size =15),        
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_deseq_control_v_moderate



# FLNF vs. Diversity ####


est_total <- c("Est = 0.017")
R_total <- c("Adj-R.sq = 0.23")
p_total <- c("p = 0.01")

model.df <- data.frame(est_total, R_total, p_total)

microcosm_flnf_div <- microcosm_data %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(color = treatment), size = 4.5) +  
  geom_smooth(method = 'lm') +
  geom_text(x=420, y = 14.5, aes(label = est_total), size = 5, data = model.df) +
  geom_text(x=420, y = 13.5, aes(label = R_total), size = 5, data = model.df) +
  geom_text(x=420, y = 12.5, aes(label = p_total), size = 5, data = model.df) +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                     labels = c(
                       "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(16, 15, 18), name = "Fumigation Exposure", 
                     labels = c("Control", "Moderate (24h)", "High (72h)")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Hill Number (q=1)") + 
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  ggtitle("Biodiversity Manipulation Experiment") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 12), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_flnf_div

# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
lm.annual <- lm(n.fix ~ shan, data = microcosm_data_annual)
summary(lm.annual)

# Perennial
lm.perennial <- lm(n.fix ~ shan, data = microcosm_data_perennial)
summary(lm.perennial)

# Forest
lm.forest <- lm(n.fix ~ shan, data = microcosm_data_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = 0.01", "Est = 0.007", "Est = 0.04" )
R <- c("Adj-R.sq = 0.087", "Adj-R.sq = 0.027", "Adj-R.sq = 0.49")
p <- c("p = 0.045", "p = 0.21", "p<0.001")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# Create Plot with system-model estimates
microcosm_data$treatment <- factor(microcosm_data$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))
microcosm_data$incubation <- as.factor(microcosm_data$incubation)


microcosm_category_flnf_div <- microcosm_data %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(shape = incubation, color = treatment), size = 4.5) +  
  facet_grid(~system) +
  geom_text(x=380, y = 15, aes(label = est), size = 5, data = model.df) +
  geom_text(x=380, y = 14, aes(label = R), size = 5, data = model.df) +
  geom_text(x=380, y = 13, aes(label = p), size = 5, data = model.df) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                             '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                    labels = c(
                      "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                      "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                      "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(16, 17, 18), name = "Fumigation Exposure", 
                     labels = c("Control", "Moderate (24hr)", "High (72hr)")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab("Diazotroph Alpha Diversity") + 
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

microcosm_category_flnf_div

# MBC vs. Diversity ####
# Generate linear models for each system: Annual, Perennial, and Forest.
# Annual 
lm.annual <- lm(n.fix ~ mbc, data = microcosm_data_annual)
summary(lm.annual)

# Perennial
lm.perennial <- lm(n.fix ~ mbc, data = microcosm_data_perennial)
summary(lm.perennial)

# Forest
lm.forest <- lm(n.fix ~ mbc, data = microcosm_data_forest)
summary(lm.forest)

# Add model information into dataframe
est <- c("Est = -0.01", "Est = 0.01", "Est = -0.001" )
R <- c("Adj-R.sq = 0.19", "Adj-R.sq = 0.38", "Adj-R.sq < 0")
p <- c("p = 0.004", "p < 0.001 ", "p = 0.61")
system <- c("annual", "perennial", "forest")

model.df <- data.frame(est, R, p, system)
model.df$system <- factor(model.df$system, levels = c("annual", "perennial", "forest"))

# Create Plot with system-model estimates
microcosm_data$treatment <- factor(microcosm_data$treatment, levels = c("T1", "T3", "T4", "T5", "T7", "SF", "CF", "DF"))
microcosm_data$system <- factor(microcosm_data$system, levels = c("annual", "perennial", "forest"))
microcosm_data$incubation <- as.factor(microcosm_data$incubation)


microcosm_category_flnf_mbc <- microcosm_data %>%
  ggplot(aes(x=mbc, y = n.fix)) +
  geom_point(aes(shape = incubation, color = treatment), size = 4.5) +  
  facet_grid(~system) +
  geom_text(x=620, y = 15, aes(label = est), size = 5, data = model.df) +
  geom_text(x=620, y = 14, aes(label = R), size = 5, data = model.df) +
  geom_text(x=620, y = 13, aes(label = p), size = 5, data = model.df) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values=c('#EEAC4E','#D7785D','#9C7C50', '#577753',
                              '#1A4C76', '#86A592', '#B0381D','#577B96'), name = "MCSE Treatment",
                     labels = c(
                       "T1: Conventional Corn/Soy/Wheat", "T3: Reduced Input Corn/Soy/Wheat",
                       "T4: Bio-Based Corn/Soy/Wheat", "T5: Perennial Poplar", "T7: Early Successional Grassland",
                       "SF: Successional Forest", "CF: Coniferous Forest", "DF: Deciduous Forest")) +
  scale_shape_manual(values =c(16, 17, 18), name = "Fumigation Exposure", 
                     labels = c("Control", "Moderate (0hr)", "High (72hr)")) +
  guides(fill = guide_legend(override.aes = list(shape=c(21)))) +
  xlab(MBC~('\u00b5g'~C~g^-1~dry~soil)) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  theme_bw() +
  theme(axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 12), 
        strip.text.x = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 12), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
microcosm_category_flnf_mbc


#Save Figures ######
# Combine Figures:
# MCSE Summer + Fall DESEQ Plot + MCSE Ordination
mcse_vs_forest_ordination_combined <- ggarrange(
  mcse_nmds_ordination,
  mcse_vs_forest_plot,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)

mcse_vs_forest_ordination_combined
# Microcosm DESEQ plots + Microcosm Ordination
microcosm_deseq_ordination_combined <- ggarrange(
  microcosm_nmds_ordination,
  microcosm_deseq, 
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = F,
  legend = "right"
)
microcosm_deseq_ordination_combined
# Microcosm diversity and FLNF with biodiversity reduction 
microcosm_div_flnf_combined <- ggarrange(
  microcosm_div,
  microcosm_flnf,
  labels = c("A", "B"),
  ncol = 1, 
  nrow = 2,
  common.legend = T,
  legend = "right"
)

# Microcosm MBC, MBN, EOC, ETN with biodiversity reduction 
microcosm_mbc_mbn_eoc_etn_combined <- ggarrange(
  microcosm_mbc,
  microcosm_mbn,
  microcosm_eoc,
  microcosm_etn,
  labels = c("A", "B", "C", "D"),
  ncol = 1, 
  nrow = 4,
  common.legend = T,
  legend = "right"
)




# FLNF v Diversity Plots 
flnf_diversity_combined <- ggarrange(
  mcse_flnf_div,
  microcosm_flnf_div,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)
flnf_diversity_combined

# FLNF v Diversity Plots: By Land-Use Category  
category_flnf_diversity_mbc_combined <- ggarrange(
  microcosm_category_flnf_div,
  microcosm_category_flnf_mbc,
  labels = c("A", "B"),
  ncol = 1,
  nrow = 2,
  common.legend = T,
  legend = "right"
)
category_flnf_diversity_mbc_combined

setwd("~/Library/CloudStorage/Dropbox/Projects/LTER_MCSE/data/Figures/")
#Windows OS
#setwd("C:/Users/Brand/Dropbox/Projects/LTER_MCSE/data/Figures") # Windows OS directory location

# Individual Figures
ggsave(filename = "mcse_div.png", plot = mcse_diversity_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "mcse_flnf.png", plot = mcse_flnf_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "mcse_mbc.png", plot = mcse_mbc_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "mcse_mbn.png", plot = mcse_mbn_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "mcse_EOC.png", plot = mcse_eoc_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "mcse_ETN.png", plot = mcse_etn_plot, w = 10, h = 5, dpi = 400)
ggsave(filename = "deseq_summer_v_fall.png", plot = mcse_summer_vs_fall_plot, w = 10, h =5, dpi =400)



ggsave(filename = "microcosm_div.png", plot = microcosm_diversity_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcosm_flnf.png", plot = microcosm_flnf_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcosm_mbc.png", plot = microcosm_mbc_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcosm_mbn.png", plot = microcosm_mbn_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcosm_EOC.png", plot = microcosm_eoc_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcosm_ETN.png", plot = microcosm_etn_plot, w = 10, h = 7, dpi = 400)
ggsave(filename = "microcom_flnf_diversity_system.png", plot = microcosm_category_flnf_div, w = 15, h =7, dpi =400)
ggsave(filename = "deseq_control_v_moderate.png", plot = microcosm_deseq_control_v_moderate, w = 8, h =7, dpi =400)


# Combined Figures
ggsave(filename = "mcse_vs_forest_ordination_combined.png", plot = mcse_vs_forest_ordination_combined, w = 10, h = 10, dpi = 400)
ggsave(filename = "microcosm_deseq_ordination_combined.png", plot = microcosm_deseq_ordination_combined, w = 10, h = 10, dpi = 500)
ggsave(filename = "flnf_diversity_combined.png", plot = flnf_diversity_combined, w = 8, h = 8, dpi = 400)
ggsave(filename = "category_flnf_diversity_mbc_combined.png", plot = category_flnf_diversity_mbc_combined, w = 14, h = 11, dpi = 400)
ggsave(filename="biodiversity_reduction_microcosm_combined.png", plot = microcosm_mbc_mbn_eoc_etn_combined, w = 11, h = 12, dpi = 600)

### SAVE MAIN TEXT FIGURES FOR JOURNAL SUBMISSION
ggsave(filename="Figure_1.png", plot = mcse_diversity_plot, w = 10, h =5, dpi = 600)
ggsave(filename="Figure_2.png", plot = mcse_vs_forest_ordination_combined, w = 10, h =10, dpi = 600)
ggsave(filename="Figure_3.png", plot = microcosm_div_flnf_combined , w = 8, h =6, dpi = 600)
ggsave(filename="Figure_4.png", plot = microcosm_deseq_ordination_combined, w = 10, h =10, dpi = 600)

ggsave(filename="Figure_5.png", plot = category_flnf_diversity_mbc_combined, w = 14, h =11, dpi = 600)




## Standard curve calibration analysis ################
setwd("~/Library/CloudStorage/Dropbox/Projects/LTER_MCSE/data/TOC_raw/")
library(readxl)
TN_curve_example <- read_excel("TN_curve_example.xlsx")


eq <- c("y = 0.0032*[TN Area] - 2.0806")
R_total <- c("R.sq = 0.997")
model.df <- data.frame(eq, R_total)


p1 <- TN_curve_example %>%
      ggplot(aes(x=`TN Area`, y = `TN [mg/l]`)) +
      geom_text(x=12500, y = 100, aes(label = eq), size = 3.5, data = model.df) +
      geom_text(x=12500, y = 95, aes(label = R_total), size = 3.5, data = model.df) +
      geom_point(aes(col = `Sample Type`), size = 2.5) +
      geom_smooth(method = 'lm')+
      theme_classic() 

TN_curve_example_filter <- TN_curve_example%>%
  filter(Sample_Name != "lo gly 5") 

eq <- c("y = 0.0024*[TN Area] + 0.8698")
R_total <- c("R.sq = 0.9987")
model.df <- data.frame(eq, R_total)


p2 <- TN_curve_example_filter %>%
      ggplot(aes(x=`TN Area`, y = `TN [mg/l]`)) +
      geom_text(x=3500, y = 20, aes(label = eq), size = 3.5, data = model.df) +
      geom_text(x=3500, y = 19, aes(label = R_total), size = 3.5, data = model.df) +
      geom_point(aes(col = `Sample Type`), size = 2.5) +
      geom_smooth(method = 'lm')+
      theme_classic() 

TN_example <- ggarrange(
  p1,
  p2,
  labels = c("Full Standard Curve", "Corrected Standard Curve"),
  ncol = 2,
  nrow = 1,
  common.legend = T,
  legend = "right"
)
ggsave("TN_standards.png", plot = TN_example, dpi = 300, w=10, h=5)

## ASE Graphical Abstract ###########
mcse_data$system <- factor(mcse_data$system, levels = c("annual", "perennial", "forest"))

mcse_data$sample_date<- factor(mcse_data$sample_date, levels = c("summer", "fall"))


category_div <-mcse_data %>%
  ggplot(aes(x=system, y = shan)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  facet_grid(~sample_date) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab("Diazotroph Diversity") +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
  

category_flnf <- mcse_data %>%
  ggplot(aes(x=system, y = n.fix)) +
  geom_boxplot(aes(fill = system), outlier.shape = NA) +
  facet_grid(~sample_date) +
  scale_fill_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

div_plot <-mcse_data %>%
  ggplot(aes(x=shan, y = n.fix)) +
  geom_point(aes(col = system), size = 2.5) +
  geom_smooth(method = 'lm', se = FALSE)+
  scale_color_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                    labels = c("Annual", "Perennial", "Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("Diazotroph Diversity") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))

mcse_data %>%
  ggplot(aes(x=soil_moisture, y = n.fix, col = system)) +
  geom_point() +
  geom_smooth(method = 'lm', se = FALSE)+
  scale_color_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                     labels = c("Annual", "Perennial", "Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("Soil Moisture") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))


div_plot_2<- microcosm_data %>%
  ggplot(aes(x=shan, y = n.fix, col = system)) +
  geom_point(size = 2.5) +
  geom_smooth(method = 'lm', se = FALSE)+
  scale_color_manual(values=c('#F0C96C', '#6E9867', '#5B65A7'), name = "Land-use Category",
                     labels = c("Annual", "Perennial", "Forest")) +
  ylab(FLNF~('\u00b5g'~N~g^-1~dry~soil~day^-1)) +
  xlab("Diazotroph Diversity") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14), 
        title = element_text(size = 14),
        axis.title.y = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"), 
        legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14), 
        plot.margin=unit(c(.5,1,.5,.5),"cm"))
  
  
  
  
  
  
