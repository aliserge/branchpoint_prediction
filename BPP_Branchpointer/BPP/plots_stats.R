library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

options(digits=10)


# BPP

mutation_impact_BPP <- function(x) {
  x_temp <- data.frame(x)
  x_temp$prediction_changed <- FALSE
  x_temp$prediction_changed[x_temp$ref_bp_pos != x_temp$alt_bp_pos |
                              x_temp$alt_zsc < x_temp$ref_zsc] <- TRUE
  x_temp$delta_score[x_temp$prediction_changed == TRUE] <- x_temp$ref_zsc[x_temp$prediction_changed == TRUE] - x_temp$alt_zsc[x_temp$prediction_changed == TRUE] 
  return(x_temp)
}

# Uploading BPP results for Clinvar and Gnomad

gnomad_BPP <- read.csv2("~/Bioinformatics/branchpoint/BPP_materials/gnomad_downsampled/gnomad_all_hg38.GRCh38.fa.bpp.combined_per_variation_supersmall.tsv", 
                        header=TRUE, sep = '\t') # individual mutations (one mutation per intron)
clinvar_BPP <- read.csv2("~/Bioinformatics/branchpoint/BPP_materials/clinvar_one_mutation_per_intron/clinvar.pathogenic.chr.GRCh38.fa.bpp.combined_per_variation.tsv", 
                         header=TRUE, sep = '\t') # individual mutations (one mutation per intron)


#### Gnomad preparation

gnomad_BPP$ref_zsc <- as.numeric(gnomad_BPP$ref_zsc)
gnomad_BPP$alt_zsc <- as.numeric(gnomad_BPP$alt_zsc)
gnomad_BPP$Group <- 'gnomad'
gnomad_BPP <- na.omit(gnomad_BPP)


# Mutation impact calculation (did the branchpoint position change? did the Z-score change?)

gnomad_BPP <- mutation_impact_BPP(gnomad_BPP)


# Filtering out the delta-scores which are too small (less than 0.001 in absolute range)

gnomad_BPP_nonzero <- filter(gnomad_BPP, abs(delta_score) > 0.001)


# How many observations do we have after filtering?

gnomad_observations <- nrow(gnomad_BPP)
cat('From all found intronic mutations', gnomad_observations, 'have non-zero impact')


# Removing extra-columns in aim of memory economy

gnomad_BPP_nonzero_small <- gnomad_BPP_nonzero %>% select(Group, ref_bp_pos, alt_bp_pos, delta_score)



#### Clinvar preparation

clinvar_BPP$ref_zsc <- as.double(clinvar_BPP$ref_zsc)
clinvar_BPP$alt_zsc <- as.double(clinvar_BPP$alt_zsc)
clinvar_BPP <- na.omit(clinvar_BPP)
clinvar_BPP$Group <- 'clinvar'


# Mutation impact calculation (did the branchpoint position change? did the Z-score change?)

clinvar_BPP <- mutation_impact_BPP(clinvar_BPP)


# Filtering out the delta-scores which are too small (less than 0.001 in absolute range)

clinvar_BPP_nonzero <- filter(clinvar_BPP, abs(delta_score) > 0.001)


# How many observations do we have after filtering?

clinvar_observations <- nrow(clinvar_BPP)
cat('From all found intronic mutations', clinvar_observations, 'have non-zero impact')


# Removing extra-columns in aim of memory economy

clinvar_BPP_nonzero_small <- clinvar_BPP_nonzero %>% select(Group, ref_bp_pos, alt_bp_pos, delta_score)



#### Combining of the Clinvar & Gnomad results and making plots

clinvar_gnomad_BPP_small <- rbind(clinvar_BPP_nonzero_small, gnomad_BPP_nonzero_small)


# Save PNG

png(filename="~/Bioinformatics/branchpoint/BPP_final_plot.png", 
    width = 1280, height = 720)
ggplot(data = clinvar_gnomad_BPP_small, aes(x = Group, y = delta_score)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.3) +
  ylab("Delta Z-score") +
  xlab("Group") +
  ggtitle("Clinvar & Gnomad via BPP") +
  scale_y_continuous(limits=c(-2.5, 2.5))
dev.off()


wilcox.test(delta_score ~ Group, data=clinvar_gnomad_BPP_small)



# Branchpointer

mutation_impact_Branchpointer <- function(x) {
  x_temp <- data.frame(x)
  x_temp$prediction_changed <- FALSE
  x_temp$prediction_changed[x_temp$dist_to_BP_REF != x_temp$dist_to_BP_ALT | x_temp$max_prob_ALT < x_temp$max_prob_REF] <- TRUE
  x_temp$delta_score[x_temp$prediction_changed == TRUE] <- x_temp$max_prob_REF[x_temp$prediction_changed == TRUE] - x_temp$max_prob_ALT[x_temp$prediction_changed == TRUE] 
  return(x_temp)
}

#### Clinvar

clinvar_Branchpointer <- read.csv("~/Bioinformatics/branchpoint/Branchpointer_materials/clinvar/hg38.branchpointer.clinvar.pathogenic.chr.summary.combined.csv", 
                                 header=TRUE, sep=",")


clinvar_Branchpointer$max_prob_REF <- as.double(clinvar_Branchpointer$max_prob_REF)
clinvar_Branchpointer$max_prob_ALT <- as.double(clinvar_Branchpointer$max_prob_ALT)
clinvar_Branchpointer$seqnames <- as.factor(clinvar_Branchpointer$seqnames)

clinvar_Branchpointer <- na.omit(clinvar_Branchpointer)
clinvar_Branchpointer$Group <- 'clinvar'

clinvar_Branchpointer <- mutation_impact_Branchpointer(clinvar_Branchpointer)

clinvar_nonzero <- filter(clinvar_Branchpointer, abs(delta_score) > 0.001)

sum(clinvar_Branchpointer$prediction_changed == TRUE) / sum(clinvar_Branchpointer$prediction_changed == FALSE)


# Save PNG

png(filename="~/Bioinformatics/branchpoint/Branchpointer_materials/clinvar_Branchpointer.png", 
    width = 1280, height = 720)
ggplot(data = clinvar_nonzero, aes(x = Group, y = delta_score)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.3) +
  ylab("Delta-score") +
  xlab("dataset") +
  ggtitle("Delta-score in ClinVar via Branchpointer")
dev.off()


#### Gnomad

gnomad_Branchpointer <- read.csv("~/Bioinformatics/branchpoint/Branchpointer_materials/summary_original/gnomad_all_hg38.branchpointer.combined.original.csv", 
                                 header=TRUE, sep=",")

gnomad_Branchpointer$max_prob_REF <- as.double(gnomad_Branchpointer$max_prob_REF)
gnomad_Branchpointer$max_prob_ALT <- as.double(gnomad_Branchpointer$max_prob_ALT)
gnomad_Branchpointer$seqnames <- as.factor(gnomad_Branchpointer$seqnames)

gnomad_Branchpointer <- na.omit(gnomad_Branchpointer)

gnomad_Branchpointer$Group <- 'gnomad'

gnomad_Branchpointer <- mutation_impact_Branchpointer(gnomad_Branchpointer)
gnomad_Branchpointer <- na.omit(gnomad_Branchpointer)

gnomad_Branchpointer_nonzero <- filter(gnomad_Branchpointer, abs(delta_score) > 0.001)


# Save PNG

png(filename="~/Bioinformatics/branchpoint/Branchpointer_materials/gnomad_Branchpointer.png", 
    width = 1280, height = 720)
ggplot(data = gnomad_Branchpointer_nonzero, aes(x = Group, y = delta_score)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.3) +
  ylab("Delta score") +
  xlab("Branchpoint") +
  ggtitle("Delta score in Gnomad via Branchpointer")
dev.off()


# Removing extra-columns in aim of memory economy

clinvar_Branchpointer_small <- clinvar_nonzero %>% select(Group, delta_score, prediction_changed)
gnomad_branchpointer_orig_small <- gnomad_Branchpointer_nonzero %>% select(Group, delta_score, prediction_changed)


#### Combining of the Clinvar & Gnomad results and making plots

clinvar_gnomad_branchpointer <- rbind(clinvar_Branchpointer_small, gnomad_branchpointer_orig_small)


# Save as PNG

png(filename="~/Bioinformatics/branchpoint/Branchpointer_materials/clinvar_gnomad_Branchpointer.png", 
    width = 1280, height = 720)

ggplot(data = clinvar_gnomad_branchpointer, aes(x = Group, y = delta_score)) + 
  geom_violin(alpha = 0.5) +
  geom_boxplot(alpha = 0.3) +
  ylab("Delta score") +
  xlab("Branchpoint") +
  ggtitle("Gnomad & Clinvar via Branchpointer")

dev.off()

wilcox.test(delta_score ~ Group, data=clinvar_gnomad_branchpointer)

