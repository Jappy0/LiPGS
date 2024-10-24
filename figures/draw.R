# Install and load necessary packages

if (!requireNamespace("qqman", quietly = TRUE)) install.packages("qqman")

if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(qqman)

library(dplyr)

 

# Assuming the data frame is already loaded and named Meta_GWAS_Summary

# and has columns: CHR, POS, P

 

# Check for missing values in the necessary columns

Meta_GWAS_Summary <- Meta_GWAS_Summary %>%

  filter(!is.na(CHR) & !is.na(POS) & !is.na(P))

 

# Calculate -log10(P) for the plot

Meta_GWAS_Summary <- Meta_GWAS_Summary %>%

  mutate(logP = -log10(P))

 

# Define the genome-wide significance threshold

genome_wide_threshold <- -log10(5e-8)

 

# Create a custom color vector

Meta_GWAS_Summary$color <- ifelse(Meta_GWAS_Summary$logP >= genome_wide_threshold, "pink",

                                  ifelse(Meta_GWAS_Summary$CHR %% 2 == 0, "blue4", "orange3"))

 

# Ensure all chromosomes from 1 to 22 are labeled

chrlabs <- as.character(1:22)

 

# Define the closest genes for annotation

closest_genes <- data.frame(

  CHR = c(1, 1, 4, 4, 4, 10, 13, 17, 17, 20),

  Gene = c("PTP4A2", "PRPF3", "LCORL", "RN7SL89P", "ANAPC10", "HK1", "DLEU1", "MAPT", "SCN4A", "STAU1")

)

 

# Identify top significant SNPs for annotation

top_snps <- Meta_GWAS_Summary %>%

  filter(logP >= genome_wide_threshold) %>%

  group_by(CHR) %>%

  slice_max(order_by = logP, n = 1) %>%

  ungroup()

 

# Merge top SNPs with closest genes

top_snps <- top_snps %>%

  left_join(closest_genes, by = "CHR")

 

# Create the Manhattan plot with adjusted axis text size and rotation

manhattan(Meta_GWAS_Summary, chr = "CHR", bp = "POS", p = "P", snp = "SNP",

          genomewideline = genome_wide_threshold, suggestiveline = -log10(1e-5),

          col = c("blue4", "orange3"), highlight = Meta_GWAS_Summary$SNP[Meta_GWAS_Summary$logP >= genome_wide_threshold],

          chrlabs = chrlabs, main = "", ylim = c(0, max(Meta_GWAS_Summary$logP) + 1),

          cex.axis = 0.7, las = 2)

 

# Annotate the top significant SNPs with gene names

with(top_snps, text(x = POS, y = logP, labels = Gene, pos = 3, cex = 0.8, col = "red"))