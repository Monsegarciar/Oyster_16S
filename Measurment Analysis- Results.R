# Measurement Analysis- Results ####
# 2022-06-29
# Author: Monse Garcia

# Loading packages
require(phyloseq)
require(ggplot2)
require(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
require(DESeq2)
library(metacoder)
library(ggtree)
library(Rcpp)
library(corrplot)

#### 1.Phyloseq Analysis with scaled volume 

physeq_highsc17 <- readRDS("Data/physeq_highsc17.rds")

taxa_scaled17 = tax_table(physeq_count18_length)

des_highsc17 <- phyloseq_to_deseq2(physeq_highsc17, ~ Volume_scale)
des_highsc17 <- DESeq(des_highsc17, test="Wald", fitType = "parametric")
resultsNames(des_highsc17)
# Intercept , volume_scale

results_highsc17 <- results(des_highsc17, name = "Volume_scale")
significant_highsc17 <- results_highsc17[which(results_highsc17$padj <0.05), ]
dim(significant_highsc17)
# 67, 6

sub_significant_highsc17 <- subset_taxa(prune_taxa(rownames(significant_highsc17), physeq_highsc17))
df_significant_highsc17 <- as.data.frame(tax_table(sub_significant_highsc17))
write.table(df_significant_highsc17, file = "Data/Significant OTU's for Scaled Volume 2017.csv", quote = FALSE, sep = ",", col.names = T)



#### Positive OTUs

significant_highsc17_pos <- significant_highsc17[significant_highsc17$log2FoldChange>0,]
dim(significant_highsc17_pos)
# 25, 6

sub_significant_highsc17_pos <- subset_taxa(prune_taxa(rownames(significant_highsc17), physeq_highsc17))
df_significant_highsc17_pos <- as.data.frame(tax_table(sub_significant_highsc17_pos))
write.table(df_significant_highsc17_pos, file = "Data/Significant Positive OTU's for Scaled Volume 2017.csv", quote = FALSE, sep = ",", col.names = T)

plot_heatmap(sub_significant_highsc17_pos, method = "NMDS", distance = "jsd", low = "#66CCFF", high = "#000033", na.value = "white", taxa.label = "Order", sample.label = "UniqueID")

#### Negative 
significant_highsc17_neg <- significant_highsc17[significant_highsc17$log2FoldChange<0, ]
dim(significant_highsc17_neg)
# 42, 6
replace_na(sub_significant_highsc17_pos)

rank_names(sub_significant_highsc17)

plot_heatmap(sub_significant_highsc17_pos, method = "NMDS", distance = "jsd", low = "#66CCFF", high = "#000033", na.value = "white", taxa.label = "Order", sample.label = "UniqueID")

?plot_heatmap
#Use the mutate_if() from dplyr/tidyverse

#Error in cmdscale(dist, k = k) : NA values not allowed in 'd'

# 2018 
physeq_highsc18 <- readRDS("Data/physeq_highsc18.rds")

taxa_scaled17 = tax_table(physeq_count18_length)

des_highsc18 <- phyloseq_to_deseq2(physeq_highsc18, ~ Volume_scale)
des_highsc18 <- DESeq(des_highsc18, test="Wald", fitType = "parametric")
resultsNames(des_highsc18)
# Intercept , volume_scale

results_highsc18 <- results(des_highsc18, name = "Volume_scale")
significant_highsc18 <- results_highsc18[which(results_highsc18$padj <0.05), ]
dim(significant_highsc18)
# 51, 6

sub_significant_highsc18 <- subset_taxa(prune_taxa(rownames(significant_highsc18), physeq_highsc18))
df_significant_highsc18 <- as.data.frame(tax_table(sub_significant_highsc18))
write.table(df_significant_highsc18, file = "Data/Significant OTU's for Scaled Volume 2018.csv", quote = FALSE, sep = ",", col.names = T)


#### Positive OTUs

significant_highsc18_pos <- significant_highsc18[significant_highsc18$log2FoldChange>0,]
dim(significant_highsc18_pos)
# 13, 6
sub_significant_highsc18_pos <- subset_taxa(prune_taxa(rownames(significant_highsc18_pos), physeq_highsc18))
df_significant_highsc18_pos <- as.data.frame(tax_table(sub_significant_highsc18_pos))
write.table(df_significant_highsc18_pos, file = "Data/Significant Positive OTU's for Scaled Volume 2018.csv", quote = FALSE, sep = ",", col.names = T)

plot_heatmap(sub_significant_highsc17_pos, method = "NMDS", distance = "jsd", low = "#66CCFF", high = "#000033", na.value = "white", taxa.label = "Order", sample.label = "UniqueID")

#### Negative OTUs

significant_highsc18_neg <- significant_highsc18[significant_highsc18$log2FoldChange<0, ]
dim(significant_highsc18_neg)

#### 2.Measurements taxa
# This section covers the significant OTUs associated with measurements; length, width, height. NOTE: These are not scaled nor are they the volume.  

# Taxa- Class
# 1. Proteobacteria 
# 2. Chloroflexi
# 3. Planctomycetota 
# 4. Desulfobacteria 
# 5. Cyanobacteria
# 6. Modulibacteria 
# 7. SAR324 clade (marine group B)

# Any papers on such taxa and relationships with oysters








