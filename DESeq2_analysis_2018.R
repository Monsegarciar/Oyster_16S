# DESeq2 Analysis on Measurments ####
# 2021-07-14
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)


# Loading data ####

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")


#Loading Physeq w/out transform_sample_counts() function ####
physeq_class18 <- readRDS("Data/physeq_class18.rds")
physeq_count18 <- readRDS("Data/physeq_count18.rds")


# Normalized Weight with DESeq2- 2018 Data
physeq_class18 = subset_samples(physeq_class18, Weight_delta != "NA")

deseq18_weight = phyloseq_to_deseq2(physeq_count18, ~ Weight_delta)
deseq18_weight = DESeq(deseq18_weight, test="Wald", fitType="parametric") 

res18_weight = results(deseq18_weight, cooksCutoff = FALSE)
alpha = 0.01
sigtab18_weight = res18_weight[which(res18_weight$padj < alpha), ]
sigtab18_weight = cbind(as(sigtab18_weight, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_weight), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab18_weight$log2FoldChange, sigtab18_weight$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_weight$Phylum = factor(as.character(sigtab18_weight$Phylum), levels=names(x))

x = tapply(sigtab18_weight$log2FoldChange, sigtab18_weight$Genus, function(x) max(x))
x = sort(x, TRUE)


sigtab18_weight$Genus = factor(as.character(sigtab18_weight$Genus), levels=names(x))
ggplot(sigtab18_weight, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

