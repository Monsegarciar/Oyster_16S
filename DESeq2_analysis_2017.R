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

meta17_data <- read.csv("Data/meta17_data_update.csv")
asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

#Loading Physeq w/out transform_sample_counts() function ####
physeq_class17 <- readRDS("Data/physeq_class17.rds")


# Normalized Weight with DESeq2- 2017 Data ####
physeq_class17 = subset_samples(physeq_class17, Weight_delta != "NA")

deseq17_weight = phyloseq_to_deseq2(physeq_class17, ~ Weight_delta)
deseq17_weight = DESeq(deseq17_weight, test="Wald", fitType="parametric")

res17_weight = results(deseq17_weight, cooksCutoff = FALSE)
alpha = 0.05 # switch to 0.05- generalized linear models
sigtab17_weight = res17_weight[which(res17_weight$padj < alpha), ]
sigtab17_weight = cbind(as(sigtab17_weight, "data.frame"), as(tax_table(physeq_class17)[rownames(sigtab17_weight), ], "matrix"))

# Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Genus 
x = tapply(sigtab17_weight$log2FoldChange, sigtab17_weight$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab17_weight$Genus = factor(as.character(sigtab17_weight$Genus), levels=names(x))

#Class
x = tapply(sigtab17_weight$log2FoldChange, sigtab17_weight$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_weight$Class = factor(as.character(sigtab17_weight$Class), levels=names(x))

ggplot(sigtab17_weight, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Height with DESeq2- 2017 Data ####
physeq_class17 = subset_samples(physeq_class17, Height_delta != "NA")

deseq17_height = phyloseq_to_deseq2(physeq_class17, ~ Height_delta)
deseq17_height = DESeq(deseq17_height, test="Wald", fitType="parametric")

res17_height = results(deseq17_height, cooksCutoff = FALSE)
alpha = 0.05 # switch to 0.05- generalized linear models
sigtab17_height = res17_height[which(res17_height$padj < alpha), ]
sigtab17_height = cbind(as(sigtab17_height, "data.frame"), as(tax_table(physeq_class17)[rownames(sigtab17_height), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab17_height$log2FoldChange, sigtab17_height$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab17_height$Genus = factor(as.character(sigtab17_height$Genus), levels=names(x))

#Class
x = tapply(sigtab17_height$log2FoldChange, sigtab17_height$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_height$Class = factor(as.character(sigtab17_height$Class), levels=names(x))

ggplot(sigtab17_height, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Length DESeq2- 2017 Data ####
physeq_class17 = subset_samples(physeq_class17, Length_delta != "NA")

deseq17_length = phyloseq_to_deseq2(physeq_class17, ~ Length_delta)
deseq17_length = DESeq(deseq17_length, test="Wald", fitType="parametric")

res17_length = results(deseq17_length, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_length = res17_length[which(res17_length$padj < alpha), ]
sigtab17_length = cbind(as(sigtab17_length, "data.frame"), as(tax_table(physeq_class17)[rownames(sigtab17_length), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab17_length$log2FoldChange, sigtab17_length$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab17_length$Genus = factor(as.character(sigtab17_length$Genus), levels=names(x))

#Class
x = tapply(sigtab17_length$log2FoldChange, sigtab17_length$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_length$Class = factor(as.character(sigtab17_length$Class), levels=names(x))

ggplot(sigtab17_length, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Width DESeq2- 2017 Data ####
physeq_class17 = subset_samples(physeq_class17, Width_delta != "NA")

deseq17_width = phyloseq_to_deseq2(physeq_class17, ~ Width_delta)
deseq17_width = DESeq(deseq17_width, test="Wald", fitType="parametric")

res17_width = results(deseq17_width, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_width = res17_length[which(res17_width$padj < alpha), ]
sigtab17_width = cbind(as(sigtab17_width, "data.frame"), as(tax_table(physeq_class17)[rownames(sigtab17_width), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab17_width$log2FoldChange, sigtab17_width$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab17_width$Genus = factor(as.character(sigtab17_width$Genus), levels=names(x))

#Class
x = tapply(sigtab17_width$log2FoldChange, sigtab17_width$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_width$Class = factor(as.character(sigtab17_width$Class), levels=names(x))

ggplot(sigtab17_width, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Taking out the OTU's ####

# Weight OTU's
otu_weight = sigtab17_weight %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus)
head(otu_weight)
?subset_samples
physeq_class17_weight = subset_taxa(physeq_class17, Class==otu_weight)
physeq_class17_weight = subset_samples(physeq_class17_weight)
physeq_class17_weight
physeq_class17

plot_bar(physeq_class17_weight)
