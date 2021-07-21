# DESeq2 Analysis on Measurments ####
# 2021-07-14
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
require(genefilter)
library(genefilter)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)

install.packages("genefilter")
?genefilter
# Loading data ####

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")


#Loading Physeq w/out transform_sample_counts() function ####
physeq_class18 <- readRDS("Data/physeq_class18.rds")
physeq_count18 <- readRDS("Data/physeq_count18.rds")


# Normalized Weight with DESeq2- 2018 Data
physeq_count18 = subset_samples(physeq_count18, Weight_delta != "NA")

deseq18_weight = phyloseq_to_deseq2(physeq_count18, ~ Weight_delta)
 
gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_weight = estimateSizeFactors(deseq18_weight, geoMeans=geoMeans, locfunc=shorth)
deseq18_weight = DESeq(deseq18_weight, test="Wald", fitType="parametric")

nrow(deseq18_weight)
nrow(OTU_count18)
ncol(OTU_count18)
View(geoMeans)

# change 1= mean across rows 
res18_weight = results(deseq18_weight, cooksCutoff = FALSE)
alpha = 0.01
sigtab18_weight = res18_weight[which(res18_weight$padj < alpha), ]
sigtab18_weight = cbind(as(sigtab18_weight, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_weight), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

# Phylum
x = tapply(sigtab18_weight$log2FoldChange, sigtab18_weight$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_weight$Phylum = factor(as.character(sigtab18_weight$Phylum), levels=names(x))

# Genus
x = tapply(sigtab18_weight$log2FoldChange, sigtab18_weight$Genus.x, function(x) max(x))
x = sort(x, TRUE)
sigtab18_weight$Genus.x = factor(as.character(sigtab18_weight$Genus.x), levels=names(x))

ggplot(sigtab18_weight, aes(x=Genus.x, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Height with DESeq2- 2018 Data ####
physeq_count18 = subset_samples(physeq_count18, Height_delta != "NA")

deseq18_height = phyloseq_to_deseq2(physeq_count18, ~ Height_delta)
gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_height = estimateSizeFactors(deseq18_height, geoMeans=geoMeans, locfunc=shorth)
deseq18_height = DESeq(deseq18_height, test="Wald", fitType="parametric")

res18_height = results(deseq18_height, cooksCutoff = FALSE)
alpha = 0.05 
sigtab18_height = res18_height[which(res18_height$padj < alpha), ]
sigtab18_height = cbind(as(sigtab18_height, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_height), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum 
x = tapply(sigtab18_height$log2FoldChange, sigtab18_height$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_height$Phylum = factor(as.character(sigtab18_height$Phylum), levels=names(x))

# Genus
x = tapply(sigtab18_height$log2FoldChange, sigtab18_height$Genus.x, function(x) max(x))
x = sort(x, TRUE)
sigtab18_height$Genus.x = factor(as.character(sigtab18_height$Genus.x), levels=names(x))

ggplot(sigtab18_height, aes(x=Genus.x, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Length DESeq2- 2018 Data ####
physeq_count18 = subset_samples(physeq_count18, Length_delta != "NA")

deseq18_length = phyloseq_to_deseq2(physeq_count18, ~ Length_delta)
gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_length = estimateSizeFactors(deseq18_length, geoMeans=geoMeans, locfunc=shorth)
deseq18_length = DESeq(deseq18_length, test="Wald", fitType="parametric")

res18_length = results(deseq18_length, cooksCutoff = FALSE)
alpha = 0.05 
sigtab18_length = res18_length[which(res18_length$padj < alpha), ]
sigtab18_length = cbind(as(sigtab18_length, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_length), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab18_length$log2FoldChange, sigtab18_length$Genus.x, function(x) max(x))
x = sort(x, TRUE)
sigtab18_length$Genus.x = factor(as.character(sigtab18_length$Genus.x), levels=names(x))

#Class
x = tapply(sigtab18_length$log2FoldChange, sigtab18_length$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab18_length$Class = factor(as.character(sigtab18_length$Class), levels=names(x))

ggplot(sigtab18_length, aes(x=Genus.x, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Normalized Width DESeq2- 2017 Data ####
physeq_count18 = subset_samples(physeq_count18, Width_delta != "NA")

deseq18_width = phyloseq_to_deseq2(physeq_count18, ~ Width_delta)
gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_width = estimateSizeFactors(deseq18_width, geoMeans=geoMeans, locfunc=shorth)
deseq18_width = DESeq(deseq18_width, test="Wald", fitType="parametric")

res18_width = results(deseq18_width, cooksCutoff = FALSE)
alpha = 0.05 
sigtab18_width = res18_length[which(res18_width$padj < alpha), ]
sigtab18_width = cbind(as(sigtab18_width, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_width), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab18_width$log2FoldChange, sigtab18_width$Genus.x, function(x) max(x))
x = sort(x, TRUE)
sigtab18_width$Genus.x = factor(as.character(sigtab18_width$Genus.x), levels=names(x))

#Class
x = tapply(sigtab18_width$log2FoldChange, sigtab18_width$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab18_width$Class = factor(as.character(sigtab18_width$Class), levels=names(x))

ggplot(sigtab18_width, aes(x=Genus.x, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# Taking significant OTU's ####

# Weight
otu_weight18 = sigtab18_weight %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Height
otu_height18 = sigtab18_height %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)
  
# Length
otu_length18 = sigtab18_length %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Width
otu_width18 = sigtab18_width %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

