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

# change 1= mean across rows 
res18_weight = results(deseq18_weight, cooksCutoff = FALSE)
alpha = 0.05
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

# Order
x = tapply(sigtab18_weight$log2FoldChange, sigtab18_weight$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab18_weight$Order = factor(as.character(sigtab18_weight$Order), levels=names(x))

ggplot(sigtab18_weight, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Weight- Phylum and Order 2018")

ggsave(filename = "Normalized Weight Phylum and Order 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 25, height = 8)  

# Taking out negative and positive OTU's

neg_otus_weight18 <- sigtab18_weight %>% 
  filter(log2FoldChange < 0)
pos_otus_weight18 <- sigtab18_weight %>% 
  filter(log2FoldChange > 0)

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

# Order
x = tapply(sigtab18_height$log2FoldChange, sigtab18_height$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab18_height$Family = factor(as.character(sigtab18_height$Family), levels=names(x))

ggplot(sigtab18_height, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Height- Phylum and Family 2018")

ggsave(filename = "Normalized Height Phylum and Family 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 30, height = 8)  

# Taking out negative and positive OTU's

neg_otus_height18 <- sigtab18_height %>% 
  filter(log2FoldChange < 0)
pos_otus_height18 <- sigtab18_height %>% 
  filter(log2FoldChange > 0)

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
# Phylum 
x = tapply(sigtab18_length$log2FoldChange, sigtab18_length$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_length$Phylum = factor(as.character(sigtab18_length$Phylum), levels=names(x))

# Family
x = tapply(sigtab18_length$log2FoldChange, sigtab18_length$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab18_length$Family = factor(as.character(sigtab18_length$Family), levels=names(x))

ggplot(sigtab18_length, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Length- Phylum and Family 2018")

ggsave(filename = "Normalized Length Phylum and Family 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 32, height = 8)  

# Taking out negative and positive OTU's

neg_otus_length18 <- sigtab18_length %>% 
  filter(log2FoldChange < 0)
pos_otus_length18 <- sigtab18_length %>% 
  filter(log2FoldChange > 0)

# Normalized Width DESeq2- 2017 Data ####
physeq_count18 = subset_samples(physeq_count18, Width_delta != "NA")

deseq18_width = phyloseq_to_deseq2(physeq_count18, ~ Width_delta)
gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_width = estimateSizeFactors(deseq18_width, geoMeans=geoMeans, locfunc=shorth)
deseq18_width = DESeq(deseq18_width, test="Wald", fitType="parametric")

res18_width = results(deseq18_width, cooksCutoff = FALSE)
alpha = 0.05 
sigtab18_width = res18_width[which(res18_width$padj < alpha), ]
sigtab18_width = cbind(as(sigtab18_width, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_width), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum
x = tapply(sigtab18_width$log2FoldChange, sigtab18_width$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_width$Phylum = factor(as.character(sigtab18_width$Phylum), levels=names(x))

# Order
x = tapply(sigtab18_width$log2FoldChange, sigtab18_width$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab18_width$Order = factor(as.character(sigtab18_width$Order), levels=names(x))

ggplot(sigtab18_width, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Width- Phylum and Order 2018")

ggsave(filename = "Normalized Width Phylum and Order 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 30, height = 8)  

# Taking out negative and positive OTU's

neg_otus_width18 <- sigtab18_width %>% 
  filter(log2FoldChange < 0)
pos_otus_width18 <- sigtab18_width %>% 
  filter(log2FoldChange > 0)

# Merging Negative OTUs

neg_otus18 = merge(neg_otus_height18, neg_otus_length18, by= "row.names")
neg_otus18 <- subset(neg_otus18, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x,
                                               baseMean.x, baseMean.y, lfcSE.x, lfcSE.y, stat.x, stat.y, 
                                               pvalue.x, pvalue.y, padj.x, padj.y, log2FoldChange.x, log2FoldChange.y))
rownames(neg_otus18)= neg_otus18$Row.names
neg_otus18 = merge(neg_otus18, neg_otus_weight18, by="row.names")
neg_otus18 <- subset(neg_otus18, select = -c(Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y,
                                             baseMean, lfcSE, stat, 
                                             pvalue, padj, log2FoldChange))
rownames(neg_otus18)= neg_otus18$Row.names
neg_otus18 <- subset(neg_otus18, select = -c(Row.names, Row.names.1))
neg_otus18 = merge(neg_otus18, neg_otus_width18, by="row.names")
neg_otus18 <- subset(neg_otus18, select = -c(Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y,
                                             baseMean, lfcSE, stat, 
                                             pvalue, padj, log2FoldChange))
colnames(neg_otus18) <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", 
                              "Family", "Genus.x","Genus.y", "Species")
# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")

#Changing rownames of new taxa table
rownames(neg_otus18)= neg_otus18$Row.names
neg_otus18$Row.names = NULL

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
taxmat18= neg_otus18
otumat18=asvtable_18

#Converting to matrix
otu_matrix18= as.matrix(otumat18, rownames = "V1")

tax_matrix18=as.matrix(taxmat18)

meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18= otu_table(otu_matrix18, taxa_are_rows = FALSE)

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

OTU_count18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))

physeq_class18_negotu = phyloseq(OTU18, TAX18, SAMP18)
physeq_class18_negotu

physeq_count18_negotu = phyloseq(OTU_count18, TAX18, SAMP18)
physeq_count18_negotu

tax_negotu18 = parse_phyloseq(physeq_count18_negotu)

set.seed(4)
tax_negotu18 %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs,
            initial_layout = "reingold-tilford",layout = "davidson-harel",
            title = "Measurement Taxa: Negative OTUs 2018",
            node_color_axis_label = "Number of OTUs")

ggsave(filename = "Heat Tree for Measurements 2018_Negative OTUs.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5) 

# Merging Positive OTUs

pos_otus18 = merge(pos_otus_width18, pos_otus_length18, by= "row.names")
pos_otus18 <- subset(pos_otus18, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x,
                                             baseMean.x, baseMean.y, lfcSE.x, lfcSE.y, stat.x, stat.y, 
                                             pvalue.x, pvalue.y, padj.x, padj.y, log2FoldChange.x, log2FoldChange.y))
rownames(pos_otus18)= pos_otus18$Row.names
pos_otus18 = merge(pos_otus18, pos_otus_height18, by="row.names")
pos_otus18 <- subset(pos_otus18, select = -c(Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y,
                                             baseMean, lfcSE, stat, 
                                             pvalue, padj, log2FoldChange,Row.names))
rownames(pos_otus18)= pos_otus18$Row.names
pos_otus18 <- subset(pos_otus18, select = -c(Row.names))
pos_otus18 = merge(pos_otus18, pos_otus_weight18, by="row.names")
pos_otus18 <- subset(pos_otus18, select = -c(Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y,
                                             baseMean, lfcSE, stat, 
                                             pvalue, padj, log2FoldChange))
colnames(pos_otus18) <- c("Row.names", "Kingdom", "Phylum", "Class", "Order", 
                          "Family", "Genus.x","Genus.y", "Species")
# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")

#Changing rownames of new taxa table
rownames(pos_otus18)= pos_otus18$Row.names
pos_otus18$Row.names = NULL

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
taxmat18= pos_otus18
otumat18=asvtable_18

#Converting to matrix
otu_matrix18= as.matrix(otumat18, rownames = "V1")

tax_matrix18=as.matrix(taxmat18)

meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18= otu_table(otu_matrix18, taxa_are_rows = FALSE)

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

OTU_count18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))

physeq_class18_posotu = phyloseq(OTU18, TAX18, SAMP18)
physeq_class18_posotu

physeq_count18_posotu = phyloseq(OTU_count18, TAX18, SAMP18)
physeq_count18_posotu

tax_posotu18 = parse_phyloseq(physeq_count18_posotu)

set.seed(4)
tax_posotu18 %>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs,
            initial_layout = "reingold-tilford",layout = "davidson-harel",
            title = "Measurement Taxa: Positive OTUs 2018",
            node_color_axis_label = "Number of OTUs")

ggsave(filename = "Heat Tree for Measurements 2018_Positive OTUs.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5) 

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

