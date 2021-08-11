# DESeq2 Analysis on Measurments ####
# 2021-07-14
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
require(data.table)
require(RColorBrewer)
require("ggpubr")
require(dplyr)
require(tidyr)
require(DESeq2)
library("ggpubr")
require(genefilter)
library(genefilter)

# Loading data ####

meta17_data <- read.csv("Data/meta17_data_update.csv")
asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

#Loading Physeq w/and w/out transform_sample_counts() function ####
physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

physeq_class17 <- readRDS("Data/physeq_class17.rds")
physeq_class17

# Normalized Weight with DESeq2- 2017 Data ####
physeq_count17 = subset_samples(physeq_count17, Weight_delta != "NA")

deseq17_weight = phyloseq_to_deseq2(physeq_count17, ~ Weight_delta)
deseq17_weight = DESeq(deseq17_weight, test="Wald", fitType="parametric")

res17_weight = results(deseq17_weight, cooksCutoff = FALSE)
alpha = 0.05 # switch to 0.05- generalized linear models
sigtab17_weight = res17_weight[which(res17_weight$padj < alpha), ]
sigtab17_weight = cbind(as(sigtab17_weight, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_weight), ], "matrix"))

weight_otus= row.names(sigtab17_weight)

physeq_count17_weight2 = prune_taxa(physeq_count17, taxa = weight_otus)
physeq_count17_weight2

physeq_count17_weight3=prune_taxa(taxa_sums(physeq_count17_weight2) > 0, physeq_count17_weight2)

View(otu_table(physeq_count17_weight3))
# Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
dim(sigtab17_weight)
# Genus 
x = tapply(sigtab17_weight$log2FoldChange, sigtab17_weight$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17_weight$Phylum = factor(as.character(sigtab17_weight$Phylum), levels=names(x))

# Class
x = tapply(sigtab17_weight$log2FoldChange, sigtab17_weight$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_weight$Class = factor(as.character(sigtab17_weight$Class), levels=names(x))

ggplot(sigtab17_weight, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Weight- Phylum and Class 2017")

ggsave(filename = "Normalized Weight Phylum and Class 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 15, height = 8)  

# Taking out negative and positive OTU's

neg_otus_weight <- sigtab17_weight %>% 
  filter(log2FoldChange < 0)
pos_otus_weight <- sigtab17_weight %>% 
  filter(log2FoldChange > 0)

# Normalized Height with DESeq2- 2017 Data ####
physeq_count17 = subset_samples(physeq_count17, Height_delta != "NA")

deseq17_height = phyloseq_to_deseq2(physeq_count17, ~ Height_delta)
deseq17_height = DESeq(deseq17_height, test="Wald", fitType="parametric")

res17_height = results(deseq17_height, cooksCutoff = FALSE)
alpha = 0.05 # switch to 0.05- generalized linear models
sigtab17_height = res17_height[which(res17_height$padj < alpha), ]
sigtab17_height = cbind(as(sigtab17_height, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_height), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum 
x = tapply(sigtab17_height$log2FoldChange, sigtab17_height$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17_height$Phylum = factor(as.character(sigtab17_height$Phylum), levels=names(x))

# Order
x = tapply(sigtab17_height$log2FoldChange, sigtab17_height$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab17_height$Order = factor(as.character(sigtab17_height$Order), levels=names(x))

ggplot(sigtab17_height, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Height- Phylum and Order 2017")

ggsave(filename = "Normalized Height Phylum and Order 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 30, height = 8)  

# Taking out negative and positive OTU's

neg_otus_height <- sigtab17_height %>% 
  filter(log2FoldChange < 0)
pos_otus_height <- sigtab17_height %>% 
  filter(log2FoldChange > 0)

# Normalized Length DESeq2- 2017 Data ####
physeq_count17 = subset_samples(physeq_count17, Length_delta != "NA")

deseq17_length = phyloseq_to_deseq2(physeq_count17, ~ Length_delta)
deseq17_length = DESeq(deseq17_length, test="Wald", fitType="parametric")

res17_length = results(deseq17_length, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_length = res17_length[which(res17_length$padj < alpha), ]
sigtab17_length = cbind(as(sigtab17_length, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_length), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Genus
x = tapply(sigtab17_length$log2FoldChange, sigtab17_length$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17_length$Phylum = factor(as.character(sigtab17_length$Phylum), levels=names(x))

#Class
x = tapply(sigtab17_length$log2FoldChange, sigtab17_length$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_length$Class = factor(as.character(sigtab17_length$Class), levels=names(x))

ggplot(sigtab17_length, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Length- Phylum and Class 2017")

ggsave(filename = "Normalized Length Phylum and Class 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 17, height = 8)  
 
# Taking out negative and positive OTU's

neg_otus_length <- sigtab17_length %>% 
  filter(log2FoldChange < 0)
pos_otus_length <- sigtab17_length %>% 
  filter(log2FoldChange > 0)

# Normalized Width DESeq2- 2017 Data ####
physeq_count17 = subset_samples(physeq_count17, Width_delta != "NA")

deseq17_width = phyloseq_to_deseq2(physeq_count17, ~ Width_delta)
deseq17_width = DESeq(deseq17_width, test="Wald", fitType="parametric")

res17_width = results(deseq17_width, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_width = res17_width[which(res17_width$padj < alpha), ]
sigtab17_width = cbind(as(sigtab17_width, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_width), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum
x = tapply(sigtab17_width$log2FoldChange, sigtab17_width$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17_width$Phylum = factor(as.character(sigtab17_width$Phylum), levels=names(x))

# Class
x = tapply(sigtab17_width$log2FoldChange, sigtab17_width$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17_width$Class = factor(as.character(sigtab17_width$Class), levels=names(x))

ggplot(sigtab17_width, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Width- Phylum and Class 2017")

ggsave(filename = "Normalized Width Phylum and Class 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 20, height = 8)  

# Taking out negative and positive OTU's

neg_otus_width <- sigtab17_width %>% 
  filter(log2FoldChange < 0)
pos_otus_width <- sigtab17_width %>% 
  filter(log2FoldChange > 0)

# Merging the positive otus ####

pos_otus_17 = merge(pos_otus_height, pos_otus_length, by= "row.names")
rownames(pos_otus_17) = pos_otus_17$Row.names
pos_otus_17 = merge(pos_otus_17, pos_otus_width, by = "row.names")
pos_otus_17 <- subset(pos_otus_17, select = -c(Row.names, Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x, 
                                                     Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y))
rownames(pos_otus_17) = pos_otus_17$Row.names
pos_otus_17 <- subset(pos_otus_17, select = -c(baseMean, baseMean.x, baseMean.y, lfcSE, lfcSE.x, lfcSE.y, stat, stat.x, stat.y, pvalue, 
                                               pvalue.x, pvalue.y, padj, padj.x, padj.y))
pos_otus_17 <- subset(pos_otus_17, select = -c(log2FoldChange, log2FoldChange.x, log2FoldChange.y))
pos_otus_17 <- subset(pos_otus_17, select = -c(Row.names))
pos_otus_17 = merge(pos_otus_17, pos_otus_weight, by= "row.names")
pos_otus_17 <- subset(pos_otus_17, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x, 
                                               baseMean, lfcSE, stat, pvalue, padj))
pos_otus_17 <- subset(pos_otus_17, select= -c(log2FoldChange))

colnames(pos_otus_17) <- c("Row.names","Kingdom", "Phylum", "Class", "Order", 
                              "Family", "Genus.x","Genus.y", "Species")
# Tax tree
meta17_data <- read.csv("Data/meta17_data_update.csv")

asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))

#Changing row names in taxa
rownames(pos_otus_17)= pos_otus_17$Row.names
pos_otus_17$Row.names = NULL

#Setting taxmat and otumat
taxmat17_weight=pos_otus_17
otumat17_weight=asvtable_17

#Converting to matrix
otu_matrix17_weight= as.matrix(otumat17_weight, rownames = rownames(asvtable_17))

tax_matrix17_weight=as.matrix(taxmat17_weight)
View(OTU17_weight)

meta17_data_weight=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_weight= otu_table(otu_matrix17_weight, taxa_are_rows = FALSE)

TAX17_weight= tax_table(tax_matrix17_weight)

SAMP17_weight= sample_data(meta17_data)

OTU_count17_weight=transform_sample_counts(OTU17_weight, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_posotu = phyloseq(OTU17_weight, TAX17_weight, SAMP17_weight)
physeq_class17_posotu

physeq_count17_posotu = phyloseq(OTU_count17_weight, TAX17_weight, SAMP17_weight)
physeq_count17_posotu

tax_otu17 = parse_phyloseq(physeq_count17_posotu)
tax_otu17%>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            initial_layout = "reingold-tilford",layout = "davidson-harel",
            title = "Measurement Taxa 2017: Positive OTUs",
            node_color_axis_label = "Number of OTUs")
ggsave(filename = "Heat Tree for Measurements 2017_Positive OTUs.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5) 

# Merging Negative otus ####

neg_otus_17 = merge(neg_otus_height, neg_otus_width, by= "row.names")
neg_otus_17 <- subset(neg_otus_17, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x,
                                              baseMean.x, baseMean.y, lfcSE.x, lfcSE.y, stat.x, stat.y, 
                                               pvalue.x, pvalue.y, padj.x, padj.y, log2FoldChange.x, log2FoldChange.y))
rownames(neg_otus_17)= neg_otus_17$Row.names
neg_otus_17= merge(neg_otus_17, neg_otus_length, by = "row.names")
neg_otus_17<- subset(neg_otus_17, select = -c(Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y, 
                                              Row.names, baseMean, lfcSE, stat, pvalue, padj, log2FoldChange))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))

#Changing row names in taxa
rownames(neg_otus_17)= neg_otus_17$Row.names
neg_otus_17$Row.names = NULL

#Setting taxmat and otumat
taxmat17=neg_otus_17
otumat17=asvtable_17

#Converting to matrix
otu_matrix17= as.matrix(otumat17, rownames = rownames(asvtable_17))

tax_matrix17=as.matrix(taxmat17)


meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= otu_table(otu_matrix17, taxa_are_rows = FALSE)

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)

OTU_count17=transform_sample_counts(OTU17, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_negotu = phyloseq(OTU17, TAX17, SAMP17)
physeq_class17_negotu

physeq_count17_negotu = phyloseq(OTU_count17, TAX17, SAMP17)
physeq_count17_negotu

set.seed(4)
tax_otu17_neg = parse_phyloseq(physeq_count17_negotu)
tax_otu17_neg%>% 
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs, 
            initial_layout = "reingold-tilford",layout = "davidson-harel",
            title = "Measurement Taxa 2017: Negative OTUs",
            node_color_axis_label = "Number of OTUs")
tax_otu17_neg
ggsave(filename = "Heat Tree for Measurements 2017_Negative OTUs.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5) 

# Taking out the significant OTU's ####

# Weight

otu_weight = sigtab17_weight %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Height

otu_height = sigtab17_height %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Length

otu_length= sigtab17_length %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)
# Width

otu_width = sigtab17_width %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

