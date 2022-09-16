# Measurement Analysis- Results ####
# 2022-06-29
# Author: Monse Garcia


# Installing Packages 
install.packages("devtools")
devtools::install_github("david-barnett/microViz@0.9.2")
install.packages("remotes")

# Loading packages
library(microViz)
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
library(vegan)
library(remotes)
library(DESeq2)
library(genefilter)
find.package("devtools")

##### 1.Phyloseq Analysis with scaled volume #####

physeq_highsc17 <- readRDS("Data/physeq_highsc17.rds")

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

#Set 0 to 1, have 1 be "white" instead of NA values

#### Positive OTUs

significant_highsc17_pos <- significant_highsc17[significant_highsc17$log2FoldChange>0,]
dim(significant_highsc17_pos)
# 25, 6

sub_significant_highsc17_pos <- subset_taxa(prune_taxa(rownames(significant_highsc17_pos), physeq_highsc17))
df_significant_highsc17_pos <- as.data.frame(tax_table(sub_significant_highsc17_pos))
write.table(df_significant_highsc17_pos, file = "Data/Significant Positive OTU's for Scaled Volume 2017.csv", quote = FALSE, sep = ",", col.names = T)

otu <- as.matrix(otu_table(sub_significant_highsc17_pos))
OTU2<- otu+1 #turn into phyloseq
View(OTU2)
# Turning into Phyloseq Object

meta17_data <- read.csv("Data/meta17_data_update.csv")

Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run23_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")
meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= OTU2

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)


OTU_count17=transform_sample_counts(OTU17, function(x) 1E6 * x/sum(x))


physeq_hi17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_hi17

physeq_hict17 = phyloseq(OTU_count17, TAX17, SAMP17)
physeq_hict17


# Saving Physeq as an RDS
saveRDS(physeq_class17, "Data/physeq_class17.rds")
physeq_class17 <- readRDS("Data/physeq_class17.rds")

saveRDS(physeq_count17, "Data/physeq_count17.rds")
physeq_count17 <- readRDS("Data/physeq_count17.rds")


plot_heatmap(physeq_hict17, method = "NMDS", distance = "bray", low = "#000033", high ="#FFFFCC"  , na.value = "white", taxa.label = "Family", sample.label = "Volume_scale", sample.order = "Volume_scale", title = "Heat Map of Postive OTUs Associated with Scaled Volume-2017" )


#### Negative 
significant_highsc17_neg <- significant_highsc17[significant_highsc17$log2FoldChange<0, ]
dim(significant_highsc17_neg)
# 42, 6

# 2018 
physeq_highsc18 <- readRDS("Data/physeq_highsc18.rds")

deseq17 = DESeq(deseq17, test = "Wald", fitType = "parametric")
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

### order by size = smallest to largest, filter out to see where they appear in a certain amount of samples (e.g. found in half of sample), if not able to work go back to measurements individually 
# https://david-barnett.github.io/microViz/reference/tax_filter.html 
# http://joey711.github.io/phyloseq/plot_heatmap-examples

# run non-scaled ones to see the difference between them (remaking them and take out OTUs only found in 1 sample)

#### Positive OTUs

significant_highsc18_pos <- significant_highsc18[significant_highsc18$log2FoldChange>0,]
dim(significant_highsc18_pos)
# 13, 6

sub_significant_highsc18_pos <- subset_taxa(prune_taxa(rownames(significant_highsc18_pos), physeq_highsc18))
df_significant_highsc18_pos <- as.data.frame(tax_table(sub_significant_highsc18_pos))
write.table(df_significant_highsc18_pos, file = "Data/Significant Positive OTU's for Scaled Volume 2018.csv", quote = FALSE, sep = ",", col.names = T)

plot_heatmap(sigtab17, method = "NMDS", distance ="bray", low = "#FFFFCC", high = "#000033", na.value = "white", taxa.label = "Order", sample.label = "UniqueID", sample.order = "Volume_scale", title = "Heat Map of Postive OTUs Associated with Scaled Volume") 

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


#### 3.Phyloseq Analysis with Volume Delta

physeq_class17 <- readRDS("Data/physeq_class17.rds")

otu <- as.data.frame(is.na(otu_table(sub_significant_highsc17_pos)))

des_vol17 <- phyloseq_to_deseq2(physeq_class17, ~ Volume_delta)
des_vol17 <- DESeq(des_vol17, test="Wald", fitType = "parametric")
resultsNames(des_vol17)
# Intercept , volume_scale

results_vol17 <- results(des_vol17, name = "Volume_delta")
significant_vol17 <- results_vol17[which(results_vol17$padj <0.05), ]
dim(significant_vol17)
# 51, 6

sub_significant_vol17 <- subset_taxa(prune_taxa(rownames(significant_vol17), physeq_class17))
df_significant_vol17 <- as.data.frame(tax_table(sub_significant_vol17))
write.table(df_significant_vol17, file = "Data/Significant OTU's for Volume Delta 2017.csv", quote = FALSE, sep = ",", col.names = T)

#phyloseq_rm_na_tax()
plot_heatmap(physeq_count17_posotuvol, method = "NMDS", distance ="bray", low = "#FFFFCC", high = "#000033", na.value = "white", taxa.label = "Order", sample.label = "UniqueID", sample.order = "Volume_delta", title = "Heat Map of Postive OTUs Associated with Volume") 
