# Phyloseq and Deseq Analysis- Measurements ####
# 2022-03-01
# Author: Monse Garcia



library(DESeq2)
require(RColorBrewer)
require(ggplot2)
require(data.table)
require("ggpubr")
require(dplyr)
require(tidyr)
library("ggpubr")
# Loading Data
physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

physeq_class17 <- readRDS("Data/physeq_class17.rds")
physeq_class17

meta17_data <- read.csv("Data/meta17_data_update.csv")

asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

# Normalized Volume with DESEQ2- 2017 Data
physeq_count17 = subset_samples(physeq_count17, Volume_delta != "NA")

deseq17_volume = phyloseq_to_deseq2(physeq_count17, ~ Volume_delta) #the design formula contains one or more numeric variables that have mean or standard deviation larger than 5 (an arbitrary threshold to trigger this message). it is generally a good idea to center and scale numeric variables in the design to improve GLM convergence.
deseq17_volume = DESeq(deseq17_volume, test="Wald", fitType="parametric")

res17_volume = results(deseq17_volume, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_volume = res17_volume[which(res17_volume$padj < alpha), ]
sigtab17_volume = cbind(as(sigtab17_volume, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_volume), ], "matrix"))

#Theme for Graph
theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum 
x = tapply(sigtab17_volume$log2FoldChange, sigtab17_volume$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17_volume$Phylum = factor(as.character(sigtab17_volume$Phylum), levels=names(x))

# Order
x = tapply(sigtab17_volume$log2FoldChange, sigtab17_volume$Order, function(x) max(x))
x = sort(x, TRUE)
sigtab17_volume$Order = factor(as.character(sigtab17_volume$Order), levels=names(x))

ggplot(sigtab17_volume, aes(x=Order, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Volume- Phylum and Order 2017")


# Taking out negative and positive OTU's

neg_otus_volume <- sigtab17_volume %>% 
  filter(log2FoldChange < 0)
pos_otus_volume <- sigtab17_volume %>% 
  filter(log2FoldChange > 0)

pos_otus_volume2 <- subset(pos_otus_volume, select = -c(baseMean, log2FoldChange, 
                                                        lfcSE, stat, pvalue, padj))
write.csv(pos_otus_volume2, file = "Data/pos_otus_volume.csv")
pos_otus17 <- read.csv("Data/pos_otus_volume.csv")


# Phyloseq ####

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$X.2=NULL
meta17_data$X=NULL
meta17_data$X.1=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))

#Changing row names in taxa
rownames(pos_otus17)= pos_otus17$X
pos_otus17$X = NULL

#Setting taxmat and otumat
taxmat17_volume=pos_otus17
otumat17_volume=asvtable_17

#Converting to matrix
otu_matrix17_volume= as.matrix(otumat17_volume, rownames = rownames(asvtable_17))

tax_matrix17_volume=as.matrix(taxmat17_volume)
View(OTU17_volume)

meta17_data_volume=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_volume= otu_table(otu_matrix17_volume, taxa_are_rows = FALSE)

TAX17_volume= tax_table(tax_matrix17_volume)

SAMP17_volume= sample_data(meta17_data)

OTU_count17_volume=transform_sample_counts(OTU17_volume, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_posotuvol = phyloseq(OTU17_volume, TAX17_volume, SAMP17_volume)
physeq_class17_posotuvol

physeq_count17_posotuvol = phyloseq(OTU_count17_volume, TAX17_volume, SAMP17_volume)
physeq_count17_posotuvol

#Saving pos otu in physeq class 2017 
saveRDS(physeq_count17_posotuvol, "Data/physeq_count17_posotuvol.rds")
saveRDS(physeq_class17_posotuvol, "Data/physeq_class17_posotuvol.rds")


physeq_class18 <- readRDS("Data/physeq_class18.rds")
physeq_count18 <- readRDS("Data/physeq_count18.rds")

# Normalized Weight with DESeq2- 2018 Data
physeq_count18 = subset_samples(physeq_count18, Volume_delta != "NA")

deseq18_volume = phyloseq_to_deseq2(physeq_count18, ~ Volume_delta)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18_volume = estimateSizeFactors(deseq18_volume, geoMeans=geoMeans, locfunc=shorth)
deseq18_volume = DESeq(deseq18_volume, test="Wald", fitType="parametric")

# change 1= mean across rows 
res18_volume = results(deseq18_volume, cooksCutoff = FALSE)
alpha = 0.05
sigtab18_volume = res18_volume[which(res18_volume$padj < alpha), ]
sigtab18_volume = cbind(as(sigtab18_volume, "data.frame"), as(tax_table(physeq_count18)[rownames(sigtab18_volume), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
# Phylum 
x = tapply(sigtab18_volume$log2FoldChange, sigtab18_volume$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_volume$Phylum = factor(as.character(sigtab18_volume$Phylum), levels=names(x))

# Order
x = tapply(sigtab18_volume$log2FoldChange, sigtab18_volume$Family, function(x) max(x))
x = sort(x, TRUE)
sigtab18_volume$Family = factor(as.character(sigtab18_volume$Family), levels=names(x))

ggplot(sigtab18_volume, aes(x=Family, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "Normalized Volume- Phylum and Family 2018")

# Taking out negative and positive OTU's

neg_otus_volume18 <- sigtab18_volume %>% 
  filter(log2FoldChange < 0)
pos_otus_volume18 <- sigtab18_volume %>% 
  filter(log2FoldChange > 0)

pos_otus_volume18 <- subset(pos_otus_volume18, select = -c(baseMean, log2FoldChange, 
                                                        lfcSE, stat, pvalue, padj))
write.csv(pos_otus_volume18, file = "Data/pos_otus_volume18.csv")
pos_otus_volume18 <- read.csv("Data/pos_otus_volume.csv")
write.csv(neg_otus_volume18, file = "Data/neg_otus_volume18.csv")

# Phyloseq Volume-2018 ####

# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")

#Changing rownames of new taxa table
rownames(pos_otus_volume18)= pos_otus_volume18$X
pos_otus_volume18$X = NULL

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
taxmat18= pos_otus_volume18
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

physeq_class18_posotuvol = phyloseq(OTU18, TAX18, SAMP18)
physeq_class18_posotuvol

physeq_count18_posotuvol = phyloseq(OTU_count18, TAX18, SAMP18)
physeq_count18_posotuvol

#Saving pos otu in physeq class 2018 
saveRDS(physeq_count18_posotuvol, "Data/physeq_count18_posotuvol.rds")
saveRDS(physeq_class18_posotuvol, "Data/physeq_class18_posotuvol.rds")



