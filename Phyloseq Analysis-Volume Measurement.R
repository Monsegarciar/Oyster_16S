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
meta17_data$UniqueID=NULL
meta17_data$X=NULL
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

