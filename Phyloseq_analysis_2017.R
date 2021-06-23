# Phyloseq Analysis 2017 Data ####
# 2021-06-21
# Monse Garcia

# Required Packages
require(phyloseq)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)

#Loading Data

meta17_data2 <- read.csv("Data/meta17_data2.csv")


asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")
Run23_taxa <- fread("Data/Run23_taxa - Copy.csv")

#Changing row names in "Run23_taxa"
rownames(Run23_taxa)= Run23_taxa$V1
Run23_taxa

#Renaming column names in "Run23_taxa"
Run23_taxa2 <- 

#Changing row names in "meta_17" data
rownames(meta17_data2)= meta17_data2$UniqueID
meta17_data2$UniqueID=NULL
meta17_data2$X=NULL
rownames(meta17_data2)

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
rownames(asvtable_17)


#Setting taxmat and otumat
taxmat=Run23_taxa
otumat=asvtable_17

?data.matrix
?as.matrix

#Converting to matrix
otu_matrix= as.matrix(otumat, rownames = "V1")

tax_matrix=as.matrix(taxmat, rownames = "V1")
colnames(tax_matrix) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                         "Genus")
meta17_data2=as.data.frame(meta17_data2)

#Setting OTU, TAX, and SAMP
OTU= otu_table(otu_matrix, taxa_are_rows = FALSE)

TAX= tax_table(tax_matrix)

SAMP= sample_data(meta17_data2)


OTU=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))

physeq_class = phyloseq(OTU, TAX, SAMP)
physeq_class


taxa_names(TAX)

# "phyloseq_class" Analysis

ntaxa(physeq_class)
  #8007 taxa

#NMDS Graph
Phy.ord <- ordinate(physeq_class, "NMDS", "bray")
p1= plot_ordination(physeq_class, Phy.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

p1 + facet_wrap(~Phylum, 6)


p3= plot_ordination(physeq_class, Phy.ord, type = "biplot", color = "Phylum", shape = "Kingdom", title = "biplot")
print(p3)

the_plot <- plot_bar(physeq_class, fill = "Family")
print(the_plot)







