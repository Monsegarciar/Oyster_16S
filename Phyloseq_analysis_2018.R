# Phyloseq Analysis 2018 Data ####
# 2021-06-27
# Author: Monse Garcia

# Loading packages 
require(phyloseq)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)

# Loading Data

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")
Run23_taxa <- fread("Data/Run23_taxa - Copy.csv")

#Changing row names in "Run23_taxa"
rownames(Run23_taxa)= Run23_taxa$V1
head(rownames(Run23_taxa))

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
taxmat18=Run23_taxa
otumat18=asvtable_18

#Converting to matrix
otu_matrix18= as.matrix(otumat18, rownames = "V1")

tax_matrix18=as.matrix(taxmat18, rownames = "V1")
colnames(tax_matrix18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                          "Genus")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18= otu_table(otu_matrix18, taxa_are_rows = FALSE)

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

OTU_count18=transform_sample_counts(OTU18, function(x) 1E6 * x/sum(x))

physeq_class18 = phyloseq(OTU18, TAX18, SAMP18)
physeq_class18

physeq_count18 = phyloseq(OTU_count18, TAX18, SAMP18)
physeq_count18

# Saving Physeq as an RDS
saveRDS(physeq_class18, "Data/physeq_class18.rds")
physeq_class18 <- readRDS("Data/physeq_class18.rds")

saveRDS(physeq_count18, "Data/physeq_count18.rds")
physeq_count18 <- readRDS("Data/physeq_count18.rds")

# Graph Analysis ####
# Ordination
Phy.ord18 <- ordinate(physeq_class18, "NMDS", "bray")
plot4= plot_ordination(physeq_class18, Phy.ord18, type = "split", 
                    color = "Phylum", shape = "Bucket2", title = "Plot Ordination: Phylum and Treatments 2018")
print(plot4)

# Alpha Diveristy####
plo_rich18= plot_richness(physeq_class18, x="Bucket2", measures=c("Chao1", "Shannon"))
print(plo_rich18)

plo_rich_species18= plot_richness(physeq_class18, x="Bucket2", color = "Species2.x",measures=c("Simpson", "Shannon"), title = "Alpha Diveristy for Treatments and Species 2018")
print(plo_rich_species18)
 






