# Phyloseq Analysis 2017 Data ####
# 2021-06-21
# Author: Monse Garcia

# Required Packages ####
require(phyloseq)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)

#Loading Data ####

meta17_data <- read.csv("Data/meta17_data2.csv")


asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")
Run23_taxa <- fread("Data/Run23_taxa - Copy.csv")

#Changing row names in "Run23_taxa"
rownames(Run23_taxa)= Run23_taxa$V1
head(rownames(Run23_taxa))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))


#Setting taxmat and otumat
taxmat17=Run23_taxa
otumat17=asvtable_17

#Converting to matrix
otu_matrix17= as.matrix(otumat17, rownames = rownames(asvtable_17))

tax_matrix17=as.matrix(taxmat17, rownames = "V1")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                         "Genus")
meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= otu_table(otu_matrix17, taxa_are_rows = FALSE)

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)


OTU_count17=transform_sample_counts(OTU, function(x) 1E6 * x/sum(x))


physeq_class17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_class17

# Saving Physeq as an RDS
saveRDS(physeq_class17, "Data/physeq_class17.rds")
physeq_class17 <- readRDS("Data/physeq_class17.rds")

# "phyloseq_class" Analysis ####

#NMDS Graph ####

#Just OTU's
Phy.ord <- ordinate(physeq_class17, "NMDS", "bray")
p1= plot_ordination(physeq_class17, Phy.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

p1 + facet_wrap(~Phylum, 6)

plot_1= plot_ordination(physeq_class17, Phy.ord, type= )

#Samples with Peacrabs and Site
p2= plot_ordination(physeq_class17, Phy.ord, type = "sample", color = "peacrabs.x", 
                    shape = "Site.x")
print(p2)

p2+ geom_polygon(aes(fill=Weight_pre)) + geom_point(size=7) + ggtitle("samples")

# Biplot
p3= plot_ordination(physeq_class17, Phy.ord, type = "biplot", color = "Treatment2", shape = "Kingdom", title = "biplot")
print(p3)

#Split Graph
p4= plot_ordination(physeq_class17, Phy.ord, type = "split", 
                    color = "Phylum", shape = "Site.x")
print(p4)









