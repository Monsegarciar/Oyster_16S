# Phyloseq Analysis for Weight 2017 Data ####
# 2021-07-18
# Author: Monse Garcia


#Loading Data ####

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


#Setting taxmat and otumat
taxmat17_weight=otu_weight
otumat17_weight=asvtable_17

#Converting to matrix
otu_matrix17_weight= as.matrix(otumat17_weight, rownames = rownames(asvtable_17))

tax_matrix17_weight=as.matrix(otu_weight)
View(OTU17_weight)

meta17_data_weight=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_weight= otu_table(otu_matrix17_weight, taxa_are_rows = FALSE)

TAX17_weight= tax_table(tax_matrix17_weight)

SAMP17_weight= sample_data(meta17_data)

# Phyloseq Class
physeq_class17_weight = phyloseq(OTU17_weight, TAX17_weight, SAMP17_weight)

physeq_class17_weight

# Plot ordination Graphs 
Phy.ord <- ordinate(physeq_class17_weight, "NMDS", "bray")
p1= plot_ordination(physeq_class17_weight, Phy.ord, type = "taxa", color = "Phylum", title = "taxa")

plot_richness(physeq_class17_weight, x="Site.x", measures=c("Simpson", "Shannon"))
plot_richness(physeq_class17_weight, x= "Site.x", color = "Treatment2", measures = c("Simpson", "Shannon"))

              