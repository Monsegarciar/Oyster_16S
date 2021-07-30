# Phyloseq Analysis for Measurements 2018 Data ####
# 2021-07-21
# Author: Monse Garcia


# Loading packages ####

require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)

#Loading Data ####

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")

# Taking out weight significant OTU's from DESeq Anlaysis####

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

# Weight Analysis ####

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Changing rownames in asvtable data
rownames(asvtable_18)= asvtable_18$V1
head(rownames(asvtable_18))

#Setting taxmat and otumat
otumat18_weight=asvtable_18
taxmat18_weight=otu_weight18

#Converting to matrix
otu_matrix18_weight= as.matrix(otumat18_weight, rownames = "V1")

tax_matrix18_weight=as.matrix(taxmat18_weight, rownames = "V2")
colnames(tax_matrix18_weight) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18_weight= otu_table(otu_matrix18_weight, taxa_are_rows = FALSE)

TAX18_weight= tax_table(tax_matrix18_weight)

SAMP18_weight= sample_data(meta_gen18_data)

OTU_count18_weight=transform_sample_counts(OTU18_weight, function(x) 1E6 * x/sum(x))

physeq_class18_weight = phyloseq(OTU18_weight, TAX18_weight, SAMP18_weight)
physeq_class18_weight

physeq_count18_weight = phyloseq(OTU_count18_weight, TAX18_weight, SAMP18_weight)
physeq_count18_weight

# Saving phyloseq's 

saveRDS(physeq_count18_weight, "Data/physeq_count18_weight.rds")

saveRDS(physeq_class18_weight, "Data/physeq_class18_weight.rds")


# Height Analysis ####

#Setting taxmat and otumat
otumat18_height=asvtable_18
taxmat18_height=otu_height18

#Converting to matrix
otu_matrix18_height= as.matrix(otumat18_height, rownames = "V1")

tax_matrix18_height=as.matrix(taxmat18_height, rownames = "V2")
colnames(tax_matrix18_height) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                   "Genus.x", "Genus.y", "Species")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18_height= otu_table(otu_matrix18_height, taxa_are_rows = FALSE)

TAX18_height= tax_table(tax_matrix18_height)

SAMP18_height= sample_data(meta_gen18_data)

OTU_count18_height=transform_sample_counts(OTU18_height, function(x) 1E6 * x/sum(x))

physeq_class18_height = phyloseq(OTU18_height, TAX18_height, SAMP18_height)
physeq_class18_height

physeq_count18_height = phyloseq(OTU_count18_height, TAX18_height, SAMP18_height)
physeq_count18_height

# Saving phyloseq's

saveRDS(physeq_count18_height, "Data/physeq_count18_height.rds")

saveRDS(physeq_class18_height, "Data/physeq_class18_height.rds")


physeq_count18_height <- readRDS("Data/physeq_count18_height.rds")
physeq_count18_height
# Length Analysis ####

#Setting taxmat and otumat
otumat18_length=asvtable_18
taxmat18_length=otu_length18

#Converting to matrix
otu_matrix18_length= as.matrix(otumat18_length, rownames = "V1")

tax_matrix18_length=as.matrix(taxmat18_length, rownames = "V2")
colnames(tax_matrix18_length) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                   "Genus.x", "Genus.y", "Species")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18_length= otu_table(otu_matrix18_length, taxa_are_rows = FALSE)

TAX18_length= tax_table(tax_matrix18_length)

SAMP18_length= sample_data(meta_gen18_data)

OTU_count18_length=transform_sample_counts(OTU18_length, function(x) 1E6 * x/sum(x))

physeq_class18_length = phyloseq(OTU18_length, TAX18_length, SAMP18_length)
physeq_class18_length

physeq_count18_length = phyloseq(OTU_count18_length, TAX18_length, SAMP18_length)
physeq_count18_length

# Saving phyloseq's

saveRDS(physeq_count18_length, "Data/physeq_count18_length.rds")

saveRDS(physeq_class18_length, "Data/physeq_class18_length.rds")

# Width Analysis ####

#Setting taxmat and otumat
otumat18_width=asvtable_18
taxmat18_width=otu_width18

#Converting to matrix
otu_matrix18_width= as.matrix(otumat18_width, rownames = "V1")

tax_matrix18_width=as.matrix(taxmat18_width, rownames = "V2")
colnames(tax_matrix18_width) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                                   "Genus.x", "Genus.y", "Species")
meta_gen18_data=as.data.frame(meta_gen18_data)

#Setting OTU, TAX, and SAMP
OTU18_width= otu_table(otu_matrix18_width, taxa_are_rows = FALSE)

TAX18_width= tax_table(tax_matrix18_width)

SAMP18_width= sample_data(meta_gen18_data)

OTU_count18_width=transform_sample_counts(OTU18_width, function(x) 1E6 * x/sum(x))

physeq_class18_width = phyloseq(OTU18_width, TAX18_width, SAMP18_width)
physeq_class18_width

physeq_count18_width = phyloseq(OTU_count18_width, TAX18_width, SAMP18_width)
physeq_count18_width

# Saving phyloseq's

saveRDS(physeq_count18_width, "Data/physeq_count18_width.rds")

saveRDS(physeq_class18_width, "Data/physeq_class18_width.rds")

# Weight Graphs ####

# Loading phyloseq
physeq_count18_weight <- readRDS("Data/physeq_count18_weight.rds")
physeq_count18_weight


# Ordination Graph 

Phy.ord18 <- ordinate(physeq_count18_weight, "NMDS", "bray")
plot_ordination(physeq_count18_weight, Phy.ord18, type = "split", 
                       color = "Phylum", shape = "Bucket2", title = "Plot Ordination for Weight: Phylum and Treatments 2018")


# Plot bar

plot_bar(physeq_count18_weight, "Bucket2", fill="Weight_delta", facet_grid=~Species2.x)
plot_bar(physeq_count18_weight, x="Bucket2", fill = "Species2.x") # shows abundance 


# Richness plots

plot_richness(physeq_count18_weight, x="Bucket2", measures=c("Shannon", "Simpson"), color = "Species2.x")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# Height Graphs ####

# Loading  phyloseq 

physeq_count18_height <- readRDS("Data/physeq_count18_height.rds")
physeq_count18_height

# Ordination Graph 

Phy.ord18 <- ordinate(physeq_count18_height, "NMDS", "bray")
plot_ordination(physeq_count18_height, Phy.ord18, type = "split", 
                color = "Phylum", shape = "Bucket2", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count18_height, "Bucket2", fill="Height_delta", facet_grid=~Species2.x)
plot_bar(physeq_count18_height, x="Bucket2", fill = "Species2.x") # shows abundance 

# Richness plots

plot_richness(physeq_count18_height, x="Bucket2", measures=c("Shannon", "Simpson"), color = "Species2.x")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))


# Length Graphs ####

# Loading  phyloseq 

physeq_count18_length <- readRDS("Data/physeq_count18_length.rds")
physeq_count18_length

# Ordination Graph 

Phy.ord18 <- ordinate(physeq_count18_length, "NMDS", "bray")
plot_ordination(physeq_count18_length, Phy.ord18, type = "split", 
                color = "Phylum", shape = "Bucket2", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count18_length, "Bucket2", fill="Height_delta", facet_grid=~Species2.x)
plot_bar(physeq_count18_length, x="Bucket2", fill = "Species2.x") # shows abundance 

# Richness plots

plot_richness(physeq_count18_length, x="Bucket2", measures=c("Shannon", "Simpson"), color = "Species2.x")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))


# Width Graphs ####

# Loading  phyloseq 

physeq_count18_width <- readRDS("Data/physeq_count18_width.rds")
physeq_count18_width

# Ordination Graph 

Phy.ord18 <- ordinate(physeq_count18_width, "NMDS", "bray")
plot_ordination(physeq_count18_width, Phy.ord18, type = "split", 
                color = "Phylum", shape = "Bucket2", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count18_width, "Bucket2", fill="Height_delta", facet_grid=~Species2.x)
plot_bar(physeq_count18_width, x="Bucket2", fill = "Species2.x") # shows abundance 

# Richness plots

plot_richness(physeq_count18_width, x="Bucket2", measures=c("Shannon", "Simpson"), color = "Species2.x")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# Merging the data to similar ones ####

# Extracting the taxa tables from each phyloseq

taxa_length18 = tax_table(physeq_count18_length)
taxa_weight18 = tax_table(physeq_count18_weight)
taxa_height18 = tax_table(physeq_count18_height)
taxa_width18 = tax_table(physeq_count18_width)

# Cleaning up and merging taxa table into new taxa table 

taxa_measure18 = merge(taxa_height18, taxa_weight18, by= "row.names")
rownames(taxa_measure18) = taxa_measure18$Row.names
taxa_measure18 = merge(taxa_measure18, taxa_length18, by= "row.names")
taxa_measure18 <- subset(taxa_measure18, select = -c(Row.names, Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x, 
                                                     Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, Species.y))
rownames(taxa_measure18)= taxa_measure18$Row.names
taxa_measure18= merge(taxa_measure18, taxa_width18, by= "row.names")
rownames(taxa_measure18) = taxa_measure18$Row.names
taxa_measure18 <- subset(taxa_measure18, select= -c(Row.names,Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x))
taxa_measure18 <- subset(taxa_measure18, select = -c(Row.names))
colnames(taxa_measure18) <- c("Kingdom", "Phylum", "Class", "Order", 
                              "Family", "Genus.x","Genus.y", "Species")

# Saving significant taxa for measurements 
write.csv(taxa_measure18, file = "Data/taxa_measure18.csv")

# Heat Tree for Measuremnts ####

