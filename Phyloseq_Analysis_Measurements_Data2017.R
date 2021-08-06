# Phyloseq Analysis for Measurements 2017 Data ####
# 2021-07-18
# Author: Monse Garcia

# 2018 Data only look at Oysters, any found in both, more starting point ####
# look for any in 2017 and 2018 for increased weight and increased length and see the similar ones####
# abundance for otu's in species and clams and save with ggplot####

install.packages("metacoder")
install.packages("vegan")
install.packages("ggtree")

# Loading packages ####

require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)
library(metacoder)
library(ggtree)

#Loading Data ####

meta17_data <- read.csv("Data/meta17_data_update.csv")

asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

# Taking out weight significant OTU's from DESeq Anlaysis####

# Weight

otu_weight = sigtab17_weight %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Height

otu_height = sigtab17_height %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Gensu.y, Species)

# Length

otu_length= sigtab17_length %>% 
  select(Kindgom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)

# Width

otu_width = sigtab17_width %>% 
  select(Kingdom, Phylum, Class, Order, Family, Genus.x, Genus.y, Species)


# Weight Analysis ####

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

OTU_count17_weight=transform_sample_counts(OTU17_weight, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_weight = phyloseq(OTU17_weight, TAX17_weight, SAMP17_weight)
physeq_class17_weight

physeq_count17_weight = phyloseq(OTU_count17_weight, TAX17_weight, SAMP17_weight)
physeq_count17_weight

# Saving Weight Phyloseq 
saveRDS(physeq_class17_weight, "Data/physeq_class17_weight.rds")
saveRDS(physeq_count17_weight, "Data/physeq_count17_weight.rds")

physeq_count17_weight <- readRDS("Data/physeq_count17_weight.rds")
physeq_class17_weight <- readRDS("Data/physeq_class17_weight.rds")

# Plot ordination Graphs 
physeq_class17_weight = subset_samples(physeq_class17_weight, sam_data(physeq_class17_weight) != "NA")
Phy.ord <- ordinate(physeq_class17_weight, "NMDS", "bray")
p1= plot_ordination(physeq_class17_weight, Phy.ord, type = "taxa", color = "Phylum", title = "taxa")
print(p1)

# Richness Plots
plot_richness(physeq_class17_weight, x="Site.x", measures=c("Simpson", "Shannon"))
plot_richness(physeq_class17_weight, x= "Site.x", color = "Treatment2", measures = c("Simpson", "Shannon"))

a_my_comparisons17_w <- list(c("NW", "OY"), c("OY", "SW"), c("NW", "SW"))
symnum.args17_w = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class17_weight, x="Site.x", measures=c("Shannon","Simpson"), color = "Site.x", title = "Boxplot of Weights and Sites 2017")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17_w, label = "p.signif", symnum.args = symnum.args17_w)

plot_ordination(physeq_count17_weight, Phy.ord, type = "split", 
                color = "Phylum", shape = "Site.x")

ggsave(filename = "Richness Plots", plot=last_plot(), path ="Data2017_plots/", width = 15, height = 8)  


# Height Analysis ####

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
taxmat17_height=otu_height
otumat17_height=asvtable_17

#Converting to matrix
otu_matrix17_height= as.matrix(otumat17_height, rownames = rownames(asvtable_17))

tax_matrix17_height=as.matrix(otu_height)
View(OTU17_height)

meta17_data_height=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_height= otu_table(otu_matrix17_height, taxa_are_rows = FALSE)

TAX17_height= tax_table(tax_matrix17_height)

SAMP17_height= sample_data(meta17_data)

OTU_count17_height=transform_sample_counts(OTU17_height, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_height = phyloseq(OTU17_height, TAX17_height, SAMP17_height)
physeq_class17_height

physeq_count17_height = phyloseq(OTU_count17_height, TAX17_height, SAMP17_height)
physeq_count17_height

# Saving Height Physeq
saveRDS(physeq_class17_height, "Data/physeq_class17_height.rds")
saveRDS(physeq_count17_height, "Data/physeq_count17_height.rds")

# Length Analysis ####

#Setting taxmat and otumat
taxmat17_length=otu_length
otumat17_length=asvtable_17

#Converting to matrix
otu_matrix17_length= as.matrix(otumat17_length, rownames = rownames(asvtable_17))

tax_matrix17_length=as.matrix(otu_length)

meta17_data_length=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_length= otu_table(otu_matrix17_length, taxa_are_rows = FALSE)

TAX17_length= tax_table(tax_matrix17_length)

SAMP17_length= sample_data(meta17_data)

OTU_count17_length=transform_sample_counts(OTU17_length, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_length = phyloseq(OTU17_length, TAX17_length, SAMP17_length)
physeq_class17_length

physeq_count17_length = phyloseq(OTU_count17_length, TAX17_length, SAMP17_length)
physeq_count17_length

saveRDS(physeq_class17_length, "Data/physeq_class17_length.rds")
saveRDS(physeq_count17_length, "Data/physeq_count17_length.rds")

# Width Analysis ####

#Setting taxmat and otumat
taxmat17_width=otu_width
otumat17_width=asvtable_17

#Converting to matrix
otu_matrix17_width= as.matrix(otumat17_width, rownames = rownames(asvtable_17))

tax_matrix17_width=as.matrix(otu_width)

meta17_data_width=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17_width= otu_table(otu_matrix17_width, taxa_are_rows = FALSE)

TAX17_width= tax_table(tax_matrix17_width)

SAMP17_width= sample_data(meta17_data)

OTU_count17_width=transform_sample_counts(OTU17_width, function(x) 1E6 * x/sum(x))

# Phyloseq Class
physeq_class17_width = phyloseq(OTU17_width, TAX17_width, SAMP17_width)
physeq_class17_width

physeq_count17_width = phyloseq(OTU_count17_width, TAX17_width, SAMP17_width)
physeq_count17_width

saveRDS(physeq_class17_width, "Data/physeq_class17_width.rds")
saveRDS(physeq_count17_width, "Data/physeq_count17_width.rds")

# Weight Graphs ####

physeq_count17_weight <- readRDS("Data/physeq_count17_weight.rds")
physeq_count17_weight

# Ordination Graph 

Phy.ord17 <- ordinate(physeq_count17_weight3, "NMDS", "bray") # does not run, but can add 1 in for the zero's 
plot_ordination(physeq_class17_weight3, Phy.ord17, type = "split", 
                       color = "Phylum", shape = "Site.x", title = "Plot Ordination for Weight: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_class17_weight, "Site.x", fill="Weight_diff", facet_grid=~Treatment2)

#mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) 
plot_bar(physeq_class17_weight, x="Site.x", fill = "Weight_diff") # shows abundance 

# Richness plots

plot_richness(physeq_count17_weight, x="Site.x", measures=c("Shannon", "Simpson"), color = "Weight_diff")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# Height Graphs ####

physeq_count17_height <- readRDS("Data/physeq_count17_height.rds")
physeq_count17_height

# Ordination Graph 

Phy.ord17 <- ordinate(physeq_count17_height, "NMDS", "bray")
plot_ordination(physeq_count17_height, Phy.ord17, type = "split", 
                color = "Phylum", shape = "Site.x", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count17_height, "Site.x", fill="Weight_diff", facet_grid=~Treatment2)

#mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) 
plot_bar(physeq_count17_height, x="Site.x", fill = "Weight_diff") # shows abundance 

# Richness plots

plot_richness(physeq_count17_height, x="Site.x", measures=c("Shannon", "Simpson"), color = "Weight_diff")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))


# Length Graphs ####

physeq_count17_length <- readRDS("Data/physeq_count17_length.rds")
physeq_count17_length

# Ordination Graph 

Phy.ord17 <- ordinate(physeq_count17_length, "NMDS", "bray")
plot_ordination(physeq_count17_length, Phy.ord17, type = "split", 
                color = "Phylum", shape = "Site.x", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count17_length, "Site.x", fill="Weight_diff", facet_grid=~Treatment2)

#mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) 
plot_bar(physeq_count17_length, x="Site.x", fill = "Weight_diff") # shows abundance 

# Richness plots

plot_richness(physeq_count17_length, x="Site.x", measures=c("Shannon", "Simpson"), color = "Weight_diff")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))


# Width Graphs ####

physeq_count17_width <- readRDS("Data/physeq_count17_width.rds")
physeq_count17_width

# Ordination Graph 

Phy.ord17 <- ordinate(physeq_count17_width, "NMDS", "bray")
plot_ordination(physeq_count17_width, Phy.ord17, type = "split", 
                color = "Phylum", shape = "Site.x", title = "Plot Ordination for Height: Phylum and Treatments 2018")

# Plot bar 

plot_bar(physeq_count17_width, "Site.x", fill="Weight_diff", facet_grid=~Treatment2)

#mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) 
plot_bar(physeq_count17_width, x="Site.x", fill = "Weight_diff") # shows abundance 

# Richness plots

plot_richness(physeq_count17_width, x="Site.x", measures=c("Shannon", "Simpson"), color = "Weight_diff")+
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

# Merging the data to see similar ones

# Extracting the taxa tables from each phyloseq

taxa_length17 = tax_table(physeq_count17_length)
taxa_weight17 = tax_table(physeq_count17_weight)
taxa_height17 = tax_table(physeq_count17_height)
taxa_width17 = tax_table(physeq_count17_width)

# Cleaning up and merging taxa table into new taxa table 

taxa_measure17 = merge(taxa_height17, taxa_weight17, by= "row.names")
rownames(taxa_measure17) = taxa_measure17$Row.names
taxa_measure17= merge(taxa_measure17, taxa_length17, by="row.names")
taxa_measure17 <- subset(taxa_measure17, select = -c(Row.names, Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, 
                                                     Species.x, Kingdom.y, Phylum.y, Class.y, Order.y, Family.y, Genus.x.y, Genus.y.y, 
                                                     Species.y))
rownames(taxa_measure17)= taxa_measure17$Row.names
taxa_measure17 = merge(taxa_measure17, taxa_width17, by= "row.names")
taxa_measure17 <- subset(taxa_measure17, select = -c(Row.names, Kingdom.x, Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x, Species.x))
rownames(taxa_measure17)= taxa_measure17$Row.names
taxa_measure17 <- subset(taxa_measure17, select = -c(Row.names))
colnames(taxa_measure17) <- c("Kingdom", "Phylum", "Class", "Order", 
                      "Family", "Genus.x","Genus.y", "Species")

# Saving significant taxa for measurements 
write.csv(taxa_measure17, file = "Data/taxa_measure17.csv")
taxa_measure17 <- read.csv("Data/taxa_measure17.csv")

# Heat Trees for Measurements ####


tax_m = parse_phyloseq(physeq_count17_weight) 
heat_tree(tax_m, node_label = taxon_names,
          node_size = n_obs(tax_m), 
          node_color = n_obs(tax_m), 
          layout = "da", initial_layout = "re", 
          title = "Taxa in Width")

heat_tree(tax_m)
tax_m %>% 
  heat_tree(node_label = taxon_names, node_size = n_obs(tax_m), 
            node_color = n_obs(tax_m), layout = "automatic", initial_layout = "fruchterman-reingold")

# Taxa for Measurements 
rownames(taxa_measure17)= taxa_measure17$X
taxa_measure17$X = NULL
tax_measure17 = parse_tax_data(taxa_measure17)
tax_measure17
heat_tree(tax_measure17, node_label = taxon_names,
          node_size = n_obs(tax_measure17), 
          node_color = n_obs(tax_measure17), 
          layout = "automatic", initial_layout = "automatic", 
          title = "Taxa in Measurements")

# geom_col()
# RFTM_dds18 <- phyloseq_to_deseq2(physeq_class, ~RFTM_score.x+Species)
# Look at taxonomy- graph specific taxa (e.g. class, genus)
# log 2 fold change and see which ones more associated with weight 
# Other variables involved include site in 2017 

