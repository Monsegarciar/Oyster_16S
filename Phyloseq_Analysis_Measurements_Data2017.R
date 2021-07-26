# Phyloseq Analysis for Measurements 2017 Data ####
# 2021-07-18
# Author: Monse Garcia

# 2018 Data only look at Oysters, any found in both, more starting point ####
# look for any in 2017 for increased weight and increased length ####

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

# abundance for otu's in species and clams and save with ggplot####

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





