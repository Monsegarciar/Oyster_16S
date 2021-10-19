# Phyloseq Analysis 2017 Data ####
# 2021-06-21
# Author: Monse Garcia

BiocManager::install("ggtree")
a# Required Packages ####
require(phyloseq)
library(data.table)
library(ggplot2)
library(plyr)
library(dplyr)
library(ggtree)
require(metacoder)
#Loading Data ####

meta17_data <- read.csv("Data/meta17_data_update.csv")


asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")
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

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))


#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]
otumat17=asvtable_17

#Converting to matrix
otu_matrix17= as.matrix(otumat17, rownames = rownames(asvtable_17))

tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                         "Genus.x", "Genus.y", "Species")
meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= otu_table(otu_matrix17, taxa_are_rows = FALSE)

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)


OTU_count17=transform_sample_counts(OTU17, function(x) 1E6 * x/sum(x))


physeq_class17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_class17

physeq_count17 = phyloseq(OTU_count17, TAX17, SAMP17)
physeq_count17


# Saving Physeq as an RDS
saveRDS(physeq_class17, "Data/physeq_class17.rds")
physeq_class17 <- readRDS("Data/physeq_class17.rds")

saveRDS(physeq_count17, "Data/physeq_count17.rds")
physeq_count17 <- readRDS("Data/physeq_count17.rds")

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


# Heat Tree ####

tax_17 = parse_phyloseq(physeq_count17)
tax_17

heat_tree(tax_17)


# 2017 Data Statistical Analysis ####

# Setting Cutpoints and Significance Values 
a_my_comparisons17 <- list(c("0", "0.5"), c("0.5", "1"), c("0", "1"), c("1","2"), c("2","3"), c("3","4"), c("4","5"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# Height
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y =Height_delta , colour = RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", caption = "2017 Data", x= "RFTM Score", y= "Normalized Height") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
ggsave(filename = "Height Growth Statistics in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)

# Weight
meta17_data_weight <- meta17_data[-c(80),]
ggplot(data = meta17_data_weight, aes(x = as.factor(RFTM_score.x), y =Weight_delta , colour = RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Weight Growth in Oysters", caption = "2017 Data", x= "RFTM Score", y= "Normalized Weight") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
ggsave(filename = "Weight Growth Statistics in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)

# Length
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y =Length_delta , colour = RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Length Growth in Oysters", caption = "2017 Data", x= "RFTM Score Type", y= "Normalized Length") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
ggsave(filename = "Length Growth Statistics in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)

# Width
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y =Width_delta , colour = RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Width Growth in Oysters", caption = "2017 Data", x= "RFTM Score", y= "Normalized Width") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
ggsave(filename = "Width Growth Statistics in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)
