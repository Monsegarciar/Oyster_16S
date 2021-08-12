# Loading packages ####

require(phyloseq)
require(ggplot2)
require(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)

# Loading Data ####
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")
Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

# Separating Oysters
metagen18_oys<- meta_gen18_data %>% 
  filter(Species2.x == "CV")
asvtable_oys18<- asvtable_18[c(1:15,21:23,28,29,34:39,45:59,61:64,70,74:79,85,91:94,98), ]

metagen18_muss <- meta_gen18_data %>% 
  filter(Species2.x == "IR")
asvtable_muss18 <- asvtable_18[c(16:20,24:27,30,40:44,65:69,80:84,95:97,99,100)]

metagen18_shell<- meta_gen18_data %>% 
  filter(Species2.x == "AM")
asvtable_shell18<- asvtable_18[c(32,33,60,73,89,90)]

metagen18_clam<- meta_gen18_data %>% 
  filter(Species2.x == "LP")
asvtable_clam18<- asvtable_18[c(31,71,72,86:88)]

#Changing row names in meta_gen18 data
rownames(metagen18_oys)= metagen18_oys$UniqueID 
head(rownames(metagen18_oys))

#Changing rownames in asvtable data
rownames(asvtable_oys18)= asvtable_oys18$V1
head(rownames(asvtable_oys18))

#Setting taxmat and otumat
otumat18_oys=asvtable_18
taxmat18_oys=Run123_taxa

#Converting to matrix
otu_matrix18_oys= as.matrix(otumat18_oys, rownames = "V1")

tax_matrix18_oys=as.matrix(taxmat18_oys, rownames = "V2")
 
metagen18_oys=as.data.frame(metagen18_oys)

#Setting OTU, TAX, and SAMP
OTU18_oys= otu_table(otu_matrix18_oys, taxa_are_rows = FALSE)

TAX18_oys= tax_table(tax_matrix18_oys)

SAMP18_oys= sample_data(metagen18_oys)

OTU_count18_oys=transform_sample_counts(OTU18_oys, function(x) 1E6 * x/sum(x))

physeq_class18_oys = phyloseq(OTU18_oys, TAX18_oys, SAMP18_oys)
physeq_class18_oys

physeq_count18_oys = phyloseq(OTU_count18_oys, TAX18_oys, SAMP18_oys)
physeq_count18_oys

# Graphs for Oysters ####

plot_bar(physeq_count18_oys, x="Bucket2", fill = "Species2.x") + geom_col()

# Weight growth in oysters
ggplot(data = metagen18_oys, aes(x = Weight_diff, y = Weight_delta, colour = Bucket2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Weight Growth in Oysters")
ggplot(data = metagen18_oys, aes(x =RFTM_score.x, y = Weight_delta)) + geom_point() #colour = "red"+ facet_wrap(~Color_Bucket)

# Weight Growth in mussels
ggplot(data = metagen18_muss, aes(x = Weight_diff, y = Weight_delta, colour = Bucket2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Weight Growth in Mussels")

# Weight Growth in clams
ggplot(data = metagen18_clam, aes(x = Weight_diff, y = Weight_delta, colour = Bucket2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Weight Growth in Clams") # no measurements for clams

# Weight Growth in shell
ggplot(data = metagen18_shell, aes(x = Weight_diff, y = Weight_delta, colour = Bucket2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Weight Growth in Shells") # no measurements for shells

# 2017 Data

meta17_data <- read.csv("Data/meta17_data_update.csv")

# Height
ggplot(data = meta17_data, aes(x = Height_diff, y = Height_delta, colour = Treatment2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Height Growth in Oysters", caption = "2017 Data", x= "Height Difference", y= "Normalized Height")
ggsave(filename = "Height Growth in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)             

# Length
ggplot(data = meta17_data, aes(x = Length_diff, y = Length_delta, colour = Treatment2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Length Growth in Oysters", caption = "2017 Data", x= "Length Difference", y= "Normalized Length")
ggsave(filename = "Length Growth in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)             

# Width
ggplot(data = meta17_data, aes(x = Width_diff, y = Width_delta, colour = Treatment2)) + geom_point() + facet_wrap(~RFTM_score.x) + labs(title = "Width Growth in Oysters", caption = "2017 Data", x= "Width Difference", y= "Normalized Width")
ggsave(filename = "Width Growth in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)             

# 2018 Data Statistical Analysis ####

# Weight 
a_my_comparisons18 <- list(c("0", "0.5"), c("0.5", "1"), c("0", "1"))
symnum.args18 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

ggplot(data = metagen18_oys, aes(x = as.factor(RFTM_score.x), y = Weight_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Weight Growth in Oysters", x= "RFTM Score", y= "Normalized Weight") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)
ggsave(filename = "Weight Growth Statistics in Oysters 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5)

# Height 
ggplot(data = metagen18_oys, aes(x = as.factor(RFTM_score.x), y = Height_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", x= "RFTM Score", y= "Normalized Height") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)
ggsave(filename = "Height Growth Statistics in Oysters 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5)

# Length
ggplot(data = metagen18_oys, aes(x = as.factor(RFTM_score.x), y = Length_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Length Growth in Oysters", x= "RFTM Score", y= "Normalized Length") 

ggplot(data = metagen18_oys, aes(x = as.factor(RFTM_score.x), y = Length_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Length Growth Statistics in Oysters", x= "RFTM Score", y= "Normalized Length") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)
ggsave(filename = "Length Growth Statistics in Oysters 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5)

# Width
ggplot(data = metagen18_oys, aes(x = as.factor(RFTM_score.x), y = Width_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Width Growth in Oysters", x= "RFTM Score", y= "Normalized Width") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)
ggsave(filename = "Width Growth Statistics in Oysters 2018.jpeg", plot=last_plot(), path ="Data2018_plots/", width = 7, height = 5)

#
ggplot(data = metagen18_oys, aes(x = Weight_diff, y = Weight_delta, colour = Bucket2)) + geom_point() + facet_wrap(~Bucket2) + labs(title = "Weight Growth in Clams") # no measurements for clams

ggplot(data = meta17_data, aes(x = Height_diff, y = Height_delta, colour = Treatment2)) + geom_point() + facet_wrap(~Treatment2) + labs(title = "Height Growth in Oysters", caption = "2017 Data", x= "Height Difference", y= "Normalized Height")

# Treatment Statistical Analysis ####
a_my_comparisons17 <- list(c("HIGH_MONO", "HIGH_POLY"), c("HIGH_POLY", "LOW_MONO"), c("LOW_MONO", "LOW_POLY"), c("HIGH_MONO", "LOW_POLY"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

# Height
ggplot(data = meta17_data, aes(x = Treatment2, y =Height_delta , colour = Treatment2)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", caption = "2017 Data", x= "Treatment", y= "Normalized Height") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

# Weight
ggplot(data = meta17_data, aes(x = Treatment2, y =Weight_delta , colour = Treatment2)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Weight Growth in Oysters", caption = "2017 Data", x= "Treatment", y= "Normalized Weight") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

# Length
ggplot(data = meta17_data, aes(x = Treatment2, y =Length_delta , colour = Treatment2)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Length Growth in Oysters", caption = "2017 Data", x= "Treatment", y= "Normalized Length") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

# Width
ggplot(data = meta17_data, aes(x = Treatment2, y =Width_delta , colour = Treatment2)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Width Growth in Oysters", caption = "2017 Data", x= "Treatment", y= "Normalized Width") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

