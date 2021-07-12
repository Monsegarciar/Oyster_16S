# Questions for Analysis ####
# 2021-06-23
# Author: Monse Garcia

# Downloading DESeq2

BiocManager::install("DESeq2")

#Packages Required
require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)
library(DESeq2)


# Loading data 

meta17_data <- read.csv("Data/meta17_data_update.csv")
asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")



#Loading Physeq w/out transform_sample_counts() function
physeq_class17 <- readRDS("Data/physeq_class17.rds")
physeq_class18 <- readRDS("Data/physeq_class18.rds")

# Example of code with ColorBrewer
mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) # How many colors
plot_bar(pp.ch,  fill="Category", x="Replicate") +
  geom_bar(aes(color=Category, fill=Category), stat="identity", position="stack")+
  facet_grid(Year~Site_Name, scales="free_x")+
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  theme(legend.position="top", legend.text=element_text(size=10), panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10))
#want a specific color with hexcodes or do red and blue or other colors instead of mycolors in the code
#search r color pallets different kinds 

#Question 1 ####
#Sites relationship with oyster species
#Example: Which sites are more prominent to bacteria or do more bacteria diversity in a specific sites. 
#Alpha diversity 
#Diversity with other factors
#Taxonomy bar to see if there are more abundance in one site than the other

#Split Graph with Site and Phylum
p4= plot_ordination(physeq_class17, Phy.ord, type = "split", 
                    color = "Phylum", shape = "Site.x")
print(p4)

##Plot bars with phylum 
table_taxa <-table(taxmat17$V3) # Used the table() function to see which phylum was most common
View(table_taxa)

#Most common phylum and site graph
pp.ch= subset_taxa(physeq_class17, Phylum=="Proteobacteria") 
plot_bar(pp.ch) #Plot bar of samples(x) and abundance (y) of Proteobacteria
bar1=plot_bar(pp.ch, x="Site.x", fill = "Genus")
print(bar1)

mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(299)
plot_bar(pp.ch,  fill="Genus", x="Treatment2") +
  geom_bar(aes(color=Genus, fill=Genus), stat="identity", position="stack")+
  facet_grid(Site.x~Species.x, scales="free_x")+
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  theme(legend.position="top", legend.text=element_text(size=10), panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10))

b2=plot_bar(pp.ch, "Site.x", fill="Genus", facet_grid=~Family) 
print(b2)
b2 + geom_point(aes(x=Family, y=Abundance), color="black", position="jitter", size=3) #came up too big to fit, might use top OTU's instead

# Top OTU's ####
#2017 Data
TopNOTUs <- names(sort(taxa_sums(physeq_class17), TRUE)[1:10])
phys10 <- prune_species(TopNOTUs, physeq_class17)


p= plot_bar(phys10, "Treatment2", fill="Site.x", facet_grid=~Genus, title= "Top NOTUs SIte and Treatment 2017") # Switching Treatment w/ site
print(p)

p + geom_bar(aes(color=Site.x, fill=Site.x), stat="identity", position="stack")

# 2018
TopNOTUs18 <- names(sort(taxa_sums(physeq_class18), TRUE)[1:20])
phys20_18 <- prune_species(TopNOTUs18, physeq_class18)

p_2= plot_bar(phys20_18, "Bucket2", fill="Species2.x", facet_grid=~Genus, 
             title= "Top NOTUs in Species and Treatment 2018") # Switching Treatment w/ site
print(p)

p_2 + geom_bar(aes(color=Species2.x, fill=Species2.x), stat="identity", position="stack", na.rm = TRUE)

# Alpha Diversity Graphics ####

#2017 Data
?plot_richness
plot_rich= plot_richness(physeq_class17, x="Site.x", measures=c("Simpson", "Shannon"))
print(plot_rich)

plot_rich2= plot_richness(physeq_class17, x= "Site.x", color = "Treatment2", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for Treatment and Species 2017")
print(plot_rich2)

#2018 Data
plo_rich18= plot_richness(physeq_class18, x="Bucket2", measures=c("Chao1", "Shannon"))
print(plo_rich18)

plo_rich_species18= plot_richness(physeq_class18, x="Bucket2", color = "Species2.x",measures=c("Simpson", "Shannon"), title = "Alpha Diveristy for Treatments and Species 2018")
print(plo_rich_species18)

#Standard deviation and mean ####

# Getting NA's in Bucket2 Column a value
meta_gen18_data$Bucket2[is.na(meta_gen18_data$Bucket2)] <- "LOW_POLY"

# 2017 Data 
richness17= estimate_richness(physeq_class17, split = TRUE, measures = c("Simpson", "Shannon"))

plot_richness(physeq_class17, x="Site.x", measures=c("Shannon", "Simpson"), color = "Site.x")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

a_my_comparisons17 <- list(c("NW", "OY"), c("OY", "SW"), c("NW", "SW"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class17, x="Site.x", measures=c("Shannon","Simpson"), color = "Site.x", title = "Boxplot of Sites 2017")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

hist(richness17$Shannon, main="Shannon index", xlab="")
hist(richness17$Simpson, main="Simpson index", xlab="")

#2018 Data
richness18= estimate_richness(physeq_class18, split= TRUE, measures= c("Simpson", "Shannon"))

plot_richness(physeq_class18, x="Bucket2", measures=c("Shannon", "Simpson"), color = "Bucket2")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

a_my_comparisons18 <- list(c("HIGH_MONO", "HIGH_POLY"), c("LOW_MONO", "LOW_POLY"), c("HIGH_POLY", "LOW_MONO"), c("HIGH_MONO", "LOW_POLY"), c("HIGH_MONO", "LOW_MONO"), c("HIGH_POLY", "LOW_POLY"))
symnum.args18 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class18, x="Bucket2", measures=c("Shannon","Simpson"), color = "Bucket2", title = "Boxplot of Treatments 2018")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)

hist(richness18$Simpson, main="Simpson index", xlab="")
hist(richness18$Shannon, main="Shannon index", xlab="")

#Question 2 ####
#Looking at pea crabs in sites or treatments 

# Plot bars 
# Plot ordination  diverse those with 

# RFTM Score and Species (Alpha Diversity)
plot_richness(physeq_class18, x= "Species2.x", color = "RFTM_score.x", measures = c("Simpson", "Shannon"), title = "RFTM Score in Species 2018")

a_my_comparisons18_Species <- list(c("AM", "CV"), c("IR", "LP"), c("CV", "IR"))
symnum.args18 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class18, x="Species2.x", measures=c("Shannon","Simpson"), color = "RFTM_score.x", title = "Boxplot of RFTM Scores and Species 2018")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18_Species, label = "p.signif", symnum.args = symnum.args18)

#Alpha Diversity of Species Significance
plot_richness(physeq_class18, x= "Species2.x", color = "Species2.x", measures = c("Simpson", "Shannon"), title = "Diveristy of Species 2018")

plot_richness(physeq_class18, x="Species2.x", measures=c("Shannon","Simpson"), color = "Species2.x", title = "Boxplot of Species 2018")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18_Species, label = "p.signif", symnum.args = symnum.args18)


# Pea crabs and Sites 
plot_richness(physeq_class17, x= "Site.x", color = "Weight_diff", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for 2017") 

a_my_comparisons17_weight <- list(c("NW", "OY"), c("OY", "SW"), c("NW", "SW"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class17, x="Site.x", measures=c("Shannon","Simpson"), color = "Weight_diff", title = "Boxplot of Weight 2018")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17_weight, label = "p.signif", symnum.args = symnum.args17)

#Question 3 ####

#Adding new column for average weight

#2017 Data
weight_pre <- as.numeric(meta17_data$Weight_pre)

meta17_data <- meta17_data %>%
  mutate(Weight_avg= (Weight_post- weight_pre)/ weight_pre)

#2018 Data
meta_gen18_data <- meta_gen18_data %>%
  mutate(Weight_avg18= (Weight_post- Weight)/ Weight)

# Adding new column for averages in length, width, height

# 2017 Data
height_pre <- as.numeric(meta17_data$Height_pre)

meta17_data <- meta17_data %>%
  mutate(Height_avg= (Height_post- height_pre)/ height_pre)

width_pre <- as.numeric(meta17_data$Width_pre)

meta17_data <- meta17_data %>%
  mutate(Width_avg= (Width_post- width_pre)/ width_pre)

length_pre <- as.numeric(meta17_data$Length_pre)

meta17_data <- meta17_data %>%
  mutate(Length_avg= (Length_post- length_pre)/ length_pre)


#2018 Data
meta_gen18_data <- meta_gen18_data %>%
  mutate(Height_avg= (Height_post- Height_pre)/ Height_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Width_avg= (Width_post- Width_pre)/ Width_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Length_avg= (Length_post- Length_pre)/ Length_pre)

# Adding Columns with weight, length, height, width differences

#2017 Data
weight_pre <- as.numeric(meta17_data$Weight_pre)

meta17_data <- meta17_data %>%
  mutate(Weight_diff= Weight_post- weight_pre)

height_pre <- as.numeric(meta17_data$Height_pre)

meta17_data <- meta17_data %>%
  mutate(Height_diff= Height_post- height_pre)

width_pre <- as.numeric(meta17_data$Width_pre)

meta17_data <- meta17_data %>%
  mutate(Width_diff= Width_post- width_pre)

length_pre <- as.numeric(meta17_data$Length_pre)

meta17_data <- meta17_data %>%
  mutate(Length_diff= Length_post- length_pre)

#2018 Data
meta_gen18_data <- meta_gen18_data %>%
  mutate(Height_diff= Height_post- Height_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Width_diff= Width_post- Width_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Length_diff= Length_post- Length_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Weight_diff= Weight_post- Weight)

meta17_data <- subset(meta17_data, select = -c(X.1, X, X.2, Number.x, Number.y))
meta_gen18_data <- subset(meta_gen18_data, select = -c(X, X.1, X.2))

#Saving new data
write.csv(meta17_data, file = "Data/meta17_data_update.csv")
write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")

# Analysis 

plot_bar(physeq_class17, "Weight_avg", fill="peacrabs.x", facet_grid=~Site.x) #add scaling to maybe see better

#Weight difference with pea crabs

#2017 Data
ggplot(meta17_data, aes(x=Site.x, y= Weight_diff))+
  geom_point(alpha =.3,aes(color= peacrabs.x)) + facet_grid(rows = vars(peacrabs.x)) # Point graph with pea crabs and weight

ggplot(meta17_data, aes(x=Site.x, y= Weight_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(peacrabs.x))) + geom_boxplot(alpha=0) + facet_grid(rows = vars(peacrabs.x), cols = vars(Site.x), scales = "free") # Jitter plot and slpit graphs with pea crabs and weight

ggplot(meta17_data, aes(x=Site.x, y= Height_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(peacrabs.x)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(peacrabs.x), cols = vars(Site.x), scales = "free")# Jitter plot with pea crabs and height

ggplot(meta17_data, aes(x=Site.x, y= Length_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(peacrabs.x)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(peacrabs.x), cols = vars(Site.x), scales = "free")# Jitter plot with pea crabs and length

ggplot(meta17_data, aes(x=Site.x, y= Width_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(peacrabs.x)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(peacrabs.x), cols = vars(Site.x), scales = "free")# Jitter plot with pea crabs and width

#2018 Data
ggplot(meta_gen18_data, aes(x=Bucket2, y= Weight_diff))+
  geom_point(alpha =.5,aes(color= as.factor(RFTM_score.y))) + geom_boxplot(alpha=0)  # Point graph with RFTM score and weight

ggplot(meta_gen18_data, aes(x=Bucket2, y= Weight_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(RFTM_score.y))) + geom_boxplot(alpha=0)+ facet_grid(rows = vars(Species2.x), cols = vars(Bucket2), scales = "free") # Jitter plot with RFTM score and weight

ggplot(meta_gen18_data, aes(x=Bucket2, y= Height_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(RFTM_score.y)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(Species2.x), cols = vars(Bucket2), scales = "free") # Jitter plot with RFTM score and height

ggplot(meta_gen18_data, aes(x=Bucket2, y= Length_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(RFTM_score.y)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(Species2.x), cols = vars(Bucket2), scales = "free")# Jitter plot with RFTM score and length

ggplot(meta_gen18_data, aes(x=Bucket2, y= Width_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(RFTM_score.y)))+ geom_boxplot(alpha=0) + facet_grid(rows = vars(Species2.x), cols = vars(Bucket2), scales = "free") #Jitter plot with RFTM score and width 

# X-axis being the Treatment type- 2017 Data
ggplot(meta17_data, aes(x=Treatment2, y= Weight_diff))+
  geom_jitter(alpha =.8,aes(color= as.factor(peacrabs.x))) + geom_boxplot(alpha=0) + facet_grid(rows = vars(peacrabs.x), cols = vars(Site.x), scales = "free") # Jitter plot and slpit graphs with pea crabs and weight

?plot_richness
# DESeq 2 ####

physeq_class17
head(sample_data(physeq_class17)$Treatment2)
physeq_class17= subset_samples(physeq_class17, )
# RFTM SCORE with DESeq2- 2017 Data
deseq17_score = phyloseq_to_deseq2(physeq_class17, ~ RFTM_score.y)
deseq17_score = DESeq(deseq17_score, test="Wald", fitType="parametric")

res17_score = results(deseq17_score, cooksCutoff = FALSE)
alpha = 0.05 # switch to 0.05- generalized linear models
sigtab17_score = res17_score[which(res17_score$padj < alpha), ]
sigtab17_score = cbind(as(sigtab17_score, "data.frame"), as(tax_table(physeq_class17)[rownames(sigtab17_score), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
  }
x = tapply(sigtab17_score$log2FoldChange, sigtab17_score$Genus, function(x) max(x))
x = sort(x, TRUE)
sigtab17_score$Genus = factor(as.character(sigtab17_score$Genus), levels=names(x))

x = tapply(sigtab17_score$log2FoldChange, sigtab17_score$Class, function(x) max(x))
x = sort(x, TRUE)
head(sigtab17_score)

sigtab17_score$Class = factor(as.character(sigtab17_score$Class), levels=names(x))
ggplot(sigtab17_score, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))

# RFTM score with DESeq2- 2018 Data
estimateSizeFactors(physeq_class18, type = 'iterate')

deseq18_score = phyloseq_to_deseq2(physeq_class18, ~ RFTM_score.x)
deseq18_score = DESeq(deseq18_score, test="Wald", fitType="parametric") # Get error: Error in estimateSizeFactorsForMatrix(counts(object), locfunc = locfunc,  : every gene contains at least one zero, cannot compute log geometric means

res18_score = results(deseq18_score, cooksCutoff = FALSE)
alpha = 0.01
sigtab18_score = res18_score[which(res18_score$padj < alpha), ]
sigtab18_score = cbind(as(sigtab18_score, "data.frame"), as(tax_table(physeq_class18)[rownames(sigtab18_score), ], "matrix"))

theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}
x = tapply(sigtab18_score$log2FoldChange, sigtab18_score$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab18_score$Phylum = factor(as.character(sigtab18_score$Phylum), levels=names(x))

x = tapply(sigtab18_score$log2FoldChange, sigtab18_score$Genus, function(x) max(x))
x = sort(x, TRUE)


sigtab18_score$Genus = factor(as.character(sigtab18_score$Genus), levels=names(x))
ggplot(sigtab18_score, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5))
