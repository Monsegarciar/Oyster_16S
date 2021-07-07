# Questions for Analysis ####
# 2021-06-23
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
library(data.table)
require(RColorBrewer)
library("ggpubr")
library(dplyr)
library(tidyr)

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

# 2017 Data 
richness17= estimate_richness(physeq_class17, split = TRUE, measures = c("Simpson", "Shannon"))

plot_richness(physeq_class17, x="Site.x", measures=c("Shannon", "Simpson"), color = "Site.x")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))

a_my_comparisons17 <- list(c("NW", "OY", "SW"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
?symnum
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

a_my_comparisons18 <- list(c("HIGH_MONO", "HIGH_POLY"), c("LOW_MONO", "LOW_POLY"), c("HIGH_POLY", "LOW_MONO"))
symnum.args18 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

plot_richness(physeq_class18, x="Bucket2", measures=c("Shannon","Simpson"), color = "Bucket2", title = "Boxplot of Treatments 2018")+
  geom_boxplot(alpha=0.6)+ 
  theme(legend.position="none", axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=12))+
  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons18, label = "p.signif", symnum.args = symnum.args18)

hist(richness18$Simpson, main="Simpson index", xlab="")
hist(richness18$Shannon, main="Shannon index", xlab="")

#Question 2 ####
#Looking at peacrabs in sites or treatments 

# Plot bars 
# Plot ordination 
# Alpha diversity- can look into this with how diverse those with peacrabs are, etc. ?

#Question 3 ####
#Weight pre and post with peacrabs

#Adding new column for average weight

#2017 Data
weight_pre <- as.numeric(meta17_data$Weight_pre)

meta17_data <- meta17_data %>%
  mutate(Weight_avg= (Weight_post- weight_pre)/ weight_pre)

#2018 Data
meta_gen18_data <- meta_gen18_data %>%
  mutate(Weight_avg18= (Weight_post- Weight)/ Weight)

#Saving new data
write.csv(meta17_data, file = "Data/meta17_data_update.csv")
write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")

# Analysis 

ggplot(meta17_data, mapping = aes(x=Weight_avg, y=peacrabs.x)) + geom_point()

plot_bar(physeq_class17, "Weight_avg", fill="peacrabs.x", facet_grid=~Site.x) #add scaling to maybe see better

plot_bar(physeq_class17, x="Site.x", fill="peacrabs.x")




