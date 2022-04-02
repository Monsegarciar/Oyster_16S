# Phyloseq Analysis 2017 Data ####
# 2021-12-17
# Author: Monse Garcia

#Packages Required 
require(RColorBrewer)
require(phyloseq)
require(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(DESeq2)
library(phyloseq)

?phyloseq

#Loading Data
pos_otus17 <- read.csv("Data/pos_otus17.csv")
physeq_class17_posotu <- readRDS("Data/physeq_class17_posotu.rds")
physeq_count17_posotu <- readRDS("Data/physeq_count17_posotu.rds")

pos_otus18 <- read.csv("Data/pos_otus18.csv")
physeq_class18_posotu <- readRDS("Data/physeq_class18_posotu.rds")
physeq_count18_posotu <- readRDS("Data/physeq_count18_posotu.rds")

#Resources 
 #https://stackoverflow.com/questions/38498684/plotting-abundance-by-food-group-in-r

 #https://github.com/joey711/phyloseq/issues/1089#issuecomment-471334036

 #https://tidyverse.tidyverse.org/


#Positive OTUs 2017 

mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(72) # How many colors 
physeq_ccount17_posotu2 = subset_taxa(physeq_count17_posotu, Genus.x != "NA")

plot_bar(physeq_ccount17_posotu2, x="UniqueID", fill= "Genus.x") +
  geom_bar(aes(color=Genus.x, fill = Genus.x), stat = "identity", position = "stack") +
  facet_grid(~Site.x, scales="free_x") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_blank(), 
        text = element_text(size=20))
  
#Positive OTUs 2018 
mycolors2= colorRampPalette(brewer.pal(8, "Dark2"))(65) # How many colors 
physeq_count18_posotu2 = subset_taxa(physeq_count18_posotu, Genus.x != "NA")

# Change x to unique id or samples 
plot_bar(physeq_count18_posotu2, x= "UniqueID", fill= "Genus.x")+
  geom_bar(aes(color=Genus.x, fill = Genus.x), stat = "identity", position = "stack") +
  facet_grid(~Bucket2, scales="free_x") +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors2) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_blank(), 
        text = element_text(size=10))



plot_bar(physeq_count18_posotu2, x= "UniqueID", fill= "Genus.x")+
  geom_bar(aes(color=Genus.x, fill = Genus.x), stat = "identity", position = "stack") +
  facet_grid(~Bucket2, scales="free_x") +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors2) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle=90, hjust=1),
        text = element_text(size=9)) #35x20Pic





