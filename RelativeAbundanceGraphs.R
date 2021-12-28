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
  
mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) # How many colors 

plot_bar(physeq_count17_posotu, x= "Sample", fill= "Genus.x")

#Positive OTUs 2018 

plot_bar(physeq_count18_posotu, x= "Sample", fill= "Genus.x")







