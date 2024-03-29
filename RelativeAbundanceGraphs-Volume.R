# Phyloseq Analysis Volume Data ####
# 2022-03-07
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


#Loading Data

physeq_class17_posotuvol <- readRDS("Data/physeq_class17_posotuvol.rds")
physeq_count17_posotuvol <- readRDS("Data/physeq_count17_posotuvol.rds")


physeq_class18_posotuvol <- readRDS("Data/physeq_class18_posotuvol.rds")
physeq_count18_posotuvol <- readRDS("Data/physeq_count18_posotuvol.rds")



#Positive OTUs 2017 

mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(78) # How many colors 
physeq_count17_posotuvol = subset_taxa(physeq_count17_posotuvol, Genus.x != "NA")

plot_bar(physeq_count17_posotuvol, x="Volume_delta", fill= "Genus.x") +
  geom_bar(aes(color=Genus.x, fill = Genus.x), stat = "identity", position = "stack") +
  facet_grid(~Site.x, scales="free_x", space = "free_x") +
  scale_fill_manual(values = mycolors) +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_text(angle=90, hjust=1), 
        text = element_text(size=9))

#Positive OTUs 2018 
mycolors2= colorRampPalette(brewer.pal(8, "Dark2"))(78) # How many colors 
physeq_count18_posotuvol = subset_taxa(physeq_count18_posotuvol, Genus.x != "NA")

# Change x to unique id or samples 
plot_bar(physeq_count18_posotuvol, x= "UniqueID", fill= "Genus.x")+
  geom_bar(aes(color=Genus.x, fill = Genus.x), stat = "identity", position = "stack") +
  facet_grid(~Bucket2, scales="free_x", space = "free_x") +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors2) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_text(angle=90, hjust=1), 
        text = element_text(size=9))#35x20Pic



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
        text = element_text(size=9)) 


####Scaling Volume Column 

# 2017
meta17_data <- read.csv("Data/meta17_data_update.csv")

meta17_data <- meta17_data %>%
  mutate(Volume_scale= scale(Volume_delta))

meta17_data <- subset(meta17_data, select = -c(X.1, X, X.2, X))

write.csv(meta17_data, file = "Data/meta17_data_update.csv")

meta17_data_lowscale <- meta17_data %>% 
  filter(Volume_scale > 0)

meta17_data_highscale <- meta17_data %>%  # Lost 8 values due to the differences in measurements not having values 
  filter(Volume_scale < 0)


#2018

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

meta_gen18_data <- meta_gen18_data %>%
  mutate(Volume_scale= scale(Volume_delta))

meta_gen18_data <- subset(meta_gen18_data, select = -c(X.1, X))

write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")

meta_gen18_data_highscale <- meta_gen18_data %>% 
  filter(Volume_scale > 0)

meta_gen18_data_lowscale <- meta_gen18_data %>%   # Lost 25 values due to the differences in measurements not having values 
  filter(Volume_scale < 0)

# Plot scale and Relative Abundance -Volume 

physeq_posscale <- prune_taxa(rownames())








