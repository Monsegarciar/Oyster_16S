# Permonova Tests  ####
# 2021-06-23
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(phyloseq)

# Permanova for 

set.seed(1)

# Calculate bray curtis distance matrix
erie_bray <- phyloseq::distance(erie_scale, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(erie))

# Adonis test
adonis(erie_bray ~ Station, data = sampledf)