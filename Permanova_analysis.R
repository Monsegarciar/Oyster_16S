# Permanova Tests  ####
# 2021-11-11
# Author: Monse Garcia

#Packages Required
require(phyloseq)
require(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(tidyr)

# Loading Data 
meta17_data <- read.csv("Data/meta17_data_update.csv")
physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17_w <- readRDS("Data/physeq_count17_w.rds")

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
physeq_count18 <- readRDS("Data/physeq_count18.rds")
physeq_class18 <- readRDS("Data/physeq_class18.rds")

# Permanova  

set.seed(1)

# Calculate bray curtis distance matrix
erie_bray <- phyloseq::distance(physeq_count17, method = "bray")

erie_bray_w <- phyloseq::distance(physeq_count17_w, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq_count17))

sampledf_w <- data.frame(sample_data(physeq_count17_w))

# Adonis test
adonis(erie_bray ~ Site.x * peacrabs.x * RFTM_score.x, data = sampledf)
adonis(erie_bray ~ Site.x + peacrabs.x + RFTM_score.x, data = sampledf)

adonis(erie_bray_w ~ Site.x * Weight_diff, data = sampledf_w)
adonis(erie_bray_w ~ Site.x + Weight_diff, data = sampledf_w)

#Dispersion test
beta <- betadisper(erie_bray, sampledf$Site.x)
permutest(beta)



#2018 

#RFTM Score and Treatment 

set.seed(1)

# Calculate bray curtis distance matrix
bray18 <- phyloseq::distance(physeq_count18, method = "bray")

# make a data frame from the sample_data
sampledf_18 <- data.frame(sample_data(physeq_count18))

# Adonis test
adonis(bray18 ~ RFTM_score.x * Bucket2, data = sampledf_18) 

# Adonis test
adonis(bray18 ~ RFTM_score.x + Bucket2, data = sampledf_18)

#Dispersion test
beta18 <- betadisper(bray18, sampledf_18$Bucket2)
permutest(beta18)


#Species and Treatments

# Calculate bray curtis distance matrix
bray18_tre <- phyloseq::distance(physeq_count18, method = "bray")

# make a data frame from the sample_data
sampledf_tre <- data.frame(sample_data(physeq_count18))

# Adonis test
adonis(formula = bray18_tre ~ Species2.x * Bucket2, data = sampledf_tre) 


#Dispersion test
beta18 <- betadisper(bray18, sampledf_18$Bucket2)
permutest(beta18)




