# 2018 Data Cleaning ####
# 2021-06-14
# Author: Monse Garcia 

#Loading Packages needed
require(phyloseq, dplyr, tidyr,stringr)
library(dplyr)
library(tidyr)
library(stringr)
# Data 

meta18<-read.csv("C:/Users/monse/OneDrive/Documents/Oyster_16S/Data/metadata_de18 - Copy.csv")
genetics_data2018 <- read.csv("Data/DE2018_alldata - Copy.csv")

#using ifelse for the new species abbreviations

meta18$Species2<- ifelse(meta18$Species== "MB", "LP", 
                         ifelse(meta18$Species== "IR", "IR", 
                                ifelse(meta18$Species=="CV","CV", "AM")))

genetics_data2018$Species2<- ifelse(genetics_data2018$Species== "MB", "LP", 
                         ifelse(meta18$Species== "IR", "IR", 
                                ifelse(meta18$Species=="CV","CV", "AM")))

#Taking out unnecessary columns for "gentics_data2018"
genetics_data2018_clean <- subset(genetics_data2018, select = -c(Species, Date_initial_measure, Mortality_Date, 
                                                                 Date_FinalMeasurement, X, RFTM_Date, Parasites))

# Taking out unnecessary columns for "meta18"
meta18_data<- subset(meta18, select = -c(X, V1, Site, Year, Species, Phase_1_DO, 
                                         Phase_1_temp, Phase_2_DO, Phase_2_Temp, 
                                         Overall_treatment))

#Creating UniqueID's

str_split_fixed(genetics_data2018_clean$Bucket, "- ", 2)

genetics_data2018_clean$Bucket2<-ifelse(data$Treatment=="HH", "HIGH_POLY", ifelse(data$Treatment=="HL", "HIGH_MONO", ifelse(data$Treatment== "LL", "LOW_MONO","LOW_POLY")))



# Adding Na's in blank columns
genetics_data2018_clean[genetics_data2018_clean == "" | genetics_data2018_clean == " "] <- NA

# Merging the datasets 
genetics_data2018_clean$Buck_color_num <- paste(genetics_data2018_clean$Bucket, genetics_data2018_clean$Color.Number, sep = "")
meta18_data$Buck_color_num <- paste(meta18_data$Color_Bucket,meta18_data$Number, sep = "")

meta_gen18_data <- merge(meta18_data, genetics_data2018_clean, by= "Buck_color_num", all.x = TRUE)

#Saving the new data
write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")

#uploading the asvtable_18
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")


#Phyloseq analysis
OTU2=asvtable_18<- fread("Data/asvtable_de18 - Copy.csv")
TAX2=Run23_taxa <- fread("Data/Run23_taxa - Copy.csv")

rownames(meta_gen18_data)= meta_gen18_data$UniqueID
physeq = phyloseq(OTU, TAX, meta_gen18_data)
physeq







