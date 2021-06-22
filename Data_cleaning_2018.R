# 2018 Data Cleaning ####
# 2021-06-14
# Author: Monse Garcia 

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

# Adding Na's in blank columns
genetics_data2018_clean[genetics_data2018_clean == "" | genetics_data2018_clean == " "] <- NA



#uploading the asvtable_18
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")
