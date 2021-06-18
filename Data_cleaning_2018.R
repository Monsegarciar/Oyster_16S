# 2018 Data Cleaning ####
# 2021-06-14
# Author: Monse Garcia 

# Changing names to new species abbreviations

meta18<-read.csv("C:/Users/monse/OneDrive/Documents/Oyster_16S/Data/metadata_de18 - Copy.csv")

select(meta18, Species)  

#using ifelse for the new species abbreviations

meta18$Species2<- ifelse(meta18$Species== "MB", "LP", 
                         ifelse(meta18$Species== "IR", "IR", 
                                ifelse(meta18$Species=="CV","CV", "AM")))

#uploading the asvtable_18
asvtable_18 <- fread("Data/asvtable_de18 - Copy.csv")
