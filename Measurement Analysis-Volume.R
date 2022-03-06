# Measurment Analysis- Volume ####
# 2022-02-17
# Author: Monse Garcia

#Loading Packages

require(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)

#Loading Data

meta17_data <- read.csv("Data/meta17_data_update.csv")

meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

### 2017 Data

# Adding new column for volume 

height_pre <- as.numeric(meta17_data$Height_pre)
length_pre <- as.numeric(meta17_data$Length_pre)
width_pre <- as.numeric(meta17_data$Width_pre)

meta17_data <- meta17_data %>%
  mutate(Volume_pre= height_pre * length_pre * width_pre)

meta17_data <- meta17_data %>%
  mutate(Volume_post= Height_post * Length_post * Width_post)

meta17_data <- meta17_data %>%
  mutate(Volume_delta= (Volume_post- Volume_pre)/ Volume_pre)

#Saving new data
write.csv(meta17_data, file = "Data/meta17_data_update.csv")

#### 2018 Data

# Adding new column for volume 

meta_gen18_data <- meta_gen18_data %>%
  mutate(Volume_pre= Length_pre * Width_pre * Height_pre)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Volume_post= Length_post * Width_post * Height_post)

meta_gen18_data <- meta_gen18_data %>%
  mutate(Volume_delta= (Volume_post - Volume_pre)/ Volume_pre)

#Saving new data
write.csv(meta_gen18_data, file = "Data/metagenetics_data18.csv")     



### Phyloseq Analysis

# 2017

Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")
asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")

#Changing row names in "Run23_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Changing row names in "asvtable_17" data
rownames(asvtable_17)= asvtable_17$V1
asvtable_17$V1=NULL
head(rownames(asvtable_17))

#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]
otumat17=asvtable_17

#Converting to matrix
otu_matrix17= as.matrix(otumat17, rownames = rownames(asvtable_17))

tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")
meta17_data=as.data.frame(meta17_data)

#Setting OTU, TAX, and SAMP
OTU17= otu_table(otu_matrix17, taxa_are_rows = FALSE)

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)


OTU_count17=transform_sample_counts(OTU17, function(x) 1E6 * x/sum(x))


physeq_class17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_class17

physeq_count17 = phyloseq(OTU_count17, TAX17, SAMP17)
physeq_count17


# Saving Physeq as an RDS
saveRDS(physeq_class17, "Data/physeq_class17.rds")
physeq_class17 <- readRDS("Data/physeq_class17.rds")

saveRDS(physeq_count17, "Data/physeq_count17.rds")
physeq_count17 <- readRDS("Data/physeq_count17.rds")




