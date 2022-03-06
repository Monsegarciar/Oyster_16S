# Phyloseq and Deseq Analysis- Measurements ####
# 2022-03-01
# Author: Monse Garcia



library(DESeq2)

# Loading Data
physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

physeq_class17 <- readRDS("Data/physeq_class17.rds")
physeq_class17

# Normalized Volume with DESEQ2- 2017 Data
physeq_count17 = subset_samples(physeq_count17, Volume_delta != "NA")

deseq17_volume = phyloseq_to_deseq2(physeq_count17, ~ Volume_delta) #the design formula contains one or more numeric variables that have mean or standard deviation larger than 5 (an arbitrary threshold to trigger this message). it is generally a good idea to center and scale numeric variables in the design to improve GLM convergence.
deseq17_volume = DESeq(deseq17_volume, test="Wald", fitType="parametric")

res17_volume = results(deseq17_volume, cooksCutoff = FALSE)
alpha = 0.05 
sigtab17_volume = res17_volume[which(res17_volume$padj < alpha), ]
sigtab17_volume = cbind(as(sigtab17_volume, "data.frame"), as(tax_table(physeq_count17)[rownames(sigtab17_volume), ], "matrix"))

