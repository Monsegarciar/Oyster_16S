
# Taxonomy Associated with Growth ####
# 2022-09-22
# Author: Monserrat Garcia


# Packages ####
require(microViz)
require(ggplot2)
require(RColorBrewer)
require(dplyr)
require(tidyr)
require(DESeq2)
require(data.table)
library(genefilter)
require(ggplot2)
require(tidyverse)
require(phyloseq)
library(metacoder)

# Filtering out OTU tables 

physeq_count18 <- readRDS("Data/physeq_count18.rds")
physeq_count18

physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

sam18 <- as.data.frame(sample_data(physeq_count18))
sam17 <- as.data.frame(sample_data(physeq_count17))


#2017 
ff <- as.factors(sample_data(physeq_count17)$Site.x)

physeq_count17 = subset_samples(physeq_count17, Weight_delta != "854")

sam17_2 <- as.data.frame(sample_data(physeq_count17))

physeq_count17 = subset_samples(physeq_count17, Volume_scale != "NA")
des_count17 <- phyloseq_to_deseq2(physeq_count17, ~ Volume_scale+Site.x)
des_count17 <- DESeq(des_count17, test="Wald", fitType = "parametric")

results17 <- results(des_count17, name = "Volume_scale")
results17
significant17 <- results17[which(results17$padj <0.05), ]
sigtab17 = cbind(as(significant17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_17_vol <- subset(sigtab17, select = -c(baseMean, lfcSE, stat, pvalue, padj))

###### Log2fold Change ####

#Phylum
x = tapply(sigtab17$log2FoldChange, sigtab17$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17$Phylum = factor(as.character(sigtab17$Phylum), levels=names(x))

# Class
x = tapply(sigtab17$log2FoldChange, sigtab17$Class, function(x) max(x))
x = sort(x, TRUE)
sigtab17$Class = factor(as.character(sigtab17$Class), levels=names(x))


mycolors3= colorRampPalette(brewer.pal(12, "Paired"))(55)
ggplot(sigtab17, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2017") +scale_fill_manual(values = mycolors3) +
  scale_color_manual(values = mycolors3)


#####Bar Plot#####
mycolors1= brewer.pal(11, "Paired")

plot_bar(physeq_v17, x= "Volume_scale", fill= "Order") +
  geom_bar(aes(color=Order, fill = Order), stat = "identity", position = "stack") +
  scale_fill_manual(values = mycolors1) +
  scale_color_manual(values = mycolors1) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_blank(), 
        text = element_text(size=10))

# Turning it into a phyloseq ####
physeq_count17
taxa <- as.matrix(sigtab_17_vol)
taxa_vol <- tax_table(taxa)
view(taxa_vol)
physeq17 = subset_taxa(prune_taxa(rownames(taxa_vol), physeq_count17))
physeq_count17
physeq17

# with site = 3871, and taxa is 10 

physeq_si17 = tax_filter(physeq17, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_si17

sf_17 <- genefilter_sample(physeq17, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq17))
sf_17
# 2193, 6

physeq17_v = prune_taxa(sf_17, physeq17)
physeq17_v
# 8 taxa 

otu <- otu_table(physeq17_v)
otu17 <- otu +1
View(otu17)

tt <- as.data.frame(tax_table(physeq17_v))
view(tt)
ss17 <- merge(tt, sigtab_17_vol, by ='row.names', all = TRUE)
view(ss17)
sss17 <- subset.data.frame(ss17, Kingdom.x != "NA")
view(sss17)
rownames(sss17)= sss17$Row.names
sss17_2 <- subset(sss17, select = -c(Row.names,Kingdom.y, Phylum.y, Class.y, Order.y, Genus.y.y, Species.y, Family.y, Genus.x.y))

log2 <- merge(sss17_2, sigtab_17_vol,by ="row.names")
log2fold <- subset(log2, select = -c(Kingdom.x, Phylum.x, Class.x, Order.x, Genus.y.x, Species.x, Family.x, Genus.x.x, log2FoldChange.y))

write.csv(log2fold, file = "Data/log2fold2017.csv")

###### Turning into Phyloseq Object #####

meta17_data <- read.csv("Data/meta17_data_update.csv")

Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run23_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in "meta_17" data
rownames(meta17_data)= meta17_data$UniqueID
meta17_data$UniqueID=NULL
meta17_data$X=NULL
head(rownames(meta17_data))

#Setting taxmat and otumat
taxmat17=Run123_taxa
taxmat17=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix17=as.matrix(taxmat17, rownames = "V2")
colnames(tax_matrix17) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")

#Setting OTU, TAX, and SAMP
OTU17= otu17

TAX17= tax_table(tax_matrix17)

SAMP17= sample_data(meta17_data)

physeq_v17 = phyloseq(OTU17, TAX17, SAMP17)
physeq_v17

saveRDS(physeq_v17, "Data/physeq_sitexvol17.rds")

taxtable_17 = as.data.frame(tax_table(physeq17_v))

taxtable17 = as.data.frame(tax_table(physeq_v17))

write.csv(taxtable_17, file = "Data/taxtable_sitexvol17.csv")

##### Heatmap #####

plot_heatmap(physeq_v17, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", taxa.label = "Phylum", sample.label = "Volume_scale", sample.order = "Volume_scale")

sub_significant17 <- subset_taxa(prune_taxa(rownames(significant17), physeq_count17))

#### Positive OTUs

significant17_pos <- significant17[significant17$log2FoldChange>0,]
dim(significant17_pos)
dd <- as.data.frame(significant17_pos)
# 25, 6

sub_significant17_pos <- subset_taxa(prune_taxa(rownames(significant17_pos), physeq_count17))
sub_significant17_pos

physeq17_v <- readRDS("Data/physeq_sitexvol17.rds")
physeq17_v
sam_v17 <- as.data.frame(sample_data(physeq17_v))
heatmap17 = parse_phyloseq(physeq17_v)

heatmap17 %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range = c(0.03, 0.04),
            node_color = n_obs,node_size_range = c(0.03, .1),
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")

#2018 Volume and Bucket ####
physeq_count18 = subset_samples(physeq_count18, Volume_scale != "NA")
physeq_count18= subset_samples(physeq_count18, Species2.x !="IR")
des_count18 <- phyloseq_to_deseq2(physeq_count18, ~ Volume_scale + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18 = estimateSizeFactors(des_count18, geoMeans=geoMeans, locfunc=shorth)
deseq18_v = DESeq(deseq18, test="Wald", fitType="parametric")

results18 <- results(deseq18_v, name = "Volume_scale")
significant18 <- results18[which(results18$padj <0.05), ]
sigtab18 = cbind(as(significant18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant18), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_18_vol <- subset(sigtab18, select = -c(baseMean, 
                                              lfcSE, stat, pvalue, padj))


# Turning into a phyloseq object####
physeq_count18
taxa18 <- as.matrix(sigtab_18_vol)
taxa_vol18 <- tax_table(taxa18)
physeq18 = subset_taxa(prune_taxa(rownames(taxa_vol18), physeq_count18)) 
physeq_count18
physeq18

# With site = 2471, 5 taxa  
physeq_si18 = tax_filter(physeq18, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_si18

#physeq_sig18 = genefilter(physeq18, filterfun_sample(function(x) x > 0, A=0.3*nsamples(physeq18)))

sf_18 <- genefilter_sample(physeq18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq18))
# 1277, 6

physeq18_v = prune_taxa(sf_18, physeq18)
physeq18_v
# 5 taxa 

df18 <- as.data.frame(tax_table(physeq18_v))
sv18 <- merge(df18, sigtab_18_vol, by ='row.names', all = TRUE)
svv18 <- subset.data.frame(sv18, Kingdom.x != "NA")
svv18_2 <- subset(svv18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(svv18_2, file = "Data/log2fold2018_volxbuck.csv")

otu18 <- otu_table(physeq_si18)
otu_v18 <- otu18 +1
View(otu_v18)
df <- as.data.frame(otu_table(physeq_count18))

# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run123_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Setting taxmat and otumat
taxmat18=table18
taxmat18=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix18=as.matrix(taxmat18, rownames = "V2")
colnames(tax_matrix18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")

#Setting OTU, TAX, and SAMP
OTU18= otu_v18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

physeq_sig18 = phyloseq(OTU18, TAX18, SAMP18)
physeq_sig18

saveRDS(physeq_sig18, "Data/physeq_volxbuck18.rds")
taxtable_18 <- as.data.frame(tax_table(physeq_si18))

taxtable18 <- as.data.frame(tax_table(physeq_sig18))

write.csv(taxtable18, file = "Data/taxtable_volxbuck18.csv")

# Graphs ####
plot_heatmap(physeq_sig18, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", taxa.label = "Family", sample.label = "Volume_scale", sample.order = "Volume_scale")

ggplot(sigtab18, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2018")

plot_richness(physeq_sig18, x= "Volume_scale", color = "Bucket2", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for Treatment and Species 2017")

# Taking out positive OTU's
neg_otus_weight <- sigtab18 %>% 
  filter(log2FoldChange < 0)
pos_otus_vol <- sigtab18 %>% 
  filter(log2FoldChange > 0)

pos_otus <- subset(pos_otus_vol, select = -c(baseMean, log2FoldChange, 
                                                        lfcSE, stat, pvalue, padj))

write.csv(pos_otus, file = "Data/po_vol18.csv")
pos_otus17 <- read.csv("Data/pos_otus17.csv")



Phy.ord <- ordinate(physeq_sig18, "NMDS", "bray")
plot_ordination(physeq_sig18, Phy.ord, type = "biplot", color = "Volume_scale", shape = "Bucket2", title = "biplot")

physeq18_v <- readRDS("Data/physeq_volxbuck18.rds")
heatmap18_vol = parse_phyloseq(physeq18_v)

heatmap18_vol %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range = c(0.05, 0.06),
            node_color = n_obs,node_size_range = c(0.05, .1), node_size_interval = range(n_obs, na.rm = TRUE, finite = TRUE),
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")

taxtable18 = as.data.frame(tax_table(physeq18_v))

# Weight #####

physeq_weight = subset_samples(physeq_count17, Weight_delta != "NA")
des_weight17 <- phyloseq_to_deseq2(physeq_weight, ~ Weight_delta + Site.x)
des_weight17 <- DESeq(des_weight17, test="Wald", fitType = "parametric")

results_weight17 <- results(des_weight17, name = "Weight_delta")
significant_weight17 <- results_weight17[which(results_weight17$padj <0.05), ]
sigtab_weight17 = cbind(as(significant_weight17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant_weight17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_weigh17 <- subset(sigtab_weight17, select = -c(baseMean, log2FoldChange, 
                                              lfcSE, stat, pvalue, padj))

physeq_count17
taxa_w <- as.matrix(sigtab_weigh17)
taxa_weight <- tax_table(taxa_w)
physeq_weight17 = subset_taxa(prune_taxa(rownames(taxa_weight), physeq_count17))
physeq_count17
physeq_weight17

physeq_filter_weight17 = tax_filter(physeq_weight17, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_filter_weight17
sfweight_17 <- genefilter_sample(physeq_weight17, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_weight17))
# 820, 6

phy_weight17 = prune_taxa(sfweight_17, physeq_weight17)
# 8 taxa 

ftax17 <- as.data.frame(tax_table(phy_weight17))
ft17 <- merge(ftax17, sigtab_weight17, by ='row.names', all = TRUE)
fft17 <- subset.data.frame(ft17, Kingdom.x != "NA")
ff17_2 <- subset(fft17, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x,Genus.y.y,
                                    padj,pvalue,stat,baseMean,lfcSE))
write.csv(ff17_2, file = "Data/log2fold2017_weight_new.csv")
#### *Note: All the taxa were filtered out when filtering for 1/3 of samples and at least appeared once, only one was present for 2017 weight ####

# 2018 

physeq_weight18 = subset_samples(physeq_count18, Weight_delta !="NA")
physeq_weight18 = subset_samples(physeq_count18, Species2.x !="IR")

des_weight18 <- phyloseq_to_deseq2(physeq_weight18, ~ Weight_delta + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_weight18 = estimateSizeFactors(des_weight18, geoMeans=geoMeans, locfunc=shorth)
des_weight18 <- DESeq(deseq_weight18, test="Wald", fitType = "parametric")

results_weight18 <- results(des_weight18, name = "Weight_delta")
significant_weight18 <- results_weight18[which(results_weight18$padj <0.05), ]
sigtab_weight18 = cbind(as(significant_weight18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_weight18), ], "matrix"))

sigtab_18_weight <- subset(sigtab_weight18, select = -c(baseMean, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_w18 <- as.matrix(sigtab_18_weight)
taxa_weight18 <- tax_table(taxa_w18)
physeq_weight18 = subset_taxa(prune_taxa(rownames(taxa_weight18), physeq_count18))
physeq_count18
physeq_weight18

df <- as.data.frame(tax_table(phy_weight18))
ss18 <- merge(df, sigtab_18_weight, by ='row.names', all = TRUE)
sss18 <- subset.data.frame(ss18, Kingdom.x != "NA")
sss18_2 <- subset(sss18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(sss18_2, file = "Data/log2fold2018_weight.csv")

# Weight and bucket= 2730, 6 taxa 
phys_weight18 = tax_filter(physeq_weight18, min_prevalence = 0.3, min_sample_abundance = 1)
phys_weight18

sfweight_18 <- genefilter_sample(physeq_weight18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_weight18))
# 1302, 6

phy_weight18 = prune_taxa(sfweight_18, physeq_weight18)
phy_weight18
#7, 8

otu_weight18 <- otu_table(phy_weight18)
otu_w18 <- otu_weight18 +1
View(otu_w18)

#Setting OTU, TAX, and SAMP
OTU18= otu_w18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

phy_weight_sig18 = phyloseq(OTU18, TAX18, SAMP18)
phy_weight_sig18

saveRDS(phy_weight_sig18, "Data/physeq_buckxweight18.rds")

physeq_weight_sig18

taxtable_weight18 <- as.data.frame(tax_table(phy_weight_sig18))

taxtable_BuckxWeigh <- as.data.frame(tax_table(phy_weight_sig18))

write.csv(taxtable_BuckxWeigh, file = "Data/taxtable_Bucket&Weight.csv")

# Plot Analysis ####

plot_heatmap(phy_weight_sig18, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", na.values = "white", taxa.label = "Family", sample.label = "Volume_scale", sample.order = "Volume_scale")

plot_richness(phy_weight_sig18, x= "Volume_scale", color = "Bucket2", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for Treatment and Species 2017")

ggplot(phy_weight_sig18, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2018")


mycolors2= colorRampPalette(brewer.pal(8, "Dark2"))(6)
plot_bar(phy_weight_sig18, x= "Volume_scale", fill= "Class")+
  geom_bar(aes(color=Class, fill = Class), stat = "identity", position = "stack") + facet_wrap(vars(Bucket2), scales = 'free_x')
  scale_fill_manual(values = mycolors2) + 
  scale_color_manual(values = mycolors2) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_text(), 
        panel.grid.minor.x = element_text(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_text(), 
        text = element_text(size=10))

physeq_buckxweight18 <- readRDS("Data/physeq_buckxweight18.rds")
physeq_buckxweight18
df = sample_data(physeq_buckxweight18)
buck <- as.data.frame(tax_table(physeq_buckxweight18))

?theme
plot_bar(phys, "Order", fill = "Phylum", facet_grid = ~Description) +
  ylab("Percentage of Sequences") 

heatmap = parse_phyloseq(physeq_buckxweight18)
samda18 <- as.data.frame(sample_data(physeq_buckxweight18))
phy_weight_sig18

ggplot(sigtab_weight18, aes(x=Phylum, y=rownames(taxtable_BuckxWeigh), color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2018")

heatmap %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range = c(0.03, 0.04),
            node_color = n_obs,node_size_range = c(0.03, .1),
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")
  

# Merging Tax Tables with significant OTUs for Weight and Volume ####

#### Loading Tax Tables ####
taxtable_BuckxWeigh <- read.csv("Data/taxtable_Bucket&Weight.csv")
taxtable_sitexvol17 <- read.csv("Data/taxtable_sitexvol17.csv")
tabtable_volxbuck18 <- read.csv("Data/taxtable_volxbuck18.csv")

#### Combining Tax Tables ####
sig_OTU_combined <- rbind(taxtable_BuckxWeigh, taxtable_sitexvol17, tabtable_volxbuck18)

write.csv(sig_OTU_combined, file = "Data/sig_OTU_combined.csv")



# Non-scaled Volume ####

physeq_count17 = subset_samples(physeq_count17, Volume_delta != "NA")
des_volume17 <- phyloseq_to_deseq2(physeq_count17, ~ Volume_delta+Site.x)
des_volume17 <- DESeq(des_volume17, test="Wald", fitType = "parametric")

results_volume17 <- results(des_volume17, name = "Volume_delta")
results_volume17
significant_volume17 <- results_volume17[which(results_volume17$padj <0.05), ]
sigtab_volume17 = cbind(as(significant_volume17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant_volume17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_Volume17 <- subset(sigtab_volume17, select = -c(baseMean, log2FoldChange, 
                                              lfcSE, stat, pvalue, padj))

# Turning it into a phyloseq ####
physeq_count17
taxa_volume17 <- as.matrix(sigtab_Volume17)
taxa_vol17 <- tax_table(taxa_volume17)
physeq_volume17 = subset_taxa(prune_taxa(rownames(taxa_vol17), physeq_count17))
physeq_count17
physeq_volume17

#Look which sites are more prominent or appear more in the significant OTUs#

physeq_sigvol17 = tax_filter(physeq_volume17, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_sigvol17

sf_vol17 <- genefilter_sample(physeq_volume17, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_volume17))
# 2193, 6

physeq17_volume = prune_taxa(sf_vol17, physeq_volume17)
physeq17_volume

# *Note: All 2017 filtered out with the filter settings above#

physeq_volume18 = subset_samples(physeq_count18, Volume_delta !="NA")

des_volume18 <- phyloseq_to_deseq2(physeq_volume18, ~ Volume_delta + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_volume18 = estimateSizeFactors(des_volume18, geoMeans=geoMeans, locfunc=shorth)
des_volume18 <- DESeq(deseq_volume18, test="Wald", fitType = "parametric")

results_volume18 <- results(des_volume18, name = "Volume_delta")
significant_volume18 <- results_volume18[which(results_volume18$padj <0.05), ]
sigtab_volume18 = cbind(as(significant_volume18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_volume18), ], "matrix"))

sigtab_18_volume <- subset(sigtab_volume18, select = -c(baseMean, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_v18 <- as.matrix(sigtab_18_volume)
taxa_volume18 <- tax_table(taxa_v18)
physeq_volume18 = subset_taxa(prune_taxa(rownames(taxa_volume18), physeq_count18))
physeq_count18
physeq_volume18


vf <- as.data.frame(tax_table(physeq18_volume))
vv18 <- merge(vf, sigtab_volume18, by ='row.names', all = TRUE)
vvv18 <- subset.data.frame(vv18, Kingdom.x != "NA")
vvv18_2 <- subset(vvv18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.y.x,Genus.x.x, Genus.y.y,Species.x,baseMean,lfcSE, stat,pvalue,padj))
write.csv(vvv18_2, file = "Data/log2fold2018_voldelta.csv")



physeq_sigvol18 = tax_filter(physeq_volume18, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_sigvol18

sf_vol18 <- genefilter_sample(physeq_volume18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_volume18))
# 3 taxa

physeq18_volume = prune_taxa(sf_vol18, physeq_volume18)
physeq18_volume

taxtable_voldelxbuck18 <- as.data.frame(tax_table(physeq_sigvol18))


heatmap18 = parse_phyloseq(physeq_sigvol18)

heatmap18 = parse_phyloseq(as.data.frame(log2fold18_weight))

heatmap18 %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range=c(0.02,0.02),
            node_color = n_obs,
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")
write.csv(taxtable_voldelxbuck18, file = "Data/taxtable_voldelxbuck18.csv")

taxtable_voldelxbuck18 <- read.csv("Data/taxtable_voldelxbuck18.csv")


# Sites in the Significant OTUs ####

physeq_buckxweight18 <- readRDS("Data/physeq_buckxweight18.rds")
site_data18 <- sample_data(physeq_buckxweight18)


physeq_sitexvol17 <- readRDS("Data/physeq_sitexvol17.rds")
site_data17 <- sample_data(physeq_sitexvol17)

physeq_volxbuck18 <- readRDS("Data/physeq_volxbuck18.rds")

##### Log2Fold Change Graphs ####
log2fold18_vol <- read.csv("Data/log2fold2018_volxbuck.csv")
log2fold18_vol

log2fold18_weight <- read.csv("Data/log2fold2018_weight.csv")
log2fold2017 <- read.csv("Data/log2fold2017.csv")


theme_set(theme_bw())
scale_fill_discrete <- function(palname = "Set1", ...) {
  scale_fill_brewer(palette = palname, ...)
}

##### 2018 Volume and Bucket ####
x = tapply(log2fold18_vol$log2FoldChange, log2fold18_vol$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold18_vol$Genus = factor(as.character(log2fold18_vol$Genus), levels=names(x))

Genus_vol18 <- length(unique(log2fold18_vol$Genus))
Genus_pal_v18 <- colorRampPalette(brewer.pal(5,"Set1"))

ggplot(log2fold18_vol, aes(x=Row.names, y=log2FoldChange, fill=Genus)) + theme_classic() +geom_hline(yintercept = 0)+
  theme(axis.text.x = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.border=element_rect(color = "black",fill = NA, size=1),
        axis.title = element_text(face = "bold"),text = element_text(size = 14) )+ 
  ylab("Log2FoldChange")+xlab("Genus")+geom_bar(stat="identity")+
  scale_x_discrete(labels= log2fold18_vol$Genus)+
  geom_text(aes(label=sprintf(log2FoldChange,fmt = "%0.2f"), fontface="bold"), vjust=-0.5, color="black",size=4.0)+
  scale_fill_manual(values = Genus_pal_v18(Genus_vol18))

##### 2018 Weight ####

x = tapply(log2fold18_weight$log2FoldChange, log2fold18_weight$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold18_weight$Genus = factor(as.character(log2fold18_weight$Genus), levels=names(x))

Genus_weigh18 <- length(unique(log2fold18_weight$Genus))
Genus_pal_w18 <- colorRampPalette(brewer.pal(7,"Set1"))

ggplot(log2fold18_weight, aes(x=Row.names, y=log2FoldChange,fill=Genus))+
  theme_classic() +geom_hline(yintercept = 0)+
  theme(axis.text.x = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border=element_rect(color = "black",fill = NA, size=0.5),
        axis.title = element_text(face = "bold"),text = element_text(size = 14))+
  ylab("Log2FoldChange")+
  scale_x_discrete(labels= log2fold18_weight$Genus)+
  geom_bar(stat="identity")+xlab("Genus")+
  geom_text(aes(label=sprintf(log2FoldChange,fmt = "%0.2f"), fontface="bold"),vjust=-0.5, color="black",size=3.5)+
  scale_fill_manual(values = Genus_pal_w18(Genus_weigh18))


##### 2017 ####

x = tapply(log2fold2017$log2FoldChange.x, log2fold2017$Genus, function(x) max(x))
x = sort(x, TRUE)
log2fold2017$Genus = factor(as.character(log2fold2017$Genus), levels=names(x))

Genus_vol17 <- length(unique(log2fold2017$Genus))
Genus_pal <- colorRampPalette(brewer.pal(7,"Set1"))

ggplot(log2fold2017, aes(x=Row.names, y=log2FoldChange.x,fill=Genus))+geom_bar(stat="identity")+  theme_classic() +
  theme(axis.text.x = element_blank(),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.border=element_rect(color = "black",fill = NA, size=1),
        axis.title = element_text(face = "bold"),text = element_text(size = 14) )+ 
  ylab("Log2FoldChange")+xlab("Genus")+
 scale_x_discrete(labels= log2fold2017$Genus) + 
  geom_bar(stat="identity")+
  geom_text(aes(label=sprintf(log2FoldChange.x,fmt = "%0.2f"), fontface="bold"), vjust=-0.5, color="black",size=4.0)+
  scale_fill_manual(values = Genus_pal(Genus_vol17))

#https://www.datanovia.com/en/blog/ggplot-theme-background-color-and-grids/#:~:text=To%20remove%20a%20particular%20panel,grid.


# Mussels for 2018 #### 

physeq_count18 = subset_samples(physeq_count18, Volume_scale != "NA")
physeq_count18= subset_samples(physeq_count18, Species2.x !="CV")
des_count18 <- phyloseq_to_deseq2(physeq_count18, ~ Volume_scale + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18 = estimateSizeFactors(des_count18, geoMeans=geoMeans, locfunc=shorth)
deseq18_v = DESeq(deseq18, test="Wald", fitType="parametric")

results18 <- results(deseq18_v, name = "Volume_scale")
significant18 <- results18[which(results18$padj <0.05), ]
sigtab18 = cbind(as(significant18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant18), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_18_muss <- subset(sigtab18, select = -c(baseMean, 
                                              lfcSE, stat, pvalue, padj))


# Turning it into a phyloseq ####
physeq_count18
taxa18 <- as.matrix(sigtab_18_muss)
taxa_muss18 <- tax_table(taxa18)
physeq18 = subset_taxa(prune_taxa(rownames(taxa_muss18), physeq_count18)) 
physeq_count18
physeq18

  
physeq_muss18 = tax_filter(physeq18, min_prevalence = 0.3, min_sample_abundance = 1)
physeq_muss18

#physeq_sig18 = genefilter(physeq18, filterfun_sample(function(x) x > 0, A=0.3*nsamples(physeq18)))

muss_18 <- genefilter_sample(physeq18, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq18))


physeq18_mu = prune_taxa(muss_18, physeq18)
physeq18_mu
# 6 taxa 

mu18 <- as.data.frame(tax_table(physeq18_mu))
mv18 <- merge(mu18, sigtab_18_muss, by ='row.names', all = TRUE)
mvv18 <- subset.data.frame(mv18, Kingdom.x != "NA")
mvv18_2 <- subset(mvv18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(mvv18_2, file = "Data/log2fold2018_volxmussels.csv")

otu_mu18 <- otu_table(physeq18_mu)
otu_vmu18 <- otu_mu18 +1
View(otu_vmu18)

# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run123_taxa"
Run123_taxa$V1=NULL
rownames(Run123_taxa)= Run123_taxa$V2
head(rownames(Run123_taxa))

#Changing row names in meta_gen18 data
rownames(meta_gen18_data)= meta_gen18_data$UniqueID 
head(rownames(meta_gen18_data))

#Setting taxmat and otumat
taxmat18=Run123_taxa
taxmat18=Run123_taxa[-c(1)]

#Converting to matrix
tax_matrix18=as.matrix(taxmat18, rownames = "V2")
colnames(tax_matrix18) <- c("Kingdom", "Phylum", "Class", "Order", "Family", 
                            "Genus.x", "Genus.y", "Species")

#Setting OTU, TAX, and SAMP
OTU18= otu_vmu18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

physeq_mus18 = phyloseq(OTU18, TAX18, SAMP18)
physeq_mus18

saveRDS(physeq_mus18, "Data/physeq_volxmussels18.rds")

taxtablemuss_18 <- as.data.frame(tax_table(physeq_mus18))

taxtablemuss18 <- as.data.frame(tax_table(physeq_mus18))

write.csv(taxtablemuss18, file = "Data/taxtable_volxmussels18.csv")


# 2018 Mussels Graph and Heatmap ####
# Loading Data #

log2fold18_volxmuss <- read.csv("Data/log2fold2018_volxmussels.csv")

Genus_volmuss18 <- length(unique(log2fold18_volxmuss$Genus))
Genus_pal_vm18 <- colorRampPalette(brewer.pal(7,"Set1"))

ggplot(log2fold18_volxmuss, aes(x=Row.names, y=log2FoldChange,fill=Genus))+theme_classic() +geom_hline(yintercept = 0)+
  theme(axis.text.x = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border=element_rect(color = "black",fill = NA, size=0.5),
        axis.title = element_text(face = "bold"),text = element_text(size = 14))+
  ylab("Log2FoldChange")+
  scale_x_discrete(labels= log2fold18_volxmuss$Genus)+
  geom_bar(stat="identity")+xlab("OTU")+
  geom_text(aes(label=sprintf(log2FoldChange,fmt = "%0.2f"), fontface="bold"),vjust=1.3, color="black",size=3.5)+
  scale_fill_manual(values = Genus_pal_vm18(Genus_volmuss18))

heatmap_muss = parse_phyloseq(physeq_mus18)

heatmap_muss %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range = c(0.03, 0.04),
            node_color = n_obs,node_size_range = c(0.03, .1),
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")

#### Weight #####
physeq_weight18 = subset_samples(physeq_count18, Weight_delta !="NA")
physeq_weight18 = subset_samples(physeq_weight18, Species2.x !="CV")


des_weight18 <- phyloseq_to_deseq2(physeq_weight18, ~ Weight_delta + Bucket2)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_weight18 = estimateSizeFactors(des_weight18, geoMeans=geoMeans, locfunc=shorth)
des_weight18 <- DESeq(deseq_weight18, test="Wald", fitType = "parametric")

results_weight18 <- results(des_weight18, name = "Weight_delta")
significant_weight18 <- results_weight18[which(results_weight18$padj <0.05), ]
sigtab_weight18 = cbind(as(significant_weight18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_weight18), ], "matrix"))

sigtab_muss_weight <- subset(sigtab_weight18, select = -c(baseMean, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_wm <- as.matrix(sigtab_muss_weight)
taxa_weightmuss <- tax_table(taxa_wm)
physeq_weightmuss = subset_taxa(prune_taxa(rownames(taxa_weightmuss), physeq_count18))
physeq_count18
physeq_weightmuss

wm <- as.data.frame(tax_table(physeq_weight_muss))
wm18 <- merge(wm, sigtab_muss_weight, by ='row.names', all = TRUE)
wmm18 <- subset.data.frame(wm18, Kingdom.x != "NA")
wmmm18_2 <- subset(wmm18, select = -c(Kingdom.x,Phylum.x, Class.x, Order.x, Family.x, Genus.x.x, Genus.y.x,Species.x))
write.csv(wmmm18_2, file = "Data/log2fold2018_mussels.csv")

# Weight and bucket= 2730, 6 taxa 
phys_weightmuss = tax_filter(physeq_weightmuss, min_prevalence = 0.3, min_sample_abundance = 1)
phys_weightmuss

muweight_18 <- genefilter_sample(physeq_weightmuss, filterfun_sample(function(x) x > 0), A=0.3*nsamples(physeq_weightmuss))
# 1302, 6

phy_weightmu = prune_taxa(muweight_18, physeq_weightmuss)
phy_weightmu
#7, 8

otu_weight18 <- otu_table(phy_weightmu)
otu_w18 <- otu_weight18 +1
View(otu_w18)

#Setting OTU, TAX, and SAMP
OTU18= otu_w18

TAX18= tax_table(tax_matrix18)

SAMP18= sample_data(meta_gen18_data)

phy_weight_muss = phyloseq(OTU18, TAX18, SAMP18)
phy_weight_muss

saveRDS(phy_weight_muss, "Data/physeq_musselsxweight18.rds")

physeq_weight_muss<- readRDS("Data/physeq_musselsxweight18.rds")

taxtable_weightmuss <- as.data.frame(tax_table(phy_weight_muss))

taxtable_mussxWeigh <- as.data.frame(tax_table(phy_weight_muss))

write.csv(taxtable_BuckxWeigh, file = "Data/taxtable_Bucket&Weight.csv")


#### Graphs for Weight #####

log2fold_mussels <- read.csv("Data/log2fold2018_mussels.csv")

Genus_weighmuss <- length(unique(log2fold_mussels$Genus))
Genus_pal_wmuss <- colorRampPalette(brewer.pal(7,"Set1"))

ggplot(log2fold_mussels, aes(x=Row.names, y=log2FoldChange,fill=Genus))+theme_classic() +geom_hline(yintercept = 0)+
  theme(axis.text.x = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),panel.border=element_rect(color = "black",fill = NA, size=0.5),
        axis.title = element_text(face = "bold"),text = element_text(size = 14))+
  ylab("Log2FoldChange")+
  scale_x_discrete(labels= log2fold_mussels$Genus)+
  geom_bar(stat="identity")+xlab("OTU")+
  geom_text(aes(label=sprintf(log2FoldChange,fmt = "%0.2f"), fontface="bold"),vjust=-0.5, color="black",size=3.5)+
  scale_fill_manual(values = Genus_pal_wmuss(Genus_weighmuss))


heatmap_muss = parse_phyloseq(physeq_weight_muss)

heatmap_muss %>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,node_label_size_range = c(0.03, 0.04),
            node_color = n_obs,node_size_range = c(0.03, .1),
            layout = "davidson-harel", initial_layout = "reingold-tilford", node_color_axis_label = "Number of Obs")



