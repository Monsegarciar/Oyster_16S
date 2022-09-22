# Filtering out OTU tables 

physeq_count18 <- readRDS("Data/physeq_count18.rds")

physeq_count17 <- readRDS("Data/physeq_count17.rds")
physeq_count17

#2017 

physeq_count17 = subset_samples(physeq_count17, Volume_scale != "NA")
des_count17 <- phyloseq_to_deseq2(physeq_count17, ~ Volume_scale)
des_count17 <- DESeq(des_count17, test="Wald", fitType = "parametric")

results17 <- results(des_count17, name = "Volume_scale")
significant17 <- results17[which(results17$padj <0.05), ]
sigtab17 = cbind(as(significant17, "data.frame"), as(tax_table(physeq_count17)[rownames(significant17), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_17_vol <- subset(sigtab17, select = -c(baseMean, log2FoldChange, 
                                                        lfcSE, stat, pvalue, padj))

#Phylum
x = tapply(sigtab17$log2FoldChange, sigtab17$Phylum, function(x) max(x))
x = sort(x, TRUE)
sigtab17$Phylum = factor(as.character(sigtab17$Phylum), levels=names(x))

ggplot(sigtab17, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2017")


#Bar Plot
mycolors2= colorRampPalette(brewer.pal(8, "Dark2"))(75)

plot_bar(physeq_sig18, x= "Volume_scale", fill= "Order")+
  geom_bar(aes(color=Order, fill = Order), stat = "identity", position = "stack") + facet_wrap(vars(Volume_scale), scales = 'free_x') +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors2) +
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
physeq17 = subset_taxa(prune_taxa(rownames(taxa_vol), physeq_count17))
physeq_count17
physeq17

physeq_si17 = tax_filter(physeq17, min_prevalence = 0.03, min_sample_abundance = 1)
physeq_si18
sf_17 <- genefilter_sample(physeq17, filterfun_sample(function(x) x > 0), A=0.5*nsamples(physeq17))
# 2193, 6

physeq17_v = prune_taxa(sf_17, physeq17)
physeq17_v
# 8 taxa 

otu <- otu_table(physeq17_v)
otu17 <- otu +1
View(otu17)

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

##### Heatmap #####

plot_heatmap(physeq_v17, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", taxa.label = "Phylum", sample.label = "Volume_scale", sample.order = "Volume_scale")


sub_significant17 <- subset_taxa(prune_taxa(rownames(significant17), physeq_count17))

#### Positive OTUs

significant17_pos <- significant17[significant17$log2FoldChange>0,]
dim(significant17_pos)
# 25, 6

sub_significant17_pos <- subset_taxa(prune_taxa(rownames(significant17_pos), physeq_count17))
sub_significant17_pos


#2018 

physeq_count18 = subset_samples(physeq_count18, Volume_scale != "NA")
des_count18 <- phyloseq_to_deseq2(physeq_count18, ~ Volume_scale)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq18 = estimateSizeFactors(des_count18, geoMeans=geoMeans, locfunc=shorth)
deseq18_v = DESeq(deseq18, test="Wald", fitType="parametric")

results18 <- results(deseq18_v, name = "Volume_scale")
significant18 <- results18[which(results18$padj <0.05), ]
sigtab18 = cbind(as(significant18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant18), ], "matrix"))

# Creating DESEq results into a tax table ####
sigtab_18_vol <- subset(sigtab18, select = -c(baseMean, log2FoldChange, 
                                              lfcSE, stat, pvalue, padj))

# Turning it into a phyloseq ####
physeq_count18
taxa18 <- as.matrix(sigtab_18_vol)
taxa_vol18 <- tax_table(taxa18)
physeq18 = subset_taxa(prune_taxa(rownames(taxa_vol18)), physeq_count18) 
physeq_count18
physeq18

physeq_si18 = tax_filter(physeq18, min_prevalence = 0.33, min_sample_abundance = 1)
physeq_si18

physeq_sig18 = genefilter(physeq18, filterfun_sample(function(x) x > 0, A=0.3*nsamples(physeq18)))

sf_18 <- genefilter_sample(physeq18, filter(function(x) x > 0), A=0.3*nsamples(physeq18))
# 1277, 6

physeq18_v = prune_taxa(sf_18, physeq18)
physeq18_v
# 4 taxa 

otu18 <- otu_table(physeq_si18)
otu_v18 <- otu18 +1
View(otu_v18)
df <- as.data.frame(otu_table(physeq_count18))

# Loading Data
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")

Run123_taxa <- fread("Data/Run123_taxa_complete - Copy.csv")

#Changing row names in "Run123_taxa"
Run123_taxa$V1=NULL
rownames(table18)= Run123_taxa$V2
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

physeq_si18

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



# Weight #####
sd <- sample_data(physeq_count17)

physeq_weight = subset_samples(physeq_count17, Weight_delta != "NA")
des_weight17 <- phyloseq_to_deseq2(physeq_weight, ~ Weight_delta)
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



# 2018 

physeq_weight18 = subset_samples(physeq_count18, Weight_delta !="NA")

des_weight18 <- phyloseq_to_deseq2(physeq_weight18, ~ Weight_delta)

gm_mean = function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0])))
geoMeans = apply(OTU_count18, 2, gm_mean)
deseq_weight18 = estimateSizeFactors(des_weight18, geoMeans=geoMeans, locfunc=shorth)
des_weight18 <- DESeq(deseq_weight18, test="Wald", fitType = "parametric")

results_weight18 <- results(des_weight18, name = "Weight_delta")
significant_weight18 <- results_weight18[which(results_weight18$padj <0.05), ]
sigtab_weight18 = cbind(as(significant_weight18, "data.frame"), as(tax_table(physeq_count18)[rownames(significant_weight18), ], "matrix"))

sigtab_18_weight <- subset(sigtab_weight18, select = -c(baseMean, log2FoldChange, 
                                                        lfcSE, stat, pvalue, padj))

physeq_count18
taxa_w18 <- as.matrix(sigtab_18_weight)
taxa_weight18 <- tax_table(taxa_w18)
physeq_weight18 = subset_taxa(prune_taxa(rownames(taxa_weight18), physeq_count18))
physeq_count18
physeq_weight18

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

physeq_weight_sig18


plot_heatmap(phy_weight_sig18, method = "NMDS", distance = "bray",low = "#FFFFFF", high ="#FF3300", na.values = "white", taxa.label = "Family", sample.label = "Volume_scale", sample.order = "Volume_scale")

plot_richness(phy_weight_sig18, x= "Volume_scale", color = "Bucket2", measures = c("Simpson", "Shannon"), title = "Alpha Diversity for Treatment and Species 2017")

ggplot(sigtab_weight18, aes(x=Class, y=log2FoldChange, color=Phylum)) + geom_point(size=6) + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5)) + labs(title = "2018")

plot_bar(phy_weight_sig18, x= "Volume_scale", fill= "Order")+
  geom_bar(aes(color=Order, fill = Order), stat = "identity", position = "stack") + facet_wrap(vars(Volume_scale), scales = 'free_x') +
  scale_fill_manual(values = mycolors2) +
  scale_color_manual(values = mycolors2) +
  theme_bw() +
  theme(legend.position = "right", panel.border = element_blank(), 
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.line = element_line(color = "black"), 
        axis.text.x = element_blank(), 
        text = element_text(size=10))
