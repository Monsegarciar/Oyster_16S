# Questions for Analysis ####
# 2021-06-23
# Author: Monse Garcia

#Packages Required
library(phyloseq)
library(ggplot2)
require(RColorBrewer)

mycolors= colorRampPalette(brewer.pal(8, "Dark2"))(2) # How many colors
plot_bar(metaSF,  fill="Category", x="Replicate") +
  geom_bar(aes(color=Category, fill=Category), stat="identity", position="stack")+
  facet_grid(Year~Site_Name, scales="free_x")+
  scale_fill_manual(values=mycolors)+
  scale_color_manual(values=mycolors)+
  theme_bw()+
  theme(legend.position="top", legend.text=element_text(size=10), panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(),
        axis.ticks.x=element_blank(), axis.line=element_line(color="black"),
        text = element_text(size=10))
#want a specific color with hexcodes or do red and blue or other colors instead of mycolors in the code
#search r color paletes different kinds 

#Question 1 ####
#Sites relationship with oyster species
#Example: Which sites are more prominent to bacteria or do more bacteria diversity in a specific sites. 
#Alpha diversity 
#Diversity with other factors
#Taxonomy bar to see if there are more abundance in one site than the other

#Split Graph with Site and Phylum
p4= plot_ordination(physeq_class, Phy.ord, type = "split", 
                    color = "Phylum", shape = "Site.x")
print(p4)
##Plot bars with phylum 
table_taxa <-table(taxmat$V3) # Used the table() function to see which phylum was most common
View(table_taxa)

#Most common phylum and site graph
pp.ch= subset_taxa(physeq_class, Phylum=="Proteobacteria") 
plot_bar(pp.ch) #Plot bar of samples(x) and abundance (y) of Proteobacteria
bar1=plot_bar(pp.ch, x="Site.x", fill = "Genus")
print(bar1)



b2=plot_bar(pp.ch, "Site.x", fill="Genus", facet_grid=~Family) 
print(b2)
b2 + geom_point(aes(x=Family, y=Abundance), color="black", position="jitter", size=3)

## Alpha Diversity Graphics

PC <- prune_species(speciesSums(physeq_class)> 0, physeq_class)
PC
plot_richness(PC) #https://github.com/joey711/phyloseq/issues/552 for troubleshooting

plot_richness(PC, measures = c("Site"))

OTU2= otu_table(otu_matrix, taxa_are_rows = FALSE)
OTU2= transform_sample_counts(OTU, as.integer)

physeq_class2 = phyloseq(OTU2, TAX, SAMP)
taxa_names(physeq_class2)
physeq_class2


#Question 2 ####
#Look at peacrabs and sites. 



#Question 3 ####
#Weight pre and post?