# Questions for Analysis ####
# 2021-06-23
# Author: Monse Garcia

#Packages Required
library(phyloseq)
library(ggplot2)

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

#b2=plot_bar(pp.ch, "Site.x", fill="Genus", facet_grid=~Family)
#print(b2)

## Alpha Diversity Graphics

PC <- prune_species(speciesSums(physeq_class)> 0, physeq_class)
PC
plot_richness(PC) #https://github.com/joey711/phyloseq/issues/552 for troubleshooting



#Question 2 ####
#Look at peacrabs and sites. 



#Question 3 ####
#Weight pre and post?