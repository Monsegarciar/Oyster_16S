# 2018 Data Cleaning ####
# 2021-06-14
# Author: Monse Garcia 


#Data

meta17<-read.csv("C:/Users/monse/OneDrive/Documents/Oyster_16S/Data/metadata_de17 - Copy.csv")
data<-read.csv("Data/DE_DATA_ForGenetics - Copy.csv")


#Adding new Column for treatment
data$Treatment2<-ifelse(data$Treatment=="HH", "HIGH_POLY", ifelse(data$Treatment=="HL", "HIGH_MONO", ifelse(data$Treatment== "LL", "LOW_MONO","LOW_POLY")))

data$Colornumber<- paste0(data$Color, data$Number)

                                                                  
#Creating unique IDs
#Example: 2017_NW_HIGH_MONO_B11_CV

data$UniqueID <- paste("2017", data$Site, data$Treatment2, data$Colornumber,data$Species, sep = "_")

meta17_data<- merge(meta17, data, by= "UniqueID", all.x= TRUE) # if you want to keep all the rows in meta17 all.x=TRUE, but for data you do all.y=TRUE

#Taking out unnecessary Columns
meta17_data2 <- subset(meta17_data, select = -c(X,V1,Phase_1_DO,Phase_1_temp,Phase_2_DO,Phase_2_Temp,
                                                Overall_treatment,Notes_pre,Notes_post,Date_post,
                                                Dry_weight_shell,POST_DEAD_ALIVE,Dry_Weight_plate,Dry_weight_final,Genetics_Weight))

#Installing "data.table"
install.packages("data.table")
library("data.table")


#Uploading asvtable_17
asvtable_17<- fread("Data/asvtable_de17 - Copy.csv")


# Loading Phyloseq 
library("phyloseq")
?"phyloseq"
?otu_table
?sample_data

#Saving meta17_data2
write.csv(meta17_data2, file = "Data/meta17_data2.csv")

rownames(meta17_data2)= meta17_data2$UniqueID

#Manipulating data with Phyloseq Example  

nrow(asvtable_17)
ncol(asvtable_17)
outmat= matrix(sample(1:100, 100, replace = TRUE), nrow = 10, ncol=10)
rownames(outmat)<- paste0("OTU", 1:nrow(outmat))
colnames(outmat)<- paste0("Sample", 1:ncol(outmat))

taxmat= matrix(sample(letters, 70, replace=TRUE), nrow=nrow(outmat), 
               ncol = 7)
rownames(taxmat) <- rownames(outmat) 
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", 
                      "Family", "Genus", "Species")

class(outmat)
class(taxmat)
OTU= otu_table(outmat, taxa_are_rows = TRUE)
TAX= tax_table(taxmat)
TAX


physeq= phyloseq(OTU, TAX)
physeq

plot_bar(physeq, fill = "Family")
a_my_comparisons17_3 <- list(c("0.5", "1"), c("0.5", "2"), c("0.5", "3"), c("0.5", "4"), c("0.5", "5"), c("2", "5"))
a_my_comparisons17_2 <- list(c("0", "1"), c("0", "2"), c("0", "3"), c("0", "4"), c("0", "5"), c("1", "5"))
a_my_comparisons17 <- list(c("0", "0.5"), c("0.5", "1"), c("1", "2"), c("2", "3"), c("3", "4"), c("4", "5"))
symnum.args17 = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))

#Statistical Analysis ####
# Height 
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Height_delta, colour= RFTM_score.x)) + geom_point() +  geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", x= "RFTM Score", y= "Normalized Height", caption = "2017 Data") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
ggsave(filename = "Height Growth Statistics in Oysters 2017.jpeg", plot=last_plot(), path ="Data2017_plots/", width = 7, height = 5)

ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Height_delta, colour= RFTM_score.x)) + geom_point() +  geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", x= "RFTM Score", y= "Normalized Height", caption = "2017 Data") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17_2, label = "p.signif", symnum.args = symnum.args17)

# 0.5 Analysis
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Height_delta, colour= RFTM_score.x)) + geom_point() +  geom_boxplot(alpha=0.3) + labs(title = "Height Growth in Oysters", x= "RFTM Score", y= "Normalized Height", caption = "2017 Data") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17_3, label = "p.signif", symnum.args = symnum.args17)

# Weight
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Weight_delta, colour= RFTM_score.x)) + geom_point()+ geom_boxplot(alpha=0.3) + labs(title = "Weight Growth in Oysters", x= "RFTM Score", y= "Normalized Weight") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)
 #+ scale_y_continuous(limits=c(0, 50)) 

# Length
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Length_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Length Growth in Oysters", x= "RFTM Score", y= "Normalized Length") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)

# Width 
ggplot(data = meta17_data, aes(x = as.factor(RFTM_score.x), y = Width_delta, colour= RFTM_score.x)) + geom_point() + geom_boxplot(alpha=0.3) + labs(title = "Width Growth in Oysters", x= "RFTM Score", y= "Normalized Width") +  stat_compare_means(method = "wilcox.test", comparisons = a_my_comparisons17, label = "p.signif", symnum.args = symnum.args17)





