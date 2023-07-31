# Taxonomy Associated with Growth ####
# 2023-04-14
# Author: Monserrat Garcia


# Packages ####

require(dplyr)
require(tidyr)
require(data.table)
require(ggplot2)
require(tidyverse)
require(RColorBrewer)
require(phyloseq)

# Loading Data####
meta_gen18_data <- read.csv("Data/metagenetics_data18.csv")
meta17_data <- read.csv("Data/meta17_data_update.csv")

# Average Weight for each year ####

oys_gen18_data <- meta_gen18_data%>%
  filter(Species2.x == "CV")

muss_gen18_data <- meta_gen18_data%>%
  filter(Species2.x == "IR")

#Excluded one sample due to it being a very great outlier 

meta17_data <- subset(meta17_data, Weight_delta<800)

data_17 <- merge(oys_gen18_data,meta17_data,all = TRUE)



data_muss <- merge(data_17, muss_gen18_data,  all=TRUE)


# Box plots
 #Oyster

 #Volume
ggplot(data_17, aes(Year, Volume_delta))+ scale_x_continuous(breaks = c(2017,2018))+
  ylab("Volume")+
  stat_boxplot(aes(Year, Volume_delta,group=Year), geom='errorbar',
                  width=0.5)+  #whiskers
  geom_boxplot(aes(Year, Volume_delta,group= Year),outlier.shape=1) +    
  stat_summary(fun.data = mean_se)

#linetype=1,, geom = "errorbar"#stat_summary(fun = mean, geom="point", size=2) +
 

  #weight
ggplot(data_17, aes(Year, Weight_delta))+ scale_x_continuous(breaks = c(2017,2018))+
  ylab("Weight")+
  stat_boxplot( aes(Year, Weight_delta,group=Year), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Year, Weight_delta,group= Year),outlier.shape=NA) +  
  stat_summary(fun.data = mean_se)


#, geom = "errorbar", stat_summary(fun = mean, geom="point", size=2) + 

 #Mussels

data_muss$Species <- paste(data_muss$Species2.y, data_muss$Species.x, data_muss$Year, sep = "_")

data_muss$Species2<- ifelse(data_muss$Species== "CV_NA_2018", "CV2018", 
                         ifelse(data_muss$Species== "NA_CV_2017", "CV2017", "IR2018"))
                                

view(data_muss$Species)
view(data_muss$Species2)

ggplot(data_muss, aes(Species2, Volume_delta))+ 
  ylab("Volume")+ xlab(" Species and Year")+
  stat_boxplot(aes(Species2, Volume_delta,group=Species2), geom='errorbar',
               width=0.5)+  #whiskers
  geom_boxplot(aes(Species2, Volume_delta,group= Species2),outlier.shape=1) +    
  stat_summary(fun.data = mean_se)


ggplot(data_muss, aes(Species2, Weight_delta))+ 
  ylab("Weight")+ xlab(" Species and Year")+
  stat_boxplot(aes(Species2, Weight_delta,group=Species2), geom='errorbar',
               width=0.5)+  #whiskers
  geom_boxplot(aes(Species2, Weight_delta,group= Species2),outlier.shape=1) +    
  stat_summary(fun.data = mean_se)





#Bar plots
# https://stackoverflow.com/questions/29768219/grouped-barplot-in-r-with-error-bars
#https://www.bing.com/ck/a?!&&p=268e35a90f39e1f9JmltdHM9MTY4NDE5NTIwMCZpZ3VpZD0zMzQ1M2QxMS0yNTkzLTYxZWUtMjI4ZS0zMDI1MjQwYTYwNmImaW5zaWQ9NTYwMA&ptn=3&hsh=3&fclid=33453d11-2593-61ee-228e-3025240a606b&u=a1L3ZpZGVvcy9zZWFyY2g_cT1iYXJwbG90cytpbityK3N0dWRpbyt3aXRoK21lYW4rYW5kK2Vycm9yYmFyJmRvY2lkPTYwMzQ4NzAxODcyNjg2NzE4MSZtaWQ9RTMwODk5OTBEQTBFREE5NjM5NTJFMzA4OTk5MERBMEVEQTk2Mzk1MiZ2aWV3PWRldGFpbCZGT1JNPVZJUkU&ntb=1


out <- subset(meta17_data, Weight_delta>400)

dat17 <- subset(data_17, !is.na(Weight_delta))

df <-dat17 %>%
  group_by(Year)%>%
  summarise(mean=mean(Weight_delta),sd=sd(Weight_delta))


ggplot(df, aes(x=Year, y=mean)) + geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.25, size=1, 
                position=position_dodge(0.9))+
  scale_x_continuous(breaks = c(2017,2018))+ylab("Weight Mean")



data17 <- subset(data_17, !is.na(Volume_delta))

df2 <-data17 %>%
  group_by(Year)%>%
  summarise(mean=mean(Volume_delta),sd=sd(Volume_delta))

ggplot(df2, aes(x=Year, y=mean)) + geom_bar(stat = "identity", position = "dodge")+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.25, size=1, 
                position=position_dodge(0.9))+
  scale_x_continuous(breaks = c(2017,2018))+ylab("Volume Mean")




####
df1  = transform(data17, mean=rowMeans(data17[Weight_delta]), sd=apply(data17[Weight_delta],1, sd))


ggplot(df1, aes(x=as.factor(Year), y=mean, fill=Year)) +
  geom_bar(position=position_dodge(), stat="identity", colour='black') +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,position=position_dodge(.9))






##
oys_gen18_data <- oys_gen18_data %>%
  mutate(weight_avg=mean(oys_gen18_data$Weight_delta, na.rm = TRUE))

muss_gen18_data <- muss_gen18_data %>%
  mutate(weight_avg=mean(muss_gen18_data$Weight_delta, na.rm = TRUE))

meta17_data <- meta17_data %>%
  mutate(weight_avg=mean(meta17_data$Weight_delta, na.rm = TRUE))

data_oys <- data.frame(oys_gen18_data$Year,oys_gen18_data$weight_avg)

data_muss <- data.frame(muss_gen18_data$Year,muss_gen18_data$weight_avg)

data17 <- data.frame(meta17_data$Year,meta17_data$weight_avg)

dataoys <-data_oys[-c(2:59), ]
datamuss <-data_muss[-c(2:30), ]
data_17 <- data17[-c(2:112), ]

names(dataoys)[names(dataoys)=='oys_gen18_data.Year'] <- 'Year'
names(dataoys)[names(dataoys)=='oys_gen18_data.weight_avg'] <- 'weight_avg'

names(datamuss)[names(datamuss)=='muss_gen18_data.Year'] <- 'Year'
names(datamuss)[names(datamuss)=='muss_gen18_data.weight_avg'] <- 'weight_avg'


names(data_17)[names(data_17)=='meta17_data.Year'] <- 'Year'
names(data_17)[names(data_17)=='meta17_data.weight_avg'] <- 'weight_avg'


data <- merge(dataoys,data_17,all = TRUE)
data2 <- merge(data, datamuss, all = TRUE)

species <- c("Oyster", "Mussel", "Oyster")
data2$Species <- species

standard <- data2 %>%
  group_by(Species) %>%
  summarize(mean=mean(weight_avg), sd=sd(weight_avg))

sd <- sd(data2$weight_avg)

ggplot(data2, aes(x=factor(Year, levels = c("2017","2018")), y=weight_avg,fill=weight_avg)) + 
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.4)+
  geom_col(position = "dodge") +
  geom_text(aes(label=sprintf(weight_avg,fmt = "%0.2f")),colour = "white", size = 4,
    vjust = 1.0, position = position_dodge(.9))+
  xlab("Year") + ylab( "Average")

#stat_summary(geom = "bar", fun = mean, position = "dodge") +
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+

# Volume ####

meta_gen18_data <- meta_gen18_data %>%
  mutate(volume_avg=mean(meta_gen18_data$Volume_delta, na.rm = TRUE))

meta17_data <- meta17_data %>%
  mutate(volume_avg=mean(meta17_data$Volume_delta, na.rm = TRUE))

data18v <- data.frame(meta_gen18_data$Year,meta_gen18_data$volume_avg)

data17v <- data.frame(meta17_data$Year,meta17_data$volume_avg)

data_18v <-data18v[-c(2:101), ]
data_17v <- data17v[-c(2:112), ]

names(data_18v)[names(data_18v)=='meta_gen18_data.Year'] <- 'Year'
names(data_18v)[names(data_18v)=='meta_gen18_data.volume_avg'] <- 'volume_avg'


names(data_17v)[names(data_17v)=='meta17_data.Year'] <- 'Year'
names(data_17v)[names(data_17v)=='meta17_data.volume_avg'] <- 'volume_avg'


data_vol <- merge(data_18v, data_17v, all = TRUE)

sd <- sd(data_vol$volume_avg)

ggplot(data_vol, aes(x=factor(Year, levels = c("2017","2018")), y=volume_avg,fill=volume_avg)) + 
  geom_col(position = "dodge") +
  geom_text(aes(label=sprintf(volume_avg,fmt = "%0.2f")),colour = "white", size = 4,
            vjust = 1.1, position = position_dodge(.9))+
  xlab("Year") + ylab( "Average") + geom_errorbar(aes(ymin = data_vol$volume_avg-sd, ymax = data_vol$volume_avg+sd), width = 0.2)

### Scaled ####

meta_gen18_data <- meta_gen18_data %>%
  mutate(volscal_avg=mean(meta_gen18_data$Volume_scale, na.rm = TRUE))

meta17_data <- meta17_data %>%
  mutate(volscal_avg=mean(meta17_data$Volume_scale, na.rm = TRUE))

data18vscal <- data.frame(meta_gen18_data$Year,meta_gen18_data$volscal_avg)

data17vscal <- data.frame(meta17_data$Year,meta17_data$volscal_avg)

data_18vscal <-data18vscal[-c(2:101), ]
data_17vscal <- data17vscal[-c(2:112), ]

names(data_18vscal)[names(data_18vscal)=='meta_gen18_data.Year'] <- 'Year'
names(data_18vscal)[names(data_18vscal)=='meta_gen18_data.volscal_avg'] <- 'volscal_avg'


names(data_17vscal)[names(data_17vscal)=='meta17_data.Year'] <- 'Year'
names(data_17vscal)[names(data_17vscal)=='meta17_data.volscal_avg'] <- 'volscal_avg'


data_vscal <- merge(data_18vscal, data_17vscal, all = TRUE)

ggplot(data_vscal, aes(x=factor(Year, levels = c("2017","2018")), y=volscal_avg,fill=volscal_avg)) + 
  stat_summary(geom = "bar", fun = mean, position = "dodge") +
  stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")+
  geom_col(position = "dodge") +
  geom_text(aes(label=sprintf(volscal_avg,fmt = "%0.2f")),colour = "white", size = 4,
            vjust = 1.1, position = position_dodge(.9))+
  xlab("Year") + ylab( "Average") 


# https://stackoverflow.com/questions/44872951/how-do-i-add-se-error-bars-to-my-barplot-in-ggplot2


