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


data_17 <- merge(oys_gen18_data,meta17_data,all = TRUE)
data_o <- merge(oys_gen18_data,meta17_data, by="Species2.x", all.y=TRUE)


# Bar plots
 #Oyster

 #Volume
ggplot(data_17, aes(Year, Volume_delta))+ scale_x_continuous(breaks = c(2017,2018))+
  ylab("Volume")+
  stat_boxplot( aes(Year, Volume_delta,group=Year), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Year, Volume_delta,group= Year),outlier.shape=1) +    
  stat_summary(fun = mean, geom="point", size=2) + 
  stat_summary(fun.data = mean_se, geom = "errorbar")


 #weight
ggplot(data_17, aes(Year, Weight_delta))+ scale_x_continuous(breaks = c(2017,2018))+
  stat_boxplot( aes(Year, Weight_delta,group=Year), 
                geom='errorbar', linetype=1, width=0.5)+  #whiskers
  geom_boxplot( aes(Year, Weight_delta,group= Year),outlier.shape=NA) + coord_cartesian(ylim = quantile(data_17$y, c(0.1, 0.9)))+  
  stat_summary(fun = mean, geom="point", size=2) + 
  stat_summary(fun.data = mean_se, geom = "errorbar")



 #Mussels















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


