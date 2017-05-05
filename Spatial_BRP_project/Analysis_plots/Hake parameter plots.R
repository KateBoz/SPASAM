#########################################################
# Visualization and adjustment of biological param
# Created by: Katelyn Bosley
# Date: 2/15/2017
########################################################


#set wd
setwd("C:\\Users\\katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Hake")


#read and data
data<-read.csv("Hakeparams.csv")
head(data)


#load required libraries
library(ggplot2)
library(tidyr)
library(reshape2)
library(maps)
library(mapdata)
library(ggmap)


#############
#plot matruity

#plot mat@ age
#melt data for ggplot
mat <- melt(data[,c(1,6,7)], id="Age")
head(mat)

ggplot(mat,aes(Age,value, col = variable ))+
  geom_line(lwd = 1)+
  geom_point(size  = 1.8)+
  scale_color_manual(values = c("red","blue"), labels=c("North","South"))+
  ylab("Maturity")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle("Maturity by Area")


###############
#plot for match/mismatch selectivity maturity

north<-melt(data[,c(1,6,8)], id="Age")
south<-melt(data[,c(1,7,8)], id="Age")

#plot NORTH
ggplot(north,aes(Age,value, col = variable ))+
  geom_line(lwd = 1)+
  geom_point(size  = 1.8)+
  scale_color_manual(values = c("red","black"), labels=c("Maturity","Fishery Sel"))+
  ylab("Maturity")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle("NORTH")

#plot SOUTH
ggplot(south,aes(Age,value, col = variable ))+
  geom_line(lwd = 1)+
  geom_point(size  = 1.8)+
  scale_color_manual(values = c("red","black"), labels=c("Maturity","Fishery Sel"))+
  ylab("Maturity")+
  theme_bw()+
  theme(legend.title = element_blank())+
  ggtitle("SOUTH")


#maps of west coast for plotting
map('state', region = c('california', 'oregon', 'washington'), xlim=c(-130,-115), ylim=c(30,55), fill=TRUE, col="gray95")
map("worldHires","Canada", xlim=c(-130,-90), ylim=c(30,60), col="gray95", fill=TRUE, add=TRUE)

######################################################################################



  





