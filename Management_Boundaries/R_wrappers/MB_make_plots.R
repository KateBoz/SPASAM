##################################################
# Quick and Dirty Make plots for CAPAM
#
# Created by:Katelyn Bosley 
# 10/3/2018
##################################################

#Brute force making plot for AB presentation

#Need to combine the .csv files from the 2 runs into a single data frame for plotting.

direct_master<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\MB_test\\Managment_Boundaries"
setwd(direct_master)

#list files in the directory
files<-list.files(direct_master) 

#pulling only the folder that you want
folder.num = 6
i=folder.num

#read in the file

data<-read.csv(paste0(direct_master,"\\",files[i],"\\ssb_table_total.csv"))

#sum from the two runs
data.est1<-subset(data,Dat=="EST"& Run==unique(data$Run)[1])
data.est2<-subset(data,Dat=="EST"& Run==unique(data$Run)[2])
data.est.tot<-data.est1 #setting dims
data.est.tot[,5:length(data.est.tot)]<-data.est1[,5:length(data.est.tot)]+data.est2[,5:length(data.est.tot)]


#calculate the bias
data.sim<-subset(data,Dat=="SIM")
data.est.bias<-data.sim
data.est.bias[,5:length(data.est.bias)]<-((data.sim[,5:length(data.est.tot)]-data.est.tot[,5:length(data.est.tot)])/data.sim[,5:length(data.est.tot)])*100
ssb.long.bias<-melt(data.est.bias[,3:dim(data.est.bias)[2]], id=c("Years","Reg"))

#bias medians
#calculate the medians
ssb.bias.meds <- ssb.long.bias%>% group_by(Reg,Years) %>%
  summarise(meds=median(value))



#combine the data frames
data.plot<-rbind(data.est.tot,subset(data,Dat=="SIM"))
ssb.long<-melt(data.plot[,2:dim(data.est.tot)[2]], id=c("Dat","Years","Reg"))
ssb.long.est<-subset(ssb.long,Dat=="EST")

#calculate the medians
ssb.meds <- ssb.long %>% group_by(Dat,Reg,Years) %>%
  summarise(meds=median(value))

# time series plot
ssb.plot.gg<-ggplot(ssb.long.est, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = subset(ssb.meds,Dat=="EST"), aes(x=Years,y=meds), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = subset(ssb.meds,Dat=="SIM"), aes(x=Years,y=meds),lty=1, lwd=0.5) + 
  geom_point(data = subset(ssb.meds,Dat=="SIM"), aes(x=Years,y=meds), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("SSB")+
  ylab("SSB")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
ssb.bias.gg<-ggplot(ssb.long.bias, aes(x=as.factor(Years), y=value)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_line(data = ssb.bias.meds, aes(x=Years,y=meds),lty=1, lwd=0.5) + 
  geom_point(data = ssb.bias.meds, aes(x=Years,y=meds), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("SSB Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(quantile(ssb.est$bias, c(.05,.95)))+
  my_theme


######################
###recruitment
#######################

rec.data<-read.csv(paste0(direct_master,"\\",files[i],"\\rec_table_total.csv"))

#sum from the two runs
rec.data.est1<-subset(rec.data,Dat=="EST"& Run==unique(rec.data$Run)[1])
rec.data.est2<-subset(rec.data,Dat=="EST"& Run==unique(rec.data$Run)[2])
rec.data.est.tot<-rec.data.est1 #setting dims
rec.data.est.tot[,5:length(rec.data.est.tot)]<-rec.data.est1[,5:length(rec.data.est.tot)]+rec.data.est2[,5:length(rec.data.est.tot)]


#calculate the bias
rec.data.sim<-subset(rec.data,Dat=="SIM")
rec.data.est.bias<-rec.data.sim
rec.data.est.bias[,5:length(rec.data.est.bias)]<-((rec.data.sim[,5:length(rec.data.est.tot)]-rec.data.est.tot[,5:length(rec.data.est.tot)])/rec.data.sim[,5:length(rec.data.est.tot)])*100
rec.long.bias<-melt(rec.data.est.bias[,3:dim(rec.data.est.bias)[2]], id=c("Years","Reg"))

#bias medians
#calculate the medians
rec.bias.meds <- rec.long.bias%>% group_by(Reg,Years) %>%
  summarise(meds=median(value))



#combine the data frames
rec.data.plot<-rbind(rec.data.est.tot,subset(rec.data,Dat=="SIM"))
rec.long<-melt(rec.data.plot[,2:dim(rec.data.est.tot)[2]], id=c("Dat","Years","Reg"))
rec.long.est<-subset(rec.long,Dat=="EST")

#calculate the medians
rec.meds <- rec.long %>% group_by(Dat,Reg,Years) %>%
  summarise(meds=median(value))

# time series plot
rec.plot.gg<-ggplot(rec.long.est, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = subset(rec.meds,Dat=="EST"), aes(x=Years,y=meds), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = subset(rec.meds,Dat=="SIM"), aes(x=Years,y=meds),lty=1, lwd=0.5) + 
  geom_point(data = subset(rec.meds,Dat=="SIM"), aes(x=Years,y=meds), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment")+
  ylab("Recruitment")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
rec.bias.gg<-ggplot(rec.long.bias, aes(x=as.factor(Years), y=value)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_line(data = rec.bias.meds, aes(x=Years,y=meds),lty=1, lwd=0.5) + 
  geom_point(data = rec.bias.meds, aes(x=Years,y=meds), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(quantile(ssb.est$bias, c(.05,.95)))+
  my_theme

grid.arrange(ncol = 1,
             top="Recruitment Estimation",
             rec.plot.gg, rec.bias.gg)

#Biomass plots
grid.arrange(ncol = 1,
             top="SSB",
             ssb.plot.gg, ssb.bias.gg)


# that is the end for now...
#################################################################





