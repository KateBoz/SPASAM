rm(list=(ls()))


##### Inputs ###############################################

runs.to.include<-c(1)  # the numbers corresponding to the runs to include in the MRE summary plots

nsim <-2
num_dir<-1

####### Where to store figures ############
Figure_dir<-"G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\movement manuscript\\runs\\Figures\\Figures"
Figure_table_dir<-"G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\movement manuscript\\runs\\Figures\\tables"

##########################################

###################### Directories that want to load (must match num_dir above) ###########################

####### BASE ############################
wd_1<<-'G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\movement manuscript\\runs\\test'
nm_1<-"Base"
################################################

####### No Movement ############################
wd_2<<-'F:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\Tagging manuscript\\Runs\\No Movement'
nm_2<-"No_Move"
################################################

##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################
##########################################################################################################################################################

       
#load required libraries
load_libraries<-function() {
  library(PBSmodelling)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(reshape)
  library(gridExtra)
  library(gplots)
  library(colorspace)
  library(RColorBrewer)
  #library(plyr)
  library(dplyr)
  library(tidyr)
  library(matrixStats) 
  library(gridExtra)
  library(grid)
  library(gtools)
  library(TeachingDemos)
  library(snowfall)
  library(parallel)
  library(snow)
  library(foreach)
  library(doSNOW)
  library(spatstat)
  library(alphahull)
  library(beanplot)
  library(png)
  library(sjPlot)
  library(xtable)
  library(forcats)
  library(ggforce)
  #library(cowplot)
}
load_libraries()
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")


#pull dimensions for building data frames for plotting
out<-readList("G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\movement manuscript\\runs\\test\\Diagnostics\\Results_converged\\run1.rep") #read in .rep file

#pull info about the model
na<-out$nages
nyrs<-out$nyrs
npops<-out$npops
nreg<-out$nregions
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)
if(npops>1){   #for metapop
  nreg=sum(nreg)}


vio.col<-"lightskyblue3"
median.col<-"gray25"

################


for (z in 1:num_dir) {
  direct_master<-eval(parse(text=paste0("wd_",z)))
######################################


diag_dir<-paste0(direct_master,"\\Diagnostics",sep="")
runs_dir<-paste0(diag_dir,"\\Results_converged",sep="")

# create a data frame for the run to hold true and estimated values from OM/EM

setwd(diag_dir)

###############################################################
# Basic Plots
###############################################################
  
 #load results
  
#load up the results
Sim_Results<-readRDS('Sim_data.RData')
Sim_Stats<-read.csv('Sim_Stats.csv')
conv.runs<-which(Sim_Stats$Converged==1)

nconv<-length(conv.runs)

#saving values
R_ave_df_sim<-matrix(NA,npops,nconv)
R_ave_df_est<-matrix(NA,npops,nconv)

R_apport_df_sim<-matrix(NA,(nyrs-1)*nreg,nconv)
R_apport_df_est<-matrix(NA,(nyrs-1)*nreg,nconv)

q_df_sim<-matrix(NA,nreg,nconv)
q_df_est<-matrix(NA,nreg,nconv)


#selectivity params
sel_beta1_df_sim<-matrix(NA,nreg,nconv)
sel_beta2_df_sim<-matrix(NA,nreg,nconv)
sel_beta3_df_sim<-matrix(NA,nreg,nconv)
sel_beta4_df_sim<-matrix(NA,nreg,nconv)
sel_beta1_df_est<-matrix(NA,nreg,nconv)
sel_beta2_df_est<-matrix(NA,nreg,nconv)
sel_beta3_df_est<-matrix(NA,nreg,nconv)
sel_beta4_df_est<-matrix(NA,nreg,nconv)


sel_beta1_surv_df_sim<-matrix(NA,nreg,nconv)
sel_beta2_surv_df_sim<-matrix(NA,nreg,nconv)
sel_beta3_surv_df_sim<-matrix(NA,nreg,nconv)
sel_beta4_surv_df_sim<-matrix(NA,nreg,nconv)
sel_beta1_surv_df_est<-matrix(NA,nreg,nconv)
sel_beta2_surv_df_est<-matrix(NA,nreg,nconv)
sel_beta3_surv_df_est<-matrix(NA,nreg,nconv)
sel_beta4_surv_df_est<-matrix(NA,nreg,nconv)


#Recruitment
rec_df_sim<-matrix(NA,nyrs*nreg,nconv)
rec_df_est<-matrix(NA,nyrs*nreg,nconv)

#Recruitment deviations
rec_devs_df_sim<-matrix(NA,nyrs*nreg,nconv)
rec_devs_df_est<-matrix(NA,(nyrs-1)*nreg,nconv)

#Init Abundance
init_abund_df_sim<-matrix(NA,nreg*nreg*na,nconv)
init_abund_df_est<-matrix(NA,nreg*nreg*na,nconv)

#SSB
ssb_df_sim<-matrix(NA,nyrs*nreg,nconv)
ssb_df_est<-matrix(NA,nyrs*nreg,nconv)

#Biomass
bio_df_sim<-matrix(NA,nyrs*nreg,nconv)
bio_df_est<-matrix(NA,nyrs*nreg,nconv)

#yield
catch_df_sim<-matrix(NA,nyrs*nreg,nconv)
catch_df_est<-matrix(NA,nyrs*nreg,nconv)

#survey bio
survey_df_sim<-matrix(NA,nyrs*nreg,nconv)
survey_df_est<-matrix(NA,nyrs*nreg,nconv)

#fmax
fmax_df_sim<-matrix(NA,nyrs*nreg,nconv)
fmax_df_est<-matrix(NA,nyrs*nreg,nconv)

#select at age
select_age_df_sim<-matrix(NA,na*nreg,nconv)
select_age_df_est<-matrix(NA,na*nreg,nconv)
select_age_survey_df_sim<-matrix(NA,na*nreg,nconv)
select_age_survey_df_est<-matrix(NA,na*nreg,nconv)

#movement
move_df_sim<-matrix(NA,nyrs*nreg*na*nreg,nconv)
move_df_est<-matrix(NA,nyrs*nreg*na*nreg,nconv)


###########################################################
# populate the matrices for plotting
for(i in 1:nconv){
R_ave_df_sim[,i]<-unlist(Sim_Results["meanR_sim",i])
R_ave_df_est[,i]<-unlist(Sim_Results["meanR_est",i])
R_apport_df_sim[,i]<-unlist(Sim_Results["apport_sim",i])
R_apport_df_est[,i]<-unlist(Sim_Results["apport_est",i])
q_df_sim[,i]<-unlist(Sim_Results["q_sim",i])
q_df_est[,i]<-unlist(Sim_Results["q_est",i])
sel_beta1_df_sim[,i]<-unlist(Sim_Results["sel_beta1_sim",i])
sel_beta2_df_sim[,i]<-unlist(Sim_Results["sel_beta2_sim",i])
sel_beta3_df_sim[,i]<-unlist(Sim_Results["sel_beta3_sim",i])
sel_beta4_df_sim[,i]<-unlist(Sim_Results["sel_beta4_sim",i])
sel_beta1_df_est[,i]<-unlist(Sim_Results["sel_beta1_est",i])
sel_beta2_df_est[,i]<-unlist(Sim_Results["sel_beta2_est",i])
sel_beta3_df_est[,i]<-unlist(Sim_Results["sel_beta3_est",i])
sel_beta4_df_est[,i]<-unlist(Sim_Results["sel_beta4_est",i])
sel_beta1_surv_df_sim[,i]<-unlist(Sim_Results["sel_beta1_surv_sim",i])
sel_beta2_surv_df_sim[,i]<-unlist(Sim_Results["sel_beta2_surv_sim",i])
sel_beta3_surv_df_sim[,i]<-unlist(Sim_Results["sel_beta3_surv_sim",i])
sel_beta4_surv_df_sim[,i]<-unlist(Sim_Results["sel_beta4_surv_sim",i])
sel_beta1_surv_df_est[,i]<-unlist(Sim_Results["sel_beta1_surv_est",i])
sel_beta2_surv_df_est[,i]<-unlist(Sim_Results["sel_beta2_surv_est",i])
sel_beta3_surv_df_est[,i]<-unlist(Sim_Results["sel_beta3_surv_est",i])
sel_beta4_surv_df_est[,i]<-unlist(Sim_Results["sel_beta4_surv_est",i])
rec_df_sim[,i]<-unlist(Sim_Results["rec_sim",i])
rec_df_est[,i]<-unlist(Sim_Results["rec_est",i])
rec_devs_df_sim[,i]=unlist(Sim_Results["rec_devs_sim",i])
rec_devs_df_est[,i]=unlist(Sim_Results["rec_devs_est",i])
init_abund_df_sim[,i]=unlist(Sim_Results["init_abund_sim",i])
init_abund_df_est[,i]=unlist(Sim_Results["init_abund_est",i])
ssb_df_sim[,i]<-unlist(Sim_Results["ssb_sim",i])
ssb_df_est[,i]<-unlist(Sim_Results["ssb_est",i])
bio_df_sim[,i]<-unlist(Sim_Results["bio_sim",i])
bio_df_est[,i]<-unlist(Sim_Results["bio_est",i])
catch_df_sim[,i]<-unlist(Sim_Results["yield_sim",i])
catch_df_est[,i]<-unlist(Sim_Results["yield_est",i])
survey_df_sim[,i]<-unlist(Sim_Results["survey_sim",i])
survey_df_est[,i]<-unlist(Sim_Results["survey_est",i])
fmax_df_sim[,i]<-unlist(Sim_Results["fmax_sim",i])
fmax_df_est[,i]<-unlist(Sim_Results["fmax_est",i])
move_df_sim[,i]<-unlist(Sim_Results["movement_sim",i])
move_df_est[,i]<-unlist(Sim_Results["movement_est",i])
select_age_df_sim[,i]<-unlist(Sim_Results["select_age_sim",i])
select_age_df_est[,i]<-unlist(Sim_Results["select_age_est",i])
select_age_survey_df_sim[,i]<-unlist(Sim_Results["select_age_survey_sim",i])
select_age_survey_df_est[,i]<-unlist(Sim_Results["select_age_survey_est",i])
}


assign(paste0("conv.rate_",z),nconv/nsim)



#########################################
############################################
# MAKE PLOTS
###########################################


  
#######################################
# Building a ggplot theme
######################################
my_theme<-
  (theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="grey20",size=8,angle=0,hjust=.5,vjust=0,face="plain"),
    plot.title=element_text(hjust=0.5,size=12),
    plot.subtitle=element_text(hjust=0.5,size=8)
  ))

######################################
# Building the plots!
######################################

#R_ave
Rave_est<-data.frame(melt(t(R_ave_df_est)))
Rave_est<-cbind(Rave_est,data.frame(melt(t(R_ave_df_sim))[3]))
names(Rave_est)<-c("Nsim","Reg","R_ave_est","R_ave_sim")
Rave_est$R_ave_bias<-((Rave_est$R_ave_est-Rave_est$R_ave_sim)/Rave_est$R_ave_sim)*100

#calc medians

#calculate the sum across areas 
rave.long<-melt(Rave_est, id=c("Reg","Nsim"))
#rave.long$Reg<-as.character(rave.long$Reg)

median_R_ave<-data.frame(rave.long %>% group_by(Reg,variable) %>% 
                           summarise(med=round(median(value),digits=2),MARE=round(median(abs(value)),digits=2), 
                                                                           Mean=round(mean(value),digits=2), 
                                                                           min=round(min(value),digits=2),
                                     max=round(max(value),digits=2)))


rave_sum_table<-c(get(paste0("nm_",z)),round(get(paste0("conv.rate_",z)),digits=2),
                  median_R_ave$MARE[median_R_ave$variable=='R_ave_bias' & median_R_ave$Reg=="1"],
                  median_R_ave$MARE[median_R_ave$variable=='R_ave_bias' & median_R_ave$Reg=="2"],
                  median_R_ave$Mean[median_R_ave$variable=='R_ave_bias' & median_R_ave$Reg=="1"],
                  median_R_ave$Mean[median_R_ave$variable=='R_ave_bias' & median_R_ave$Reg=="2"])
names(rave_sum_table)<-c("Run","Conv_Rate","MARE_Reg_1","MARE_Reg_2",
                        "MRE_Reg_1","MRE_Reg_2")

#plot
rave.plot.gg<-ggplot(Rave_est, aes(x=as.factor(Reg), y=R_ave_est, group=Reg)) +
  #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  geom_point(data=subset(median_R_ave,variable=="R_ave_est"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
  geom_point(data=subset(median_R_ave,variable=="R_ave_sim"), aes(x=Reg,y=med), col="black", shape=16,cex=1.0) + 
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+
  ylab("Mean Recruitment")+
  xlab("Area")+
  my_theme



rave.bias.gg<-ggplot(Rave_est, aes(x=as.factor(Reg), y=R_ave_bias, group=Reg)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'black',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0",alpha=0.6)+
  geom_point(data=subset(median_R_ave,variable=="R_ave_bias"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
  #geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+
  ylab("Relative % Difference")+
  xlab("Area")+
  #ylim(quantile(Rave_est$R_ave_bias, c(.05,.95)))+
  my_theme

#####################################
# Recruitment plot
#####################################

#build data.frame
rec.data<-data.frame(Dat=c(rep("SIM",nrow(rec_df_sim)),rep("EST",nrow(rec_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

rec.data<-cbind(rec.data,rbind(rec_df_sim,rec_df_est))
rec.long<-melt(rec.data, id=c("Dat","Years","Reg"))
rec.long$Reg<-as.factor(as.character(rec.data$Reg))

#calculate the sum across areas 
total.rec<-data.frame(rec.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))


total.rec$Reg<-rep("System",nrow(total.rec))
rec.long<-rbind(rec.long,total.rec)


#separate again for plotting
rec.est<-rec.long[rec.long$Dat=="EST",]
rec.sim<-rec.long[rec.long$Dat=="SIM",]


#calculate the percent bias
rec.est$val.true<-rec.sim$value
rec.est$bias=((rec.est$value-rec.est$val.true)/rec.est$val.true)*100
#rec.est$bias=(rec.est$val.true-rec.est$value)

#calc medians table
rec.est.med <- rec.est %>% group_by(Reg,Years) %>%
  summarise(med.est=round(median(value),digits=2),med.sim=round(median(val.true),digits=2), med.bias = round(median(bias),digits=2),
            MARE=round(median(abs(bias)),digits=2),
            mean.est=round(mean(value),digits=2),mean.sim=round(mean(val.true),digits=2), MRE = round(mean(bias),digits=2),
            max=round(max(bias),digits=2),min=round(min(bias),digits=2))


rec_sum_table<-c(get(paste0("nm_",z)),round(get(paste0("conv.rate_",z)),digits=2),rec.est.med$MARE[2],rec.est.med$MARE[32],rec.est.med$MARE[62],rec.est.med$MARE[30],rec.est.med$MARE[60],rec.est.med$MARE[90],
                 round(median(abs(rec.est$bias[rec.est$Reg=='1'])),digits=2),round(median(abs(rec.est$bias[rec.est$Reg=='2'])),digits=2),
                 round(median(abs(rec.est$bias[rec.est$Reg=='System'])),digits=2),
                 round(median(abs(rec.est$bias)),digits=2),
                 rec.est.med$MRE[2],rec.est.med$MRE[32],rec.est.med$MRE[62],rec.est.med$MRE[30],rec.est.med$MRE[60],rec.est.med$MRE[90],
                 round(mean(rec.est$bias[rec.est$Reg=='1']),digits=2),round(mean(rec.est$bias[rec.est$Reg=='2']),digits=2),
                 round(mean(rec.est$bias[rec.est$Reg=='System']),digits=2),
                 round(mean(rec.est$bias),digits=2))
names(rec_sum_table)<-c("Run","Conv_Rate","MARE_Reg_1_YR_2","MARE_Reg_2_YR_2","MARE_Tot_YR_2","MARE_Reg_1_YR_30","MARE_Reg_2_YR_30","MARE_Tot_YR_30",
                        "MARE_YR_Reg_1","MARE_YR_Reg_2","MARE_YR_Total","MARE_YR+Reg",
                        "MRE_Reg_1_YR_2","MRE_Reg_2_YR_2","MRE_Tot_YR_2","MRE_Reg_1_YR_30","MRE_Reg_2_YR_30","MRE_Tot_YR_30","MRE_YR_Reg_1","MRE_YR_Reg_2","MRE_YR_Total",
                        "MRE_YR+Reg")
#generate Rec Plot
rec.plot.gg<-ggplot(rec.est, aes(x=as.factor(Years), y=value)) +
  #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  geom_point(data = rec.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = rec.est.med, aes(x=Years,y=med.sim),lty=1,lwd=0.5) + 
  geom_point(data = rec.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Recruitment")+
  xlab("Year")+
  ylim(0,25000)+
  facet_grid(Reg~.)+
  my_theme


my_theme.fig4<-
  (#theme_bw()+
    theme(
      plot.margin = unit(c(0, .15, 0, .15), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="black",size=6,angle=0,hjust=0.45,vjust=1,face="plain"),
      axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=6),
      axis.title.y = element_text(size=6),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,vjust=0,size=8,margin=margin(2,2,3,2)),
      plot.subtitle=element_text(hjust=0.5,size=8),
      strip.text.x = element_text(margin = margin(0,0,.15,0, "cm"),
        size = 5, color = "black"),
      strip.text.y = element_text(
        size = 5, color = "black", face = "italic"),
      strip.background = element_rect(
        fill=NA),
      legend.position = c(0.5, 0.5),
      legend.text=element_text(size=5),
      legend.title=element_text(size=5),
      legend.key.size = unit(.15,"cm"),
      panel.background = element_blank(), #element_rect(colour = "black", size=.75, fill=NA),
      #axis.line.x = element_line(),
      #axis.line.y = element_line(),
      legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid',size=.25),
      #legend.title=element_text(),
      legend.margin=margin(c(1,2,1,2)),
      panel.border = element_rect(colour = "black", fill=NA, size=.75)
    ))



rec.name<-c('1'="Population 1 Recruitment")
rec.plot.gg.fig4<-ggplot(subset(rec.est, Reg %in% '1'), aes(x=as.factor(Years), y=value)) +
  geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
  geom_point(data = subset(rec.est.med, Reg %in% '1'), aes(x=Years,y=med.est), fill='black', colour='black',shape=20,size=.75) + 
  geom_line(data =subset(rec.est.med, Reg %in% '1'), aes(x=Years,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  ylab('Recruits (1000s Fish)')+
  xlab(NULL)+
  #facet_grid(.~Reg,labeller=as_labeller(rec.name))+
  ylim(0,25000)+
  my_theme.fig4


#generate Rec bias plot
rec.bias.gg<-ggplot(rec.est, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'black',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  geom_line(data = rec.est.med, aes(x=Years,y=med.bias),lty=1,lwd=0.5) + 
  geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) +
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.,scales="free")+
  ylim(-100,100)+
  my_theme

######################################
# Biomass

#build data.frame
bio.data<-data.frame(Dat=c(rep("SIM",nrow(bio_df_sim)),rep("EST",nrow(bio_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

bio.data<-cbind(bio.data,rbind(bio_df_sim,bio_df_est))
bio.long<-melt(bio.data, id=c("Dat","Years","Reg"))
bio.long$Reg<-as.character(bio.data$Reg)

#calculate the sum across areas 
total.bio<-data.frame(bio.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.bio$Reg<-rep("System",nrow(total.bio))
bio.long<-rbind(total.bio,bio.long)


#separate again for plotting
bio.est<-bio.long[bio.long$Dat=="EST",]
bio.sim<-bio.long[bio.long$Dat=="SIM",]


#calculate the percent bias
bio.est$val.true<-bio.sim$value
bio.est$bias=((bio.est$value-bio.est$val.true)/bio.est$val.true)*100

#calc medians table
bio.est.med <- bio.est %>% group_by(Reg,Years) %>%
  summarise(med.est=round(median(value),digits=2),med.sim=round(median(val.true),digits=2), med.bias = round(median(bias),digits=2),
            MARE=round(median(abs(bias)),digits=2),
            mean.est=round(mean(value),digits=2),mean.sim=round(mean(val.true),digits=2), MRE = round(mean(bias),digits=2),
            max=round(max(bias),digits=2),min=round(min(bias),digits=2))
bio_sum_table<-c(get(paste0("nm_",z)),get(paste0("conv.rate_",z)),bio.est.med$MARE[2],bio.est.med$MARE[32],bio.est.med$MARE[62],bio.est.med$MARE[30],bio.est.med$MARE[60],bio.est.med$MARE[90],
                 round(median(abs(bio.est$bias[bio.est$Reg=='1'])),digits=2),round(median(abs(bio.est$bias[bio.est$Reg=='2'])),digits=2),
                 round(median(abs(bio.est$bias[bio.est$Reg=='System'])),digits=2),
                 round(median(abs(bio.est$bias)),digits=2),
                 bio.est.med$MRE[2],bio.est.med$MRE[32],bio.est.med$MRE[62],bio.est.med$MRE[30],bio.est.med$MRE[60],bio.est.med$MRE[90],
                 round(mean(bio.est$bias[bio.est$Reg=='1']),digits=2),round(mean(bio.est$bias[bio.est$Reg=='2']),digits=2),
                 round(mean(bio.est$bias[bio.est$Reg=='System']),digits=2),
                 round(mean(bio.est$bias),digits=2))
names(bio_sum_table)<-c("Run","Conv_Rate","MARE_Reg_1_YR_2","MARE_Reg_2_YR_2","MARE_Tot_YR_2","MARE_Reg_1_YR_30","MARE_Reg_2_YR_30","MARE_Tot_YR_30",
                        "MARE_YR_Reg_1","MARE_YR_Reg_2","MARE_YR_Total","MARE_YR+Reg",
                        "MRE_Reg_1_YR_2","MRE_Reg_2_YR_2","MRE_Tot_YR_2","MRE_Reg_1_YR_30","MRE_Reg_2_YR_30","MRE_Tot_YR_30","MRE_YR_Reg_1","MRE_YR_Reg_2","MRE_YR_Total",
                        "MRE_YR+Reg")



#generate bio Plot
bio.plot.gg<-ggplot(bio.est, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  geom_point(data = bio.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = bio.est.med, aes(x=Years,y=med.sim),lty=1, lwd=0.5) + 
  geom_point(data = bio.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Biomass")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(0,150000)+
  my_theme


bio.name<-c('Syste,'="Metapopulation Biomass")
bio.plot.gg.fig5a<-ggplot(subset(bio.est, Reg %in% 'System'), aes(x=as.factor(Years), y=value)) +
  geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
  geom_point(data = subset(bio.est.med, Reg %in% 'System'), aes(x=Years,y=med.est), fill='black', colour='black',shape=20,size=.75) + 
  geom_line(data =subset(bio.est.med, Reg %in% 'System'), aes(x=Years,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
  ylab('Biomass (mt)')+
  xlab(NULL)+
  #facet_grid(.~Reg,labeller=as_labeller(f.name))+
  ylim(0,160000)+
  my_theme.fig4


#generate Rec bias plot
bio.bias.gg<-ggplot(bio.est, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'black',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0",alpha=0.6)+
  geom_line(data = bio.est.med, aes(x=Years,y=med.bias),lty=1, lwd=0.5) + 
  geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(-100,100)+
  my_theme

##############################################################
###F Max

#build data.frame
fmax.data.sim<-data.frame(Dat=rep("SIM",nrow(fmax_df_sim)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))
fmax.data.est<-data.frame(Dat=rep("EST",nrow(fmax_df_est)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))

fmax.data<-rbind(fmax.data.sim,fmax.data.est)

fmax.data<-cbind(fmax.data,rbind(fmax_df_sim,fmax_df_est))
fmax.long<-melt(fmax.data, id=c("Dat","Years","Reg"))
fmax.long$Reg<-as.character(fmax.data$Reg)

#calculate the sum across areas 
#total.fmax<-data.frame(fmax.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

#total.fmax$Reg<-rep("System",nrow(total.fmax))
#fmax.long<-rbind(total.fmax,fmax.long)


#separate again for plotting
fmax.est<-fmax.long[fmax.long$Dat=="EST",]
fmax.sim<-fmax.long[fmax.long$Dat=="SIM",]


#calculate the percent bias
fmax.est$val.true<-fmax.sim$value
fmax.est$bias=((fmax.est$value-fmax.est$val.true)/fmax.est$val.true)*100
#fmax.est$bias=((fmax.est$val.true-fmax.est$value))

#calc medians table
fmax.est.med <- fmax.est %>% group_by(Reg,Years) %>%
  summarise(med.est=round(median(value),digits=2),med.sim=round(median(val.true),digits=2), med.bias = round(median(bias),digits=2),
            MARE=round(median(abs(bias)),digits=2),
            mean.est=round(mean(value),digits=2),mean.sim=round(mean(val.true),digits=2), MRE = round(mean(bias),digits=2),
            max=round(max(bias),digits=2),min=round(min(bias),digits=2))

fmax_sum_table<-c(get(paste0("nm_",z)),round(get(paste0("conv.rate_",z)),digits=2),fmax.est.med$MARE[2],fmax.est.med$MARE[32],fmax.est.med$MARE[30],fmax.est.med$MARE[60],
                 round(median(abs(fmax.est$bias[fmax.est$Reg=='1'])),digits=2),round(median(abs(fmax.est$bias[fmax.est$Reg=='2'])),digits=2),
                 round(median(abs(fmax.est$bias)),digits=2),
                 fmax.est.med$MRE[2],fmax.est.med$MRE[32],fmax.est.med$MRE[30],fmax.est.med$MRE[60],
                 round(mean(fmax.est$bias[fmax.est$Reg=='1']),digits=2),round(mean(fmax.est$bias[fmax.est$Reg=='2']),digits=2),
                 round(mean(fmax.est$bias),digits=2))
names(fmax_sum_table)<-c("Run","Conv_Rate","MARE_Reg_1_YR_2","MARE_Reg_2_YR_2","MARE_Reg_1_YR_30","MARE_Reg_2_YR_30",
                        "MARE_YR_Reg_1","MARE_YR_Reg_2","MARE_YR+Reg",
                        "MRE_Reg_1_YR_2","MRE_Reg_2_YR_2","MRE_Reg_1_YR_30","MRE_Reg_2_YR_30","MRE_YR_Reg_1","MRE_YR_Reg_2",
                        "MRE_YR+Reg")
#generate fmax Plot
fmax.plot.gg<-ggplot(fmax.est, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  geom_point(data = fmax.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_point(data = fmax.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
  geom_line(data = fmax.est.med, aes(x=Years,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("F")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(0,2.)+
  my_theme

f.name<-c('1'="Population 1 Fishing Mortality")
fmax.plot.gg.fig4<-ggplot(subset(fmax.est, Reg %in% '1'), aes(x=as.factor(Years), y=value)) +
  geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
  geom_point(data = subset(fmax.est.med, Reg %in% '1'), aes(x=Years,y=med.est), fill='black', colour='black',shape=20,size=.75) + 
  #geom_point(data = subset(rec.est.med, Reg %in% '1'), aes(x=Years,y=med.sim), fill=median.col, shape=20,size=.5) + 
  geom_line(data =subset(fmax.est.med, Reg %in% '1'), aes(x=Years,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
  ylab('F')+
  xlab(NULL)+
  #facet_grid(.~Reg,labeller=as_labeller(f.name))+
  ylim(0,1.75)+
  my_theme.fig4


fmax.plot.gg.fig5a<-ggplot(subset(fmax.est, Reg %in% '2'), aes(x=as.factor(Years), y=value)) +
  geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
  geom_point(data = subset(fmax.est.med, Reg %in% '2'), aes(x=Years,y=med.est), fill='black', colour='black',shape=20,size=.75) + 
  geom_line(data =subset(fmax.est.med, Reg %in% '2'), aes(x=Years,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
  ylab('F')+
  xlab(NULL)+
  #facet_grid(.~Reg,labeller=as_labeller(f.name))+
  ylim(0,1.75)+
  my_theme.fig4

#bias plot
fmax.bias.gg<-ggplot(fmax.est, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'black',size=0.5,lty=2,lwd=2)+
  geom_violin(fill=vio.col,trim=T,bw="nrd0")+
  geom_line(data = fmax.est.med, aes(x=Years,y=med.bias),lty=1, lwd=0.5) + 
  geom_point(data = fmax.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(-100,100)+
  my_theme

#Movement

movement.data.sim<-data.frame(Dat=rep("SIM",nrow(move_df_sim)),Year=rep(1:nyrs,each=na*nreg, times=nreg),Reg_from=rep(1:nreg,each=nyrs*nreg*na),Reg_to=rep(1:nreg,times=nyrs*nreg*na),
                              Age=rep(1:na,each=nreg,times=nyrs*nreg))

movement.data.est<-data.frame(Dat=rep("EST",nrow(move_df_est)),Year=rep(1:nyrs,each=na*nreg, times=nreg),Reg_from=rep(1:nreg,each=nyrs*nreg*na),Reg_to=rep(1:nreg,times=nyrs*nreg*na),
                              Age=rep(1:na,each=nreg,times=nyrs*nreg))

movement.data<-rbind(movement.data.sim,movement.data.est)
movement.data<-cbind(movement.data,rbind(move_df_sim,move_df_est))

move.long<-melt(movement.data, id=c("Dat","Year","Reg_from","Reg_to","Age"))
move.long$Reg_from<-as.character(move.long$Reg_from)
move.long$Reg_to<-as.character(move.long$Reg_to)
move.long$Age<-as.character(move.long$Age)

#separate again for plotting
move.est<-move.long[move.long$Dat=="EST",]
move.sim<-move.long[move.long$Dat=="SIM",]

#calculate the percent bias
move.est$val.true<-move.sim$value
move.est$bias=((move.est$value-move.est$val.true)/move.est$val.true)*100
#move.est$bias=(move.est$val.true-move.est$value)


#calc medians table
move.est.med <- move.est %>% group_by(Reg_from,Reg_to,Age,Year) %>%
  summarise(med.est=round(median(value),digits=2),med.sim=round(median(val.true),digits=2), med.bias = round(median(bias),digits=2),
            MARE=round(median(abs(bias)),digits=2),
            mean.est=round(mean(value),digits=2),mean.sim=round(mean(val.true),digits=2), MRE = round(mean(bias),digits=2),
            max=round(max(bias),digits=2),min=round(min(bias),digits=2))

move_sum_table_year_age<-c(get(paste0("nm_",z)),round(get(paste0("conv.rate_",z)),digits=2),
                  round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2'])),digits=2),
                  round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1'])),digits=2),
                  round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1'])),digits=2),
                  round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2'])),digits=2),
                  round(median(abs(move.est$bias)),digits=2),
                  round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2']),digits=2),
                  round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1']),digits=2),
                  round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1']),digits=2),
                  round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2']),digits=2),
                  round(mean(move.est$bias),digits=2)
                 )
names(move_sum_table_year_age)<-c("Run","Conv_Rate","MARE_YR+Age_Reg_1_Reg_2","MARE_YR+Age_Reg_2_Reg_1","MARE_YR+Age_Reg_1_Reg_1","MARE_YR+Age_Reg_2_Reg_2","MARE_YR+Reg+Age",
                         "MRE_YR+Age_Reg_1_Reg_2","MRE_YR+Age_Reg_2_Reg_1","MRE_YR+Age_Reg_1_Reg_1","MRE_YR+Age_Reg_2_Reg_2","MRE_YR+Reg+Age"
)

move_sum_table_age<-c(get(paste0("nm_",z)),round(get(paste0("conv.rate_",z)),digits=2),
                       round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='1'])),digits=2),
                       round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='1'])),digits=2),
                       round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='1'])),digits=2),
                       round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='1'])),digits=2),
                       round(median(abs(move.est$bias[move.est$Age=='1'])),digits=2),
                       round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='1']),digits=2),
                       round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='1']),digits=2),
                       round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='1']),digits=2),
                       round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='1']),digits=2),
                       round(mean(move.est$bias[move.est$Age=='1']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='2'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='2'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='2'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='2'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='2'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='2']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='2']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='2']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='2']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='2']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='3'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='3'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='3'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='3'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='3'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='3']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='3']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='3']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='3']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='3']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='4'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='4'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='4'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='4'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='4'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='4']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='4']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='4']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='4']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='4']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='5'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='5'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='5'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='5'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='5'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='5']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='5']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='5']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='5']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='5']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='6'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='6'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='6'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='6'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='6'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='6']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='6']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='6']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='6']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='6']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='7'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='7'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='7'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='7'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='7'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='7']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='7']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='7']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='7']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='7']),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='8'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='8'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='8'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='8'])),digits=2),
                      round(median(abs(move.est$bias[move.est$Age=='8'])),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='2' & move.est$Age=='8']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='1' & move.est$Age=='8']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='1' & move.est$Reg_to=='1' & move.est$Age=='8']),digits=2),
                      round(mean(move.est$bias[move.est$Reg_from=='2' & move.est$Reg_to=='2' & move.est$Age=='8']),digits=2),
                      round(mean(move.est$bias[move.est$Age=='8']),digits=2)
)
names(move_sum_table_age)<-c("Run","Conv_Rate","MARE_Age_1_Reg_1_Reg_2","MARE_Age_1_Reg_2_Reg_1","MARE_Age_1_Reg_1_Reg_1","MARE_Age_1_Reg_2_Reg_2","MARE_YR+Reg_Age_1",
                         "MRE_Age_1_Reg_1_Reg_2","MRE_Age_1_Reg_2_Reg_1","MRE_Age_1_Reg_1_Reg_1","MRE_Age_1_Reg_2_Reg_2","MRE_YR+RegAge_1",
                         "MARE_Age_2_Reg_1_Reg_2","MARE_Age_2_Reg_2_Reg_1","MARE_Age_2_Reg_1_Reg_1","MARE_Age_2_Reg_2_Reg_2","MARE_YR+Reg_Age_2",
                         "MRE_Age_2_Reg_1_Reg_2","MRE_Age_2_Reg_2_Reg_1","MRE_Age_2_Reg_1_Reg_1","MRE_Age_2_Reg_2_Reg_2","MRE_YR+RegAge_2"
                         ,"MARE_Age_3_Reg_1_Reg_2","MARE_Age_3_Reg_2_Reg_1","MARE_Age_3_Reg_1_Reg_1","MARE_Age_3_Reg_2_Reg_2","MARE_YR+Reg_Age_3",
                         "MRE_Age_3_Reg_1_Reg_2","MRE_Age_3_Reg_2_Reg_1","MRE_Age_3_Reg_1_Reg_1","MRE_Age_3_Reg_2_Reg_2","MRE_YR+RegAge_3"
                         ,"MARE_Age_4_Reg_1_Reg_2","MARE_Age_4_Reg_2_Reg_1","MARE_Age_4_Reg_1_Reg_1","MARE_Age_4_Reg_2_Reg_2","MARE_YR+Reg_Age_4",
                         "MRE_Age_4_Reg_1_Reg_2","MRE_Age_4_Reg_2_Reg_1","MRE_Age_4_Reg_1_Reg_1","MRE_Age_4_Reg_2_Reg_2","MRE_YR+RegAge_4"
                         ,"MARE_Age_5_Reg_1_Reg_2","MARE_Age_5_Reg_2_Reg_1","MARE_Age_5_Reg_1_Reg_1","MARE_Age_5_Reg_2_Reg_2","MARE_YR+Reg_Age_5",
                         "MRE_Age_5_Reg_1_Reg_2","MRE_Age_5_Reg_2_Reg_1","MRE_Age_5_Reg_1_Reg_1","MRE_Age_5_Reg_2_Reg_2","MRE_YR+RegAge_5"
                         ,"MARE_Age_6_Reg_1_Reg_2","MARE_Age_6_Reg_2_Reg_1","MARE_Age_6_Reg_1_Reg_1","MARE_Age_6_Reg_2_Reg_2","MARE_YR+Reg_Age_6",
                         "MRE_Age_6_Reg_1_Reg_2","MRE_Age_6_Reg_2_Reg_1","MRE_Age_6_Reg_1_Reg_1","MRE_Age_6_Reg_2_Reg_2","MRE_YR+RegAge_6"
                         ,"MARE_Age_7_Reg_1_Reg_2","MARE_Age_7_Reg_2_Reg_1","MARE_Age_7_Reg_1_Reg_1","MARE_Age_7_Reg_2_Reg_2","MARE_YR+Reg_Age_7",
                         "MRE_Age_7_Reg_1_Reg_2","MRE_Age_7_Reg_2_Reg_1","MRE_Age_7_Reg_1_Reg_1","MRE_Age_7_Reg_2_Reg_2","MRE_YR+RegAge_7"
                         ,"MARE_Age_8_Reg_1_Reg_2","MARE_Age_8_Reg_2_Reg_1","MARE_Age_8_Reg_1_Reg_1","MARE_Age_8_Reg_2_Reg_2","MARE_YR+Reg_Age_8",
                         "MRE_Age_8_Reg_1_Reg_2","MRE_Age_8_Reg_2_Reg_1","MRE_Age_8_Reg_1_Reg_1","MRE_Age_8_Reg_2_Reg_2","MRE_YR+RegAge_8"
)



move.plot.gg<-ggplot(move.est, aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
  geom_violin(aes(fill=Reg_to),trim=T,bw='nrd0',position = position_dodge(width=0.8), alpha=0.2)+
  geom_line(data = move.est.med, aes(x=Year,y=med.est, group=Reg_to))+
  geom_point(data = move.est.med, aes(x=Year,y=med.est, group=Reg_to),position = position_dodge(width=0.8), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = move.est.med, aes(x=Year,y=med.sim, group=Reg_to),lty=2)+
  geom_point(data = move.est.med, aes(x=Year,y=med.sim, group=Reg_to, fill=Reg_to),position = position_dodge(width=0.8),shape=21,size=1.0) + 
  scale_color_grey()+ #"Move To",palette = "Set1")+
  scale_fill_grey()+ #scale_fill_brewer("Move To",palette = "Set1")+
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Movement Rate")+
  xlab("Year")+
  facet_grid(Reg_from~Age)+
  #ylim(-5,5)+
  my_theme

move.name<-c('1'="Population 1 Movement")
move.plot.gg.fig4<-ggplot(subset(move.est, Reg_from %in% '1'), aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
  geom_violin(aes(fill=Reg_to),trim=T,bw='nrd0',position = position_dodge(width=0.6), alpha=0.8,scale='width')+
  geom_point(data = subset(move.est.med, Reg_from %in% '1'), aes(x=Year,y=med.est, group=Reg_to),
             position = position_dodge(width=0.8), 
             fill='black',colour='black', shape=20,size=.75) + 
  geom_line(data = subset(move.est.med, Reg_from %in% '1'), aes(x=Year,y=med.sim, group=Reg_to),lty=1)+
  scale_color_manual(values=c("black", median.col))+
  scale_fill_manual(values=c("black", median.col))+
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  ylab('Proportion Moving')+
  xlab(NULL)+
  #facet_grid(.~Reg_from,labeller=as_labeller(move.name))+
  my_theme.fig4

move.bias.gg<-ggplot(move.est, aes(x=as.factor(Year), y=bias, col=Reg_to), group=Reg_to) +
  geom_violin(aes(fill=Reg_to),trim=T,bw="nrd0",position = position_dodge(width=0.8), alpha=0.2)+
  geom_hline(aes(yintercept = 0, group = Reg_to), colour = 'black',size=0.5,lty=2)+
  geom_line(data = move.est.med, aes(x=Year,y=med.bias, group=Reg_to))+
  geom_point(data = move.est.med, aes(x=Year,y=med.bias,group=Reg_to),position = position_dodge(width=0.8), fill=median.col, shape=21,size=1.5) + 
  scale_color_grey()+ #"Move To",palette = "Set1")+
  scale_fill_grey()+ #scale_fill_brewer("Move To",palette = "Set1")+
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle(get(paste0("nm_",z)))+
  labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg_from~.)+
  ylim(-100,100)+
  my_theme



assign(paste0("rave.plot.gg_",z),rave.plot.gg)
assign(paste0("rave.bias.gg_",z),rave.bias.gg)
assign(paste0("rave.est.med_",z),median_R_ave)
assign(paste0("rave.est_",z),Rave_est)
assign(paste0("rave.sum.table_",z),as.data.frame(t(rave_sum_table)))

assign(paste0("rec.plot.gg_",z),rec.plot.gg)
assign(paste0("rec.plot.gg.fig4_",z),rec.plot.gg.fig4)
assign(paste0("rec.bias.gg_",z),rec.bias.gg)
assign(paste0("rec.est.med_",z),rec.est.med)
assign(paste0("rec.est_",z),rec.est)
assign(paste0("rec.sum.table_",z),as.data.frame(t(rec_sum_table)))

assign(paste0("bio.plot.gg_",z),bio.plot.gg)
assign(paste0("bio.plot.gg.fig5a_",z),bio.plot.gg.fig5a)
assign(paste0("bio.bias.gg_",z),bio.bias.gg)
assign(paste0("bio.est.med_",z),bio.est.med)
assign(paste0("bio.est_",z),bio.est)
assign(paste0("bio.sum.table_",z),as.data.frame(t(bio_sum_table)))

assign(paste0("fmax.plot.gg_",z),fmax.plot.gg)
assign(paste0("fmax.plot.gg.fig4_",z),fmax.plot.gg.fig4)
assign(paste0("fmax.plot.gg.fig5a_",z),fmax.plot.gg.fig5a)
assign(paste0("fmax.bias.gg_",z),fmax.bias.gg)
assign(paste0("fmax.est.med_",z),fmax.est.med)
assign(paste0("fmax.est_",z),fmax.est)
assign(paste0("fmax.sum.table_",z),as.data.frame(t(fmax_sum_table)))

assign(paste0("move.plot.gg_",z),move.plot.gg)
assign(paste0("move.plot.gg.fig4_",z),move.plot.gg.fig4)
assign(paste0("move.bias.gg_",z),move.bias.gg)
assign(paste0("move.est.med_",z),move.est.med)
assign(paste0("move.est_",z),move.est)
assign(paste0("move.sum.table.age_",z),as.data.frame(t(move_sum_table_age)))
assign(paste0("move.sum.table.year.age_",z),as.data.frame(t(move_sum_table_year_age)))


bio.table<-spread(bio.long, key = variable, value = value)
rec.table<-spread(rec.long, key = variable, value = value)
fmax.table<-spread(fmax.long, key = variable, value = value)
move.table<-spread(move.long, key = variable, value = value)

setwd(Figure_table_dir)

write.csv(bio.table,paste0("bio_table_",z,".csv"))
write.csv(rec.table,paste0("rec_table_",z,".csv"))
write.csv(fmax.table,paste0("fmax_table_",z,".csv"))
write.csv(move.table,paste0("move_table_",z,".csv"))

}

library(plyr) #need this for below calcs, but does piping issues with above calcs 
               #(can't have it loaded when doing %>% statements)

setwd(Figure_table_dir)

rave.summary<-sapply(1:num_dir, function(i) list(get(paste0("rave.sum.table_",i))))
rave.summary<-rbindlist(rave.summary)
write.csv(rave.summary,"RAVE MARE and MRE.csv")

rec.summary<-sapply(1:num_dir, function(i) list(get(paste0("rec.sum.table_",i))))
rec.summary<-rbindlist(rec.summary)
write.csv(rec.summary,"Recruitment MARE and MRE.csv")

bio.summary<-sapply(1:num_dir, function(i) list(get(paste0("bio.sum.table_",i))))
bio.summary<-rbindlist(bio.summary)
write.csv(bio.summary,"Biomass MARE and MRE.csv")

fmax.summary<-sapply(1:num_dir, function(i) list(get(paste0("fmax.sum.table_",i))))
fmax.summary<-rbindlist(fmax.summary)
write.csv(fmax.summary,"Fishing Mortality MARE and MRE.csv")

move.summary.age<-sapply(1:num_dir, function(i) list(get(paste0("move.sum.table.age_",i))))
move.summary.age<-rbindlist(move.summary.age)
write.csv(move.summary.age,"Movement MARE and MRE by Age.csv")

move.summary.year.age<-sapply(1:num_dir, function(i) list(get(paste0("move.sum.table.year.age_",i))))
move.summary.year.age<-rbindlist(move.summary.year.age)
write.csv(move.summary.year.age,"Movement MARE and MRE Age+Year Combined.csv")

MARE.table<-(cbind(bio.summary[,c(1:2,9:11)],fmax.summary[,c(7:8)],rec.summary[,c(9:11)],move.summary.year.age[,c(3:4)]))
MRE.table<-(cbind(bio.summary[,c(1:2,19:21)],fmax.summary[,c(14:15)],rec.summary[,c(19:21)],move.summary.year.age[,c(8:9)]))
write.csv(MARE.table,"MARE.csv")
write.csv(MRE.table,"MRE.csv")

#################################################################################################
setwd(Figure_dir)



MARE.graph<-MARE.table[,-c(1:2)]
MRE.graph<-MRE.table[,-c(1:2)]
MARE.graph<-data.frame(sapply(MARE.graph, function(x) as.numeric(as.character(x))))
MRE.graph<-data.frame(sapply(MRE.graph, function(x) as.numeric(as.character(x))))


MRE.bio.graph<-sapply(1:num_dir, function(i) list(get(paste0("bio.est_",i))))
names(MRE.bio.graph)<-sapply(1:num_dir, function(i) get(paste0("nm_",i)))
MRE.bio.graph<-ldply(MRE.bio.graph,data.frame)
MRE.bio.graph$parameter<-'B'

MRE.rec.graph<-sapply(1:num_dir, function(i) list(get(paste0("rec.est_",i))))
names(MRE.rec.graph)<-sapply(1:num_dir, function(i) get(paste0("nm_",i)))
MRE.rec.graph<-ldply(MRE.rec.graph,data.frame)
MRE.rec.graph$parameter<-'R'

MRE.f.graph<-sapply(1:num_dir, function(i) list(get(paste0("fmax.est_",i))))
names(MRE.f.graph)<-sapply(1:num_dir, function(i) get(paste0("nm_",i)))
MRE.f.graph<-ldply(MRE.f.graph,data.frame)
MRE.f.graph$parameter<-'Z_F'

MRE.T.graph<-sapply(1:num_dir, function(i) list(get(paste0("move.est_",i))))
names(MRE.T.graph)<-sapply(1:num_dir, function(i) get(paste0("nm_",i)))
MRE.T.graph<-ldply(MRE.T.graph,data.frame)
MRE.T.graph$parameter<-'Z_T'
MRE.T.graph<-MRE.T.graph[-which(MRE.T.graph$Reg_from==MRE.T.graph$Reg_to),]
MRE.T.graph<-MRE.T.graph[,-5] #remove 'reg_to' column
MRE.T.graph<-dplyr::rename(MRE.T.graph,"Reg"="Reg_from")
MRE.T.graph<-dplyr::rename(MRE.T.graph,"Years"="Year")
MRE.T.graph<-MRE.T.graph[,c(1:3,5:6,4,7:10)]

region.names<-c('1'="Population 1", '2'="Population 2",'System'="Metapopulation",'B'="Biomass",
                'Z_F'="Fishing Mortality",'R'="Recruitment",'Z_T'="Movement")

########################################################################################################
########################################################################################################
########################################################################################################


my_theme.hor<-
  (#theme_bw()+
    theme(
      #panel.grid.major = element_blank(), 
      #panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="grey20",size=7,angle=50,hjust=1,vjust=1,face="plain"),
      axis.text.y = element_text(colour="grey20",size=7,angle=0,hjust=.5,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=9),
      axis.title.y = element_text(size=9),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,size=10),
      plot.subtitle=element_text(hjust=0.5,size=8)
    ))

my_theme.ver<-
  (#theme_bw()+
    theme(
      #panel.grid.major = element_blank(), 
      #panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="grey20",size=7,angle=50,hjust=1,vjust=1,face="plain"),
      axis.text.y = element_text(colour="grey20",size=7,angle=0,hjust=1,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=9),
      axis.title.y = element_text(size=9),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,size=10),
      plot.subtitle=element_text(hjust=0.5,size=8)
    ))

my_theme.ver.alt<-
  (#theme_bw()+
    theme(
      #panel.grid.major = element_blank(), 
      #panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="grey20",size=8,angle=0,hjust=0.45,vjust=1,face="plain"),
      axis.text.y = element_text(colour="grey20",size=8,angle=0,hjust=1,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=9),
      axis.title.y = element_text(size=9),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,size=10),
      plot.subtitle=element_text(hjust=0.5,size=8),
      strip.text.x = element_text(
        size = 10, color = "black", face = "bold"),
      strip.text.y = element_text(
        size = 10, color = "black", face = "italic"),
      strip.background = element_rect(
        fill="grey80" )
    ))
if(any(grepl("package:plyr", search()))) detach("package:plyr") else message("plyr not loaded")

########################################################################################################
####### FINAL VIOLIN PLOTS MANUSCRIPT #################################################################################################
########################################################################################################
names.keep<-c(sapply(runs.to.include, function(i) get(paste0("nm_",i))))
MRE.bio.graph2<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep,]
MRE.f.graph2<-MRE.f.graph[MRE.f.graph$.id %in% names.keep,]
MRE.rec.graph2<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep,]
MRE.T.graph2<-MRE.T.graph[MRE.T.graph$.id %in% names.keep,]
MRE.T.graph2$Reg[which(MRE.T.graph2$Reg=="1")]<-"Population 1"
MRE.T.graph2$Reg[which(MRE.T.graph2$Reg=="2")]<-"Population 2"

MRE.graph.final<-rbind.data.frame(MRE.bio.graph2,MRE.rec.graph2,MRE.f.graph2) #,MRE.T.graph2)
MRE.graph.final<-MRE.graph.final[-which(MRE.graph.final$Reg=='System'),]

MRE.bio.graph.term<-MRE.bio.graph2[MRE.bio.graph2$Years ==nyrs,]
MRE.f.graph.term<-MRE.f.graph2[MRE.f.graph2$Years ==nyrs,]
MRE.rec.graph.term<-MRE.rec.graph2[MRE.rec.graph2$Years ==nyrs,]
MRE.T.graph.term<-MRE.T.graph2[MRE.T.graph2$Years ==nyrs,]

MRE.graph.final.term<-rbind.data.frame(MRE.bio.graph.term,MRE.rec.graph.term,MRE.f.graph.term) #,MRE.T.graph.term)
MRE.graph.final.term<-MRE.graph.final.term[-which(MRE.graph.final.term$Reg=='System'),]

MRE.graph.final$.id<-ordered(MRE.graph.final$.id,levels=names.keep)
MRE.graph.final.term$.id<-ordered(MRE.graph.final.term$.id,levels=names.keep)

MRE.T.graph2$.id<-ordered(MRE.T.graph2$.id,levels=names.keep)
MRE.T.graph.term$.id<-ordered(MRE.T.graph.term$.id,levels=names.keep)

#############################################
MRE.graph.final<-MRE.graph.final %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot<-ggplot(MRE.graph.final, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("Fig 3_Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot)
dev.off()
#############################################

#############################################
MRE.graph.final.term<-MRE.graph.final.term %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term<-ggplot(MRE.graph.final.term, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("Fig 4_Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term)
dev.off()







#############################################
########################################################################################################
########################################################################################################
######## Movement Bias by Age ################################################################################################
########################################################################################################
#############################

#############################################
MRE.T.graph2<-MRE.T.graph2 %>% group_by(Reg,parameter,.id,Age) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

age.names<-c('1'="Age-1", '2'="Age-2",'3'="Age-3", '4'="Age-4",'5'="Age-5", '6'="Age-6",'7'="Age-7", '8'="Age-8",'Population 1'="Population 1",'Population 2'="Population 2")

sinlge.plot.age.T<-ggplot(MRE.T.graph2, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~Age,scales="free_x",labeller= as_labeller(age.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.T.graph2$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("Fig 5_Bias Plots Movement by Age All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.age.T)
dev.off()
#############################################

#############################################
MRE.graph.final.term<-MRE.graph.final.term %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term<-ggplot(MRE.graph.final.term, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("Fig 5_Bias Plots Movement by Age Terminal Year Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term)
dev.off()












#############################################
########################################################################################################
########################################################################################################
######## Figure 5 Timeseries for Various Tag Timeseries################################################################################################
########################################################################################################
########################################################################################################
tiff("Fig 5_Timeseries of R,F,T Tag Timeseries Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 3,top=textGrob("Timeseries Estimates",gp=gpar(fontsize=15,font=3)), 
             #as.table=FALSE,
             bottom=textGrob("Year",gp=gpar(fontsize=12,font=1)), 
             arrangeGrob(move.plot.gg.fig4_1,move.plot.gg.fig4_10,move.plot.gg.fig4_15,
                         move.plot.gg.fig4_29,move.plot.gg.fig4_3, top=textGrob("Movement",gp=gpar(fontsize=9,font=2),hjust=0.21),nrow=5), 
             arrangeGrob(rec.plot.gg.fig4_1,rec.plot.gg.fig4_10,rec.plot.gg.fig4_15,
                         rec.plot.gg.fig4_29,rec.plot.gg.fig4_3, top=textGrob("Recruitment",gp=gpar(fontsize=9,font=2),hjust=.225),nrow=5), 
             arrangeGrob(fmax.plot.gg.fig4_1,fmax.plot.gg.fig4_10,fmax.plot.gg.fig4_15,
                        fmax.plot.gg.fig4_29,fmax.plot.gg.fig4_3,top=textGrob("Fishing Mortality",gp=gpar(fontsize=9,font=2),hjust=.37),nrow=5)
             )
dev.off()


########################################################################################################
########################################################################################################
######## Figure 6 Timeseries for Tag Mixing ################################################################################################
########################################################################################################
########################################################################################################
tiff("Fig 6_Timeseries of F,B Tag Mixing Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 3,top=textGrob("Timeseries Estimates",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Year",gp=gpar(fontsize=12,font=1)), 
             arrangeGrob(fmax.plot.gg.fig4_51,fmax.plot.gg.fig4_53, top=textGrob("Population 1 Fishing Mortality",gp=gpar(fontsize=9,font=2),hjust=0.4),nrow=2), 
             arrangeGrob(fmax.plot.gg.fig5a_51,fmax.plot.gg.fig5a_53, top=textGrob("Population 2 Fishing Mortality",gp=gpar(fontsize=9,font=2),hjust=.4),nrow=2), 
             arrangeGrob(bio.plot.gg.fig5a_51,bio.plot.gg.fig5a_53,top=textGrob("Metapopulation Biomass",gp=gpar(fontsize=9,font=2),hjust=.38),nrow=2)
)
dev.off()


########################################################################################################
########################################################################################################
######## Figure 7 Tag Plot ################################################################################################
########################################################################################################
########################################################################################################

runs.to.include.tag.plot<-c(10,6,16,19,13,1,25,26,15,51,35,12,8,18,7,29,22,14,21,11,17,3)
names.keep.tag.plot<-c(sapply(runs.to.include.tag.plot, function(i) get(paste0("nm_",i))))
tag.plot.move<-move.summary[move.summary$Run %in% names.keep.tag.plot,]
tag.plot.fmax<-fmax.summary[fmax.summary$Run %in% names.keep.tag.plot,]
tag.plot.bio<-bio.summary[bio.summary$Run %in% names.keep.tag.plot,]
tag.plot.rec<-rec.summary[rec.summary$Run %in% names.keep.tag.plot,]
tag.plot.cost<-c(14,0,30,8,8,60,2,10,20,2,10,20,2,10,20,7,7,17.5,10.5,3,0,14)
tag.plot.move<-tag.plot.move[,c(1,7)]
tag.plot.move$parameter<-'T'
tag.plot.move$cost<-tag.plot.cost
tag.plot.fmax<-tag.plot.fmax[,c(1,9)]
tag.plot.fmax$parameter<-'F'
tag.plot.fmax$cost<-tag.plot.cost
tag.plot.rec<-tag.plot.rec[,c(1,12)]
tag.plot.rec$parameter<-'R'
tag.plot.rec$cost<-tag.plot.cost
tag.plot.bio<-tag.plot.bio[,c(1,12)]
tag.plot.bio$parameter<-'B'
tag.plot.bio$cost<-tag.plot.cost
tag.plot<-rbind.data.frame(tag.plot.move,tag.plot.bio,tag.plot.fmax,tag.plot.rec)
names(tag.plot)<-c('Run','MARE','Par','cost')
tag.plot$MARE<-data.frame(sapply(tag.plot$MARE, function(x) as.numeric(as.character(x))))
tag.plot$cost<-data.frame(sapply(tag.plot$cost, function(x) as.numeric(as.character(x))))
tag.plot$color<-FALSE
tag.plot[tag.plot$Run=="Base",5]<-TRUE
tag.plot[tag.plot$Run=="No_Tag",5]<-TRUE


require(ggrepel)
my_theme.fig5<-
  (
    theme(
      plot.margin = unit(c(0, .15, 0, .15), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="black",size=7,angle=0,hjust=0.45,vjust=1,face="plain"),
      axis.text.y = element_text(colour="black",size=7,angle=0,hjust=1,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=9),
      axis.title.y = element_text(size=9),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,vjust=0,size=12,margin=margin(2,2,7,2)),
      plot.subtitle=element_text(hjust=0.5,size=8),
      strip.text.x = element_text(
        size = 9, color = "black", face = "bold"),
      strip.text.y = element_text(
        size = 9, color = "black", face = "italic"),
      strip.background = element_rect(
        fill="grey80",linetype='solid',colour='black',size=1. ),
      legend.position = c(0.5, 0.5),
      legend.text=element_text(size=5),
      legend.title=element_text(size=5),
      legend.key.size = unit(.15,"cm"),
      panel.background = element_blank(), #element_rect(colour = "black", size=.75, fill=NA),
      #axis.line.x = element_line(),
      #axis.line.y = element_line(),
      legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid',size=.25),
      #legend.title=element_text(),
      legend.margin=margin(c(1,2,1,2)),
      panel.border = element_rect(colour = "black", fill=NA, size=.75)
    ))

tag.name<-c('B'="Biomass", 'F'="Fishing Mortality",'R'="Recruitment",'T'="Movement")

tiff("Fig 7_MARE vs. Tag Cost Final_ALT.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
ggplot(tag.plot, aes(x=MARE, y=cost)) +
  geom_point()+ #fill='black', colour='black',shape=20,size=.75) + 
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_label_repel(aes(fill=tag.plot$color,label=tag.plot$Run),label.padding=.14, #vjust=1,
                  size = 2.2,force=25, max.iter=50000, point.padding=.2,segment.size = .3,min.segment.length = unit(.1, 'lines'),
                  show.legend = FALSE,fontface='bold',colour='white',segment.colour='black')+
  ggtitle('Tagging Program Relative Cost and Associated MARE')+
  ylab('Relative Cost')+
  xlab('MARE')+
  facet_wrap(Par~.,scales="free",labeller=as_labeller(tag.name))+
  ylim(-7,64)+
  scale_fill_manual(values = c('grey50','black'))+

  my_theme.fig5
dev.off()





###################################################################################################
########################################################################################################
######## SM Figures ################################################################################################
################################################################################################################################################################################################################
########################################################################################################
########################################################################################################
######### Tag Timeseries ###############################################################################################
########################################################################################################
########################################################################################################

runs.to.include.alt1<-c(1,3,10,6:8,11:19,29)  # the numbers corresponding to the runs to include in the MRE summary plots

names.keep.alt1<-c(sapply(runs.to.include.alt1, function(i) get(paste0("nm_",i))))
MRE.bio.graph2.alt1<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep.alt1,]
MRE.f.graph2.alt1<-MRE.f.graph[MRE.f.graph$.id %in% names.keep.alt1,]
MRE.rec.graph2.alt1<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep.alt1,]
MRE.T.graph2.alt1<-MRE.T.graph[MRE.T.graph$.id %in% names.keep.alt1,]

MRE.graph.final.alt1<-rbind.data.frame(MRE.bio.graph2.alt1,MRE.rec.graph2.alt1,MRE.f.graph2.alt1,MRE.T.graph2.alt1)
MRE.graph.final.alt1<-MRE.graph.final.alt1[-which(MRE.graph.final.alt1$Reg=='System'),]

MRE.bio.graph.term.alt1<-MRE.bio.graph2.alt1[MRE.bio.graph2.alt1$Years ==nyrs,]
MRE.f.graph.term.alt1<-MRE.f.graph2.alt1[MRE.f.graph2.alt1$Years ==nyrs,]
MRE.rec.graph.term.alt1<-MRE.rec.graph2.alt1[MRE.rec.graph2.alt1$Years ==nyrs,]
MRE.T.graph.term.alt1<-MRE.T.graph2.alt1[MRE.T.graph2.alt1$Years ==nyrs,]

MRE.graph.final.term.alt1<-rbind.data.frame(MRE.bio.graph.term.alt1,MRE.rec.graph.term.alt1,MRE.f.graph.term.alt1,MRE.T.graph.term.alt1)
MRE.graph.final.term.alt1<-MRE.graph.final.term.alt1[-which(MRE.graph.final.term.alt1$Reg=='System'),]

MRE.graph.final.alt1$.id<-ordered(MRE.graph.final.alt1$.id,levels=names.keep.alt1)
MRE.graph.final.term.alt1$.id<-ordered(MRE.graph.final.term.alt1$.id,levels=names.keep.alt1)

#############################################
MRE.graph.final.alt1<-MRE.graph.final.alt1 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.alt1<-ggplot(MRE.graph.final.alt1, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.alt1$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_1_Tag_Timeseries Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.alt1)
dev.off()
#############################################

#############################################
MRE.graph.final.term.alt1<-MRE.graph.final.term.alt1 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term.alt1<-ggplot(MRE.graph.final.term.alt1, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term.alt1$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_2_Tag_Timeseries Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term.alt1)
dev.off()
#############################################

########################################################################################################
########################################################################################################
######### Tag Deployment ###############################################################################################
########################################################################################################
########################################################################################################

runs.to.include.alt2<-c(1,3,20:29)  # the numbers corresponding to the runs to include in the MRE summary plots

names.keep.alt2<-c(sapply(runs.to.include.alt2, function(i) get(paste0("nm_",i))))
MRE.bio.graph2.alt2<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep.alt2,]
MRE.f.graph2.alt2<-MRE.f.graph[MRE.f.graph$.id %in% names.keep.alt2,]
MRE.rec.graph2.alt2<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep.alt2,]
MRE.T.graph2.alt2<-MRE.T.graph[MRE.T.graph$.id %in% names.keep.alt2,]

MRE.graph.final.alt2<-rbind.data.frame(MRE.bio.graph2.alt2,MRE.rec.graph2.alt2,MRE.f.graph2.alt2,MRE.T.graph2.alt2)
MRE.graph.final.alt2<-MRE.graph.final.alt2[-which(MRE.graph.final.alt2$Reg=='System'),]

MRE.bio.graph.term.alt2<-MRE.bio.graph2.alt2[MRE.bio.graph2.alt2$Years ==nyrs,]
MRE.f.graph.term.alt2<-MRE.f.graph2.alt2[MRE.f.graph2.alt2$Years ==nyrs,]
MRE.rec.graph.term.alt2<-MRE.rec.graph2.alt2[MRE.rec.graph2.alt2$Years ==nyrs,]
MRE.T.graph.term.alt2<-MRE.T.graph2.alt2[MRE.T.graph2.alt2$Years ==nyrs,]

MRE.graph.final.term.alt2<-rbind.data.frame(MRE.bio.graph.term.alt2,MRE.rec.graph.term.alt2,MRE.f.graph.term.alt2,MRE.T.graph.term.alt2)
MRE.graph.final.term.alt2<-MRE.graph.final.term.alt2[-which(MRE.graph.final.term.alt2$Reg=='System'),]

MRE.graph.final.alt2$.id<-ordered(MRE.graph.final.alt2$.id,levels=names.keep.alt2)
MRE.graph.final.term.alt2$.id<-ordered(MRE.graph.final.term.alt2$.id,levels=names.keep.alt2)

#############################################
MRE.graph.final.alt2<-MRE.graph.final.alt2 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.alt2<-ggplot(MRE.graph.final.alt2, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.alt2$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_3_Tag_Deployment Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.alt2)
dev.off()
#############################################

#############################################
MRE.graph.final.term.alt2<-MRE.graph.final.term.alt2 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term.alt2<-ggplot(MRE.graph.final.term.alt2, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term.alt2$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_4_Tag_Deployment Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term.alt2)
dev.off()
#############################################




########################################################################################################
########################################################################################################
######### Data Assumptions ###############################################################################################
########################################################################################################
########################################################################################################

runs.to.include.alt3<-c(1,3,30,63,40:43,46:55)  # the numbers corresponding to the runs to include in the MRE summary plots

names.keep.alt3<-c(sapply(runs.to.include.alt3, function(i) get(paste0("nm_",i))))
MRE.bio.graph2.alt3<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep.alt3,]
MRE.f.graph2.alt3<-MRE.f.graph[MRE.f.graph$.id %in% names.keep.alt3,]
MRE.rec.graph2.alt3<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep.alt3,]
MRE.T.graph2.alt3<-MRE.T.graph[MRE.T.graph$.id %in% names.keep.alt3,]

MRE.graph.final.alt3<-rbind.data.frame(MRE.bio.graph2.alt3,MRE.rec.graph2.alt3,MRE.f.graph2.alt3,MRE.T.graph2.alt3)
MRE.graph.final.alt3<-MRE.graph.final.alt3[-which(MRE.graph.final.alt3$Reg=='System'),]

MRE.bio.graph.term.alt3<-MRE.bio.graph2.alt3[MRE.bio.graph2.alt3$Years ==nyrs,]
MRE.f.graph.term.alt3<-MRE.f.graph2.alt3[MRE.f.graph2.alt3$Years ==nyrs,]
MRE.rec.graph.term.alt3<-MRE.rec.graph2.alt3[MRE.rec.graph2.alt3$Years ==nyrs,]
MRE.T.graph.term.alt3<-MRE.T.graph2.alt3[MRE.T.graph2.alt3$Years ==nyrs,]

MRE.graph.final.term.alt3<-rbind.data.frame(MRE.bio.graph.term.alt3,MRE.rec.graph.term.alt3,MRE.f.graph.term.alt3,MRE.T.graph.term.alt3)
MRE.graph.final.term.alt3<-MRE.graph.final.term.alt3[-which(MRE.graph.final.term.alt3$Reg=='System'),]

MRE.graph.final.alt3$.id<-ordered(MRE.graph.final.alt3$.id,levels=names.keep.alt3)
MRE.graph.final.term.alt3$.id<-ordered(MRE.graph.final.term.alt3$.id,levels=names.keep.alt3)

#############################################
MRE.graph.final.alt3<-MRE.graph.final.alt3 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.alt3<-ggplot(MRE.graph.final.alt3, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.alt3$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_5_Tag Assumptions Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.alt3)
dev.off()
#############################################

#############################################
MRE.graph.final.term.alt3<-MRE.graph.final.term.alt3 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term.alt3<-ggplot(MRE.graph.final.term.alt3, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term.alt3$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_6_Tag Assumptions Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term.alt3)
dev.off()
#############################################



########################################################################################################
########################################################################################################
######### Data Quality ###############################################################################################
########################################################################################################
########################################################################################################

runs.to.include.alt4<-c(1,3:5,31:32,44:45,33:35)  # the numbers corresponding to the runs to include in the MRE summary plots

names.keep.alt4<-c(sapply(runs.to.include.alt4, function(i) get(paste0("nm_",i))))
MRE.bio.graph2.alt4<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep.alt4,]
MRE.f.graph2.alt4<-MRE.f.graph[MRE.f.graph$.id %in% names.keep.alt4,]
MRE.rec.graph2.alt4<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep.alt4,]
MRE.T.graph2.alt4<-MRE.T.graph[MRE.T.graph$.id %in% names.keep.alt4,]

MRE.graph.final.alt4<-rbind.data.frame(MRE.bio.graph2.alt4,MRE.rec.graph2.alt4,MRE.f.graph2.alt4,MRE.T.graph2.alt4)
MRE.graph.final.alt4<-MRE.graph.final.alt4[-which(MRE.graph.final.alt4$Reg=='System'),]

MRE.bio.graph.term.alt4<-MRE.bio.graph2.alt4[MRE.bio.graph2.alt4$Years ==nyrs,]
MRE.f.graph.term.alt4<-MRE.f.graph2.alt4[MRE.f.graph2.alt4$Years ==nyrs,]
MRE.rec.graph.term.alt4<-MRE.rec.graph2.alt4[MRE.rec.graph2.alt4$Years ==nyrs,]
MRE.T.graph.term.alt4<-MRE.T.graph2.alt4[MRE.T.graph2.alt4$Years ==nyrs,]

MRE.graph.final.term.alt4<-rbind.data.frame(MRE.bio.graph.term.alt4,MRE.rec.graph.term.alt4,MRE.f.graph.term.alt4,MRE.T.graph.term.alt4)
MRE.graph.final.term.alt4<-MRE.graph.final.term.alt4[-which(MRE.graph.final.term.alt4$Reg=='System'),]

MRE.graph.final.alt4$.id<-ordered(MRE.graph.final.alt4$.id,levels=names.keep.alt4)
MRE.graph.final.term.alt4$.id<-ordered(MRE.graph.final.term.alt4$.id,levels=names.keep.alt4)

#############################################
MRE.graph.final.alt4<-MRE.graph.final.alt4 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.alt4<-ggplot(MRE.graph.final.alt4, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.alt4$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_7_Data Quality Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.alt4)
dev.off()
#############################################

#############################################
MRE.graph.final.term.alt4<-MRE.graph.final.term.alt4 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term.alt4<-ggplot(MRE.graph.final.term.alt4, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term.alt4$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_8_Data Quality Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term.alt4)
dev.off()
#############################################



########################################################################################################
########################################################################################################
######### Parametrization ###############################################################################################
########################################################################################################
########################################################################################################

runs.to.include.alt5<-c(1,62,2,37,36,39,57,56,38,58,59)  # the numbers corresponding to the runs to include in the MRE summary plots

names.keep.alt5<-c(sapply(runs.to.include.alt5, function(i) get(paste0("nm_",i))))
MRE.bio.graph2.alt5<-MRE.bio.graph[MRE.bio.graph$.id %in% names.keep.alt5,]
MRE.f.graph2.alt5<-MRE.f.graph[MRE.f.graph$.id %in% names.keep.alt5,]
MRE.rec.graph2.alt5<-MRE.rec.graph[MRE.rec.graph$.id %in% names.keep.alt5,]
MRE.T.graph2.alt5<-MRE.T.graph[MRE.T.graph$.id %in% names.keep.alt5,]

MRE.graph.final.alt5<-rbind.data.frame(MRE.bio.graph2.alt5,MRE.rec.graph2.alt5,MRE.f.graph2.alt5,MRE.T.graph2.alt5)
MRE.graph.final.alt5<-MRE.graph.final.alt5[-which(MRE.graph.final.alt5$Reg=='System'),]

MRE.bio.graph.term.alt5<-MRE.bio.graph2.alt5[MRE.bio.graph2.alt5$Years ==nyrs,]
MRE.f.graph.term.alt5<-MRE.f.graph2.alt5[MRE.f.graph2.alt5$Years ==nyrs,]
MRE.rec.graph.term.alt5<-MRE.rec.graph2.alt5[MRE.rec.graph2.alt5$Years ==nyrs,]
MRE.T.graph.term.alt5<-MRE.T.graph2.alt5[MRE.T.graph2.alt5$Years ==nyrs,]

MRE.graph.final.term.alt5<-rbind.data.frame(MRE.bio.graph.term.alt5,MRE.rec.graph.term.alt5,MRE.f.graph.term.alt5,MRE.T.graph.term.alt5)
MRE.graph.final.term.alt5<-MRE.graph.final.term.alt5[-which(MRE.graph.final.term.alt5$Reg=='System'),]

MRE.graph.final.alt5$.id<-ordered(MRE.graph.final.alt5$.id,levels=names.keep.alt5)
MRE.graph.final.term.alt5$.id<-ordered(MRE.graph.final.term.alt5$.id,levels=names.keep.alt5)

#############################################
MRE.graph.final.alt5<-MRE.graph.final.alt5 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.alt5<-ggplot(MRE.graph.final.alt5, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6, scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA,coef=0)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.alt5$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_9_Parametrization Bias Plots All Years Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Aggregate Bias Across Years",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.alt5)
dev.off()
#############################################

#############################################
MRE.graph.final.term.alt5<-MRE.graph.final.term.alt5 %>% group_by(Reg,parameter,.id) %>% 
  filter(bias<quantile(bias,.999), bias>quantile(bias,.001))  %>% as.data.frame() # 

sinlge.plot.term.alt5<-ggplot(MRE.graph.final.term.alt5, aes(x=.id, y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=0.5,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6,scale="width")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=0.5, color='black')+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(vars(parameter,Reg),dir='v',scales="free_x",labeller="label_both")+ # as_labeller(region.names),)+ #switch="both")+ #,strip.position="top")+
  facet_grid(Reg~parameter,scales="free_x",labeller= as_labeller(region.names))+ #, switch="y")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.graph.final.term.alt5$.id))))+
  #ylim(-100,100)+
  geom_blank()+
  coord_flip()+
  my_theme.ver.alt

##################################
tiff("SM_Fig_10_Parametrization Bias Plots Terminal Final.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 1,top=textGrob("Terminal Year Bias",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Percent Relative Error",gp=gpar(fontsize=12,font=1)), 
             left=textGrob("Scenario",gp=gpar(fontsize=12,font=1),rot=90),
             sinlge.plot.term.alt5)
dev.off()
#############################################












#################################################################

pdf("Simulation_Comparison_Base_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_2,move.plot.gg_3,move.plot.gg_29)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_2,rec.plot.gg_3,rec.plot.gg_29)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_2,bio.plot.gg_3,bio.plot.gg_29)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_2,fmax.plot.gg_3,fmax.plot.gg_29)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Base_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_2,move.bias.gg_3,move.bias.gg_29)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_2,rec.bias.gg_3,rec.bias.gg_29)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_2,bio.bias.gg_3,bio.bias.gg_29)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_2,fmax.bias.gg_3,fmax.bias.gg_29)

dev.off()








#################################################################

pdf("Simulation_Comparison_Base_No_Move_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_2,move.plot.gg_3,move.plot.gg_62)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_2,rec.plot.gg_3,rec.plot.gg_62)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_2,bio.plot.gg_3,bio.plot.gg_62)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_2,fmax.plot.gg_3,fmax.plot.gg_62)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Base_No_Move_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_2,move.bias.gg_3,move.bias.gg_62)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_2,rec.bias.gg_3,rec.bias.gg_62)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_2,bio.bias.gg_3,bio.bias.gg_62)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_2,fmax.bias.gg_3,fmax.bias.gg_62)

dev.off()









setwd(Figure_dir)
pdf("Simulation_Comparison_Move_Par_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_3,move.plot.gg_57,move.plot.gg_36,move.plot.gg_37,move.plot.gg_39)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_3,rec.plot.gg_57,rec.plot.gg_36,rec.plot.gg_37,rec.plot.gg_39)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_3,bio.plot.gg_57,bio.plot.gg_36,bio.plot.gg_37,bio.plot.gg_39)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_3,fmax.plot.gg_57,fmax.plot.gg_36,fmax.plot.gg_37,fmax.plot.gg_39)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Move_Par_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_3,move.bias.gg_57,move.bias.gg_36,move.bias.gg_37,move.bias.gg_39)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_3,rec.bias.gg_57,rec.bias.gg_36,rec.bias.gg_37,rec.bias.gg_39)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_3,bio.bias.gg_57,bio.bias.gg_36,bio.bias.gg_37,bio.bias.gg_39)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_3,fmax.bias.gg_57,fmax.bias.gg_36,fmax.bias.gg_37,fmax.bias.gg_39)

dev.off()








setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timeseries_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_6,move.plot.gg_7,move.plot.gg_8,move.plot.gg_9,move.plot.gg_10)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_6,rec.plot.gg_7,rec.plot.gg_8,rec.plot.gg_9,rec.plot.gg_10)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_6,bio.plot.gg_7,bio.plot.gg_8,bio.plot.gg_9,bio.plot.gg_10)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_6,fmax.plot.gg_7,fmax.plot.gg_8,fmax.plot.gg_9,fmax.plot.gg_10)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timeseries_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_6,move.bias.gg_7,move.bias.gg_8,move.bias.gg_9,move.bias.gg_10)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_6,rec.bias.gg_7,rec.bias.gg_8,rec.bias.gg_9,rec.bias.gg_10)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_6,bio.bias.gg_7,bio.bias.gg_8,bio.bias.gg_9,bio.bias.gg_10)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_6,fmax.bias.gg_7,fmax.bias.gg_8,fmax.bias.gg_9,fmax.bias.gg_10)

dev.off()










setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_Begin_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_11,move.plot.gg_12,move.plot.gg_13)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_11,rec.plot.gg_12,rec.plot.gg_13)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_11,bio.plot.gg_12,bio.plot.gg_13)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_11,fmax.plot.gg_12,fmax.plot.gg_13)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_Begin_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_11,move.bias.gg_12,move.bias.gg_13)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_11,rec.bias.gg_12,rec.bias.gg_13)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_11,bio.bias.gg_12,bio.bias.gg_13)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_11,fmax.bias.gg_12,fmax.bias.gg_13)

dev.off()








setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_Middle_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_14,move.plot.gg_15,move.plot.gg_16)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_14,rec.plot.gg_15,rec.plot.gg_16)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_14,bio.plot.gg_15,bio.plot.gg_16)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_14,fmax.plot.gg_15,fmax.plot.gg_16)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_Middle_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_14,move.bias.gg_15,move.bias.gg_16)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_14,rec.bias.gg_15,rec.bias.gg_16)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_14,bio.bias.gg_15,bio.bias.gg_16)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_14,fmax.bias.gg_15,fmax.bias.gg_16)

dev.off()










setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_End_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_17,move.plot.gg_18,move.plot.gg_19)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_17,rec.plot.gg_18,rec.plot.gg_19)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_17,bio.plot.gg_18,bio.plot.gg_19)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_17,fmax.plot.gg_18,fmax.plot.gg_19)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_End_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_17,move.bias.gg_18,move.bias.gg_19)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_17,rec.bias.gg_18,rec.bias.gg_19)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_17,bio.bias.gg_18,bio.bias.gg_19)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_17,fmax.bias.gg_18,fmax.bias.gg_19)

dev.off()








setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_5_Years_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_12,move.plot.gg_15,move.plot.gg_18)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_12,rec.plot.gg_15,rec.plot.gg_18)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_12,bio.plot.gg_15,bio.plot.gg_18)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_12,fmax.plot.gg_15,fmax.plot.gg_18)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Timing_5_Years_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_12,move.bias.gg_15,move.bias.gg_18)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_12,rec.bias.gg_15,rec.bias.gg_18)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_12,bio.bias.gg_15,bio.bias.gg_18)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_12,fmax.bias.gg_15,fmax.bias.gg_18)

dev.off()









setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Distribution_Reg_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_20,move.plot.gg_21,move.plot.gg_22,move.plot.gg_23,move.plot.gg_29)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_20,rec.plot.gg_21,rec.plot.gg_22,rec.plot.gg_23,rec.plot.gg_29)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_20,bio.plot.gg_21,bio.plot.gg_22,bio.plot.gg_23,bio.plot.gg_29)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_20,fmax.plot.gg_21,fmax.plot.gg_22,fmax.plot.gg_23,fmax.plot.gg_29)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Distribution_Reg_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_20,move.bias.gg_21,move.bias.gg_2,move.bias.gg_23,move.bias.gg_29)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_20,rec.bias.gg_21,rec.bias.gg_22,rec.bias.gg_23,rec.bias.gg_29)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_20,bio.bias.gg_21,bio.bias.gg_22,bio.bias.gg_23,bio.bias.gg_29)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_20,fmax.bias.gg_21,fmax.bias.gg_22,fmax.bias.gg_23,fmax.bias.gg_29)

dev.off()









setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Age_Dist_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_24,move.plot.gg_30)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_24,rec.plot.gg_30)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_24,bio.plot.gg_30)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_24,fmax.plot.gg_30)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Age_Dist_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_24,move.bias.gg_30)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_24,rec.bias.gg_30)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_24,bio.bias.gg_30)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_24,fmax.bias.gg_30)

dev.off()














setwd(Figure_dir)
pdf("Simulation_Comparison_Number_Tags_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_3,move.plot.gg_25,move.plot.gg_26,move.plot.gg_27,move.plot.gg_28)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_3,rec.plot.gg_25,rec.plot.gg_26,rec.plot.gg_27,rec.plot.gg_28)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_3,bio.plot.gg_25,bio.plot.gg_26,bio.plot.gg_27,bio.plot.gg_28)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_3,fmax.plot.gg_25,fmax.plot.gg_26,fmax.plot.gg_27,fmax.plot.gg_28)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Number_Tags_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_3,move.bias.gg_25,move.bias.gg_26,move.bias.gg_27,move.bias.gg_28)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_3,rec.bias.gg_25,rec.bias.gg_26,rec.bias.gg_27,rec.bias.gg_28)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_3,bio.bias.gg_25,bio.bias.gg_26,bio.bias.gg_27,bio.bias.gg_28)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_3,fmax.bias.gg_25,fmax.bias.gg_26,fmax.bias.gg_27,fmax.bias.gg_28)

dev.off()









setwd(Figure_dir)
pdf("Simulation_Comparison_Data_Quality_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_3,move.plot.gg_4,move.plot.gg_5,move.plot.gg_33,move.plot.gg_34)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_3,rec.plot.gg_4,rec.plot.gg_5,rec.plot.gg_33,rec.plot.gg_34)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_3,bio.plot.gg_4,bio.plot.gg_5,bio.plot.gg_33,bio.plot.gg_34)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_3,fmax.plot.gg_4,fmax.plot.gg_5,fmax.plot.gg_33,fmax.plot.gg_34)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Data_Quality_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_3,move.bias.gg_4,move.bias.gg_5,move.bias.gg_33,move.bias.gg_34)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_3,rec.bias.gg_4,rec.bias.gg_5,rec.bias.gg_33,rec.bias.gg_34)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_3,bio.bias.gg_4,bio.bias.gg_5,bio.bias.gg_33,bio.bias.gg_34)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_3,fmax.bias.gg_4,fmax.bias.gg_5,fmax.bias.gg_33,fmax.bias.gg_34)

dev.off()










setwd(Figure_dir)
pdf("Simulation_Comparison_Data_Weighting_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1, move.plot.gg_35,move.plot.gg_44,move.plot.gg_45)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_35,rec.plot.gg_44,rec.plot.gg_45)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_35,bio.plot.gg_44,bio.plot.gg_45)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_35,fmax.plot.gg_44,fmax.plot.gg_45)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Data_Weighting_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1, move.bias.gg_35,move.bias.gg_44,move.bias.gg_45)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_35,rec.bias.gg_44,rec.bias.gg_45)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_35,bio.bias.gg_44,bio.bias.gg_45)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_35,fmax.bias.gg_44,fmax.bias.gg_45)

dev.off()






setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Life_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_40,move.plot.gg_41,move.plot.gg_42,move.plot.gg_43)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_40,rec.plot.gg_41,rec.plot.gg_42,rec.plot.gg_43)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_40,bio.plot.gg_41,bio.plot.gg_42,bio.plot.gg_43)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_40,fmax.plot.gg_41,fmax.plot.gg_42,fmax.plot.gg_43)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Life_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_40,move.bias.gg_41,move.bias.gg_42,move.bias.gg_43)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_40,rec.bias.gg_41,rec.bias.gg_42,rec.bias.gg_43)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_40,bio.bias.gg_41,bio.bias.gg_42,bio.bias.gg_43)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_40,fmax.bias.gg_41,fmax.bias.gg_42,fmax.bias.gg_43)

dev.off()






setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.75_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_46,move.plot.gg_47,move.plot.gg_48,move.plot.gg_49)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_46,rec.plot.gg_47,rec.plot.gg_48,rec.plot.gg_49)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_46,bio.plot.gg_47,bio.plot.gg_48,bio.plot.gg_49)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_46,fmax.plot.gg_47,fmax.plot.gg_48,fmax.plot.gg_49)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.75_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_46,move.bias.gg_47,move.bias.gg_48,move.bias.gg_49)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_46,rec.bias.gg_47,rec.bias.gg_48,rec.bias.gg_49)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_46,bio.bias.gg_47,bio.bias.gg_48,bio.bias.gg_49)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_46,fmax.bias.gg_47,fmax.bias.gg_48,fmax.bias.gg_49)

dev.off()










setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.5_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_50,move.plot.gg_51,move.plot.gg_52,move.plot.gg_53)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_50,rec.plot.gg_51,rec.plot.gg_52,rec.plot.gg_53)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_50,bio.plot.gg_51,bio.plot.gg_52,bio.plot.gg_53)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_50,fmax.plot.gg_51,fmax.plot.gg_52,fmax.plot.gg_53)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.5_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_50,move.bias.gg_51,move.bias.gg_52,move.bias.gg_53)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_50,rec.bias.gg_51,rec.bias.gg_52,rec.bias.gg_53)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_50,bio.bias.gg_51,bio.bias.gg_52,bio.bias.gg_53)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_50,fmax.bias.gg_51,fmax.bias.gg_52,fmax.bias.gg_53)

dev.off()








setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.75_Opp_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_47,move.plot.gg_49,move.plot.gg_54,move.plot.gg_55)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_47,rec.plot.gg_49,rec.plot.gg_54,rec.plot.gg_55)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_47,bio.plot.gg_49,bio.plot.gg_54,bio.plot.gg_55)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_47,fmax.plot.gg_49,fmax.plot.gg_54,fmax.plot.gg_55)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Tag_Mixing_F_0.75_Opp_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_47,move.bias.gg_49,move.bias.gg_54,move.bias.gg_55)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_47,rec.bias.gg_49,rec.bias.gg_54,rec.bias.gg_55)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_47,bio.bias.gg_49,bio.bias.gg_54,bio.bias.gg_55)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_47,fmax.bias.gg_49,fmax.bias.gg_54,fmax.bias.gg_55)

dev.off()








setwd(Figure_dir)
pdf("Simulation_Comparison_Parametrization_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.plot.gg_1, move.plot.gg_56,move.plot.gg_58,move.plot.gg_59,move.plot.gg_60,move.plot.gg_38)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.plot.gg_1, rec.plot.gg_56,rec.plot.gg_58,rec.plot.gg_59,rec.plot.gg_60,rec.plot.gg_38)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.plot.gg_1, bio.plot.gg_56,bio.plot.gg_58,bio.plot.gg_59,bio.plot.gg_60,bio.plot.gg_38)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.plot.gg_1, fmax.plot.gg_56,fmax.plot.gg_58,fmax.plot.gg_59,fmax.plot.gg_60,fmax.plot.gg_38)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Parametrization_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,
             top="Movement",
             move.bias.gg_1, move.bias.gg_56,move.bias.gg_58,move.bias.gg_59,move.bias.gg_60,move.bias.gg_38)


grid.arrange(ncol = 3,
             top="Recruitment",
             rec.bias.gg_1, rec.bias.gg_56,rec.bias.gg_58,rec.bias.gg_59,rec.bias.gg_60,rec.bias.gg_38)


grid.arrange(ncol = 3,
             top="Biomass",
             bio.bias.gg_1, bio.bias.gg_56,bio.bias.gg_58,bio.bias.gg_59,bio.bias.gg_60,bio.bias.gg_38)


grid.arrange(ncol = 3,
             top="Fishing Mortality",
             fmax.bias.gg_1, fmax.bias.gg_56,fmax.bias.gg_58,fmax.bias.gg_59,fmax.bias.gg_60,fmax.bias.gg_38)

dev.off()






setwd(Figure_dir)
pdf("Simulation_Comparison_3_Areas_Values.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.plot.gg_1,move.plot.gg_3, move.plot.gg_62,move.plot.gg_63)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.plot.gg_1,rec.plot.gg_3, rec.plot.gg_62,rec.plot.gg_63)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.plot.gg_1,bio.plot.gg_3, bio.plot.gg_62,bio.plot.gg_63)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.plot.gg_1,fmax.plot.gg_3, fmax.plot.gg_62,fmax.plot.gg_63)

dev.off()



setwd(Figure_dir)
pdf("Simulation_Comparison_Parametrization_Bias.pdf",paper="letter",height = 11, width=8)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 2,
             top="Movement",
             move.bias.gg_1,move.bias.gg_3, move.bias.gg_62,move.bias.gg_63)


grid.arrange(ncol = 2,
             top="Recruitment",
             rec.bias.gg_1,rec.bias.gg_3, rec.bias.gg_62,rec.bias.gg_63)


grid.arrange(ncol = 2,
             top="Biomass",
             bio.bias.gg_1,bio.bias.gg_3, bio.bias.gg_62,bio.bias.gg_63)


grid.arrange(ncol = 2,
             top="Fishing Mortality",
             fmax.bias.gg_1,fmax.bias.gg_3, fmax.bias.gg_62,fmax.bias.gg_63)

dev.off()






bio.plot<-ggplot(MRE.bio.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  facet_wrap(.~Reg,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor

f.plot<-ggplot(MRE.f.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Fishing Mortality")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor

rec.plot<-ggplot(MRE.rec.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_grid(Reg~.,labeller=as_labeller(region.names))+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor

##############################
tiff("Bias Plots Horizontal.tif",width=190,height=240, unit='mm',res=300)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,top="Aggregate Bias Across Years",left="Percent Relative Error", bottom="Scenario",
             bio.plot,rec.plot,f.plot)
dev.off()
##################################################



bio.plot.vert<-ggplot(MRE.bio.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  facet_wrap(.~Reg,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  coord_flip()+
  my_theme.ver

f.plot.vert<-ggplot(MRE.f.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Fishing Mortality")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  coord_flip()+
  my_theme.ver

rec.plot.vert<-ggplot(MRE.rec.graph2, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_grid(Reg~.,labeller=as_labeller(region.names))+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  coord_flip()+
  my_theme.ver+
  theme(axis.text.y = element_blank())

##################################
tiff("Bias Plots Vertical.tif",width=190,height=240, unit='mm',res=300)
par(mar=c(6,6,6,6))
#plot_grid( bio.plot.vert,rec.plot.vert,f.plot.vert, align = "hv", ncol = 3) #, rel_heights = c(1/4, 1/4, 1/2))

grid.arrange(ncol = 3,top="Aggregate Bias Across Years", bottom="Percent Relative Error", left="Scenario",
             bio.plot.vert,rec.plot.vert,f.plot.vert)
dev.off()

######################################################################################################


bio.plot.term<-ggplot(MRE.bio.graph.term, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  facet_wrap(.~Reg,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor

f.plot.term<-ggplot(MRE.f.graph.term, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0",alpha=0.6)+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Fishing Mortality")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top")+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor

rec.plot.term<-ggplot(MRE.rec.graph.term, aes(x=as.factor(.id), y=bias)) +
  geom_hline(yintercept = 0, colour = 'black',size=1,lty=2)+
  geom_violin(trim=T,bw="nrd0")+  #fill=vio.col,
  geom_boxplot(width=.35,outlier.shape=NA)+
  #geom_boxplot(width=0.1,outlier.shape=NA)+
  stat_summary(fun.y=mean, geom="point", size=1, color='black')+
  #geom_line(MRE.bio.graph2, aes(x=as.factor(.id), y=mean(bias)),lty=1, lwd=0.5) + 
  #geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  #scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment")+
  #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Relative % Difference")+
  xlab(NULL)+
  ylab(NULL)+
  #facet_grid(Reg~.,labeller=as_labeller(region.names))+
  facet_wrap(Reg~.,labeller=as_labeller(region.names),strip.position="top",ncol=1)+
  #scale_x_discrete(limits=rev(levels(as.factor(MRE.bio.graph2$.id))))+
  ylim(-100,100)+
  #coord_flip()+
  my_theme.hor


###########################
tiff("Bias Plots Terminal.tif",width=190,height=240, unit='mm',res=300)
par(mar=c(6,6,6,6))


grid.arrange(ncol = 3,top="Terminal Year Bias",
             bio.plot.term,rec.plot.term,f.plot.term)
dev.off()
########################################

