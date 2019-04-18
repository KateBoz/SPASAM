rm(list=(ls()))



####### Where to store figures ############
Figure_dir<-"G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\Tagging manuscript\\manuscript\\Figures\\figures"

##########################################

###################### Directories that want to load (must match num_dir above) ###########################

####### BASE ############################
wd_1<<-'G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\Tagging manuscript\\Runs\\Exploratory Runs\\Self Consistency Runs\\Est w True Data'
nm_1<-"Base"
################################################


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
  library(plyr)
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
  
  #library(cowplot)
}
load_libraries()


nsim <-5
num_dir<-1

out<-readList("G:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\Tagging manuscript\\Runs\\Exploratory Runs\\Self Consistency Runs\\Est w True Data\\Diagnostics\\Results_converged\\Run1.rep") #read in .rep file

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
  move_df_sim<-matrix(NA,nyrs*nreg*nreg,nconv)
  move_df_est<-matrix(NA,nyrs*nreg*nreg,nconv)
  
  
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
    move_df_sim[,i]<-unlist(Sim_Results["movement_year_sim",i])
    move_df_est[,i]<-unlist(Sim_Results["movement_year_est",i])
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
  #rec.long<-rbind(rec.long,total.rec)
  
  
  #separate again for plotting
  rec.est<-rec.long[rec.long$Dat=="EST",]
  rec.sim<-rec.long[rec.long$Dat=="SIM",]
  
  
  #calculate the percent bias
  rec.est$val.true<-rec.sim$value
  rec.est$bias=((rec.est$value-rec.est$val.true)/rec.est$val.true)*100
  #rec.est$bias=(rec.est$val.true-rec.est$value)
  
  #calc medians table

  #generate Rec Plot
  rec.plot.gg<-ggplot(rec.est, aes(x=as.factor(Years), y=value)) +
    #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
    #geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
    geom_point( aes(x=Years,y=value), fill=median.col, shape=21,size=2.0) + 
    geom_line(aes(x=Years,y=val.true),lty=1,lwd=0.5) + 
    #geom_point(data = rec.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
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
   # geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
    geom_point(data = subset(rec.est, Reg %in% '1'), aes(x=Years,y=val.true), fill='black', colour='black',shape=19, size=2) + 
    geom_line(data =subset(rec.est, Reg %in% '1'), aes(x=Years,y=value),lty=1) + 
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    ylab('Recruits (1000s Fish)')+
    xlab(NULL)+
    ggtitle("Population 1")+
    
    #facet_grid(.~Reg,labeller=as_labeller(rec.name))+
    ylim(0,20000)+
    my_theme.fig4
  
  
  rec.plot.gg.fig5a<-ggplot(subset(rec.est, Reg %in% '2'), aes(x=as.factor(Years), y=value)) +
    #geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
    geom_point(data = subset(rec.est, Reg %in% '2'), aes(x=Years,y=val.true), fill='black', colour='black',shape=19, size=2) + 
    geom_line(data =subset(rec.est, Reg %in% '2'), aes(x=Years,y=value),lty=1) + 
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
    ylab('Recruits (1000s Fish)')+
    xlab(NULL)+
    ggtitle("Population 2")+
    
    #facet_grid(.~Reg,labeller=as_labeller(rec.name))+
    ylim(0,20000)+
    my_theme.fig4
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
 # bio.long<-rbind(total.bio,bio.long)
  
  
  #separate again for plotting
  bio.est<-bio.long[bio.long$Dat=="EST",]
  bio.sim<-bio.long[bio.long$Dat=="SIM",]
  
  
  #calculate the percent bias
  bio.est$val.true<-bio.sim$value
  bio.est$bias=((bio.est$value-bio.est$val.true)/bio.est$val.true)*100
 # bio.est$Years<-factor(years,levels=years)

  
  #generate bio Plot
  bio.plot.gg<-ggplot(bio.est, aes(x=as.factor(Years), y=value)) +
  #  geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
    geom_point(data = bio.est, aes(x=Years,y=value), fill=median.col, shape=21,size=2.0) + 
    geom_line(data = bio.est, aes(x=Years,y=val.true),lty=1, lwd=0.5) + 
   # geom_point(data = bio.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
    #scale_x_discrete(breaks=seq(0,nyrs,5))+
    ggtitle(get(paste0("nm_",z)))+
    labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Biomass")+
    xlab("Year")+
    facet_grid(Reg~.)+
    ylim(0,100000)+
    my_theme
  
  bio.name<-c('Syste,'="Metapopulation Biomass")
  
  bio.plot.gg.fig4<-ggplot(subset(bio.est, Reg %in% '1'), aes(x=Years, y=value)) +
    #geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
    geom_point(data = subset(bio.est, Reg %in% '1'), aes(x=Years,y=val.true), fill='black', colour='black',shape=19, size=2) + 
    geom_line(data =subset(bio.est, Reg %in% '1'), aes(x=Years,y=value),lty=1) + 
   scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
    ylab('Biomass (mt)')+
    xlab(NULL)+
    ggtitle("Population 1")+
    
    #facet_grid(.~Reg,labeller=as_labeller(f.name))+
    ylim(0,100000)+
    #xlim(1,30)+
    my_theme.fig4
  
  
  bio.plot.gg.fig5a<-ggplot(subset(bio.est, Reg %in% '2'), aes(x=(Years), y=value)) +
    #geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
    geom_point(data = subset(bio.est, Reg %in% '2'), aes(x=Years,y=val.true), fill='black', colour='black',shape=19, size=2) + 
    geom_line(data =subset(bio.est, Reg %in% '2'), aes(x=Years,y=value),lty=1) + 
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
    ylab('Biomass (mt)')+
    xlab(NULL)+
    ggtitle("Population 2")+
    
    #facet_grid(.~Reg,labeller=as_labeller(f.name))+
    ylim(0,100000)+
    my_theme.fig4
  
  
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
  
   #generate fmax Plot
  fmax.plot.gg<-ggplot(fmax.est, aes(x=as.factor(Years), y=value)) +
    #geom_violin(fill=vio.col,trim=T,bw="nrd0", alpha=0.6)+
  #  geom_point(data = fmax.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
    geom_point(data = fmax.est, aes(x=Years,y=value), fill="black", shape=16,size=1.0) + 
    geom_line(data = fmax.est, aes(x=Years,y=val.true),lty=1) + 
    scale_x_discrete(breaks=seq(0,nyrs,5))+
    ggtitle(get(paste0("nm_",z)))+
    labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("F")+
    xlab("Year")+
    facet_grid(Reg~.)+
    ylim(0,2.)+
    my_theme
  
  
  
  
  f.name<-c('1'="Population 1 Fishing Mortality")
  fmax.plot.gg.fig4<-ggplot(subset(fmax.est, Reg %in% '1'), aes(x=as.factor(Years), y=value)) +
    geom_line(data =subset(fmax.est, Reg %in% '1'), aes(x=Years,y=value,color="Estimated"),lty=1) + 
    geom_point(data = subset(fmax.est, Reg %in% '1'), aes(x=Years,y=val.true,color="True"),shape=19, size=2) + 

    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
    ylab('F')+
    xlab(NULL)+
    ggtitle("Population 1")+
    #facet_grid(.~Reg,labeller=as_labeller(f.name))+
    ylim(0,1.25)+
    my_theme.fig4+ 
    labs(colour="Model")+
    scale_colour_manual(values = c("black","black"))+
    theme(legend.position=c(.5,.9),legend.title=element_blank(),legend.margin=margin(c(4,4,4,4)),legend.text=element_text(size=8),
          legend.key.size = unit(.3,"cm"),legend.spacing.x = unit(.2, 'cm'))+
    guides (colour = guide_legend (ncol=2,reverse=TRUE,override.aes = list(linetype = c("blank","solid"), shape = c(19, NA))))
  
  
  fmax.plot.gg.fig5a<-ggplot(subset(fmax.est, Reg %in% '2'), aes(x=as.factor(Years), y=value)) +
    #geom_violin(fill='gray',trim=T,bw="nrd0", alpha=0.6,scale='width')+
    geom_point(data = subset(fmax.est, Reg %in% '2'), aes(x=Years,y=val.true), fill='black', colour='black',shape=19, size=2) + 
    geom_line(data =subset(fmax.est, Reg %in% '2'), aes(x=Years,y=value),lty=1) + 
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  
    ylab('F')+
    xlab(NULL)+
    ggtitle("Population 2")+
    
    #facet_grid(.~Reg,labeller=as_labeller(f.name))+
    ylim(0,1.25)+
    my_theme.fig4
 
  #Movement
  
  movement.data.sim<-data.frame(Dat=rep("SIM",nrow(move_df_sim)),Year=rep(1:nyrs,each=nreg),Reg_from=rep(1:nreg,each=nyrs*nreg),Reg_to=rep(1:nreg,nyrs*nreg))
  
  movement.data.est<-data.frame(Dat=rep("EST",nrow(move_df_est)),Year=rep(1:nyrs,each=nreg),Reg_from=rep(1:nreg,each=nyrs*nreg),Reg_to=rep(1:nreg,nyrs*nreg))                                
  
  
  
  movement.data<-rbind(movement.data.sim,movement.data.est)
  movement.data<-cbind(movement.data,rbind(move_df_sim,move_df_est))
  
  move.long<-melt(movement.data, id=c("Dat","Year","Reg_from","Reg_to"))
  move.long$Reg_from<-as.character(move.long$Reg_from)
  move.long$Reg_to<-as.character(move.long$Reg_to)
  
  
  
  #separate again for plotting
  move.est<-move.long[move.long$Dat=="EST",]
  move.sim<-move.long[move.long$Dat=="SIM",]
  
  
  #calculate the percent bias
  move.est$val.true<-move.sim$value
  move.est$bias=((move.est$value-move.est$val.true)/move.est$val.true)*100
  #move.est$bias=(move.est$val.true-move.est$value)
  

  
  move.plot.gg<-ggplot(move.est, aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
    geom_violin(aes(fill=Reg_to),trim=T,bw='nrd0',position = position_dodge(width=0.8), alpha=0.2)+
    geom_line(data = move.est, aes(x=Year,y=value, group=Reg_to))+
    geom_point(data = move.est, aes(x=Year,y=value, group=Reg_to),position = position_dodge(width=0.8), fill=median.col, shape=21,size=2.0) + 
    geom_line(data = move.est, aes(x=Year,y=val.true, group=Reg_to),lty=2)+
    geom_point(data = move.est, aes(x=Year,y=val.true, group=Reg_to, fill=Reg_to),position = position_dodge(width=0.8),shape=21,size=1.0) + 
    scale_color_grey()+ #"Move To",palette = "Set1")+
    scale_fill_grey()+ #scale_fill_brewer("Move To",palette = "Set1")+
    scale_x_discrete(breaks=seq(0,nyrs,5))+
    ggtitle(get(paste0("nm_",z)))+
    labs(subtitle=paste0("Convergence Rate = ", get(paste0("conv.rate_",z)),sep=""))+  ylab("Movement Rate")+
    xlab("Year")+
    facet_grid(Reg_from~.)+
    #ylim(-5,5)+
    my_theme
  
  move.name<-c('1'="Population 1 Movement")
  move.plot.gg.fig4<-ggplot(subset(move.est, Reg_from %in% '1'), aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
  #  geom_violin(aes(fill=Reg_to),trim=T,bw='nrd0',position = position_dodge(width=0.6), alpha=0.8,scale='width')+
    geom_point(data = subset(move.est, Reg_from %in% '1'), aes(x=Year,y=val.true, group=Reg_to,fill=Reg_to,shape=Reg_to)
              ,size=2 ) + 
    geom_line(data = subset(move.est, Reg_from %in% '1'), aes(x=Year,y=value, group=Reg_to,col=Reg_to,lty=Reg_to))+
    scale_color_manual(values=c("black", "grey"))+
    scale_fill_manual(values=c("black", "grey"))+
    ggtitle('Population 1')+
    ylab('Proportion Moving')+
    xlab(NULL)+
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #facet_grid(.~Reg_from,labeller=as_labeller(move.name))+
    my_theme.fig4
  
  move.plot.gg.fig5a<-ggplot(subset(move.est, Reg_from %in% '2'), aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
    #  geom_violin(aes(fill=Reg_to),trim=T,bw='nrd0',position = position_dodge(width=0.6), alpha=0.8,scale='width')+
    geom_point(data = subset(move.est, Reg_from %in% '2'), aes(x=Year,y=val.true, group=Reg_to,fill=Reg_to,shape=Reg_to)
    ,size=2) + 
    geom_line(data = subset(move.est, Reg_from %in% '2'), aes(x=Year,y=value, group=Reg_to,col=Reg_to,lty=Reg_to))+
    scale_color_manual(values=c("black", "grey"))+
    scale_fill_manual(values=c("black", "grey"))+
    ggtitle('Population 2')+
    ylab('Proportion Moving')+
    xlab(NULL)+
    scale_x_continuous(breaks=c(1,5,10,15,20,25,30),labels=c("1","5","10","15","20","25","30"))+
    #facet_grid(.~Reg_from,labeller=as_labeller(move.name))+
    my_theme.fig4
  
  assign(paste0("rec.plot.gg_",z),rec.plot.gg)
  assign(paste0("rec.plot.gg.fig4_",z),rec.plot.gg.fig4)
  assign(paste0("rec.plot.gg.fig5a_",z),rec.plot.gg.fig5a)
  assign(paste0("rec.est_",z),rec.est)

  assign(paste0("bio.plot.gg_",z),bio.plot.gg)
  assign(paste0("bio.plot.gg.fig5a_",z),bio.plot.gg.fig5a)
  assign(paste0("bio.plot.gg.fig4_",z),bio.plot.gg.fig4)
  assign(paste0("bio.est_",z),bio.est)

  assign(paste0("fmax.plot.gg_",z),fmax.plot.gg)
  assign(paste0("fmax.plot.gg.fig4_",z),fmax.plot.gg.fig4)
  assign(paste0("fmax.plot.gg.fig5a_",z),fmax.plot.gg.fig5a)
  assign(paste0("fmax.est_",z),fmax.est)

  assign(paste0("move.plot.gg_",z),move.plot.gg)
  assign(paste0("move.plot.gg.fig4_",z),move.plot.gg.fig4)
  assign(paste0("move.plot.gg.fig5a_",z),move.plot.gg.fig5a)
  assign(paste0("move.est_",z),move.est)


}



########################################################################################################
########################################################################################################
######## Figure 2 Timeseries for operating model################################################################################################
########################################################################################################
########################################################################################################

setwd(Figure_dir)


tiff("Fig 2_Timeseries of Operating Model Values.tif",width=240,height=190, unit='mm',res=500)
par(mar=c(6,6,6,6))

grid.arrange(ncol = 4,top=textGrob("Operating Model Timeseries Values",gp=gpar(fontsize=15,font=3)), 
             bottom=textGrob("Year",gp=gpar(fontsize=12,font=1)), 
             arrangeGrob(move.plot.gg.fig4_1, move.plot.gg.fig5a_1,top=textGrob("Movement",gp=gpar(fontsize=9,font=2),hjust=0.22),nrow=2), 
             arrangeGrob(fmax.plot.gg.fig4_1, fmax.plot.gg.fig5a_1,top=textGrob("Fishing Mortality",gp=gpar(fontsize=9,font=2),hjust=0.37),nrow=2), 
             arrangeGrob(bio.plot.gg.fig4_1,bio.plot.gg.fig5a_1,top=textGrob("Biomass",gp=gpar(fontsize=9,font=2),hjust=.1),nrow=2),
             arrangeGrob(rec.plot.gg.fig4_1,rec.plot.gg.fig5a_1,top=textGrob("Recruitment",gp=gpar(fontsize=9,font=2),hjust=.26),nrow=2)
)
dev.off()
