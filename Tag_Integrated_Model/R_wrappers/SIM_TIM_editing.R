###################################################
# Code to run many simulations of TIM OM/EM
# Created by: Katelyn Bosley
###################################################


#remove junk from workspace
rm(list=(ls()))


#load required libraries
load_libraries<-function() {
  library(PBSmodelling)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(gplots)
  library(colorspace)
  library(RColorBrewer)
  library(dplyr)
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
}
load_libraries()




#######################################################################
#Setting up the simulations
##################################


diag.run<-1
# ==0 NO Do NOT run the diagnostics plots for a single run
# ==1 YES run the diagnostics plots for a single run


# Set number of simulations to perform
nsim <-24  

######################
#plot parameters
#####################

# Select color for violins
vio.col<-"lightskyblue3"

#select color for median points
median.col<-"grey95"


###########################################################
### setting up the directories
###############################

# To run this simulation a folder for each scenario will need to be placed in the master directory. Each folder will need separate folders named 'Operating_Model' and 'EM_model' with the exe files and .DAT for the OM only. The Code will do the rest. Be sure that the TIM_diagnostics.R code is also in the master directory. 
#

# set master file with holding the runs 
direct_master<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SIM_TEST"
setwd(direct_master)

#list files in the directory
files<-list.files(direct_master) # these folders in the master will be the individual scenarios 

#select the file with the scenario you want to run
#if only running 1 folder set i to the number corresponding to the folder you want to run
folder.num=1

i=folder.num

#if running the several folders use the loop - This will come later
#for(i in 1:3){

##############################################################

#OM Location
OM_direct<-paste0(direct_master,"\\",files[i],"\\Operating_Model",sep="")
OM_name<-"TIM_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-paste0(direct_master,"\\",files[i],"\\Estimation_Model",sep="") #location of run(s)
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.

#Build diagnostics folder
dir.create(paste0(direct_master,"\\",files[i],"\\Diagnostics",sep=""))
diag_direct<-paste0(direct_master,"\\",files[i],"\\Diagnostics",sep="")


#
# SKIP TO PLOTTING SECTION IF SIMS ARE COMPLETED!
#

{ #run the hold dang thing...

###############################################################
# RUN TIM DIAGNOSTICS
###############################################################

# generates a single run diagnotic of the OM and EM and
# moves file to the diagnostics folder

if(diag.run==1){
source("TIM_Diagnostics_master.R") # make sure this code is in the direct_master folder.

invisible(file.rename(from=paste0(EM_direct,"\\Model_Diagnostics.pdf",sep=""),to=paste0(diag_direct,"\\Model_Diagnostics.pdf",sep="")))

invisible(file.rename(from=paste0(EM_direct,"\\Tag_Residuals.pdf",sep=""),to=paste0(diag_direct,"\\Tag_Residuals.pdf",sep="")))
}

############################################################
############## Run simulations #############################
############################################################

{ # run the simulations section

#Set up convergence record to holds likelihood components and eventually gradient 

Sim.Stats<-data.frame(SimN=seq(1:nsim), Converged=rep(NA,nsim), Tag_like=rep(NA,nsim), Catch_like=rep(NA,nsim), Survey_like=rep(NA,nsim), Fish_age_like=rep(NA,nsim), survey_age_like = rep(NA,nsim), Rec_like=rep(NA,nsim))

#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#set up text progress bar
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


#j=1

ls=foreach(j=1:nsim,.options.snow = opts,.combine='rbind',.packages =c('PBSmodelling')) %dopar% {
  
 temp<-as.vector(rep(NA,(ncol(Sim.Stats)-1)))
  
  dir.create(paste0(diag_direct,"\\Run",j,sep=""))
  dir.create(paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep=""))
  dir.create(paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep=""))
  
#Move folders and .dat files  
  invisible(file.copy(from=paste0(OM_direct,"\\TIM_OM.exe",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(OM_direct,"\\TIM_OM.dat",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(EM_direct,"\\TIM_EM.exe",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep="")))
  

  
#change seed for OM
setwd(paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep=""))
      
SIM.DAT=readLines("TIM_OM.dat",n=-1)
SIM.DAT[(grep("myseed_yield",SIM.DAT)+1)]=j
SIM.DAT[(grep("myseed_survey",SIM.DAT)+1)]=100000+j
SIM.DAT[(grep("myseed_F",SIM.DAT)+1)]=200000+nsim
SIM.DAT[(grep("myseed_rec_devs",SIM.DAT)+1)]=300000+j
SIM.DAT[(grep("myseed_rec_apport",SIM.DAT)+1)]=400000+j
SIM.DAT[(grep("myseed_rec_index",SIM.DAT)+1)]=500000+j
SIM.DAT[(grep("myseed_survey_age",SIM.DAT)+1)]=600000+j
SIM.DAT[(grep("myseed_catch_age",SIM.DAT)+1)]=700000+j
SIM.DAT[(grep("myseed_tag",SIM.DAT)+1)]=800000+j

writeLines(SIM.DAT,"TIM_OM.dat")

#run the OM
invisible(shell(paste0(OM_name," -nox -nohess"),wait=T))

# move the .rep from the OM_direct and change name. 
from<-paste0(paste0(diag_direct,"\\Run",j,"\\Operating_Model\\",OM_name,".rep",sep=""))
to<-paste0(paste0(diag_direct,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""))
file.copy(from = from,  to = to)

#Delete OM folder to save space
unlink(paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep=""),recursive = T)


################################
#run the EM

setwd(paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep=""))
invisible(shell(paste0(EM_name," -nox -ind"),wait=T))


if(file.exists("TIM_EM.cor")==TRUE)
   
{
  temp[1]=1 #set convergence
  
#Things to save values from .REP  
  out<-readList("TIM_EM.rep")
  temp[2:length(temp)]=c(out$tag_like,out$catch_like,out$survey_like,out$fish_age_like,out$survey_age_like,out$rec_like)
  
 #grad<-readLines("gradient.dat") #need to work on this
  
}


if(file.exists("TIM_EM.cor")==FALSE)
{
  
  temp[1]=0 #set convergence
  
  #Things to save values from .REP  
  out<-readList("TIM_EM.rep")
  temp[2:length(temp)]=c(out$tag_like,out$catch_like,out$survey_like,out$fish_age_like,out$survey_age_like,out$rec_like)

}
  
  #then delete the junk
  invisible(shell("del TIM_EM.exe")) 
  invisible(shell("del admodel.cov"))
  invisible(shell("del admodel.dep"))
  invisible(shell("del admodel.hes")) 
  invisible(shell("del TIM_EM.b01"))
  invisible(shell("del TIM_EM.P01"))
  invisible(shell("del TIM_EM.R01"))
  invisible(shell("del TIM_EM.b02"))
  invisible(shell("del TIM_EM.P02"))
  invisible(shell("del TIM_EM.R02"))
  invisible(shell("del TIM_EM.b03"))
  invisible(shell("del TIM_EM.P03"))
  invisible(shell("del TIM_EM.R03"))
  invisible(shell("del TIM_EM.bar"))
  invisible(shell("del TIM_EM.eva"))
  invisible(shell("del TIM_EM.log"))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
 
return(temp)

}

close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 


#Save Sim Statistics from runs  - We can add more stuff here eventually
Sim.Stats[,2:ncol(Sim.Stats)]<-ls[,1:ncol(ls)]
Sim.Stats$Total_like=rowSums(Sim.Stats)
 
#ending
unlink(getwd())
setwd(diag_direct)

#save convergence file
write.csv(file="Sim_Stats.csv",Sim.Stats)

} # end of simulations and value save

##############################################################
#Grabbing Values from EM runs
##############################################################

{ # run the value grab

# create a data frame for the run to hold true and estimated values from OM/EM

#reload directories if we are skipping the sims
setwd(direct_master)
files<-list.files(direct_master)

#i=folder.num # if you are working only one folder or if you are running this section separate from sims

OM_direct<-paste0(direct_master,"\\",files[i],"\\Operating_Model",sep="")
OM_name<-"TIM_OM" #name of the OM you are wanting to run
EM_direct<-paste0(direct_master,"\\",files[i],"\\Estimation_Model",sep="") 
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.
diag_direct<-paste0(direct_master,"\\",files[i],"\\Diagnostics",sep="")

#pull dimensions for building data frames for plotting

#get dimensions for scenario
out<-readList(paste0(EM_direct,"\\",EM_name,".rep")) #read in .rep file

#pull info about the model
na<-out$nages
nyrs<-out$nyrs
npops<-out$npops
nreg<-out$nregions
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)

#for metapop
if(npops>1){
  nreg=sum(nreg)}

########################################################################
# Building Run in parallel to save the values we want to keep
#######################################################################

#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)


#set up text progress bar
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

vg=foreach(k=1:nsim,.options.snow = opts,.combine='cbind',.packages = c('PBSmodelling','data.table','matrixStats')) %dopar% {

#read in rep
out<-readList(paste0(diag_direct,"\\Run",k,"\\Estimation_Model\\",EM_name,".rep",sep=""))
  
#Save R_ave
meanR_sim_temp<-out$R_ave_TRUE
meanR_est_temp<-out$R_ave
  
#Save_apport
apport_sim_temp<-data.table(Apport_Sim = c(t(out$Rec_Prop_TRUE)))
apport_est_temp<-data.table(Apport_Est = c(t(out$Rec_Prop)))  

#Save Rec
rec_sim_temp<-data.table(Rec_Sim = c(t(out$recruits_BM_TRUE)))
rec_est_temp<-data.table(Rec_Est= c(t(out$recruits_BM)))
  
#Save SSB
ssb_sim_temp<-data.table(SSB_Sim = c(t(out$SSB_region_TRUE)))
ssb_est_temp<-data.table(SSB_Est = c(t(out$SSB_region)))
  
#Save Bio
bio_sim_temp<-data.table(Bio_Sim = c(t(out$biomass_AM_TRUE)))
bio_est_temp<-data.table(Bio_Est = c(t(out$biomass_AM)))

#Save Survey
survey_sim_temp<-data.table(Surv_Sim = c(t(out$OBS_survey_fleet_bio)))
survey_est_temp<-data.table(Survey_Est = c(t(out$survey_fleet_bio)))

#Save Yield
yield_sim_temp<-data.table(Catch_Sim = c(t(out$OBS_yield_fleet)))
yield_est_temp<-data.table(Catch_Est = c(t(out$yield_fleet)))


#Save F
fmax_sim_temp<-data.table(Fmax_Sim=c(t(rowMaxs(out$F))))
fmax_est_temp<-data.table(Fmax_Est=c(t(rowMaxs(out$F_TRUE))))


#we can add other stuff here 

#return all the matrices as a list
  return(list(meanR_sim=meanR_sim_temp,
              meanR_est=meanR_est_temp,
              apport_sim=apport_sim_temp,
              apport_est=apport_est_temp,
              rec_sim=rec_sim_temp,
              rec_est=rec_est_temp,
              ssb_sim=ssb_sim_temp,
              ssb_est=ssb_est_temp,
              bio_sim=bio_sim_temp,
              bio_est=bio_est_temp,
              survey_sim=survey_sim_temp,
              survey_est=survey_est_temp,
              yield_sim=yield_sim_temp,
              yield_est=yield_est_temp,
              fmax_sim=fmax_sim_temp,
              fmax_est=fmax_est_temp))
              
} #end parallel
  

close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 

#saves all the results in a list for later plotting
setwd(diag_direct)
saveRDS(vg, file="Sim_run.RData")

}

################################################################################
# Load Sim results for plotting
#################################################################################
#reload directory information for plotting
 
{ #start here if you have already completed the SIMS and saved the data

####################################
# ONLY NEED THIS IF NOT RUNNING SIMS
####################################
#if we are just working with a sim that has already been done, we can start here...
#setwd(direct_master)
#files<-list.files(direct_master)
#i=1

#OM_direct<-paste0(direct_master,"\\",files[i],"\\Operating_Model",sep="")
#OM_name<-"TIM_OM" #name of the OM you are wanting to run
#EM_direct<-paste0(direct_master,"\\",files[i],"\\Estimation_Model",sep="") 
#EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.
#diag_direct<-paste0(direct_master,"\\",files[i],"\\Diagnostics",sep="")

#pull dimensions for building data frames for plotting

#get dimensions for scenario
#out<-readList(paste0(EM_direct,"\\",EM_name,".rep")) #read in .rep file

#pull info about the model
#na<-out$nages
#nyrs<-out$nyrs
#npops<-out$npops
#nreg<-out$nregions
#years<-seq(1:out$nyrs)
#ages<-seq(1:out$nages)

#for metapop
#if(npops>1){
#  nreg=sum(nreg)}
################################################################

  

###############################################################
# Start plotting
###############################################################
  
#load up the results
setwd(diag_direct)
Sim_Results<-readRDS('Sim_run.RData')
Sim.Stats<-read.csv('Sim_Stats.csv')

#Recruitment
rec_df_sim<-matrix(NA,nyrs*nreg,nsim)
rec_df_est<-matrix(NA,nyrs*nreg,nsim)

#SSB
ssb_df_sim<-matrix(NA,nyrs*nreg,nsim)
ssb_df_est<-matrix(NA,nyrs*nreg,nsim)

#Biomass
bio_df_sim<-matrix(NA,nyrs*nreg,nsim)
bio_df_est<-matrix(NA,nyrs*nreg,nsim)

#yield
catch_df_sim<-matrix(NA,nyrs*nreg,nsim)
catch_df_est<-matrix(NA,nyrs*nreg,nsim)

#survey
survey_df_sim<-matrix(NA,nyrs*nreg,nsim)
survey_df_est<-matrix(NA,nyrs*nreg,nsim)

#fmax
fmax_df_sim<-matrix(NA,nyrs*nreg,nsim)
fmax_df_est<-matrix(NA,nyrs*nreg,nsim)



# F and T to come later

# populate the matrices for plotting
for(i in 1:nsim){
rec_df_sim[,i]<-unlist(Sim_Results[5,i])
rec_df_est[,i]<-unlist(Sim_Results[6,i])
ssb_df_sim[,i]<-unlist(Sim_Results[7,i])
ssb_df_est[,i]<-unlist(Sim_Results[8,i])
bio_df_sim[,i]<-unlist(Sim_Results[9,i])
bio_df_est[,i]<-unlist(Sim_Results[10,i])
catch_df_sim[,i]<-unlist(Sim_Results[11,i])
catch_df_est[,i]<-unlist(Sim_Results[12,i])
survey_df_sim[,i]<-unlist(Sim_Results[13,i])
survey_df_est[,i]<-unlist(Sim_Results[14,i])
fmax_df_sim[,i]<-unlist(Sim_Results[15,i])
fmax_df_est[,i]<-unlist(Sim_Results[16,i])
# add in survey and catch... KB

}


#######################################
# Building a ggplot theme
######################################
my_theme<-
  (theme_bw()+
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(colour="grey20",size=8,angle=0,hjust=.5,vjust=0,face="plain")
  ))


#####################################################
# Likelihood plots
#####################################################

#plots panel of likelihoods

conv<-sum(Sim.Stats$Converged)/nsim

#starts with col 3

bins<-nsim/2

Like.hist<-function(j=5){
breaks<-(max(Sim.Stats[j])-min(Sim.Stats[j]))/bins

ggplot(Sim.Stats, aes(x=Sim.Stats[j])) + 
  geom_histogram(binwidth = breaks, col="black", fill="grey80")+
  ylab("Frequency")+
  xlab("Likelihood")+
  ggtitle(colnames(Sim.Stats[j]))+
  my_theme

}


#######################################
# Extracting only the "good" converged runs 
#######################################

#soon




#####################################
# Recruitment plot
#####################################

#build data.frame
rec.data<-data.frame(Dat=c(rep("SIM",nrow(rec_df_sim)),rep("EST",nrow(rec_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

rec.data<-cbind(rec.data,rbind(rec_df_sim,rec_df_est))
rec.long<-melt(rec.data, id=c("Dat","Years","Reg"))
rec.long$Reg<-as.character(rec.data$Reg)

#calculate the sum across areas 
total.rec<-data.frame(rec.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.rec$Reg<-rep("System",nrow(total.rec))
rec.long<-rbind(total.rec,rec.long)


#separate again for plotting
rec.est<-rec.long[rec.long$Dat=="EST",]
rec.sim<-rec.long[rec.long$Dat=="SIM",]


#calculate the percent bias
rec.est$val.true<-rec.sim$value
rec.est$bias=((rec.est$val.true-rec.est$value)/rec.est$val.true)*100

#calc medians table
rec.est.med <- rec.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
rec.meds <- rec.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate Rec Plot
rec.plot.gg<-ggplot(rec.meds, aes(x=as.factor(Years), y=log(value))) +
  geom_hline(aes(yintercept = log(med), group = Reg), colour = 'darkred', size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T)+
  geom_line(data = rec.est.med, aes(x=Years,y=log(med.sim)),lty=1) + 
  geom_point(data = rec.est.med, aes(x=Years,y=log(med.sim)), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Total Recruitment")+
  ylab("log(Recruitment)")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(-5,5)+
  my_theme


#generate Rec bias plot
rec.bias.gg<-ggplot(rec.meds, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'darkred',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T)+
  geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Recruitment Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  ylim(-500,500)+
  my_theme


#save rec bias calcs
write.csv(rec.est.med,"Rec_Bias.csv")

#####################################
# SSB plot


#build data.frame
ssb.data<-data.frame(Dat=c(rep("SIM",nrow(ssb_df_sim)),rep("EST",nrow(ssb_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

ssb.data<-cbind(ssb.data,rbind(ssb_df_sim,ssb_df_est))
ssb.long<-melt(ssb.data, id=c("Dat","Years","Reg"))
ssb.long$Reg<-as.character(ssb.data$Reg)

#calculate the sum across areas 
total.ssb<-data.frame(ssb.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.ssb$Reg<-rep("System",nrow(total.ssb))
ssb.long<-rbind(total.ssb,ssb.long)


#separate again for plotting
ssb.est<-ssb.long[ssb.long$Dat=="EST",]
ssb.sim<-ssb.long[ssb.long$Dat=="SIM",]


#calculate the percent bias
ssb.est$val.true<-ssb.sim$value
ssb.est$bias=((ssb.est$val.true-ssb.est$value)/ssb.est$val.true)*100

#calc medians table
ssb.est.med <- ssb.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
ssb.meds <- ssb.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate ssb Plot
ssb.plot.gg<-ggplot(ssb.meds, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=T)+
  geom_line(data = ssb.est.med, aes(x=Years,y=med.sim),lty=1) + 
  geom_point(data = ssb.est.med, aes(x=Years,y=med.sim), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("SSB")+
  ylab("SSB")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
ssb.bias.gg<-ggplot(ssb.meds, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'darkred',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T)+
  geom_point(data = ssb.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("SSB Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-500,500)+
  my_theme


#save rec bias calcs
write.csv(ssb.est.med,"SSB_Bias.csv")


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
bio.est$bias=((bio.est$val.true-bio.est$value)/bio.est$val.true)*100

#calc medians table
bio.est.med <- bio.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
bio.meds <- bio.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate bio Plot
bio.plot.gg<-ggplot(bio.meds, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=T)+
  geom_line(data = bio.est.med, aes(x=Years,y=med.sim),lty=1) + 
  geom_point(data = bio.est.med, aes(x=Years,y=med.sim), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass")+
  ylab("Biomass")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
bio.bias.gg<-ggplot(bio.meds, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'darkred',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T)+
  geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-500,500)+
  my_theme


#save rec bias calcs
write.csv(bio.est.med,"Biomass_Bias.csv")


##############################################################
###F Max

#build data.frame
fmax.data<-data.frame(Dat=c(rep("SIM",nrow(fmax_df_sim)),rep("EST",nrow(fmax_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

fmax.data<-cbind(fmax.data,rbind(fmax_df_sim,fmax_df_est))
fmax.long<-melt(fmax.data, id=c("Dat","Years","Reg"))
fmax.long$Reg<-as.character(fmax.data$Reg)

#calculate the sum across areas 
total.fmax<-data.frame(fmax.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.fmax$Reg<-rep("System",nrow(total.fmax))
fmax.long<-rbind(total.fmax,fmax.long)


#separate again for plotting
fmax.est<-fmax.long[fmax.long$Dat=="EST",]
fmax.sim<-fmax.long[fmax.long$Dat=="SIM",]


#calculate the percent bias
fmax.est$val.true<-fmax.sim$value
fmax.est$bias=((fmax.est$val.true-fmax.est$value)/fmax.est$val.true)*100

#calc medians table
fmax.est.med <- fmax.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
fmax.meds <- fmax.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate fmax Plot
fmax.plot.gg<-ggplot(fmax.meds, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=T)+
  geom_line(data = fmax.est.med, aes(x=Years,y=med.sim),lty=1) + 
  geom_point(data = fmax.est.med, aes(x=Years,y=med.sim), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Fully Selected F")+
  ylab("F")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
fmax.bias.gg<-ggplot(fmax.meds, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'darkred',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=T)+
  geom_point(data = fmax.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("F Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-500,500)+
  my_theme


#save rec bias calcs
write.csv(fmax.est.med,"Fmax_Bias.csv")




#######################################
# Create a summary PDF
#######################################
{
pdf("Simulation_Summary.pdf",paper="letter",height = 10, width=7)
par(mar=c(6,6,6,6))


# Some important info about the run

if(out$npops_OM>1 & out$npops>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Metapopulation",sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Spatial Heterogeneity",sep ="\n")
}
if(out$npops_OM==1 && out$nregions_OM==1 && out$npops==1 && out$nregions==1){
text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Panmictic", sep = "\n")
}

text2<-paste0("\nTotal # of Sims: ",nsim)

text3<-paste0("Proporton of Sims Converged: ", sum(Sim.Stats$Converged)," of ", nsim, " (",conv,")")

time<-paste0("Report Date:  ", Sys.time())

text.all<-paste(text1,text2,text3,sep = "\n")


# Create a text grob
tgrob <- textGrob(text.all,just = "centre")



grid.arrange(ncol=1,nrow=4,
             bottom=time,
             center=textGrob("TIM Simulation Summary", gp=gpar(fontsize=18,font=3)),
             center=tgrob)



#add rec plots
grid.arrange(ncol = 1,
             top="Recruitment Estimation",
             rec.plot.gg, rec.bias.gg)


#add ssb plots
grid.arrange(ncol = 1,
             top="SSB",
             ssb.plot.gg, ssb.bias.gg)

#add bio plots
grid.arrange(ncol = 1,
             top="Biomass",
             bio.plot.gg, bio.bias.gg)

#add fmax plots
grid.arrange(ncol = 1,
             top="Fully Selected F",
             fmax.plot.gg, fmax.bias.gg)



##############################
# Plot Likelihood components
##############################

grid.arrange(ncol = 2,nrow=4,
             top="Likelihood Components",
             ggplotGrob(Like.hist(4)),
             ggplotGrob(Like.hist(5)),
             ggplotGrob(Like.hist(6)),
             ggplotGrob(Like.hist(7)),
             ggplotGrob(Like.hist(8)),             
             ggplotGrob(Like.hist(9)),
             ggplotGrob(Like.hist(10)))



dev.off()

} #end PDF code

} #end plotting/data section

} # the whole kit and caboodle
###############################################
#some notes for things to add

# CALCULATIONS
#Calculate summary stats
# - Parameter bias (% relative error for single parameters)
# - MARE (or whatever Andre has determined is best measure) for timeseries values (e.g., recruitment, F, biomass, #etc.) - MARE is good
# - Correlation (read in correlation matrix, have R code, and report all correlations above X% in table)
# - Convergence rates
# - Gradients
# - Likelihood components


#PLOTS
#Timeseries Trajectories compared to true (bio, ssb, recruitment, F, T.) (bean/violin plots overlayed on true point estimates)
#Fit to data
#Timeseries (yield, survey)
#Age comp (bubble plots?)
#Residuals (bubble plots? Timeseries)
#Parameter Bias (violin plots of %RE or MARE)
#Likelihood components (report out max gradient, convergence, etc.)
#Table of correlations



