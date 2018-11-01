##################################################
###################################################
# SIMPLE SIMS CODE FOR SPASAM MODELS
# Created by: Katelyn Bosley
# Updated 10/2/2018
###################################################


##############################################################################
#
# This code will run a simulation experiment for the SPASAM models. 
# This "simple" version is one that only runs the sims and saves a 
# overall report of the runs (convergence rates/etc).
# It will also run the diagnostics code which is a one-off of the model
# to see 
# if there was major failure and provide clues of the cause. 
# Good Luck and Enjoy!
#
############################################################################


#remove junk from workspace at beginning of experiment
rm(list=(ls()))


#load required libraries - only need the ones for running the sims
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
  library(tidyr)
}
load_libraries()



###################################
#Setting up the simulations
##################################

#1) what are the model names
OM_name<-"MB_OM" #name of the OM you are wanting to run
EM_name<-"MB_EM" ###name of .dat, .tpl., .rep, etc.


#2) Do you want to run the diagnostics code before running the sims?

# if diagnostic code is run before the sims, the par file from the diagnostic run will be used create a pin for simulation runs.

diag.run<-1
   # ==0 NO Do NOT run the diagnostics plots for a single run before the sims
   # ==1 YES run the diagnostics plots for a single run

# if running the diagnostic run set residual switch for different values plotted

resid.switch=2
   # ==1 straight residual (TRUE-ESTIMATED; not % of true)
   # ==2 Relative percent difference ((TRUE/ESTIMATED)/TRUE *100)


#3) Do you want to run simulations?
run.sims<-1
# ==0 DO NOT run simulation experiment. This was completed already at an earlier time
# ==1 Run simulation experiement

#4) Set number of simulations to perform
nsim <-50

#5) Do you want to save important values from the runs?
save.values<-1
# ==0 DO NOT save values from the sims run
# ==1 SAVE values from the sims run

#6) Do you want to make plots with summary of the run?
make.plots<-1
# ==0 DO NOT create summary document
# ==1 Create summary document with graphs and run statistics

# If making plots, Select color for violins
vio.col<-"lightskyblue3"

# select color for median points
median.col<-"grey95"


#7) Do you want to plot non-converged runs? This will be needed/useful for MISMATCH
keep.all<-1
# ==0 DO NOT plot non converged runs
# ==1 DO plot non-converged runs


#########################################
### setting up the directories
#########################################

# To run this simulation a folder for each scenario will need to be placed in the master directory. Each scenario folder will need separate folders named 'Operating_Model' and 'Estimation_model' with the .exe files and .dat for the OM only. The code will do the rest. 

#If the diagnostic switch ==1 the TIM_diagnostics_SIM.R code will also have to be in the master directory.


#8) set master file with holding the runs 

direct_master<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\MB_test\\Management_Boundaries"
setwd(direct_master)

#list files in the directory
files<-list.files(direct_master) # these folders in the master will be the individual scenarios 

#select the file with the scenario you want to run
#if only running 1 folder set i to the number corresponding to the folder you want to run


folder.num=1

i=folder.num

##run.sims<-function(folder.num=i){
#run the whole code
{ #GO!

#if running the several folders use the loop - This will come later
#for(i in 1:3){

##############################################################
#OM Location
OM_direct<-paste0(direct_master,"\\",files[i],"\\Operating_Model",sep="")
#OM_name<-"MB_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-paste0(direct_master,"\\",files[i],"\\Estimation_Model",sep="") #location of run(s)
#EM_name<-"MB_EM" ###name of .dat, .tpl., .rep, etc.

#Build diagnostics/results folder
dir.create(paste0(direct_master,"\\",files[i],"\\Diagnostics",sep=""))
diag_direct<-paste0(direct_master,"\\",files[i],"\\Diagnostics",sep="")


#Build directories and folders to keep results from run

#Runs directory for individual runs
dir.create(paste0(diag_direct,"\\Runs",sep=""))
runs_dir<-paste0(diag_direct,"\\Runs",sep="")

#save good results
dir.create(paste0(diag_direct,"\\Results_converged",sep=""))
results_good<-paste0(diag_direct,"\\Results_converged",sep="")

#save bad results
dir.create(paste0(diag_direct,"\\Results_NOT_conv",sep=""))
results_bad<-paste0(diag_direct,"\\Results_NOT_conv",sep="")


# Running the simulations...Here we go!

#run simulations
{

###############################################################
# RUN TIM DIAGNOSTICS
###############################################################

# generates a single run diagnotic of the OM and EM and
# moves files from diagnostic run to the diagnostics folder
# this run generates the .pin file for the simulations

if(diag.run==1){
source("TIM_Diagnostics_master.R") # make sure this code is in the direct_master folder.

invisible(file.rename(from=paste0(EM_direct,"\\Model_Diagnostics.pdf",sep=""),to=paste0(diag_direct,"\\Model_Diagnostics.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Tag_Residuals.pdf",sep=""),to=paste0(diag_direct,"\\Tag_Residuals.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Correlation_Matrix.pdf",sep=""),to=paste0(diag_direct,"\\Correlation_Matrix.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Standard_Error_table.pdf",sep=""),to=paste0(diag_direct,"\\Standard_Error_table.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\output.txt",sep=""),to=paste0(diag_direct,"\\diag_output.txt",sep="")))
}


if(run.sims==1){
############################################################
############## RUN SIMULATIONS #############################
############################################################

# run the simulations section

#Set up convergence record to holds likelihood components and gradient 
Sim.Stats<-data.frame(SimN=seq(1:nsim), Converged=rep(NA,nsim),Max_Gradient=rep(NA,nsim),Obj_fun=rep(NA,nsim), Tag_like=rep(NA,nsim), Catch_like=rep(NA,nsim), Survey_like=rep(NA,nsim), Fish_age_like=rep(NA,nsim), survey_age_like = rep(NA,nsim), Rec_like=rep(NA,nsim))


#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#set up text progress bar
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


#j=1 # for debugging

ls=foreach(j=1:nsim,.options.snow = opts,.combine='rbind',.packages =c('PBSmodelling')) %dopar% {
  
 temp<-as.vector(rep(NA,(ncol(Sim.Stats)-1)))
  
  dir.create(paste0(runs_dir,"\\Run",j,sep=""))
  dir.create(paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep=""))
  dir.create(paste0(runs_dir,"\\Run",j,"\\Estimation_Model",sep=""))
  
#Move files to sims folders  
  invisible(file.copy(from=paste0(OM_direct,"\\",OM_name,".exe",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(OM_direct,"\\",OM_name,".dat",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(EM_direct,"\\", EM_name,".exe",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Estimation_Model",sep="")))
  invisible(file.copy(from=paste0(EM_direct,"\\", EM_name,".par",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".pin",sep="")))
  

#change seed for OM
setwd(paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep=""))
      
SIM.DAT=readLines(paste0(OM_name,".dat"),n=-1)
SIM.DAT[(grep("myseed_yield",SIM.DAT)+1)]=411+j
SIM.DAT[(grep("myseed_survey",SIM.DAT)+1)]=1110+j
SIM.DAT[(grep("myseed_rec_index",SIM.DAT)+1)]=5610+j
SIM.DAT[(grep("myseed_survey_age",SIM.DAT)+1)]=6831+j
SIM.DAT[(grep("myseed_catch_age",SIM.DAT)+1)]=7157+j
SIM.DAT[(grep("myseed_tag",SIM.DAT)+1)]=10009+j

writeLines(SIM.DAT,paste0(OM_name,".dat"))

#run the OM
invisible(shell(paste0(OM_name," -nox -nohess"),wait=T))

# move the .rep from the OM_direct and change name. 
from<-paste0(paste0(runs_dir,"\\Run",j,"\\Operating_Model\\",OM_name,".rep",sep=""))
to<-paste0(paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""))
file.copy(from = from,  to = to)

#set WD to estimation model folder
setwd(paste0(runs_dir,"\\Run",j,"\\Estimation_Model",sep=""))

#Delete OM folder to save space
#unlink(paste(runs_dir,"\\Run",j,"\\Operating_Model",sep=""),recursive = T)


################################
#run the EM

invisible(shell(paste0(EM_name," -nox -ind"),wait=T))


# if run converged...save results and move to results_good folder
if(file.exists(paste0(EM_name,".cor"))==T){
  
  #set convergence
  temp[1]=1 

 #save report  
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".rep",sep=""),  
            to = paste0(results_good,"\\Run",j,".rep",sep=""))
  #save std
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".std",sep=""),  
            to = paste0(results_good,"\\Run",j,".std",sep=""))
  #save cor
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".cor",sep=""),  
            to = paste0(results_good,"\\Run",j,".cor",sep=""))
  #save par
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".par",sep=""),  
            to = paste0(results_good,"\\Run",j,".par",sep=""))
  #save dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""),  
            to = paste0(results_good,"\\Run",j,".dat",sep=""))
  #save gradient file
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\gradient.dat",sep=""),  
            to = paste0(results_good,"\\gradient",j,".dat",sep=""))
  }
  
# if not converged..save results to results_bad
if(file.exists(paste0(EM_name,".cor"))==FALSE){
  #set convergence
  temp[1]=0 

  #save report  
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".rep",sep=""),  
            to = paste0(results_bad,"\\Run",j,".rep",sep=""))
  #save std
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".std",sep=""),  
            to = paste0(results_bad,"\\Run",j,".std",sep=""))
  #save cor
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".cor",sep=""),  
            to = paste0(results_bad,"\\Run",j,".cor",sep=""))
  #save par
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".par",sep=""),  
            to = paste0(results_bad,"\\Run",j,".par",sep=""))
  #save dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""),  
            to = paste0(results_bad,"\\Run",j,".dat",sep=""))
  
  #save gradient file
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\gradient.dat",sep=""),  
            to = paste0(results_bad,"\\gradient",j,".dat",sep=""))
  }
  
  
#save max gradient
  grad<-read.table("gradient.dat",header=T) #need to work on this
  temp[2]=max(grad$Gradient)
  
#Things to save values from .REP  
  out<-readList(paste0(EM_name,".rep"))
  temp[3:length(temp)]=c(out$f,out$tag_like,out$catch_like,out$survey_like,out$fish_age_like,out$survey_age_like,out$rec_like)

return(temp)
  
}


# end parallel job
close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 

#delete run folder
setwd(diag_direct)
unlink(paste(runs_dir),recursive = T)


#Save Sim Statistics from runs  - We can add more stuff here eventually
Sim.Stats[,2:ncol(Sim.Stats)]<-ls[,1:ncol(ls)]
Sim.Stats$Total_like=rowSums(Sim.Stats)
 
#save convergence file
write.csv(file="Sim_Stats.csv",Sim.Stats)

} 
  
}# end of simulations and summary save

# This is a simple value grab for the MB presentation at CAPAM
##################################################################

if(save.values==1){
  
  # create a data frame for the run to hold true and estimated values from OM/EM
  
  #setdirectory
  setwd(direct_master)
  
  Sim_Stats<-read.csv(paste0(diag_direct,'\\Sim_Stats.csv'))
  
  
  #copy non-converged run files into converged run folder if plotting all.
  if(keep.all==1){
  conv.runs<-seq(1:nsim)
    copy<-list.files(paste0(diag_direct,"\\Results_NOT_conv",sep=""))
    
    for(z in 1:length(copy)){
      file.copy(from=paste0(diag_direct,"\\Results_NOT_conv\\",copy[z],sep=""),to=paste0(diag_direct,"\\Results_converged\\",copy[z],sep=""))
    }
  }
  
  if(keep.all==0){
  conv.runs<-which(Sim_Stats$Converged==1)
  }
  
  #pull dimensions for building data frames for plotting
  out<-readList(paste0(EM_direct,"\\",EM_name,".rep")) #read in .rep file
  
  ########################################################################
  # Building value save in parallel
  #######################################################################
  
  #set up parallel 
  no_cores <- detectCores()
  cl<-makeCluster(no_cores)
  registerDoSNOW(cl)
  
  
  #set up text progress bar
  pb <- txtProgressBar(max = length(conv.runs), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  
  #k=1 #for debugging
  
  vg=foreach(k=conv.runs,.options.snow = opts,.combine='cbind',.packages = c('PBSmodelling','data.table','matrixStats')) %dopar% {
    
    
    #read in rep
    out<-readList(paste0(results_good,"\\Run",k,".rep",sep=""))
    
    #Save R_ave
    meanR_sim_temp<-data.table(meanR_sim = c(t(out$R_ave_TRUE)))
    meanR_est_temp<-data.table(meanR_est = c(t(out$R_ave)))
    
    #Save_apport
    apport_sim_temp<-data.table(Apport_Sim = c(t(out$Rec_Prop_TRUE)))
    apport_est_temp<-data.table(Apport_Est = c(t(out$Rec_Prop)))  
    
    #Save Rec Total
    rec_sim_temp<-data.table(Rec_Sim = c(t(out$recruits_BM_TRUE)))
    rec_est_temp<-data.table(Rec_Est= c(t(out$recruits_BM)))
    
    #Save Rec deviations
    rec_devs_sim_temp<-data.table(Rec_Devs_Sim = c(t(out$rec_devs_TRUE)))
    rec_devs_est_temp<-data.table(Rec_Devs_Est= c(t(out$rec_devs)))
    
    #save Init_abund
    init_abund_sim_temp<-data.table(Init_Abund_Sim = c(t(out$Init_Abund_TRUE)))
    init_abund_est_temp<-data.table(Init_Abund_Est= c(t(out$Init_Abund)))
    
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
    fmax_sim_temp<-data.table(Fmax_Sim=c(t(rowMaxs(out$F_TRUE))))
    fmax_est_temp<-data.table(Fmax_Est=c(t(rowMaxs(out$F))))
    
    #Save M
    #M_sim_temp<-data.table(c(t(out$M_TRUE)))
    #M_est_temp<-data.table(c(t(out$M)))
    
    #save reporting rate
    #rr_sim_temp<-data.table(c(t(out$report_rate_TRUE)))
    #rr_est_temp<-data.table(c(t(out$report_rate)))
    
    #Save q
    q_sim_temp<-data.table(q_Sim=c(t(out$q_survey_TRUE)))
    q_est_temp<-data.table(q_Est=c(t(out$q_survey)))
    
    #Save selectivity
    
    #fishery 
    sel_beta1_sim_temp<-data.table(sel_beta1_Sim=c(t(out$sel_beta1_TRUE)))
    sel_beta2_sim_temp<-data.table(sel_beta2_Sim=c(t(out$sel_beta2_TRUE)))
    sel_beta3_sim_temp<-data.table(sel_beta3_Sim=c(t(out$sel_beta3_TRUE)))
    sel_beta4_sim_temp<-data.table(sel_beta4_Sim=c(t(out$sel_beta4_TRUE)))
    
    sel_beta1_est_temp<-data.table(sel_beta1_Est=c(t(out$sel_beta1)))
    sel_beta2_est_temp<-data.table(sel_beta2_Est=c(t(out$sel_beta2)))
    sel_beta3_est_temp<-data.table(sel_beta3_Est=c(t(out$sel_beta3)))
    sel_beta4_est_temp<-data.table(sel_beta4_Est=c(t(out$sel_beta4)))
    
    #survey 
    sel_beta1_surv_sim_temp<-data.table(sel_beta1_Sim=c(t(out$sel_beta1_survey_TRUE)))
    sel_beta2_surv_sim_temp<-data.table(sel_beta2_Sim=c(t(out$sel_beta2_survey_TRUE)))
    sel_beta3_surv_sim_temp<-data.table(sel_beta3_Sim=c(t(out$sel_beta3_survey_TRUE)))
    sel_beta4_surv_sim_temp<-data.table(sel_beta4_Sim=c(t(out$sel_beta4_survey_TRUE)))
    
    sel_beta1_surv_est_temp<-data.table(sel_beta1_Est=c(t(out$sel_beta1_survey)))
    sel_beta2_surv_est_temp<-data.table(sel_beta2_Est=c(t(out$sel_beta2_survey)))
    sel_beta3_surv_est_temp<-data.table(sel_beta3_Est=c(t(out$sel_beta3_survey)))
    sel_beta4_surv_est_temp<-data.table(sel_beta4_Est=c(t(out$sel_beta4_survey)))
    
    #selectivity age for plotting
    select_age_sim_temp<-data.table(out$selectivity_age_TRUE)
    select_age_est_temp<-data.table(out$selectivity_age)
    
    select_age_survey_sim_temp<-data.table(out$survey_selectivity_age_TRUE)
    select_age_survey_est_temp<-data.table(out$survey_selectivity_age)
    
    
    #Save Movement
    ##pull.t<-unlist(out[grep("T_true", names(out), value = TRUE)])
    ##pull<-unlist(out[grep("T_est", names(out), value = TRUE)])
    
    ##movement_sim_temp<-data.table(movement_SIM=c(t(pull.t)))
    ##movement_est_temp<-data.table(movement_EST=c(t(pull)))
    
    ###################
    #Save age comps
    #################
    
    #survey
    pull.surv.sim<-out[grep("OBS_survey_prop", names(out), value = TRUE)]
    pull.surv.est<-out[grep("EST_survey_age_prop", names(out), value = TRUE)]
    
    surv_sim_temp<-data.table(do.call("rbind",pull.surv.sim))
    surv_est_temp<-data.table(do.call("rbind",pull.surv.est))
    
    #fishery
    pull.catch.sim<-out[grep("OBS_catch_prop", names(out), value = TRUE)]
    pull.catch.est<-out[grep("EST_catch_age_fleet_prop", names(out), value = TRUE)]
    
    catch_sim_temp<-data.table(do.call("rbind",pull.catch.sim))
    catch_est_temp<-data.table(do.call("rbind",pull.catch.est))
    
    
    
    #return all the matrices as a list
    return(list(meanR_sim=meanR_sim_temp,
                meanR_est=meanR_est_temp,
                apport_sim=apport_sim_temp,
                apport_est=apport_est_temp,
                rec_sim=rec_sim_temp,
                rec_est=rec_est_temp,
                rec_devs_sim=rec_devs_sim_temp,
                rec_devs_est=rec_devs_est_temp,
                init_abund_sim=init_abund_sim_temp,
                init_abund_est=init_abund_est_temp,
                ssb_sim=ssb_sim_temp,
                ssb_est=ssb_est_temp,
                bio_sim=bio_sim_temp,
                bio_est=bio_est_temp,
                survey_sim=survey_sim_temp,
                survey_est=survey_est_temp,
                yield_sim=yield_sim_temp,
                yield_est=yield_est_temp,
                fmax_sim=fmax_sim_temp,
                fmax_est=fmax_est_temp,
                ##M_sim=M_sim_temp,
                ##M_est=M_est_temp,
                ##rr_sim=rr_sim_temp,
                ##rr_est=rr_est_temp,
                sel_beta1_sim=sel_beta1_sim_temp,
                sel_beta2_sim=sel_beta2_sim_temp,
                sel_beta3_sim=sel_beta3_sim_temp,
                sel_beta4_sim=sel_beta4_sim_temp,
                sel_beta1_est=sel_beta1_est_temp,
                sel_beta2_est=sel_beta2_est_temp,
                sel_beta3_est=sel_beta3_est_temp,
                sel_beta4_est=sel_beta4_est_temp,
                sel_beta1_surv_sim=sel_beta1_surv_sim_temp,
                sel_beta2_surv_sim=sel_beta2_surv_sim_temp,
                sel_beta3_surv_sim=sel_beta3_surv_sim_temp,
                sel_beta4_surv_sim=sel_beta4_surv_sim_temp,
                sel_beta1_surv_est=sel_beta1_surv_est_temp,
                sel_beta2_surv_est=sel_beta2_surv_est_temp,
                sel_beta3_surv_est=sel_beta3_surv_est_temp,
                sel_beta4_surv_est=sel_beta4_surv_est_temp,
                select_age_sim=select_age_sim_temp,
                select_age_est=select_age_est_temp,
                select_age_survey_sim=select_age_survey_sim_temp,
                select_age_survey_est=select_age_survey_est_temp,
                q_sim=q_sim_temp,
                q_est=q_est_temp,
                #movement_sim=movement_sim_temp,
                #movement_est=movement_est_temp,
                surv_comp_sim=surv_sim_temp,
                surv_comp_est=surv_est_temp,
                catch_comp_sim=catch_sim_temp,
                catch_comp_est=catch_est_temp
                
    ))
    
  } #end parallel
  
  
  close(pb)
  stopCluster(cl) #end the cluster for parallel processing
  closeAllConnections() 
  
  #saves all the results in a list for later plotting
  setwd(diag_direct)
  saveRDS(vg, file="Sim_data.RData")
  
} #end save values
  

##################################################
### Build Data frames for plotting
###################################################
# 

if(make.plots==1){
  
  #pull dimensions for building data frames for plots
  out<-readList(paste0(EM_direct,"\\",EM_name,".rep")) #read in .rep file
  
  #pull info about the model
  na<-out$nages
  nyrs<-out$nyrs
  npops_OM<-out$npops_OM
  npops<-out$npops
  nreg_OM<-out$nregions_OM
  nreg<-out$nregions
  years<-seq(1:out$nyrs)
  ages<-seq(1:out$nages)
  #nrel<-out$nyrs_release
  
  
#for running the meta pop example. Might need fixing if more complex
if(npops_OM>1){
  nreg_OM=sum(nreg_OM)}
  
if(npops>1){
 nreg=sum(nreg)}
  

#load results form the value grab generating plots
  setwd(diag_direct)
  Sim_Results<-readRDS('Sim_data.RData')
  Sim_Stats<-read.csv('Sim_Stats.csv')
  

#get convergence information
nconv<-length(which(Sim_Stats$Converged==1))
pconv<-round(nconv/nsim,3)

#keep and plot all the runs converged and not converged
if(keep.all==1){
      nconv=nsim}
    
##############################################
# Unpacking Results
#############################################
    
    #saving values
    R_ave_df_sim<-matrix(NA,npops_OM,nconv)
    R_ave_df_est<-matrix(NA,npops,nconv)
    
    R_apport_df_sim<-matrix(NA,(nyrs)*nreg_OM,nconv)
    R_apport_df_est<-matrix(NA,(nyrs)*nreg,nconv)
    
    q_df_sim<-matrix(NA,nreg_OM,nconv)
    q_df_est<-matrix(NA,nreg,nconv)
    
    
    #selectivity params
    sel_beta1_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta2_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta3_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta4_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta1_df_est<-matrix(NA,nreg,nconv)
    sel_beta2_df_est<-matrix(NA,nreg,nconv)
    sel_beta3_df_est<-matrix(NA,nreg,nconv)
    sel_beta4_df_est<-matrix(NA,nreg,nconv)
    
    
    sel_beta1_surv_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta2_surv_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta3_surv_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta4_surv_df_sim<-matrix(NA,nreg_OM,nconv)
    sel_beta1_surv_df_est<-matrix(NA,nreg,nconv)
    sel_beta2_surv_df_est<-matrix(NA,nreg,nconv)
    sel_beta3_surv_df_est<-matrix(NA,nreg,nconv)
    sel_beta4_surv_df_est<-matrix(NA,nreg,nconv)
    
    #mortality
    #M_df_sim<-matrix(NA,na*npops_OM,nconv)
    #M_df_est<-matrix(NA,nyrs*npops*nreg*na,nconv)
    
    #reporting rate
    #rr_df_sim<-matrix(NA,nrel*nreg_OM,nconv)
    #rr_df_est<-matrix(NA,nrel*nreg,nconv)
    
    #Recruitment
    rec_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    rec_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #Recruitment deviations
    rec_devs_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    rec_devs_df_est<-matrix(NA,(nyrs-1)*nreg,nconv)
    
    #Init Abundance
    init_abund_df_sim<-matrix(NA,nreg_OM*nreg_OM*na,nconv)
    init_abund_df_est<-matrix(NA,nreg*nreg*na,nconv)
    
    #SSB
    ssb_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    ssb_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #Biomass
    bio_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    bio_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #yield
    catch_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    catch_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #survey bio
    survey_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    survey_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #fmax
    fmax_df_sim<-matrix(NA,nyrs*nreg_OM,nconv)
    fmax_df_est<-matrix(NA,nyrs*nreg,nconv)
    
    #select at age
    select_age_df_sim<-matrix(NA,na*nreg_OM,nconv)
    select_age_df_est<-matrix(NA,na*nreg,nconv)
    select_age_survey_df_sim<-matrix(NA,na*nreg_OM,nconv)
    select_age_survey_df_est<-matrix(NA,na*nreg,nconv)
    
    
    #movement
    #move_df_sim<-matrix(NA,nyrs*nreg_OM*nreg_OM*na,nconv)
    #move_df_est<-matrix(NA,nyrs*nreg*nreg*na,nconv)
    
###########################################################
# populate the matrices for plotting
    
#i=1 #for debugging
    
for(i in 1:nconv){
      R_ave_df_sim[,i]<-unlist(Sim_Results["meanR_sim",i])
      R_ave_df_est[,i]<-unlist(Sim_Results["meanR_est",i])
      #R_apport_df_sim[,i]<-unlist(Sim_Results["apport_sim",i])
      #R_apport_df_est[,i]<-unlist(Sim_Results["apport_est",i])
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
      #M_df_sim[,i]<-unlist(Sim_Results["M_sim",i])
      #M_df_est[,i]<-unlist(Sim_Results["M_est",i])
      #rr_df_sim[,i]<-unlist(Sim_Results["rr_sim",i])
      #rr_df_est[,i]<-unlist(Sim_Results["rr_est",i])
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
      
      #move_df_sim[,i]<-unlist(Sim_Results["movement_sim",i])
      #move_df_est[,i]<-unlist(Sim_Results["movement_est",i])
      
      select_age_df_sim[,i]<-unlist(Sim_Results["select_age_sim",i])
      select_age_df_est[,i]<-unlist(Sim_Results["select_age_est",i])
      select_age_survey_df_sim[,i]<-unlist(Sim_Results["select_age_survey_sim",i])
      select_age_survey_df_est[,i]<-unlist(Sim_Results["select_age_survey_est",i])
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
#Plot histograms of the likelihoods
    
#calculate the percent convervged
conv<-round(sum(Sim_Stats$Converged)/nsim,3)
    
#how many bins for the number of simulations
#bins<-round(nsim*0.75,0)
    
    Like.hist<-function(j=5){
      #breaks<-(max(Sim_Stats[j])-min(Sim_Stats[j]))/2
      
      ggplot(Sim_Stats, aes(x=Sim_Stats[j])) + 
        geom_histogram(bins = 15, col="black", fill="grey80")+
        ylab("Frequency")+
        xlab("Likelihood")+
        ggtitle(colnames(Sim_Stats[j]))+
        my_theme
      
    }
    
################################################
#plot from runs...ONLY the important ones for now
################################################
    
#Q plots
    
#Q_ave
    
#for matching population types
if(nreg_OM==nreg){
      q_est<-data.frame(melt(t(q_df_est)))
      q_est<-cbind(q_est,data.frame(melt(t(q_df_sim))[3]))
    }
    
# Spatial to panmictic
    if(nreg_OM>nreg){
      q_est<-data.frame(melt(t(q_df_est)))
      q_est<-cbind(q_est,data.frame(melt(t(colMeans(q_df_sim)))[3]))
    }
    
    names(q_est)<-c("Nsim","Reg","q_est","q_sim")
    q_est$q_bias<-((q_est$q_sim-q_est$q_est)/q_est$q_sim)*100
    
#calc medians
    
    #calculate the sum across areas 
    q.long<-melt(q_est, id=c("Reg","Nsim"))
    q.long$Reg<-as.character(q.long$Reg)
    
    median_q<-data.frame(q.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))
    
    
    
#plot
q.plot.gg<-ggplot(q_est, aes(x=as.factor(Reg), y=q_est, group=Reg)) +
      #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F, alpha=0.6,bw="SJ")+
      geom_point(data=subset(median_q,variable=="q_est"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
      geom_point(data=subset(median_q,variable=="q_sim"), aes(x=Reg,y=med), col="black", shape=16,cex=1.0) + 
      ggtitle("Survey Catchability")+
      ylab("Survey Catchability")+
      xlab("Area")+
      my_theme
    

q.bias.gg<-ggplot(q_est, aes(x=as.factor(Reg), y=q_bias, group=Reg)) +
      geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=T,bw="SJ",alpha=0.6)+
      geom_point(data=subset(median_q,variable=="q_bias"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
      #geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
      scale_x_discrete(breaks=seq(0,nyrs,1))+
      ggtitle("Survey Catchability Bias")+
      ylab("Relative % Difference")+
      xlab("Year")+
      #ylim(quantile(q_est$q_bias, c(.05,.95)))+
      my_theme
    
write.csv(median_q,"q_medians.csv")
    
    
    
    
###############################
#Recruitment Params
###########################
#R_ave
    
if(npops==npops_OM){
      Rave_est<-data.frame(melt(t(R_ave_df_est)))
      Rave_est<-cbind(Rave_est,data.frame(melt(t(R_ave_df_sim))[3]))
      names(Rave_est)<-c("Nsim","Reg","R_ave_est","R_ave_sim")
      Rave_est$R_ave_bias<-((Rave_est$R_ave_sim-Rave_est$R_ave_est)/Rave_est$R_ave_sim)*100
      
#calc medians
#calculate the sum across areas 
      rave.long<-melt(Rave_est, id=c("Reg","Nsim"))
      rave.long$Reg<-as.character(rave.long$Reg)
      
      median_R_ave<-data.frame(rave.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))
      
    }
    
    if(npops>npops_OM){
      Rave_est<-data.frame(melt(t(colSums(R_ave_df_est))))
      Rave_est<-cbind(Rave_est,data.frame(melt(t(R_ave_df_sim))[3]))
      names(Rave_est)<-c("Reg","Nsim","R_ave_est","R_ave_sim")
      Rave_est$R_ave_bias<-((Rave_est$R_ave_sim-Rave_est$R_ave_est)/Rave_est$R_ave_sim)*100
      
      #calc medians
      #calculate the sum across areas 
      rave.long<-melt(Rave_est, id=c("Reg","Nsim"))
      rave.long$Reg<-as.character(rave.long$Reg)
      
      median_R_ave<-data.frame(rave.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))
      
    }
    
    
    if(npops<npops_OM){
      Rave_est<-data.frame(melt(t(R_ave_df_est)))
      Rave_est<-cbind(Rave_est,data.frame(melt(t(colSums(R_ave_df_sim))[3])))
      names(Rave_est)<-c("Nsim","Reg","R_ave_est","R_ave_sim")
      Rave_est$R_ave_bias<-((Rave_est$R_ave_sim-Rave_est$R_ave_est)/Rave_est$R_ave_sim)*100
      
      #calc medians
      #calculate the sum across areas 
      rave.long<-melt(Rave_est, id=c("Reg","Nsim"))
      rave.long$Reg<-as.character(rave.long$Reg)
      
      median_R_ave<-data.frame(rave.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))
      
    }
    
    #plot
    rave.plot.gg<-ggplot(Rave_est, aes(x=as.factor(Reg), y=R_ave_est, group=Reg)) +
      #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
      geom_point(data=subset(median_R_ave,variable=="R_ave_est"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
      geom_point(data=subset(median_R_ave,variable=="R_ave_sim"), aes(x=Reg,y=med), col="black", shape=16,cex=1.0) + 
      ggtitle("Mean Recruitment")+
      ylab("Mean Recruitment")+
      xlab("Area")+
      my_theme
    
    
    rave.bias.gg<-ggplot(Rave_est, aes(x=as.factor(Reg), y=R_ave_bias, group=Reg)) +
      geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
      geom_point(data=subset(median_R_ave,variable=="R_ave_bias"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
      #geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
      scale_x_discrete(breaks=seq(0,nyrs,1))+
      ggtitle("Mean Recruitment Bias")+
      ylab("Relative % Difference")+
      xlab("Area")+
      #ylim(quantile(Rave_est$R_ave_bias, c(.05,.95)))+
      my_theme
    
    write.csv(median_R_ave,"Rave_medians.csv")
    
########################
#Rec timeseries

#build data.frame
    
    #for matching population types
    if((npops_OM==npops&&nreg_OM==nreg)||(npops_OM>npops&&nreg>1)){
      rec.data<-data.frame(Dat=c(rep("SIM",nrow(rec_df_sim)),rep("EST",nrow(rec_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))
      rec.data<-cbind(rec.data,rbind(rec_df_sim,rec_df_est))
      rec.long<-melt(rec.data, id=c("Dat","Years","Reg"))
      rec.long$Reg<-as.factor(as.character(rec.data$Reg))
      
      #calculate the sum across areas 
      total.rec<-data.frame(rec.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))
      total.rec$Reg<-rep("System",nrow(total.rec))
      rec.long<-rbind(rec.long,total.rec)
    }
    
    #For metamictic/Metapop to panmictic
    if((npops_OM==npops&&nreg_OM>nreg)||(npops_OM>npops&&npops==1&&nreg==1)){
      
      #aggregating simulation recruitment
      rec.data<-data.frame(Dat=rep("SIM",nrow(rec_df_sim)),Years=rep(years,nreg_OM),Reg=rep(1:nreg_OM,each=nyrs))
      rec.data<-cbind(rec.data,rbind(rec_df_sim))
      rec.long.sim<-melt(rec.data, id=c("Dat","Years","Reg"))
      total.rec.sim<-data.frame(rec.long.sim %>% group_by(Dat,variable, Years) %>% summarise(value=sum(value)))
      
      #adding Estimated recruitment  
      rec.data.est<-data.frame(Dat=rep("EST",nrow(rec_df_est)),Years=rep(years,nreg))
      rec.data.est<-cbind(rec.data.est,rbind(rec_df_est))
      rec.long.est<-melt(rec.data.est, id=c("Dat","Years"))
      
      rec.long<-rbind(rec.long.est,total.rec.sim)
      rec.long$Reg<-"System"
      rec.long$Reg<-as.factor(as.character(rec.long$Reg))
    }
    
    
#separate again for plotting
rec.est<-rec.long[rec.long$Dat=="EST",]
rec.sim<-rec.long[rec.long$Dat=="SIM",]
    
    
#calculate the percent bias
rec.est$val.true<-rec.sim$value
rec.est$bias=((rec.est$val.true-rec.est$value)/rec.est$val.true)*100
#rec.est$bias=(rec.est$val.true-rec.est$value)
    
#calc medians table
rec.est.med <- rec.est %>% group_by(Reg,Years) %>%
      summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))
    
#calc vals for plots
rec.meds <- rec.est %>% group_by(Reg) %>%
      mutate(med.est = median(value),med.true=median(val.true))
    
    
#generate Rec Plot
    rec.plot.gg<-ggplot(rec.est, aes(x=as.factor(Years), y=value)) +
      #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
      geom_point(data = rec.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
      geom_line(data = rec.est.med, aes(x=Years,y=med.sim),lty=1,lwd=0.5) + 
      geom_point(data = rec.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("Total Recruitment")+
      ylab("Recruitment")+
      xlab("Year")+
      facet_grid(Reg~.)+
      my_theme
    
    
    #generate Rec bias plot
    rec.bias.gg<-ggplot(rec.est, aes(x=as.factor(Years), y=bias)) +
      geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
      geom_line(data = rec.est.med, aes(x=Years,y=med.bias),lty=1,lwd=0.5) + 
      geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) +
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("Recruitment Bias")+
      ylab("Relative % Difference")+
      xlab("Year")+
      facet_grid(Reg~.,scales="free")+
      #ylim(quantile(rec.est$bias, c(.05,.95)))+
      my_theme
    
    
    #save rec bias calcs
    write.csv(rec.est.med,"Rec_Bias.csv")
    
#####################################################
#SSB plot
    
#build data.frame
    
#for matching population types
    if((npops_OM==npops&&nreg_OM==nreg)||(npops_OM>npops&&nreg>1)){
      ssb.data<-data.frame(Dat=c(rep("SIM",nrow(ssb_df_sim)),rep("EST",nrow(ssb_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))
      
      ssb.data<-cbind(ssb.data,rbind(ssb_df_sim,ssb_df_est))
      ssb.long<-melt(ssb.data, id=c("Dat","Years","Reg"))
      ssb.long$Reg<-as.character(ssb.data$Reg)
      
      #calculate the sum across areas 
      total.ssb<-data.frame(ssb.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))
      
      total.ssb$Reg<-rep("System",nrow(total.ssb))
      ssb.long<-rbind(total.ssb,ssb.long)
    }
    
#Metamictic/Metapop to Panmictic
if((npops_OM==npops&&nreg_OM>nreg)||(npops_OM>npops&&npops==1&&nreg==1)){
ssb.data<-data.frame(Dat=rep("SIM",nrow(ssb_df_sim)),Years=rep(years,nreg_OM),Reg=rep(1:nreg_OM,each=nyrs))
ssb.data<-cbind(ssb.data,rbind(ssb_df_sim))
ssb.long.sim<-melt(ssb.data, id=c("Dat","Years","Reg"))
total.ssb.sim<-data.frame(ssb.long.sim %>% group_by(Dat,variable, Years) %>% summarise(value=sum(value)))
      
#adding Estimated ssb  
ssb.data.est<-data.frame(Dat=rep("EST",nrow(ssb_df_est)),Years=rep(years,nreg))
ssb.data.est<-cbind(ssb.data.est,rbind(ssb_df_est))
ssb.long.est<-melt(ssb.data.est, id=c("Dat","Years"))
      
ssb.long<-rbind(ssb.long.est,total.ssb.sim)
ssb.long$Reg<-"System"
ssb.long$Reg<-as.factor(as.character(ssb.long$Reg))
}
    
    
    
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
    ssb.plot.gg<-ggplot(ssb.est, aes(x=as.factor(Years), y=value)) +
      geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
      geom_point(data = ssb.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
      geom_line(data = ssb.est.med, aes(x=Years,y=med.sim),lty=1, lwd=0.5) + 
      geom_point(data = ssb.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("SSB")+
      ylab("SSB")+
      xlab("Year")+
      facet_grid(Reg~.)+
      #ylim(-5,5)+
      my_theme
    
    
    #generate Rec bias plot
    ssb.bias.gg<-ggplot(ssb.est, aes(x=as.factor(Years), y=bias)) +
      geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
      geom_line(data = ssb.est.med, aes(x=Years,y=med.bias),lty=1, lwd=0.5) + 
      geom_point(data = ssb.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("SSB Bias")+
      ylab("Relative % Difference")+
      xlab("Year")+
      facet_grid(Reg~.)+
      #ylim(quantile(ssb.est$bias, c(.05,.95)))+
      my_theme
    
    
    #save rec bias calcs
    write.csv(ssb.est.med,"SSB_Bias.csv")
    
    
    ####################################
    #F plots
    
    #matching population types
    if((npops_OM==npops&&nreg_OM==nreg)||(npops_OM>npops&&nreg>1)){
      #build data.frame
      fmax.data.sim<-data.frame(Dat=rep("SIM",nrow(fmax_df_sim)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))
      fmax.data.est<-data.frame(Dat=rep("EST",nrow(fmax_df_est)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))
      
      fmax.data<-rbind(fmax.data.sim,fmax.data.est)
      
      fmax.data<-cbind(fmax.data,rbind(fmax_df_sim,fmax_df_est))
      fmax.long<-melt(fmax.data, id=c("Dat","Years","Reg"))
      fmax.long$Reg<-as.character(fmax.data$Reg)
      
      #calculate the sum across areas 
      total.fmax<-data.frame(fmax.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))
      
      total.fmax$Reg<-rep("System",nrow(total.fmax))
      fmax.long<-rbind(total.fmax,fmax.long)
    }
    
    #Spatial to panmictic
    if((npops_OM==npops&&nreg_OM>nreg)||(npops_OM>npops&&npops==1&&nreg==1)){
      
      fmax.data.sim.t<-data.frame(Dat=rep("SIM",nrow(fmax_df_sim)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))
      fmax.data.est.t<-data.frame(Dat=rep("EST",nrow(fmax_df_est)),Years=rep(years,nreg),Reg=rep(1:nreg,each=nyrs))
      
      #aggregating the values for panmictic
      fmax.data.sim<-cbind(fmax.data.sim.t,rbind(fmax_df_sim))
      fmax.long.sim<-melt(fmax.data.sim, id=c("Dat","Years","Reg"))
      total.fmax.sim<-data.frame(fmax.long.sim %>% group_by(Dat,variable, Years) %>% summarise(value=mean(value)))
      total.fmax.sim$Reg<-"System"
      
      #setting up the estimatated values
      fmax.data.est<-cbind(fmax.data.est.t,rbind(fmax_df_est))
      fmax.long.est<-melt(fmax.data.est, id=c("Dat","Years","Reg"))
      total.fmax.est<-data.frame(fmax.long.est %>% group_by(Dat,variable, Years) %>% summarise(value=mean(value)))
      total.fmax.est$Reg<-"System"
      
      fmax.long<-rbind(total.fmax.est,total.fmax.sim)
    }
    
    #separate again for plotting
    fmax.est<-fmax.long[fmax.long$Dat=="EST",]
    fmax.sim<-fmax.long[fmax.long$Dat=="SIM",]
    
    
    #calculate the percent bias
    fmax.est$val.true<-fmax.sim$value
    fmax.est$bias=((fmax.est$val.true-fmax.est$value)/fmax.est$val.true)*100
    #fmax.est$bias=((fmax.est$val.true-fmax.est$value))
    
    #removing the top and bottom 1% for the plots...
#    fmax.est<-subset(fmax.est, bias>(quantile(fmax.est$bias, .01) && bias<(quantile(fmax.est$bias, .99))))
    
    #calc medians table
    fmax.est.med <- fmax.est %>% group_by(Reg,Years) %>%
      summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))
    
    #calc vals for plots
    fmax.meds <- fmax.est %>% group_by(Reg) %>%
      mutate(med = median(value),med.true=median(val.true))
    
    #incase there are crazy outlier values
    #fmax.est<-subset(fmax.est, bias>(quantile(fmax.est$bias, .01) && bias<(quantile(fmax.est$bias, .99))))
    
    #generate fmax Plot
    fmax.plot.gg<-ggplot(fmax.est, aes(x=as.factor(Years), y=value)) +
      geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
      geom_point(data = fmax.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
      geom_point(data = fmax.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
      geom_line(data = fmax.est.med, aes(x=Years,y=med.sim),lty=1) + 
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("Fully Selected F")+
      ylab("F")+
      xlab("Year")+
      facet_grid(Reg~.)+
      #ylim(-5,5)+
      #ylim(quantile(fmax.est$bias, c(.01,.99)))+
      my_theme
    
    
    #bias plot
    fmax.bias.gg<-ggplot(fmax.est, aes(x=as.factor(Years), y=bias)) +
      geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
      geom_violin(fill=vio.col,trim=F,bw="SJ")+
      geom_line(data = fmax.est.med, aes(x=Years,y=med.bias),lty=1, lwd=0.5) + 
      geom_point(data = fmax.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
      scale_x_discrete(breaks=seq(0,nyrs,5))+
      ggtitle("F Bias")+
      ylab("Relative % Difference")+
      xlab("Year")+
      facet_grid(Reg~.)+
      my_theme
    
    
    #save rec bias calcs
    write.csv(fmax.est.med,"Fmax_Bias.csv")
    
    
    ################################
    #Build summary document
    ###############################
    
    pdf("Simulation_Summary.pdf",paper="letter",height = 11, width=8)
    par(mar=c(6,6,6,6))
    
    
    #What are the population assumptions of OM and EM?
    if(out$npops_OM>1 & out$npops>1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Metapopulation",sep = "\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions>1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Spatial Heterogeneity",sep ="\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM==1 && out$npops==1 && out$nfleets==1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Panmictic", sep = "\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM==1 && out$nfleets>1 && out$nfleets>1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","OM Population Structure: Fleets-as-Areas", "EM Population Structure: Fleets-as-Areas", sep = "\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions==1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Panmictic", sep = "\n")
    }
    
    if(out$npops_OM>1 && out$npops==1 && out$nregions==1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Metapopulation", "EM Population Structure: Panmictic", sep = "\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions==1 && out$nfleets>1 ){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Fleets-as-Areas", sep = "\n")
    }
    
    if(out$npops_OM==1 && out$nregions_OM>1 && out$npops>1 && out$nregions==1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Metapopulation", sep = "\n")
    }
    
    if(out$npops_OM>1 && sum(out$nregions_OM)>1 && out$npops==1 && out$nregions>1){
      text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Metapopulation", "EM Population Structure: Spatial Heterogeneity", sep = "\n")
    }
    
    
    text2<-paste0("\nTotal # of Sims: ",nsim)
    
    text3<-paste0("Proportion of Sims Converged: ", sum(Sim_Stats$Converged)," of ", nsim, " (",pconv,")")
    
    time<-paste0("Report Date:  ", Sys.time())
    
    text.all<-paste(text1,text2,text3,sep = "\n")
    
    
    # Create a text grob
    tgrob <- textGrob(text.all,just = "centre")
    
    
    
    grid.arrange(ncol=1,nrow=4,
                 bottom=time,
                 center=textGrob("TIM Simulation Summary", gp=gpar(fontsize=18,font=3)),
                 center=tgrob)
    
    
    #Recruitment plots
    grid.arrange(ncol = 2,
                 top="Input Parameters",
                 q.plot.gg,
                 q.bias.gg,
                 rave.plot.gg,
                 rave.bias.gg)
    
    #grid.arrange(ncol = 1,
    #             top="Recruitment Apportionment",
    #             r.apport.plot.gg, r.apport.bias.gg)
    
    
    grid.arrange(ncol = 1,
                 top="Recruitment Estimation",
                 rec.plot.gg, rec.bias.gg)
    
    #Biomass plots
    grid.arrange(ncol = 1,
                 top="SSB",
                 ssb.plot.gg, ssb.bias.gg)
    
    #Add fmax plots
    grid.arrange(ncol = 1,
                 top="Fully Selected F",
                 fmax.plot.gg, fmax.bias.gg)
    
    
    ##############################
    # Plot Likelihood components
    
    
    grid.arrange(ncol = 2,nrow=4,
                 top="Likelihood Components",
                 ggplotGrob(Like.hist(4)),
                 ggplotGrob(Like.hist(5)),
                 ggplotGrob(Like.hist(6)),
                 ggplotGrob(Like.hist(7)),
                 ggplotGrob(Like.hist(8)),             
                 ggplotGrob(Like.hist(9)),
                 ggplotGrob(Like.hist(10)),
                 #ggplotGrob(Like.hist(11)), #rec penalty
                 ggplotGrob(Like.hist(12))
    )
    
    
    
    dev.off()
    
  } # end of the code
  

#Pulling out the SSB values and creating a csv
ssb.table<-spread(ssb.long, key = variable, value = value)
rec.table<-spread(rec.long, key = variable, value = value)
fmax.table<-spread(fmax.long, key = variable, value = value)

write.csv(ssb.table,"ssb_table.csv")
write.csv(rec.table,"rec_table.csv")
write.csv(fmax.table,"fmax_table.csv")


} #end of plots code 
  
#######################


