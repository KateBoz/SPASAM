###################################################
# Code to run many simulations of TIM OM/EM
# Created by: Katelyn Bosley
###################################################

# This code will run a simulation experiment for the SPASAM model
# 

###################################
#Setting up the simulations
##################################


#1) Do you want to run the diagnostics code before running the sims?

diag.run<-1
# ==0 NO Do NOT run the diagnostics plots for a single run
# ==1 YES run the diagnostics plots for a single run


#2) Set number of simulations to perform
nsim <-20

######################
#plot parameters
#####################

#3) Select color for violins
vio.col<-"lightskyblue3"

#4) select color for median points
median.col<-"grey95"


#########################################
### setting up the directories
#########################################

# To run this simulation a folder for each scenario will need to be placed in the master directory. Each folder will need separate folders named 'Operating_Model' and 'EM_model' with the exe files and .DAT for the OM only. The Code will do the rest. Be sure that the TIM_diagnostics.R code is also in the master directory if diag_switch==1 
#

#5) set master file with holding the runs 
direct_master<-wd  #"F:\\NOAA FILES\\Research\\SPASAM\\CAPAM Runs\\simple example"
#setwd(direct_master)

#list files in the directory
#files<-list.files(direct_master) # these folders in the master will be the individual scenarios 

#select the file with the scenario you want to run
#if only running 1 folder set i to the number corresponding to the folder you want to run
#folder.num=1

#i=folder.num



{ #run the whole thing

  

#if running the several folders use the loop - This will come later
#for(i in 1:3){

##############################################################
#OM Location
OM_direct<-paste0(direct_master,"\\Operating_Model",sep="")
OM_name<-"TIM_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-paste0(direct_master,"\\Estimation_Model",sep="") #location of run(s)
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.

#Build diagnostics/results folder
dir.create(paste0(direct_master,"\\Diagnostics",sep=""))
diag_direct<-paste0(direct_master,"\\Diagnostics",sep="")

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

#save the primary output files needed to correct file
#need .rep, .std

##################################################
# SKIP TO PLOTTING SECTION IF SIMS ARE COMPLETED!#

#Simulations...

#run simulations
{
###############################################################
# RUN TIM DIAGNOSTICS
###############################################################


# generates a single run diagnotic of the OM and EM and
# moves files from diagnostic run to the diagnostics folder

if(diag.run==1){
source("TIM_Diagnostics_master.R") # make sure this code is in the direct_master folder.

invisible(file.rename(from=paste0(EM_direct,"\\Model_Diagnostics.pdf",sep=""),to=paste0(diag_direct,"\\Model_Diagnostics.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Tag_Residuals.pdf",sep=""),to=paste0(diag_direct,"\\Tag_Residuals.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Correlation_Matrix.pdf",sep=""),to=paste0(diag_direct,"\\Correlation_Matrix.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Standard_Error_table.pdf",sep=""),to=paste0(diag_direct,"\\Standard_Error_table.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\output.txt",sep=""),to=paste0(diag_direct,"\\diag_output.txt",sep="")))
}


############################################################
############## Run simulations #############################
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


#j=1

ls=foreach(j=1:nsim,.options.snow = opts,.combine='rbind',.packages =c('PBSmodelling')) %dopar% {
  
 temp<-as.vector(rep(NA,(ncol(Sim.Stats)-1)))
  
  dir.create(paste0(runs_dir,"\\Run",j,sep=""))
  dir.create(paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep=""))
  dir.create(paste0(runs_dir,"\\Run",j,"\\Estimation_Model",sep=""))
  
#Move folders and .dat files  
  invisible(file.copy(from=paste0(OM_direct,"\\",OM_name,".exe",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(OM_direct,"\\",OM_name,".dat",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(EM_direct,"\\", EM_name,".exe",sep=""),to=paste0(runs_dir,"\\Run",j,"\\Estimation_Model",sep="")))
  

  
#change seed for OM
setwd(paste0(runs_dir,"\\Run",j,"\\Operating_Model",sep=""))
      
SIM.DAT=readLines(paste0(OM_name,".dat"),n=-1)
SIM.DAT[(grep("myseed_yield",SIM.DAT)+1)]=411+j
SIM.DAT[(grep("myseed_survey",SIM.DAT)+1)]=1110+j
SIM.DAT[(grep("myseed_F",SIM.DAT)+1)]=45651+j
SIM.DAT[(grep("myseed_T",SIM.DAT)+1)]=32+j
SIM.DAT[(grep("myseed_rec_devs",SIM.DAT)+1)]=23440+j
#SIM.DAT[(grep("myseed_rec_apport",SIM.DAT)+1)]=4199+j
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
if(file.exists("TIM_EM.cor")==TRUE)
   
{
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
  
  #save OM dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Operating_Model\\",OM_name,".dat",sep=""),  
            to = paste0(results_good,"\\Run",j,"_OM.dat",sep=""))
  
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Operating_Model\\",OM_name,".rep",sep=""),  
            to = paste0(results_good,"\\Run",j,"_OM.rep",sep=""))
  
  #save EM dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""),  
            to = paste0(results_good,"\\Run",j,".dat",sep=""))
  
  #save cov
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\gradient.dat",sep=""),  
            to = paste0(results_good,"\\Run",j,".gradient",sep=""))

  }
  
# if not converged..save results to results_bad
if(file.exists("TIM_EM.cor")==FALSE)
{
  #set convergence
  temp[1]=0 

  #save report  
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".rep",sep=""),  
            to = paste0(results_bad,"\\Run",j,".rep",sep=""))
  
  #save OM dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Operating_Model\\",OM_name,".dat",sep=""),  
            to = paste0(results_bad,"\\Run",j,"_OM.dat",sep=""))
  
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Operating_Model\\",OM_name,".rep",sep=""),  
            to = paste0(results_bad,"\\Run",j,"_OM.rep",sep=""))
  
  #save EM dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""),  
            to = paste0(results_bad,"\\Run",j,".dat",sep=""))
  
  #save par
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".par",sep=""),  
            to = paste0(results_bad,"\\Run",j,".par",sep=""))
  }
  
  
#save max gradient
  grad<-read.table("gradient.dat",header=T) #need to work on this
  temp[2]=max(grad$Gradient)
  
#Things to save values from .REP  
  out<-readList("TIM_EM.rep")
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


}  # end of simulations and value save


##############################################################
##############################################################
#Grabbing Values from converged report files
##############################################################
##############################################################

# run the value grab
{
  
# create a data frame for the run to hold true and estimated values from OM/EM

#reload directories if we are skipping the sims
setwd(direct_master)
files<-list.files(direct_master)

#i=folder.num # if you are working only one folder or if you are running this section separate from sims

OM_direct<-paste0(direct_master,"\\Operating_Model",sep="")
OM_name<-OM_name #name of the OM you are wanting to run

EM_direct<-paste0(direct_master,"\\Estimation_Model",sep="") 
EM_name<-EM_name ###name of .dat, .tpl., .rep, etc.

diag_direct<-paste0(direct_master,"\\Diagnostics",sep="")


Sim_Stats<-read.csv(paste0(diag_direct,'\\Sim_Stats.csv'))

#create a vector for converged runs to calculate bias and create plots
conv.runs<-which(Sim_Stats$Converged==1)



#pull dimensions for building data frames for plotting
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
pb <- txtProgressBar(max = length(conv.runs), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


#k=1

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

pull<-unlist(out[grep("T_est", names(out), value = TRUE)])
movement_est_temp<-data.table(movement_Est=c(t(pull)))

pull.t<-unlist(out[grep("T_true", names(out), value = TRUE)])
movement_sim_temp<-data.table(movement_Sim=c(t(pull.t)))

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
              movement_sim=movement_sim_temp,
              movement_est=movement_est_temp,
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

}

################################################################################
# Load Sim results for plotting
#################################################################################
#reload directory information for plotting
 
#start here if you have already completed the SIMS and saved the data

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
# Basic Plots
###############################################################
  
{ #load results
  
#load up the results
setwd(diag_direct)
Sim_Results<-readRDS('Sim_data.RData')
Sim_Stats<-read.csv('Sim_Stats.csv')

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


#age comps - eventually?


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


}


#########################################
############################################
# MAKE PLOTS
###########################################

{
  
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
bins<-round(nsim*0.8,0)

Like.hist<-function(j=5){
breaks<-round((max(Sim_Stats[j])-min(Sim_Stats[j]))/bins)

ggplot(Sim_Stats, aes(x=Sim_Stats[j])) + 
  geom_histogram(binwidth = breaks, col="black", fill="grey80")+
  ylab("Frequency")+
  xlab("Likelihood")+
  ggtitle(colnames(Sim_Stats[j]))+
  my_theme

}



######################################
# Building the plots!
######################################

#R_ave
q_est<-data.frame(melt(t(q_df_est)))
q_est<-cbind(q_est,data.frame(melt(t(q_df_sim))[3]))
names(q_est)<-c("Nsim","Reg","q_est","q_sim")
q_est$q_bias<-((q_est$q_est-q_est$q_sim)/q_est$q_sim)*100

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



#R_ave
Rave_est<-data.frame(melt(t(R_ave_df_est)))
Rave_est<-cbind(Rave_est,data.frame(melt(t(R_ave_df_sim))[3]))
names(Rave_est)<-c("Nsim","Reg","R_ave_est","R_ave_sim")
Rave_est$R_ave_bias<-((Rave_est$R_ave_est-Rave_est$R_ave_sim)/Rave_est$R_ave_sim)*100

#calc medians

#calculate the sum across areas 
rave.long<-melt(Rave_est, id=c("Reg","Nsim"))
rave.long$Reg<-as.character(rave.long$Reg)

median_R_ave<-data.frame(rave.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))



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

#selectivity coeffs here eventually

######################################################################
#R_apport
R_apport_est<-data.frame(Year=rep(2:nyrs,nreg),Reg=rep(1:nreg,each=nyrs-1))
R_apport_est<-cbind(R_apport_est,data.frame(R_apport_df_est))
R_apport_est_melt<-melt(R_apport_est,id=c("Year","Reg"))

R_apport_sim<-data.frame(Year=rep(2:nyrs,nreg),Reg=rep(1:nreg,each=nyrs-1))
R_apport_sim<-cbind(R_apport_sim,data.frame(R_apport_df_sim))
R_apport_sim_melt<-melt(R_apport_sim,id=c("Year","Reg"))

names(R_apport_est_melt)<-c("Year","Reg","Nsim","R_apport_est")
R_apport_est_melt$R_apport_sim<-R_apport_sim_melt[,4]
R_apport_est_melt$R_apport_bias<-((R_apport_est_melt$R_apport_est-R_apport_est_melt$R_apport_sim)/R_apport_est_melt$R_apport_sim)*100


r.apport.long<-melt(R_apport_est_melt, id=c("Year","Reg","Nsim"))
r.apport.long$Reg<-as.character(r.apport.long$Reg)

median_R_apport<-data.frame(r.apport.long%>% group_by(Reg,variable) %>% summarise(med=median(value),min=min(value),max=max(value)))


r.apport.plot.gg<-ggplot(R_apport_est_melt, aes(x=as.factor(Reg), y=R_apport_est, group=Reg)) +
  #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data=subset(median_R_apport,variable=="R_apport_est"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
  geom_point(data=subset(median_R_apport,variable=="R_apport_sim"), aes(x=Reg,y=med), col="black", shape=16,cex=2.0) + 
  ggtitle("Recruitment Apportionment")+
  ylab("Recruitment Apportionment")+
  xlab("Area")+
  ylim(0,1)+
  my_theme

r.apport.bias.gg<-ggplot(R_apport_est_melt, aes(x=as.factor(Reg), y=R_apport_bias, group=Reg)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_point(data=subset(median_R_apport,variable=="R_apport_bias"), aes(x=Reg,y=med),fill=median.col, shape=21,size=2.0) + 
  #geom_point(data = rec.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Recruitment Apportionment Bias")+
  ylab("Relative % Difference")+
  xlab("Area")+
  #ylim(quantile(R_apport_est_melt$R_apport_bias, c(.05,.95)))+
  my_theme

write.csv(median_R_apport,"R_apport_medians.csv")

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

#################################################
# Rec deviations

#build data.frame
rec.dev.data.sim<-data.frame(Dat=rep("SIM",nrow(rec_devs_df_sim)),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))
rec.dev.data.sim<-cbind(rec.dev.data.sim,rec_devs_df_sim)
rec.dev.sim.temp<-subset(rec.dev.data.sim,!Years==1)

rec.dev.data<-data.frame(Dat=rep("EST",nrow(rec_devs_df_est)),Years=rep((2:nyrs),nreg*2),Reg=rep(1:nreg,each=nyrs-1))
rec.dev.data<-cbind(rec.dev.data,rec_devs_df_est)
rec.dev.data<-rbind(rec.dev.data,rec.dev.sim.temp)


rec.dev.long<-melt(rec.dev.data, id=c("Dat","Years","Reg"))
rec.dev.long$Reg<-as.character(rec.dev.data$Reg)

#calculate the sum across areas 
mean.rec.dev<-data.frame(rec.dev.long %>% group_by(Dat, Years, variable) %>% summarise(value=mean(value)))

mean.rec.dev$Reg<-rep("System Mean",nrow(mean.rec.dev))
rec.dev.long<-rbind(mean.rec.dev,rec.dev.long)


#separate again for plotting
rec.dev.est<-rec.dev.long[rec.dev.long$Dat=="EST",]
rec.dev.sim<-rec.dev.long[rec.dev.long$Dat=="SIM",]


#calculate the percent bias
rec.dev.est$val.true<-rec.dev.sim$value
rec.dev.est$bias=((rec.dev.est$value-rec.dev.est$val.true)/rec.dev.est$val.true)*100

#calc medians table
rec.dev.est.med <- rec.dev.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
rec.dev.meds <- rec.dev.est %>% group_by(Reg) %>%
  mutate(med.est = median(value),med.true=median(val.true))


#generate Rec Plot
rec.dev.plot.gg<-ggplot(rec.dev.est, aes(x=as.factor(Years), y=value)) +
  #geom_hline(aes(yintercept = med.est, group = Reg), colour = 'darkred', size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_point(data = rec.dev.est.med, aes(x=as.factor(Years),y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = rec.dev.est.med, aes(x=Years-1,y=med.sim),lty=1,lwd=0.5) + 
  geom_point(data = rec.dev.est.med, aes(x=as.factor(Years),y=med.sim), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(2,nyrs+1,4))+
  ggtitle("Recruitment Deviations")+
  ylab("Recruitment Deviation")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
rec.dev.bias.gg<-ggplot(rec.dev.est, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_line(data = rec.dev.est.med, aes(x=Years-1,y=med.bias),lty=1,lwd=0.5) + 
  geom_point(data = rec.dev.est.med, aes(x=as.factor(Years),y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(2,nyrs+1,4))+
  ggtitle("Recruitment Deviation Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.,scales="free")+
  #ylim(quantile(rec.dev.est$bias, c(.05,.95)))+
  my_theme


#save rec bias calcs
write.csv(rec.dev.est.med,"Rec_Dev_Bias.csv")



############################################
# Initial abundance


init.abund.data<-data.frame(Dat=c(rep("SIM",nrow(init_abund_df_est)),rep("EST",nrow(init_abund_df_est))),Age=rep(1:na,(nreg*nreg)*2),Reg=rep(1:nreg,each=na*nreg))

init.abund.data<-cbind(init.abund.data,rbind(init_abund_df_sim,init_abund_df_est))
init.abund.melt<-melt(init.abund.data, id=c("Dat","Age","Reg"))


#sum across region
#calculate the sum across areas 
total.init.abund<-data.frame(init.abund.melt %>% group_by(Dat, Age, Reg,variable) %>% summarise(value=sum(value)))

#add the system
total.init.abund2<-data.frame(total.init.abund %>% group_by(Dat, Age,variable) %>% summarise(value=sum(value)))
total.init.abund2$Reg<-"System"

total.init.abund<-rbind(total.init.abund,total.init.abund2)

init.abund.est<-subset(total.init.abund,Dat=="EST")
init.sim<-subset(total.init.abund,Dat=="SIM")


init.abund.est$value.true<-init.sim$value
init.abund.est$bias<-((init.abund.est$value-init.abund.est$value.true)/init.abund.est$value.true)*100


init.abund.meds<-data.frame(init.abund.est %>% group_by(Age, Reg) %>% summarise(median.est=median(value),median.sim=median(value.true),median.bias=median(bias),min.bias=min(bias),max.bias=max(bias)))


init.abund.plot.gg<-ggplot(init.abund.est, aes(x=as.factor(Age), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = init.abund.meds, aes(x=Age,y=median.est), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = init.abund.meds, aes(x=Age,y=median.sim),lty=1, lwd=0.5) + 
  geom_point(data = init.abund.meds, aes(x=Age,y=median.sim), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,na,1))+
  ggtitle("Initial Abundance")+
  ylab("Initial Abundance")+
  xlab("Age")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#bias
init.abund.bias.gg<-ggplot(init.abund.est, aes(x=as.factor(Age), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_line(data = init.abund.meds, aes(x=Age,y=median.bias),lty=1,lwd=0.5) + 
  geom_point(data = init.abund.meds, aes(x=as.factor(Age),y=median.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(2,nyrs+1,4))+
  ggtitle("Initial Abundance Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(quantile(init.abund.est$bias, c(.05,.95)))+
  my_theme


write.csv(init.abund.meds,"Initial_abund_medians.csv")


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
ssb.est$bias=((ssb.est$value-ssb.est$val.true)/ssb.est$val.true)*100

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
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
bio.meds <- bio.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate bio Plot
bio.plot.gg<-ggplot(bio.est, aes(x=as.factor(Years), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = bio.est.med, aes(x=Years,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = bio.est.med, aes(x=Years,y=med.sim),lty=1, lwd=0.5) + 
  geom_point(data = bio.est.med, aes(x=Years,y=med.sim), fill="black", shape=16,size=1.0) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass")+
  ylab("Biomass")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
bio.bias.gg<-ggplot(bio.est, aes(x=as.factor(Years), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_line(data = bio.est.med, aes(x=Years,y=med.bias),lty=1, lwd=0.5) + 
  geom_point(data = bio.est.med, aes(x=Years,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,5))+
  ggtitle("Biomass Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg~.)+
  #ylim(quantile(bio.est$bias, c(.05,.95)))+
  my_theme


#save rec bias calcs
write.csv(bio.est.med,"Biomass_Bias.csv")


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
total.fmax<-data.frame(fmax.long %>% group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.fmax$Reg<-rep("System",nrow(total.fmax))
fmax.long<-rbind(total.fmax,fmax.long)


#separate again for plotting
fmax.est<-fmax.long[fmax.long$Dat=="EST",]
fmax.sim<-fmax.long[fmax.long$Dat=="SIM",]


#calculate the percent bias
fmax.est$val.true<-fmax.sim$value
fmax.est$bias=((fmax.est$value-fmax.est$val.true)/fmax.est$val.true)*100
#fmax.est$bias=((fmax.est$val.true-fmax.est$value))

#calc medians table
fmax.est.med <- fmax.est %>% group_by(Reg,Years) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
fmax.meds <- fmax.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


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
  #ylim(quantile(fmax.est$bias, c(.05,.95)))+
  my_theme


#save rec bias calcs
write.csv(fmax.est.med,"Fmax_Bias.csv")



# Selectivity at age
#build data.frame
f.select.data<-data.frame(Dat=c(rep("SIM",nrow(select_age_df_sim)),rep("EST",nrow(select_age_df_est))),Age=rep(ages,nreg*2),Reg=rep(1:nreg,each=na))

f.select.data<-cbind(f.select.data,rbind(select_age_df_sim,select_age_df_est))
f.select.long<-melt(f.select.data, id=c("Dat","Age","Reg"))
f.select.long$Reg<-as.character(f.select.long$Reg)



#separate again for plotting
f.select.est<-f.select.long[f.select.long$Dat=="EST",]
f.select.sim<-f.select.long[f.select.long$Dat=="SIM",]


#calculate the percent bias
f.select.est$val.true<-f.select.sim$value
f.select.est$bias=((f.select.est$value-f.select.est$val.true)/f.select.est$val.true)*100

#calc medians table
f.select.est.med <- f.select.est %>% group_by(Reg,Age) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
f.select.meds <- f.select.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate fmax Plot
f.select.plot.gg<-ggplot(f.select.est, aes(x=as.factor(Age), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = f.select.est.med, aes(x=Age,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_point(data = f.select.est.med, aes(x=Age,y=med.sim), fill="black", shape=16,size=1.0) + 
  geom_line(data = f.select.est.med, aes(x=Age,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Fishery Selectivity")+
  ylab("Selectivity")+
  xlab("Age")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
f.select.bias.gg<-ggplot(f.select.est, aes(x=as.factor(Age), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ",alpha=0.6)+
  geom_line(data = f.select.est.med, aes(x=Age,y=med.bias),lty=1, lwd=0.5) + 
  geom_point(data = f.select.est.med, aes(x=Age,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Fishery Selectivity Bias")+
  ylab("Relative % Difference")+
  xlab("Age")+
  facet_grid(Reg~.)+
  my_theme


#save rec bias calcs
write.csv(f.select.est,"F_select_age_bias.csv")



####################################
#survey selectivity
s.select.data<-data.frame(Dat=c(rep("SIM",nrow(select_age_survey_df_sim)),rep("EST",nrow(select_age_survey_df_est))),Age=rep(ages,nreg*2),Reg=rep(1:nreg,each=na))

s.select.data<-cbind(s.select.data,rbind(select_age_survey_df_sim,select_age_survey_df_est))
s.select.long<-melt(s.select.data, id=c("Dat","Age","Reg"))
s.select.long$Reg<-as.character(s.select.long$Reg)



#separate again for plotting
s.select.est<-s.select.long[s.select.long$Dat=="EST",]
s.select.sim<-s.select.long[s.select.long$Dat=="SIM",]


#calculate the percent bias
s.select.est$val.true<-s.select.sim$value
s.select.est$bias=((s.select.est$value-s.select.est$val.true)/s.select.est$val.true)*100

#calc medians table
s.select.est.med <- s.select.est %>% group_by(Reg,Age) %>%
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))

#calc vals for plots
s.select.meds <- s.select.est %>% group_by(Reg) %>%
  mutate(med = median(value),med.true=median(val.true))


#generate fmax Plot
s.select.plot.gg<-ggplot(s.select.est, aes(x=as.factor(Age), y=value)) +
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_point(data = s.select.est.med, aes(x=Age,y=med.est), fill=median.col, shape=21,size=2.0) + 
  geom_point(data = s.select.est.med, aes(x=Age,y=med.sim), fill="black", shape=16,size=1.0) + 
  geom_line(data = s.select.est.med, aes(x=Age,y=med.sim),lty=1) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Survey Selectivity")+
  ylab("Selectivity")+
  xlab("Age")+
  facet_grid(Reg~.)+
  #ylim(-5,5)+
  my_theme


#generate Rec bias plot
s.select.bias.gg<-ggplot(s.select.est, aes(x=as.factor(Age), y=bias)) +
  geom_hline(aes(yintercept = 0, group = Reg), colour = 'red',size=0.5,lty=2)+
  geom_violin(fill=vio.col,trim=F,bw="SJ", alpha=0.6)+
  geom_line(data = s.select.est.med, aes(x=Age,y=med.bias),lty=1, lwd=0.5) + 
  geom_point(data = s.select.est.med, aes(x=Age,y=med.bias), fill=median.col, shape=21,size=1.5) + 
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Survey Selectivity Bias")+
  ylab("Relative % Difference")+
  xlab("Age")+
  facet_grid(Reg~.)+
  #ylim(-500,500)+
  my_theme


#save rec bias calcs
write.csv(s.select.est.med,"Surv_select_age_bias.csv")



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
  summarise(med.est=median(value),med.sim=median(val.true), med.bias = median(bias),max=max(bias),min=min(bias))


#calc vals for plots
move.meds <- move.est %>% group_by(Reg_from,Reg_to) %>%
  mutate(med = median(value),med.true=median(val.true))


move.plot.gg.alt<-ggplot(move.est, aes(x=as.factor(Year), y=value, col=Reg_to), group=Reg_to) +
  geom_violin(aes(fill=Reg_to),trim=F,bw="SJ",position = position_dodge(width=0.8), alpha=0.2)+
  geom_line(data = move.est.med, aes(x=Year,y=med.est, group=Reg_to))+
  geom_point(data = move.est.med, aes(x=Year,y=med.est, group=Reg_to),position = position_dodge(width=0.8), fill=median.col, shape=21,size=2.0) + 
  geom_line(data = move.est.med, aes(x=Year,y=med.sim, group=Reg_to),lty=2)+
  geom_point(data = move.est.med, aes(x=Year,y=med.sim, group=Reg_to, fill=Reg_to),position = position_dodge(width=0.8),shape=21,size=1.0) + 
  scale_color_brewer("Move To",palette = "Set1")+
  scale_fill_brewer("Move To",palette = "Set1")+
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Movement")+
  ylab("Movement Rate")+
  xlab("Year")+
  #facet_grid_paginate(Reg_from~Age, ncol = 2, nrow = 2, page = i)+
  facet_grid(Reg_from~Age)+
  #ylim(-5,5)+
  my_theme


move.bias.gg.alt<-ggplot(move.est, aes(x=as.factor(Year), y=bias, col=Reg_to), group=Reg_to) +
  geom_violin(aes(fill=Reg_to),trim=F,bw="SJ",position = position_dodge(width=0.8), alpha=0.2)+
  geom_hline(aes(yintercept = 0, group = Reg_to), colour = 'black',size=0.5,lty=2)+
  geom_line(data = move.est.med, aes(x=Year,y=med.bias, group=Reg_to))+
  geom_point(data = move.est.med, aes(x=Year,y=med.bias,group=Reg_to),position = position_dodge(width=0.8), fill=median.col, shape=21,size=1.5) + 
  scale_color_brewer("Move To",palette = "Set1")+
  scale_fill_brewer("Move To",palette = "Set1")+
  scale_x_discrete(breaks=seq(0,nyrs,1))+
  ggtitle("Movement Bias")+
  ylab("Relative % Difference")+
  xlab("Year")+
  facet_grid(Reg_from~Age)+
  #ylim(quantile(move.est$bias, c(.05,.95)))+
  my_theme


write.csv(move.est.med,"Move_Bias.csv")

ssb.table<-spread(ssb.long, key = variable, value = value)
rec.table<-spread(rec.long, key = variable, value = value)
fmax.table<-spread(fmax.long, key = variable, value = value)
move.table<-spread(move.long, key = variable, value = value)

write.csv(ssb.table,"ssb_table.csv")
write.csv(rec.table,"rec_table.csv")
write.csv(fmax.table,"fmax_table.csv")
write.csv(move.table,"move_table.csv")

}

#######################################
# Create a summary PDF
#######################################
{
pdf("Simulation_Summary.pdf",paper="letter",height = 11, width=8)
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

text3<-paste0("Proportion of Sims Converged: ", sum(Sim_Stats$Converged)," of ", nsim, " (",conv,")")

time<-paste0("Report Date:  ", Sys.time())

text.all<-paste(text1,text2,text3,sep = "\n")


# Create a text grob
tgrob <- textGrob(text.all,just = "centre")



grid.arrange(ncol=1,nrow=4,
             bottom=time,
             center=textGrob("TIM Simulation Summary", gp=gpar(fontsize=18,font=3)),
             center=tgrob)



#add rec plots
grid.arrange(ncol = 2,
             top="Input Parameters",
             q.plot.gg,rave.plot.gg,
             q.bias.gg,rave.bias.gg)


grid.arrange(ncol = 1,
             top="Recruitment Apportionment",
             r.apport.plot.gg, r.apport.bias.gg)


grid.arrange(ncol = 1,
             top="Recruitment Estimation",
             rec.plot.gg, rec.bias.gg)

grid.arrange(ncol = 1,
             top="Recruitment Deviations",
             rec.dev.plot.gg, rec.dev.bias.gg)

grid.arrange(ncol = 1,
             top="Initial Abundance",
             init.abund.plot.gg, init.abund.bias.gg)

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

grid.arrange(ncol = 1,
             top="Fishery Selectivity",
             f.select.plot.gg, f.select.bias.gg)

grid.arrange(ncol = 1,
             top="Fishery Selectivity",
             s.select.plot.gg, s.select.bias.gg)

#grid.arrange(ncol = 1,
#             top="Movement",
#             move.plot.gg, move.bias.gg)


#pdf("Simulation_Summary.pdf",paper="letter",height = 11, width=8)


for(i in 1:(na/2)){
 print( move.plot.gg.alt+facet_grid_paginate(Reg_from~Age, ncol = 2, nrow = 2, page = i))
}
for(i in 1:(na/2)){
  print( move.bias.gg.alt+facet_grid_paginate(Reg_from~Age, ncol = 2, nrow = 2, page = i))
}

##############################
# Plot Likelihood components
##############################

#grid.arrange(ncol = 2,nrow=4,
#             top="Likelihood Components",
#             ggplotGrob(Like.hist(4)),
#             ggplotGrob(Like.hist(5)),
#             ggplotGrob(Like.hist(6)),
#             ggplotGrob(Like.hist(7)),
#             ggplotGrob(Like.hist(8)),             
#             ggplotGrob(Like.hist(9)),
#             ggplotGrob(Like.hist(10)),
             #ggplotGrob(Like.hist(11)), #rec penalty
#             ggplotGrob(Like.hist(12))
#             )



dev.off()


}#end PDF code

}#end of full run - fingers crossed



###############################################
#some notes for things to add


#PLOTS
#Timeseries Trajectories compared to true (bio, ssb, recruitment, F, T.) (bean/violin plots overlayed on true point estimates)
#Fit to data
#Timeseries (yield, survey)
#Age comp (bubble plots?)
#Residuals (bubble plots? Timeseries)
#Parameter Bias (violin plots of %RE or MARE)
#Likelihood components (report out max gradient, convergence, etc.)
#Table of correlations



