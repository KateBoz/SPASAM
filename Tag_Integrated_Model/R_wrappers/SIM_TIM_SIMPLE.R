##################################################
###################################################
# SIMPLE SIMS CODE FOR SPASAM MODELS
# Created by: Katelyn Bosley
# Updated 9/21/2018
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
}
load_libraries()



###################################
#Setting up the simulations
##################################


#1) Do you want to run the diagnostics code before running the sims?

# if diagnostic code is run before the sims, the par file from the diagnostic run will be used create a pin for simulation runs.

diag.run<-1
   # ==0 NO Do NOT run the diagnostics plots for a single run before the sims
   # ==1 YES run the diagnostics plots for a single run

# if running the diagnostic run set residual switch for different values plotted

resid.switch=2
   # ==1 straight residual (TRUE-ESTIMATED; not % of true)
   # ==2 Relative percent difference ((TRUE/ESTIMATED)/TRUE *100)



#2) Set number of simulations to perform
nsim <-50

#########################################
### setting up the directories
#########################################

# To run this simulation a folder for each scenario will need to be placed in the master directory. Each scenario folder will need separate folders named 'Operating_Model' and 'Estimation_model' with the .exe files and .dat for the OM only. The code will do the rest. 

#If the diagnostic switch ==1 the TIM_diagnostics_SIM.R code will also have to be in the master directory.


#3) set master file with holding the runs 
direct_master<-"C:\\Users\\katelyn.bosley\\Desktop\\Mismatch\\CAPAM_RUNS\\"
setwd(direct_master)

#list files in the directory to choose the correct one
files<-list.files(direct_master) # these folders in the master will be the individual scenarios 

#select the file with the scenario you want to run
#if only running 1 folder set i to the number corresponding to the folder you want to run

folder.num=9
i=folder.num

#run.sims<-function(i=folder.num,nsim=nsim,direct_master=direct_master){
#run the whole code
{ 
  
#if running the several folders use the loop - This will come later
#for(i in 1:3){

##############################################################
#OM Location
OM_direct<-paste0(direct_master,files[i],"\\Operating_Model",sep="")
OM_name<-"TIM_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-paste0(direct_master,files[i],"\\Estimation_Model",sep="") #location of run(s)
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.

#Build diagnostics/results folder
dir.create(paste0(direct_master,files[i],"\\Diagnostics",sep=""))
diag_direct<-paste0(direct_master,files[i],"\\Diagnostics",sep="")


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
source("TIM_Diagnostics_SIM.R") # make sure this code is in the direct_master folder.

invisible(file.rename(from=paste0(EM_direct,"\\Model_Diagnostics.pdf",sep=""),to=paste0(diag_direct,"\\Model_Diagnostics.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Tag_Residuals.pdf",sep=""),to=paste0(diag_direct,"\\Tag_Residuals.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Correlation_Matrix.pdf",sep=""),to=paste0(diag_direct,"\\Correlation_Matrix.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\Standard_Error_table.pdf",sep=""),to=paste0(diag_direct,"\\Standard_Error_table.pdf",sep="")))
invisible(file.rename(from=paste0(EM_direct,"\\output.txt",sep=""),to=paste0(diag_direct,"\\diag_output.txt",sep="")))
}



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
  #save dat
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""),  
            to = paste0(results_good,"\\Run",j,".dat",sep=""))
  #save gradient file
  file.copy(from = paste0(runs_dir,"\\Run",j,"\\Estimation_Model\\gradient.dat",sep=""),  
            to = paste0(results_good,"\\gradient",j,".dat",sep=""))
  }
  
# if not converged..save results to results_bad
if(file.exists("TIM_EM.cor")==FALSE)
{
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
            to = paste0(results_good,"\\gradient",j,".dat",sep=""))
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


}  # end of simulations and summary save

}
  

