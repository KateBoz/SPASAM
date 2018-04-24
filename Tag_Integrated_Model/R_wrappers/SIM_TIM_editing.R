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




nsim <-10   # How many simulations to perform



##################################
### setting up the directories
##################################

# master file with holding the runs 
direct_master<-"F:\\TIM_editing\\"

#list files in the directory
files<-list.files(direct_master)

#select the file you want to run
#if only running 1 folder set i to the number corresponding to the folder you want to run
i=1


##############################################################
#if running the whole folder - This will come later
#for(i in 1:length(files)){

#OM Location
OM_direct<-paste0(direct_master,files[i],"\\Operating_Model",sep="")
OM_name<-"TIM_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-paste0(direct_master,files[i],"\\Estimation_Model",sep="") #location of run(s)
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.

#Build diagnostics folder
dir.create(paste0(direct_master,files[i],"\\Diagnostics",sep=""))
diag_direct<-paste0(direct_master,files[i],"\\Diagnostics",sep="")

############## Run simulations #############################


#Set up convergence record

conv<-data.frame(SimN=seq(1:nsim),Converged=rep(NA,nsim))


#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#set up text progress bar
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


ls=foreach(j=1:nsim,.options.snow = opts) %dopar% {
  
 temp<-NA #to keep track of convergence
  
  dir.create(paste0(diag_direct,"\\Run",j,sep=""))
  dir.create(paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep=""))
  dir.create(paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep=""))

  invisible(file.copy(from=paste0(OM_direct,"\\TIM_OM.exe",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(OM_direct,"\\TIM_OM.tpl",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep="")))
  invisible(file.copy(from=paste0(OM_direct,"\\TIM_OM.dat",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Operating_Model",sep="")))
  
  invisible(file.copy(from=paste0(EM_direct,"\\TIM_EM.exe",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep="")))
  invisible(file.copy(from=paste0(EM_direct,"\\TIM_EM.tpl",sep=""),to=paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep="")))

  
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

#del junk
#invisible(shell(paste("del TIM_OM.par", sep="")))
#invisible(shell(paste("del TIM_OM.rep", sep="")))
invisible(shell("del TIM_OM.exe")) 
invisible(shell("del TIM_OM.bar"))
invisible(shell("del TIM_OM.log"))
invisible(shell("del fmin.log"))
invisible(shell("del variance"))


# move the .rep from the OM_direct and change name. 
from<-paste0(paste0(diag_direct,"\\Run",j,"\\Operating_Model\\",OM_name,".rep",sep=""))
to<-paste0(paste0(diag_direct,"\\Run",j,"\\Estimation_Model\\",EM_name,".dat",sep=""))

file.copy(from = from,  to = to)


#Run the estimation model

#run the EM
setwd(paste0(diag_direct,"\\Run",j,"\\Estimation_Model",sep=""))
invisible(shell(paste0(EM_name," -nox -ind"),wait=T))


if(file.exists("TIM_EM.cor")==TRUE)
   
{
  temp=1
  
  invisible(shell("del admodel.dep"))
  invisible(shell("del admodel.hes")) 
  invisible(shell("del admodel.cov"))
  invisible(shell(paste("del TIM_EM.b01",sep="")))
  invisible(shell(paste("del TIM_EM.P01",sep="")))
  invisible(shell(paste("del TIM_EM.R01",sep="")))
  invisible(shell(paste("del TIM_EM.b02",sep="")))
  invisible(shell(paste("del TIM_EM.P02",sep="")))
  invisible(shell(paste("del TIM_EM.R02",sep="")))
  invisible(shell(paste("del TIM_EM.b03",sep="")))
  invisible(shell(paste("del TIM_EM.P03",sep="")))
  invisible(shell(paste("del TIM_EM.R03",sep="")))
  invisible(shell(paste("del TIM_EM.bar",sep="")))
  invisible(shell(paste("del TIM_EM.eva",sep="")))
  invisible(shell(paste("del TIM_EM.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))

}


if(file.exists("TIM_EM.cor")==FALSE)
{
  temp=0
  
  #then delete the junk
  invisible(shell("del TIM_EM.exe")) 
  invisible(shell("del admodel.dep"))
  invisible(shell("del admodel.hes")) 
  invisible(shell("del admodel.cov"))
  invisible(shell(paste("del TIM_EM.b01",sep="")))
  invisible(shell(paste("del TIM_EM.P01",sep="")))
  invisible(shell(paste("del TIM_EM.R01",sep="")))
  invisible(shell(paste("del TIM_EM.b02",sep="")))
  invisible(shell(paste("del TIM_EM.P02",sep="")))
  invisible(shell(paste("del TIM_EM.R02",sep="")))
  invisible(shell(paste("del TIM_EM.b03",sep="")))
  invisible(shell(paste("del TIM_EM.P03",sep="")))
  invisible(shell(paste("del TIM_EM.R03",sep="")))
  invisible(shell(paste("del TIM_EM.bar",sep="")))
  invisible(shell(paste("del TIM_EM.eva",sep="")))
  invisible(shell(paste("del TIM_EM.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
 
}

return(temp)

}


close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 

#Convergence data from run
for(i in 1:nsim){conv[i,2]<-ls[[i]]}
 
#ending
unlink(getwd())
setwd(direct_master)


##############################################################
#Grabbing Values
##############################################################

# create a data frame for the run to hold true and estimated values from OM/EM

#pull dimensions for building data frames for plotting
df_build<-(paste0(diag_direct,"\\Run",1,"\\Estimation_Model"))

#get dimensions for scenario
out<-readList(paste0(df_build,"\\",EM_name,".rep")) #read in .rep file

#pull info about the model
na<-out$nages
nyrs<-out$nyrs
npops<-out$npops
nreg<-out$nregions
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)


if(npops>1){
  nreg=sum(nreg)}


########################################################################
# Building Run parallel to save the values we want to keep


#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)


#set up text progress bar
pb <- txtProgressBar(max = nsim, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


vg=foreach(i=1:nsim,.options.snow = opts,.combine='cbind',.packages = c('PBSmodelling','data.table')) %dopar% {

#read in rep
  out<-readList(paste0(diag_direct,"\\Run",i,"\\Estimation_Model\\",EM_name,".rep",sep=""))
  
#Save rec
  rec_sim_temp<-data.table(Rec_Sim = c(t(out$recruits_BM_TRUE)))
  rec_est_temp<-data.table(Rec_Est= c(t(out$recruits_BM)))
  
#Save SSB
  ssb_sim_temp<-data.table(SSB_Sim = c(t(out$SSB_region_TRUE)))
  ssb_est_temp<-data.table(SSB_Est = c(t(out$SSB_region)))
  
#Save_apport
  apport_sim_temp<-data.table(Apport_Sim = c(t(out$Rec_apport_TRUE)))
  apport_est_temp<-data.table(Apport_Est = c(t(out$Rec_apport)))

#

#return all the matrices as a list
  return(list(rec_sim=rec_sim_temp,rec_est=rec_est_temp,ssb_sim=ssb_sim_temp,ssb_est=ssb_est_temp,apport_sim=apport_sim_temp,apport_est=apport_est_temp))

} #end parallel
  

close(pb)
stopCluster(cl) #end the cluster for parallel processing
closeAllConnections() 

#save this for later
setwd(diag_direct)
saveRDS(vg, file="Sim_run.RData")

################################################################################
# Load Sim results for plotting
#################################################################################
#reload for later plots
setwd(diag_direct)
Sim_Results<-readRDS('Sim_run.RData')



#Recruitment
rec_df_sim<-matrix(NA,nyrs*nreg,nsim)
rec_df_est<-matrix(NA,nyrs*nreg,nsim)

#SSB
ssb_df_sim<-matrix(NA,nyrs*nreg,nsim)
ssb_df_est<-matrix(NA,nyrs*nreg,nsim)

#Rec Apport
app_df_sim<-matrix(NA,nyrs*nreg,nsim)
app_df_est<-matrix(NA,nyrs*nreg,nsim)



# populate the matrices for plotting
for(i in 1:nsim){
rec_df_sim[,i]<-unlist(Sim_Results[1,i])
rec_df_est[,i]<-unlist(Sim_Results[2,i])
ssb_df_sim[,i]<-unlist(Sim_Results[3,i])
ssb_df_est[,i]<-unlist(Sim_Results[4,i])

#app_df_sim[,i]<-unlist(Sim_Results[5,i])
#app_df_est[,i]<-unlist(Sim_Results[6,i])
}


#####################################
# Recruitment plot

#build data.frame
rec.data<-data.frame(Dat=c(rep("SIM",nrow(rec_df_sim)),rep("EST",nrow(rec_df_est))),Years=rep(years,nreg*2),Reg=rep(1:nreg,each=nyrs))

rec.data<-cbind(rec.data,rbind(rec_df_sim,rec_df_est))
rec.long<-melt(rec.data, id=c("Dat","Years","Reg"))
rec.long$Reg<-as.character(rec.data$Reg)

#calculate the sum across areas 
total.rec<-data.frame(rec.long %>%
  group_by(Dat, Years, variable) %>% summarise(value=sum(value)))

total.rec$Reg<-rep("System",nrow(total.rec))
rec.long<-rbind(total.rec,rec.long)

#split for multipanel plot
split.by.reg<-split(rec.long,rec.long$Reg)
split.by.reg[[1]]

#no scientific notation
options(scipen=1000)


#build plot function
rec.plot<-function(k=1){

#test plot
beanplot(jitter(value,0.01) ~ interaction(Dat,factor(Years)),split.by.reg[[k]], ll = 0.04,main = split.by.reg[[k]]$Reg[k], side = "both", xlab="Years",ylab="Total Recruitment",col = list("grey40", c("grey90", "black")),axes=T, overallline = "median",beanlinewd = 1)

axis(2)

legend("topright", fill = c("grey40", "grey90"),
       legend = c("True", "Estimated"), box.lty=0, bty="n")

}


#####################################
# SSB plot

#Create matrices to get data for recruitment




########################################################
# Create a print out of plots for each sim
########################################################

#Print rec plots

# print Rec plots
RecPlots <- lapply (1:(nreg+1), function(i) {
  
  png("rec_plot.png",height=4.5,width=8, units = "in", res=800)
  
  rec.plot(i)
  
  dev.off()
  rasterGrob(readPNG("rec_plot.png", native = FALSE),
             interpolate = FALSE)
})



#Print ssb plots

# will go here next


#######################################
# Create a summary PDF
#######################################

pdf("Simulation_Summary.pdf",paper="special",width=8,height=11)

plot(0:10, type = "n", xaxt="n", yaxt="n", bty="n", xlab = "", ylab = "")
text(6, 8, "Spatial Model Simulation Summary")
text(6, 7.5, paste0("Total number of Sims: ",nsim))

do.call(grid.arrange, c(RecPlots, ncol = 1))

dev.off()





