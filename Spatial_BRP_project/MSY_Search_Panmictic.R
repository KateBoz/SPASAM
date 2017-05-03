######################################################
# Another attempt to run the MSY search on the panmictic population
# Created by Katelyn Bosley
# Date: 3/21/2017
###################################################
#need to reset the working directory in for run

# remove previous objects from workspace
rm(list = ls())

#Sablefish
WD<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\MENHADEN\\Panmictic"

setwd(WD)
#wd<-WD

#install libraries
load_libraries<-function() {
suppressWarnings(suppressMessages(require(PBSmodelling)))
suppressWarnings(suppressMessages(require(matrixStats)))
suppressWarnings(suppressMessages(require(TeachingDemos)))
suppressWarnings(suppressMessages(require(snowfall)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library(snow)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doSNOW)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(spatstat)))
suppressWarnings(suppressMessages(library(alphahull)))
}
load_libraries()


#Setting up the F values to iterate over
F.name<-"input_F"
F.start<-0
F.end<-5
it<-0.025

F.Test<-seq(F.start,F.end,it) #F values to cycle through
Ntrial<-length(F.Test) # count of total number of trials

#read original .dat
update=readLines("Spatial_BRP_panmictic.dat",n=-1)
#pull important values - for generalized version  - crossover from previous
nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("npopulations",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))

#setting up the files to run through
dir.create(paste0(WD,"\\MSY Results",sep=""))
dir.create(paste0(WD,"\\MSY Results\\Figures",sep=""))
dir.create(paste0(WD,"\\MSY Results\\Report Files",sep=""))


###########################################
# THE SEARCH FUNCTION
##########################################

#MSY search function
MSY_search<-function(wd=WD,F.test=F.Test,ntrial=Ntrial) {
  
#do parallel processing
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)
  
#set up text progress bar
pb <- txtProgressBar(max = ntrial , style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
  
stime <- system.time({
  
  #run parallel  
  ls<-foreach(i=1:ntrial,.options.snow = opts) %dopar% {

#for(i in 1:ntrial){
  dir.create(paste0(wd,"\\MSY Results\\Run",i,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.exe",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.dat",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.tpl",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.tpl",sep="")))
  
  
  #set new directory for running the model
  setwd(paste0(wd,"\\MSY Results\\Run",i,sep="")) # now set the working direcctory as the run# file
  update=readLines("Spatial_BRP_panmictic.dat",n=-1)
  update[grep("input_F",update)+1]<-F.test[i]
  
  
  #update[(grep(F.grab,update)+1):(grep(F.grab,update)+sum(nregions))]=F.region
  
  writeLines(update,"Spatial_BRP_panmictic.dat")
  
  ### Run ADMB with updated F
  invisible(shell("Spatial_BRP_panmictic -nohess",wait=T)) #show.output.on.console=FALSE))  
  
  #clean non-needed files and move to results folder
  invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))
  
  
   } #end of code for MSY search
})

stopCluster(cl)  #end parallel

stime 

#setting up the values from the runs to plot
#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = ntrial,ncol=10))
names(msy_results)<-c("trial","F","biomass_total_start","biomass_total_end","yield_total","harvest_rate_total_bio","depletion_total","SSB_total_start","SSB_total_end","Bratio_total")


#set wd to report files
wd_results<-paste0(wd,"\\MSY Results\\Report Files",sep="")
wd_figs<-paste0(wd,"\\MSY Results\\Figures",sep="")

#set to results for grabbing values
setwd(wd_results)


# function to pull in results and grab MSY
get_msy<-function() {
#run  loop to get all the values into msy_results into a csv for plotting
for(i in 1:ntrial){
out=readList(paste0("Report",i,".rep",sep=""))

#store results to a full spreadsheet
temp<-c(i,F.test[i],out$biomass_total[1],out$biomass_total[nyrs],out$yield_total[nyrs],out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],out$SSB_total[1],out$SSB_total[nyrs], out$Bratio_total[nyrs])


msy_results[i,]<-temp

#results
write.csv(msy_results,"MSY_results.csv")

#MSY vals
t<-which(msy_results$yield_total==max(msy_results$yield_total))
msy_pan<-msy_results[t,]
write.csv(msy_pan,"MSY_true.csv")

} #end value grab

#move results files to figs folder
invisible(file.copy(from=paste0(wd_results,"\\MSY_results.csv",sep=""),to=paste0(wd_figs,"\\MSY_results.csv",sep="")))
invisible(file.remove(paste0(wd_results,"\\MSY_results.csv",sep="")))
  
invisible(file.copy(from=paste0(wd_results,"\\MSY_true.csv",sep=""),to=paste0(wd_figs,"\\MSY_true.csv",sep="")))
invisible(file.remove(paste0(wd_results,"\\MSY_true.csv",sep="")))
  
#move over full report also for fun
invisible(file.copy(from=paste0(wd_results,"\\Report",t,".rep",sep=""),to=paste0(wd_figs,"\\Report",t,".rep",sep="")))

}

get_msy()


############################################
#plotting code
############################################

# Make some plots

##############
#MSY PLOTS
#############

#plotting function
MSY_plots<-function() {
  msy_results<-read.csv("msy_results.csv") #read in the data again if there were changes above

  par(mfrow = c(2,2))
  
#SSB Ratio vs. Yield
options(scipen = 50, digits = 4) # fix sci notation
plot(msy_results$Bratio_total,msy_results$yield_total, type = 'b',ylab='Yield',xlab = " Equilibruim SSB Ratio", lwd = 2)

#Harvest Rate vs Yield
plot(msy_results$harvest_rate_total_bio,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Harvest Rate", lwd = 2)

#Equilibruim Biomass vs Yield
plot(msy_results$biomass_total_end,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)

#Harvest rate vs biomass
plot(msy_results$biomass_total_end,msy_results$harvest_rate_total_bio, type = 'b',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)
}


#save the plot
setwd(wd_figs)
pdf("panmictic_plots.pdf")
MSY_plots()
dev.off()


#remove the old files from run
for(i in 1:ntrial){
  unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
  unlink(paste0(wd,"\\MSY Results\\Report Files",sep = ""),recursive = T)
}


} #end panmictic search functikon

MSY_search()

########################################################












