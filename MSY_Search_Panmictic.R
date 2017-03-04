######################################################
# Another attempt to run the MSY search on the panmictic population
# Created by Katelyn Bosley
# Date: 2/15/2017
###################################################
#need to reset the working directory in for run


# remove previous objects from workspace
rm(list = ls())

# set the working directory
#wd<-getwd()

#

wd<-"G:\\SPASAM CODING\\MS_1_CODE\\Hake\\Panmictic"

#wd<-"C:\\Users\\Katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Hake\\Panmictic"
#wd<-"C:\\Users\\katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Sablefish\\Panmictic"
setwd(wd)
wd<<-wd


#install libraries
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


#Setting up the F values to iterate over
F.name<-"input_F"
F.start<-0.05
F.end<-0.8
it<-0.05

F.test<-seq(F.start,F.end,it) #F values to cycle through
ntrial<-length(F.test) # count of total number of trials


# if using lots of different structures
#F.region1<-matrix(perm[ntrial,],ncol=nfleets,nrow=sum(nregions),byrow=T)  ####NEEDS TO BE GENERALIZED FOR UNEVEN FLEETS BY REGION (i.e. ncol)
#F.region<-apply(matrix(as.character(perm[ntrial,]),ncol=nfleets,nrow=sum(nregions),byrow=T),1,paste0,sep=" ",collapse="")


#setting up the files to run through
dir.create(paste0(wd,"\\MSY Results",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Figures",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Report Files",sep=""))


#MSY search function
MSY_search<-function() {
#set up the files for doing the iterations
for(i in 1:ntrial){
  dir.create(paste0(wd,"\\MSY Results\\Run",i,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.exe",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.dat",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP_panmictic.tpl",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.tpl",sep="")))
  
  
  #set new directory for running the model
  setwd(paste0(wd,"\\MSY Results\\Run",i,sep="")) # now set the working direcctory as the run# file
  update=readLines("Spatial_BRP_panmictic.dat",n=-1)
  
  
  #pull important values - for generalized version 
  nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
  nstocks<-as.numeric(update[(grep("nstocks",update)+1)])
  nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
  nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))
  
  update[grep("input_F",update)+1]<-F.test[i]
  
  
  #update[(grep(F.grab,update)+1):(grep(F.grab,update)+sum(nregions))]=F.region
  
  writeLines(update,"Spatial_BRP_panmictic.dat")
  
  ### Run ADMB with updated F
  invisible(shell("Spatial_BRP_panmictic -nohess",wait=T)) #show.output.on.console=FALSE))  
  
  #clean non-needed files and move to results folder
  invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP_panmictic.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))
  
  
  } #end of code for MSY search
  
}

#run the function
MSY_search()



#setting up the values from the runs to plot

#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = ntrial,ncol=8))

names(msy_results)<-c("trial","F","biomass_total","yield_total","harvest_rate_total_bio","depletion_total","SSB_total","Bratio_total")


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
temp<-c(i,F.test[i],out$biomass_total[nyrs],out$yield_total[nyrs],out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],out$SSB_total[nyrs], out$Bratio_total[nyrs])

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
plot(msy_results$biomass_total,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)

#Harvest rate vs biomass
plot(msy_results$biomass_total,msy_results$harvest_rate_total_bio, type = 'b',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)
}


#save the plot
setwd(wd_figs)
pdf("panmictic_plots.pdf")
MSY_plots()
dev.off()


















