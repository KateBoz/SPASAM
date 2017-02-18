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

#DIRECTORIES
#
#HAKE runs
#wd<-"C:\\Users\\Katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Hake\\Maturity_mismatch"
#wd<-"C:\\Users\\Katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Hake\\Panmictic"

#SABLEFISH runs
wd<-"C:\\Users\\katelyn.bosley\\Desktop\\SPASAM CODING\\MS_1_CODE\\Sablefish\\Match_selectivity"

#
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
F.start<-0
F.end<-0.6
it<-0.1


#read in .dat file to get values for setting up the runs 

update=readLines("Spatial_BRP.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("nstocks",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))


#set up the combinations of permutations per region
F.test<-seq(F.start,F.end,it) #F values to cycle through
ntrial<-length(F.test) # count of total number of trials

#set up the permutations of F by pop/fleet/region
permutation<-permutations(ntrial,sum(nfleets),F.test,repeats.allowed=TRUE)  # determine all permutations of F in each stock
n_perm<-nrow(permutation)

#setting up the files to run through
dir.create(paste0(wd,"\\MSY Results",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Figures",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Report Files",sep=""))

#set up a progress bar
pb <- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Simulation Run 0 of ",n_perm,sep=""),max=100)

#set up the files for doing the iterations
for(i in 1:n_perm){
  
  #progress bar
  setWinProgressBar(pb,(i/n_perm*100),label=paste("Simulation Run", i,"of", n_perm,"Completed"))
  
  dir.create(paste0(wd,"\\MSY Results\\Run",i,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.exe",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.dat",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.tpl",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.tpl",sep="")))
  
  
  #set new directory for running the model
  setwd(paste0(wd,"\\MSY Results\\Run",i,sep="")) # now set the working direcctory as the run# file
  update=readLines("Spatial_BRP.dat",n=-1)
  
  
  #pull important values - for generalized version 
 
  update[(grep("input_F",update)+1):(grep("input_F",update)+sum(nregions))]=permutation[i,]
  
  writeLines(update,"Spatial_BRP.dat")
  
  ### Run ADMB with updated F
  invisible(shell("Spatial_BRP -nohess",wait=T)) #show.output.on.console=FALSE))  
  
  #clean non-needed files and move to results folder
  invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))
  
} #end of code for MSY search

close(pb) #close the progress bar


  
#########################################################
#setting up and grabbing the values from the runs to plot
#########################################################

#need to generalize this port for plotting results
#par_names<-names(out)

N_par_reg<-5 # number of parameters with regional values # need to fix this up..
N_par_pop<-7 # number of parameters for stock

#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = nrow(permutation),ncol=((N_par_reg*nregions)+N_par_pop))) # number of parameters with multiple regions+number of total pop values
 
#fill in the table slowly for error checking
names(msy_results)[1]<-c("perm")


#fill in the rest will loops so can be changed as needed 

for(i in 1:nregions) 
  {
  names(msy_results)[1+i]<-paste0("F.",i)
  names(msy_results)[2+nregions]<-"biomass_total"
  names(msy_results)[2+nregions+i]<-paste("yield_region.",i,sep = "")
  names(msy_results)[3+nregions*2]<-"yield_total"
  names(msy_results)[3+nregions*2+i]<-paste("u_region.",i,sep = "")
  names(msy_results)[4+nregions*3]<-"u_region_total"
  names(msy_results)[4+(nregions*3)+i]<-paste("depletion_region",i,sep="")
  names(msy_results)[5+(nregions*4)]<-"depletion_total"
  names(msy_results)[5+(nregions*4)+i]<-paste("SSB_region",i,sep="")
  names(msy_results)[6+nregions*5]<-"SSB_total"
  names(msy_results)[7+nregions*5]<-"Bratio_total"
}

#names(msy_results)


#set wd to report files
wd_results<-paste0(wd,"\\MSY Results\\Report Files",sep="")
wd_figs<-paste0(wd,"\\MSY Results\\Figures",sep="")

#set to results for grabbing values
setwd(wd_results)

#progress bar for finding MSY
pb<- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Grabbing values from 0 of ",n_perm," simulations",sep=""),max=100)

#run loop to get all the values into msy_results into a csv for plotting
for(i in 1:n_perm){
  
#start progress  
setWinProgressBar(pb,(i/n_perm*100),label=paste("Grabbing values from ", i,"of ", n_perm," simulations",sep=""))
  
#read report 
out=readList(paste0("Report",i,".rep",sep=""))

#store results to a full spreadsheet add them in slowly for easy changes- check to see if it matches the .csv made above

temp<-c(i,permutation[i,],
        out$biomass_total[nyrs],
        out$yield_region[,nyrs],
        out$yield_total[nyrs],
        out$harvest_rate_region_bio[,nyrs],
        out$harvest_rate_total_bio[nyrs],
        out$depletion_region[,nyrs],
        out$depletion_total[nyrs],
        out$SSB_region[,nyrs],
        out$SSB_total[nyrs], 
        out$Bratio_total[nyrs])

msy_results[i,]<-temp

#results
write.csv(msy_results,"MSY_results.csv")

} #end value grab


# Getting the MSY vals
t<-which(msy_results$yield_total==max(msy_results$yield_total))
msy_pan<-msy_results[t,]
write.csv(msy_pan,"MSY_true.csv")

close(pb)


#move results files to figs folder
invisible(file.copy(from=paste0(wd_results,"\\MSY_results.csv",sep=""),to=paste0(wd_figs,"\\MSY_results.csv",sep="")))
invisible(file.remove(paste0(wd_results,"\\MSY_results.csv",sep="")))

invisible(file.copy(from=paste0(wd_results,"\\MSY_true.csv",sep=""),to=paste0(wd_figs,"\\MSY_true.csv",sep="")))
invisible(file.remove(paste0(wd_results,"\\MSY_true.csv",sep="")))

#move over full report also for fun
invisible(file.copy(from=paste0(wd_results,"\\Report",t,".rep",sep=""),to=paste0(wd_figs,"\\Report",t,".rep",sep="")))



############################################
#plotting code
############################################


##############
#MSY PLOTS


#plotting function
MSY_plots<-function() {
  msy_results<-read.csv("MSY_results.csv") #read in the data again if there were changes above

  par(mfrow = c(2,2))
  
#SSB Ratio vs. Yield
options(scipen = 50, digits = 4) # fix sci notation
plot(msy_results$Bratio_total,msy_results$yield_total, type = 'p',ylab='Yield',xlab = " Equilibruim SSB Ratio", lwd = 2)

#Harvest Rate vs Yield
plot(msy_results$u_region_total,msy_results$yield_total, type = 'p',ylab='Yield',xlab = "Harvest Rate", lwd = 2)

#Equilibruim Biomass vs Yield
plot(msy_results$biomass_total,msy_results$yield_total, type = 'p',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)

#Harvest rate vs biomass
plot(msy_results$biomass_total,msy_results$u_region_total, type = 'p',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)
}

#save the plot
setwd(wd_figs)
pdf("spatial_1_plots.pdf")
MSY_plots()
dev.off()






######################
# MODEL PLOTS - for fun

out=readList(paste0("Report",t,".rep",sep=""))
par_names<-names(out)
par_names


#create a seq of years = super boring!
#plot
years<-seq(1,200,1)

plot(years,out$yield_total)


##########################################################
# RUN THE FULL MODEL WITH NEW harvest rate
##########################################################

#pull the report file from the MSY run and make a couple plots


















