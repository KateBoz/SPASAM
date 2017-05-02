###############################################################
# Conducting the MSY Search using the Spatial_BRP model for 
# the "one population/multiple areas" spatial structure
# Created by Katelyn Bosley
# Date: 3/27/2017
##############################################################

#need to reset the working directory to the folder with the MSY__search files

# remove previous objects from workspace
rm(list = ls())


# load libraries function
load_libraries<-function() {
  suppressWarnings(suppressMessages(library(PBSmodelling)))
  suppressWarnings(suppressMessages(library(matrixStats)))
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


#Setting up the F values to iterate over for all the runs
F.name<-"input_F"
F.start<-0.005
F.end<-5
it<-0.25


############################################################################
# set the working directory to location where the folders are to iterate through

#HAKE runs
#folder<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SABLEFISH\\test"

#create a list of working directorys to iterate over
#runs<-vector("list",1)

#for(i in 1:length(runs)) {
#  runs[[i]]<-paste0(folder,"\\",i,sep="")
#}


#loop over the folders for each run
#for(i in 1:length(runs)) {
#   
#   WD<-runs[[i]] #setting the new WD
#   setwd(WD)
#   WD<<-WD
   
   
   
####################################################################
#if not looping over the different folders
WD<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\HAKE\\BASE_main"
setwd(WD)
WD<<-WD  
###################################################################   
  
#to bypass MSY_search for edting
#wd<-WD

# the search function  - run the below function to make things easier for completing runs. Still working out the kinks
   
MSY_search<-function(wd=WD) {

#read in .dat file to get values for setting up the runs-carryover from DG code
update=readLines("Spatial_BRP.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("npopulations",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])


#need to adjust by number of populations
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nstocks))],ncol=nstocks))

#
#nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))


#set up the combinations of permutations per region
F.test<-seq(F.start,F.end,it) #F values to cycle through
ntrial<-length(F.test) # count of total number of trials

#set up the permutations of F by pop/fleet/region
#had to make adjustment to nregions... not fleets
permutation<-permutations(ntrial,sum(nregions),F.test,repeats.allowed=TRUE)  # determine all permutations of F in each stock
n_perm<-nrow(permutation)

#setting up the files to put results
dir.create(paste0(wd,"\\MSY Results",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Figures",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Report Files",sep=""))

#set up the parallel processing
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)
  
#set up text progress bar
pb <- txtProgressBar(max = n_perm, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#actual MSY_Search function
#do the parallel processing over the loops
  
#test if parallel is working
stime <- system.time({
  
#run parallel  
ls<-foreach(i=1:n_perm,.options.snow = opts) %dopar% {
  dir.create(paste0(wd,"\\MSY Results\\Run",i,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.exe",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.dat",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.tpl",sep=""),to=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.tpl",sep="")))
  
  #set new directory for running the model
  setwd(paste0(wd,"\\MSY Results\\Run",i,sep="")) # now set the working directory as the run# file
  
  update=readLines("Spatial_BRP.dat",n=-1)
  
  #pull important values - for generalized version 
  update[(grep("input_F",update)+1):(grep("input_F",update)+sum(nregions))]=permutation[i,]
  
  writeLines(update,"Spatial_BRP.dat")
  
  ### Run ADMB with updated F
  invisible(shell("Spatial_BRP -nohess",wait=T))
  
  #clean non-needed files and move to results folder
  invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))

} #end of code for MSY search

})

stime

stopCluster(cl) #end the cluster for parallel processing


#########################################################
#setting up and grabbing the values from the runs to plot
#########################################################

#need to generalize this port for plotting results
#par_names<-names(out)

#set up wd to report files
wd_results<-paste0(wd,"\\MSY Results\\Report Files",sep="")
wd_figs<-paste0(wd,"\\MSY Results\\Figures",sep="")


N_par_reg<-6 # number of parameters with regional values # need to fix this up..
N_par_pop<-9 # number of parameters for stock

#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = nrow(permutation),ncol=((N_par_reg*nregions)+N_par_pop))) # number of parameters with multiple regions+number of total pop values
 
#fill in the table slowly for error checking
names(msy_results)[1]<-c("perm")


#fill in the rest will loops so can be changed as needed 

for(i in 1:nregions) 
  {
  names(msy_results)[1+i]<-paste0("F.",i)
  names(msy_results)[2+nregions]<-"biomass_total_start"
  names(msy_results)[3+nregions]<-"biomass_total_end"
  names(msy_results)[3+nregions+i]<-paste("yield_region.",i,sep = "")
  names(msy_results)[4+nregions*2]<-"yield_total"
  names(msy_results)[4+nregions*2+i]<-paste("u_region.",i,sep = "")
  names(msy_results)[5+nregions*3]<-"u_region_total"
  names(msy_results)[5+(nregions*3)+i]<-paste("depletion_region.",i,sep="")
  names(msy_results)[6+(nregions*4)]<-"depletion_total"
  names(msy_results)[6+(nregions*4)+i]<-paste("SSB_start_region.",i,sep="")
  names(msy_results)[7+nregions*5]<-"SSB_total_start"
  names(msy_results)[7+(nregions*5)+i]<-paste("SSB_end_region.",i,sep="")
  names(msy_results)[8+nregions*6]<-"SSB_total_end"
  names(msy_results)[9+nregions*6]<-"Bratio_total"
}

#names(msy_results)

#set to results for grabbing values
setwd(wd_results)

#set up parellel
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#setting up another progress bar
pb <- txtProgressBar(max = n_perm, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)


stime2 <- system.time({
#run loop to get all the values into msy_results into a csv for plotting
grab<-foreach(i=1:n_perm,.options.snow = opts) %dopar% { 
  
#read report 
out=PBSmodelling::readList(paste0("Report",i,".rep",sep=""))

#store results to a full spreadsheet add them in slowly for easy changes- check to see if it matches the .csv made above

temp<-c(i,permutation[i,],
        out$biomass_total[1],
        out$biomass_total[nyrs],
        out$yield_region[,nyrs],
        out$yield_total[nyrs],
        out$harvest_rate_region_bio[,nyrs],
        out$harvest_rate_total_bio[nyrs],
        out$depletion_region[,nyrs],
        out$depletion_total[nyrs],
        out$SSB_region[,1],
        out$SSB_total[1],
        out$SSB_region[,nyrs],
        out$SSB_total[nyrs], 
        out$Bratio_total[nyrs])
}


for(i in 1:n_perm){
msy_results[i,]<-grab[[i]]

}
}) #end value grab

close(pb)
stopCluster(cl) #end the cluster for parallel processing

stime2 

#results
write.csv(msy_results,"MSY_results.csv")


# Getting the MSY vals
t<-which(msy_results$yield_total==max(msy_results$yield_total))
msy_true<-msy_results[t,]
write.csv(msy_true,"MSY_true.csv")


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
plot(msy_results$biomass_total_end,msy_results$yield_total, type = 'p',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)

#Harvest rate vs biomass
plot(msy_results$u_region_total, msy_results$biomass_total_end,type = 'p',xlab='Harvest Rate',ylab = "Equilibrium Biomass", lwd = 2)

}

#save the plot
setwd(wd_figs)
pdf("spatial_1_plots.pdf")
MSY_plots()
dev.off()

#remove the old files from run
for(i in 1:n_perm){
  unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
  unlink(paste0(wd,"\\MSY Results\\Report Files",sep = ""),recursive = T)
  }

} #end MSY_search function

MSY_search()

#} #end of full folder loops






