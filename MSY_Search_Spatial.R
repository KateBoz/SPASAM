######################################################
# Another attempt to run the MSY search on the panmictic population
# Created by Katelyn Bosley
# Date: 2/15/2017
###################################################
#need to reset the working directory in for run and run the run_MSY function


# remove previous objects from workspace
rm(list = ls())

# set the working directory

#DIRECTORIES

#HAKE runs

#SABLEFISH runs
WD<-"G:\\SPASAM CODING\\MS_1_CODE\\Sablefish\\Base_model"

#MENHADEN runs

#
setwd(WD)
WD<<-WD

#install libraries - don't need all these but carry over
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
F.start<-0.1
F.end<-0.8
it<-0.1

# the function  - run the below function to make things easier for completing runs. Still working out the kinks
MSY_search<-function(wd=WD) {

#read in .dat file to get values for setting up the runs-carryover from DG code
update=readLines("Spatial_BRP.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("npopulations",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))


#set up the combinations of permutations per region
F.test<-seq(F.start,F.end,it) #F values to cycle through
ntrial<-length(F.test) # count of total number of trials

#set up the permutations of F by pop/fleet/region
permutation<-permutations(ntrial,sum(nfleets),F.test,repeats.allowed=TRUE)  # determine all permutations of F in each stock
n_perm<-nrow(permutation)

#setting up the files to put results
dir.create(paste0(wd,"\\MSY Results",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Figures",sep=""))
dir.create(paste0(wd,"\\MSY Results\\Report Files",sep=""))

#set up the parallel processing
no_cores <- detectCores() - 1
cl<-makeCluster(no_cores)
registerDoSNOW(cl)
  
#set up text progress bar
pb <- txtProgressBar(max = n_perm, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#do the parallel processing over the loops
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

stopCluster(cl) #end the cluster for parallel processing



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
setWinProgressBar(pb,(i/n_perm*100),label=paste("Grabbing values from ", i," of ", n_perm," simulations",sep=""))
  
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

close(pb)


# Getting the MSY vals
t<-which(msy_results$yield_total==max(msy_results$yield_total))
msy_pan<-msy_results[t,]
write.csv(msy_pan,"MSY_true.csv")



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
plot(msy_results$biomass_total,msy_results$yield_total, type = 'p',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)

#Harvest rate vs biomass
plot(msy_results$biomass_total,msy_results$u_region_total, type = 'p',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)

}

#save the plot
setwd(wd_figs)
pdf("spatial_1_plots.pdf")
MSY_plots()
dev.off()

} #end MSY_search function

MSY_search()






########## FOR THE FUTURE!!######################




#################################################################
# Setting up the function to run the program with u-msy or TAC msy
#################################################################

#a work in progress! - can set up as a function also

setwd(wd)

#index some other folders
wd_results<-paste0(wd,"\\MSY Results\\Report Files",sep="")
wd_figs<-paste0(wd,"\\MSY Results\\Figures",sep="")


#create new directory for real run
dir.create(paste0(wd,"\\ApplyTAC",sep=""))
dir.create(paste0(wd,"\\ApplyTAC\\TACresults",sep=""))

#move over the files for running the 
invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.exe",sep=""),to=paste0(wd,"\\ApplyTAC\\TACresults\\Spatial_BRP.exe",sep="")))
invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.dat",sep=""),to=paste0(wd,"\\ApplyTAC\\TACresults\\Spatial_BRP.dat",sep="")))
invisible(file.copy(from=paste0(wd,"\\Spatial_BRP.tpl",sep=""),to=paste0(wd,"\\ApplyTAC\\TACresults\\Spatial_BRP.tpl",sep="")))

#
#set new wd
wd_TACresults<-paste0(wd,"\\ApplyTAC\\TACresults",sep="")
setwd(wd_TACresults)


#pull in TAC and u values for running model
#pull in MSY_true from MSY search run
MSY_rep<-read.csv(paste0(wd_figs,"\\MSY_true.csv",sep=""))

#pull pop_values
update=readLines("Spatial_BRP.dat",n = -1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("nstocks",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))

#full in TAC values
TAC_cols<-rep(NA,nregions)
TAC<-rep(NA,nregions)

u_cols<-rep(NA,nregions)
u<-rep(NA,nregions)

for (i in 1:nregions){

 TAC_cols[i]<-paste0("yield_region.",i,sep = "")
 u_cols[i]<-paste0("u_region.",i,sep = "")
}

#filling values
TAC<-MSY_rep[TAC_cols]
u<-MSY_rep[u_cols]


#set values for changing switches
model_switch<-2 #1 for TAC 2 for umsy

#updating values in .dat with TAC and model switch
update[(grep("model_type_swtich",update)+4)]= as.character(model_switch) # there is a misspelling here I need to fix...
update[(grep("input_TAC",update)+1):(grep("input_TAC",update)+sum(nregions))]=as.character(TAC)
update[(grep("input_u",update)+1):(grep("input_u",update)+sum(nregions))]=as.character(u)

writeLines(update,"Spatial_BRP.dat")

###############################################################

#run the new model

invisible(shell("Spatial_BRP -nohess",wait=T)) #show.output.on.console=FALSE))  



##############################################
# Make some figures
#############################################


#create a directory for making figures

dir.create(paste0(wd_TACresults,"\\TAC_figs",sep=""))
invisible(file.copy(from=paste0(wd_TACresults,"\\Spatial_BRP.rep",sep=""),to=paste0(wd_TACresults,"\\TAC_figs\\Spatial_BRP.rep",sep="")))

TAC_figs<-paste0(wd_TACresults,"\\TAC_figs",sep = "")

setwd(TAC_figs)


######################
# MODEL PLOTS - for fun

out=readList("Spatial_BRP.rep")
par_names<-names(out)
par_names


#create a seq of years = super boring!
#plot
years<-seq(1,200,1)

#ssb
plot(years,out$SSB_total, ylim = c(0,300), ylab = "SSB MSY", type = 'l',lwd = 2)
points(years,out$SSB_region[1,], type = "l", col = "red", lwd = 2)
points(years,out$SSB_region[2,], type = "l", col = "blue", lwd = 2)
points(years,out$SSB_region[3,], type = "l", col = "green", lwd = 2)
text(0,300,"Match Sel, Match M, vary by area", pos = 4)


#depletion
plot(years,out$depletion_total, ylim = c(0,1.2), ylab = "SSB MSY", type = 'l',lwd = 2)
points(years,out$depletion_region[1,], type = "l", col = "red", lwd = 2)
points(years,out$depletion_region[2,], type = "l", col = "blue", lwd = 2)
points(years,out$depletion_region[3,], type = "l", col = "green", lwd = 2)
text(150,1.2,"Match Sel, Match M, vary by area", pos = 1)

###################################################################
# There is much more to do!
####################################################################


















