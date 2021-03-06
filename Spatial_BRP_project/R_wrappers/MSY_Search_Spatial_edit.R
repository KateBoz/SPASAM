###############################################################
# Code for conducting the MSY Search using the Spatial_BRP model for 
# multiple spatial structures
# Created by Katelyn Bosley
# Date: 5/8/2017
##############################################################

#remove junk
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
  #suppressWarnings(suppressMessages(library(spatstat)))
  #suppressWarnings(suppressMessages(library(alphahull)))
}
load_libraries()


#need to reset the working directory to the folder with the MSY__search files
####################################################################
#if not looping over the different folders
WD<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\MY DOCS\\SPASAM STUFF_MS1\\MS1_results\\HAKE\\BASE_Alt_Done\\MSY Results_grid"
setwd(WD)
WD<<-WD  

#to bypass MSY_search to edit
wd<-WD

###################################################################  
#Setting up the F values to iterate over for all the runs

#hake
F.name<-"input_F"
F.start<-0
F.end<-2.3
it<-0.025

#sablefish
#F.name<-"input_F"
#F.start<-0
#F.end<-2.3
#it<-0.025


# select the population type
# 1 - panmictic
# 2 - multiple area
# 3 - metapop
# 4 - natal homing
pop.type<-2


# select the Recruitment calculation type
# 1 - R_ave
# 2 - BH
# 3 - Cyclical
rec.type<-2

################################################################

#run the whole thing!
MSY_search<-function(wd=WD) {

#set pop.type in the function
pop.type = pop.type

#read in .dat file to get values for setting up the runs-carryover from DG code
update=readLines("Spatial_BRP.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("npopulations",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])

#need to adjust by number of populations
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nstocks))],ncol=nstocks))
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

#set up wd to report files
wd_results<-paste0(wd,"\\MSY Results\\Report Files",sep="")
wd_figs<-paste0(wd,"\\MSY Results\\Figures",sep="")

##############################################################################
#set up the .csv to hold all the runs

#panmictic
if(pop.type==1){
#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = ntrial,ncol=10))
names(msy_results)<-c("trial","F","biomass_total_start","biomass_total_end","yield_total","harvest_rate_total_bio","depletion_total","SSB_total_start","SSB_total_end","Bratio_total")
}


# multi-area
if(pop.type==2){

N_par_reg<-8 # number of parameters with regional values # need to fix this up..
N_par_pop<-9 # number of parameters for stock

#picking out only the values I want.
msy_results<-data.frame(matrix(NA,nrow = nrow(permutation),ncol=((N_par_reg*nregions)+N_par_pop))) # number of parameters with multiple regions+number of total pop values

#fill in the table slowly for error checking
names(msy_results)[1]<-c("perm")
#fill in the rest will loops so can be changed as needed 

for(i in 1:nregions) 
{
  names(msy_results)[1+i]<-paste0("F.",i)
  names(msy_results)[1+nregions+i]<-paste("biomass_start_region.",i,sep = "")
  names(msy_results)[2+nregions*2]<-"biomass_total_start"
  names(msy_results)[2+nregions*2+i]<-paste("biomass_end_region.",i,sep = "")
  names(msy_results)[3+nregions*3]<-"biomass_total_end"
  names(msy_results)[3+nregions*3+i]<-paste("yield_region.",i,sep = "")
  names(msy_results)[4+nregions*4]<-"yield_total"
  names(msy_results)[4+nregions*4+i]<-paste("u_region.",i,sep = "")
  names(msy_results)[5+nregions*5]<-"u_region_total"
  names(msy_results)[5+(nregions*5)+i]<-paste("depletion_region.",i,sep="")
  names(msy_results)[6+(nregions*6)]<-"depletion_total"
  names(msy_results)[6+(nregions*6)+i]<-paste("SSB_start_region.",i,sep="")
  names(msy_results)[7+nregions*7]<-"SSB_total_start"
  names(msy_results)[7+(nregions*7)+i]<-paste("SSB_end_region.",i,sep="")
  names(msy_results)[8+nregions*8]<-"SSB_total_end"
  names(msy_results)[9+nregions*8]<-"Bratio_total"
 }
}


#metapopulation
if(pop.type==3) {
  
  N_par_reg<-9 # number of parameters with regional values
  N_par_pop<-9 # number of parameters for totals
  
  #picking out only the values I want.
  msy_results<-data.frame(matrix(NA,nrow = nrow(permutation),ncol=((N_par_reg*nstocks)+N_par_pop))) # number of parameters with multiple regions+number of total pop values
  
  #fill in the table slowly for error checking
  names(msy_results)[1]<-c("perm")
  #fill in the rest will loops so can be changed as needed  - stocks instead of area
  
  for(i in 1:nstocks) 
  {
    names(msy_results)[1+i]<-paste0("F.",i)
    names(msy_results)[1+nstocks+i]<-paste("biomass_start_population.",i,sep = "")
    names(msy_results)[2+nstocks*2]<-"biomass_total_start"
    names(msy_results)[2+nstocks*2+i]<-paste("biomass_end_population.",i,sep = "")
    names(msy_results)[3+nstocks*3]<-"biomass_total_end"
    names(msy_results)[3+nstocks*3+i]<-paste("yield_population.",i,sep = "")
    names(msy_results)[4+nstocks*4]<-"yield_total"
    names(msy_results)[4+nstocks*4+i]<-paste("u_population.",i,sep = "")
    names(msy_results)[5+nstocks*5]<-"u_region_total"
    names(msy_results)[5+(nstocks*5)+i]<-paste("depletion_population.",i,sep="")
    names(msy_results)[6+(nstocks*6)]<-"depletion_total"
    names(msy_results)[6+(nstocks*6)+i]<-paste("SSB_start_population.",i,sep="")
    names(msy_results)[7+nstocks*7]<-"SSB_total_start"
    names(msy_results)[7+(nstocks*7)+i]<-paste("SSB_end_population.",i,sep="")
    names(msy_results)[8+nstocks*8]<-"SSB_total_end"
    names(msy_results)[8+(nstocks*8)+i]<-paste("Bratio_population.",i,sep="")
    names(msy_results)[9+nstocks*9]<-"Bratio_total"
  }
}

#names(msy_results)



# set up natal homing results here 



#set up parallel 
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#set up text progress bar
pb <- txtProgressBar(max = n_perm, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)



#########################################################################
# MSY Search in loops
#########################################################################

# simple search for panmictic population structure
if (pop.type==1){
#run  loop to get all the values into msy_results into a csv for plotting
  
  #do the parallel processing over the loops
stime <- system.time({
    #run parallel  
ls=foreach(i=1:ntrial,.options.snow = opts) %dopar% {

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
 
  #read the dat
  out=PBSmodelling::readList("Spatial_BRP.rep")
  
  #store results to a full spreadsheet
  temp<-c(i,F.test[i],out$biomass_total[1],
          out$biomass_total[nyrs],out$yield_total[nyrs],
          out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],
          out$SSB_total[1],out$SSB_total[nyrs], out$Bratio_total[nyrs])

  
  #clean non-needed files and move to results folder
  invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))
  
  return(temp)
  
} #end of code for MSY search


for(i in 1:ntrial){
  msy_results[i,]<-ls[[i]]
  #unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
}

  }) #end parallel
  
} 



# run parallel for other spatial structures
if(pop.type>1) {
  
#do the parallel processing over the loops
stime <- system.time({
  #run parallel  
  ls=foreach(i=1:n_perm,.options.snow = opts) %dopar% {
    
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
    
    #read the dat
    out=PBSmodelling::readList("Spatial_BRP.rep")
    
    #store results to a full spreadsheet add them in slowly for easy changes- check to see if it matches the .csv made above
  
    
  # multi-area  
    if(pop.type==2) {
    temp<-c(i,permutation[i,],
            out$biomass_AM[,1],
            out$biomass_total[1],
            out$biomass_AM[,nyrs],
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
    
  # meta population
    if(pop.type==3) {
      temp<-c(i,permutation[i,],
              out$biomass_population[,1],
              out$biomass_total[1],
              out$biomass_population[,nyrs],
              out$biomass_total[nyrs],
              out$yield_population[,nyrs],
              out$yield_total[nyrs],
              out$harvest_rate_population_bio[,nyrs],
              out$harvest_rate_total_bio[nyrs],
              out$depletion_population[,nyrs],
              out$depletion_total[nyrs],
              out$SSB_population[,1],
              out$SSB_total[1],
              out$SSB_population[,nyrs],
              out$SSB_total[nyrs], 
              out$Bratio_population[,nyrs],
              out$Bratio_total[nyrs])
      
    }
      
    
  # add natal homing here
    
    #clean non-needed files and move to results folder
    invisible(file.remove(paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.exe",sep="")))
    invisible(file.copy(from=paste0(wd,"\\MSY Results\\Run",i,"\\Spatial_BRP.rep",sep=""),to=paste0(wd,"\\MSY Results\\Report Files\\Report",i,".rep",sep="")))
    
    return(temp)
    
  } #end of code for MSY search
  
  
  for(i in 1:n_perm){
    msy_results[i,]<-ls[[i]]
    #unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
  }
  
 }) #end parallel

} # end if statement for pop.type>1


#stime
close(pb)
stopCluster(cl) #end the cluster for parallel processing

#save results
setwd(wd_figs)
write.csv(msy_results,"MSY_results.csv")

# Getting the MSY vals
t<-which(msy_results$yield_total==max(msy_results$yield_total))
msy_true<-msy_results[t,]

msy_true$SSB_ratio<-msy_true$SSB_total_start/msy_true$SSB_total_end

write.csv(msy_true,"MSY_true.csv")

#move over full report also for fun
invisible(file.copy(from=paste0(wd_results,"\\Report",t,".rep",sep=""),to=paste0(wd_figs,"\\Report",t,".rep",sep="")))


############################################
#plotting code
############################################

#simple plotting function
if(rec.type==1 || rec.type ==3){
MSY_plots<-function() {
  
  #add some more if statements to fix the R_ave and BH _ plots...
  
  if(pop.type==1){
    msy_results<-read.csv("msy_results.csv") #read in the data again if there were changes above
    
    par(mfrow = c(2,2))
    
    #SSB Ratio vs. Yield
    options(scipen = 50, digits = 4) # fix sci notation
    plot(msy_results$SSB_ratio,msy_results$yield_total, type = 'b',ylab='Yield',xlab = " Equilibruim SSB Ratio", lwd = 2)
    
    #Harvest Rate vs Yield
    plot(msy_results$harvest_rate_total_bio,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Harvest Rate", lwd = 2)
    
    #Equilibruim Biomass vs Yield
    plot(msy_results$biomass_total_end,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)
    
    #Harvest rate vs biomass
    plot(msy_results$biomass_total_end,msy_results$harvest_rate_total_bio, type = 'b',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)
  }
  
  
  if(pop.type>1) {
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
  plot(msy_results$u_region_total, msy_results$biomass_total_end,type = 'p',xlab='Harvest Rate',ylab = "Equilibrium Biomass",     lwd = 2)
  }
  
 } 
}


if(rec.type==2) {
  MSY_plots<-function() {
    
    #add some more if statements to fix the R_ave and BH _ plots...
    
    if(pop.type==1){
      msy_results<-read.csv("msy_results.csv") #read in the data again if there were changes above
      
      par(mfrow = c(2,2))
      
      #SSB Ratio vs. Yield
      options(scipen = 50, digits = 4) # fix sci notation
      plot(msy_results$Bratio,msy_results$yield_total, type = 'b',ylab='Yield',xlab = " Equilibruim SSB Ratio", lwd = 2)
      
      #Harvest Rate vs Yield
      plot(msy_results$harvest_rate_total_bio,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Harvest Rate", lwd = 2)
      
      #Equilibruim Biomass vs Yield
      plot(msy_results$biomass_total_end,msy_results$yield_total, type = 'b',ylab='Yield',xlab = "Equilibrium Biomass", lwd = 2)
      
      #Harvest rate vs biomass
      plot(msy_results$biomass_total_end,msy_results$harvest_rate_total_bio, type = 'b',ylab='Harvest Rate',xlab = "Equilibrium Biomass", lwd = 2)
    }
    
    
    if(pop.type>1) {
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
      plot(msy_results$u_region_total, msy_results$biomass_total_end,type = 'p',xlab='Harvest Rate',ylab = "Equilibrium Biomass",     lwd = 2)
    }
  } 
}


#save the plot
setwd(wd_figs)
pdf("spatial_1_plots.pdf")
MSY_plots()
dev.off()

#remove the old files from run
if(pop.type==1){
  for(i in 1:ntrial){
    #delete folders
    unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
    unlink(paste0(wd,"\\MSY Results\\Report Files",sep = ""),recursive = T)
  }}

if(pop.type>1){
for(i in 1:n_perm){
  #delete folders
  unlink(paste0(wd,"\\MSY Results\\","Run",i,sep = ""),recursive = T)
  unlink(paste0(wd,"\\MSY Results\\Report Files",sep = ""),recursive = T)
}}

} #end MSY_search function

MSY_search()

#} #end of full folder loops

