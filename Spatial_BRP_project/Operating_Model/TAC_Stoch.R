###########################################################
# Stochastic runs for TAC 
# Created by: Katelyn Bosley
# Date: 7/17/2017
###########################################################

#remove junk
rm(list = ls())

#debug(utils:::unpackPkgZip)
#install.packages("PBSmodelling")

#load libraries
# load libraries function
load_libraries<-function() {
  suppressWarnings(suppressMessages(library(PBSmodelling)))
  suppressWarnings(suppressMessages(library(matrixStats)))
  suppressWarnings(suppressMessages(require(snowfall)))
  suppressWarnings(suppressMessages(library(parallel)))
  suppressWarnings(suppressMessages(library(snow)))
  suppressWarnings(suppressMessages(library(foreach)))
  suppressWarnings(suppressMessages(library(doSNOW)))
}
load_libraries()



directory<-"C:\\Users\\katelyn.bosley\\Desktop\\_TAC_temp_runs\\stoch_runs_menhaden"
folders<-list.files(path = directory)


#Looping over the folders

for (j in 1:length(folders)){
  
#setwd
wd<-paste0(directory,"\\",folders[j],sep = "")
setwd(wd)


##############################
#For saving the results
#
# select the population type
# 1 - panmictic
# 2 - multiple area
# 3 - metapop
# 4 - natal homing
pop.type<-2

################################
# Load info for making the data frame

#set the number of runs
nr<-100
run.index<-seq(1:nr)

#which stochastic element of TAC allocation
# 1 - Survey Biomass
# 2 - Rec Index
allo.type<-1

#################################

################################
# Load get info from .dat for saving values

{ #run the whole code
  
dat<-readLines("Pop_TAC.dat")
#msy<-readList("Metapop.rep") # this is here so I can grab the meta pop values
nregions<-as.numeric(dat[grep("nregions",dat)+1])
npops<-as.numeric(dat[grep("npopulations",dat)+1])
nyrs<-as.numeric(dat[grep("nyrs",dat)+1])
nages<-as.numeric(dat[grep("nages",dat)+1])
#  

#build data frame for saving values


  #panmictic
  if(pop.type==1){
    #picking msy only the values I want.
    msy_results<-data.frame(matrix(NA,nrow = nr,ncol=10))
    names(msy_results)<-c("Model","F","biomass_total_start","biomass_total_end",
                          "yield_total","harvest_rate_total_bio","depletion_total",
                          "SSB_total_start","SSB_total_end","Bratio_total")
  }
  
  
  # if using the multiple area model
  if(pop.type==2){
    
    N_par_reg<-6 # number of parameters with regional values # need to fix this up..
    N_par_pop<-9 # number of parameters for stock
    
    #picking msy only the values I want.
    msy_results<-data.frame(matrix(NA,nrow = nr,ncol=(N_par_reg*nregions)+N_par_pop)) # number of parameters with multiple regions+number of total pop values
    
    #fill in the table slowly for error checking
    names(msy_results)[1]<-c("Model")
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
  }
  
  #metapopulation
  if(pop.type==3) {
    
    N_par_reg<-9 # number of parameters with regional values
    N_par_pop<-9 # number of parameters for totals
    
    #picking msy only the values I want.
    msy_results<-data.frame(matrix(NA,nrow = nr,ncol=((N_par_reg*npops)+N_par_pop))) # number of parameters with multiple regions+number of total pop values
    
    #fill in the table slowly for error checking
    names(msy_results)[1]<-c("Model")
    #fill in the rest will loops so can be changed as needed  - stocks instead of area
    
    for(i in 1:npops) 
    {
      names(msy_results)[1+i]<-paste0("F.",i)
      names(msy_results)[1+npops+i]<-paste("biomass_start_population.",i,sep = "")
      names(msy_results)[2+npops*2]<-"biomass_total_start"
      names(msy_results)[2+npops*2+i]<-paste("biomass_end_population.",i,sep = "")
      names(msy_results)[3+npops*3]<-"biomass_total_end"
      names(msy_results)[3+npops*3+i]<-paste("yield_population.",i,sep = "")
      names(msy_results)[4+npops*4]<-"yield_total"
      names(msy_results)[4+npops*4+i]<-paste("u_population.",i,sep = "")
      names(msy_results)[5+npops*5]<-"u_region_total"
      names(msy_results)[5+(npops*5)+i]<-paste("depletion_population.",i,sep="")
      names(msy_results)[6+(npops*6)]<-"depletion_total"
      names(msy_results)[6+(npops*6)+i]<-paste("SSB_start_population.",i,sep="")
      names(msy_results)[7+npops*7]<-"SSB_total_start"
      names(msy_results)[7+(npops*7)+i]<-paste("SSB_end_population.",i,sep="")
      names(msy_results)[8+npops*8]<-"SSB_total_end"
      names(msy_results)[8+(npops*8)+i]<-paste("Bratio_population.",i,sep="")
      names(msy_results)[9+npops*9]<-"Bratio_total"
    }
  }
  

#msy_results


##################################################
# run the model and save the results
##################################################

#set up temp file for simulations
dir.create(paste0(wd,"\\sim",sep="")) #create the results directory in the WD with run number

#setting up parallel
no_cores <- detectCores()
cl<-makeCluster(no_cores)
registerDoSNOW(cl)

#set up text progress bar
pb <- txtProgressBar(max = nr, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

#run the stochastic simulation

loop.ls=foreach(i=1:nr,.options.snow = opts) %dopar% {
  
  dir.create(paste0(wd,"\\sim\\Run",i,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(wd,"\\Pop_TAC.exe",sep=""),to=paste0(wd,"\\sim\\Run",i,"\\Pop_TAC.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Pop_TAC.dat",sep=""),to=paste0(wd,"\\sim\\Run",i,"\\Pop_TAC.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"\\Pop_TAC.tpl",sep=""),to=paste0(wd,"\\sim\\Run",i,"\\Pop_TAC.tpl",sep="")))
  
  #set new directory for running the model
  setwd(paste0(wd,"\\sim\\Run",i,sep="")) # now set the working directory as the run# file
  

#update myseed for random run
  
#for survey
if(allo.type == 1) { 
  
new.rand<-readLines("Pop_TAC.dat", n=-1)
new.rand[(grep("myseed_survey_rand",new.rand)+1)]<-run.index[i]
writeLines(new.rand, "Pop_TAC.dat")
}

#for rec
if(allo.type == 2) {
  new.rand<-readLines("Pop_TAC.dat", n=-1)
  new.rand[(grep("myseed_rec_index",new.rand)+1)]<-run.index[i]+test
  writeLines(new.rand, "Pop_TAC.dat")
}

#run the model

invisible(shell("Pop_TAC -nohess",wait=T))

out=PBSmodelling::readList("Pop_TAC.rep")


#if there are only 2 regions

  if(nregions==2){
    
          #panmictic
      if(pop.type==1){
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),out$biomass_total[1],
                out$biomass_total[nyrs],out$yield_total[nyrs],
                out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],
                out$SSB_total[1],out$SSB_total[nyrs], out$Bratio_total[nyrs])
      }
      
      
      # multi-area  
      if(pop.type==2) {
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),
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
      
      
      # meta population
      if(pop.type==3) {
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),
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
      
      return(temp)
      #msy_results[i,]<-as.numeric(temp)
      
    }
    

  
#for the 3 area models
  
  if(nregions==3){
    
      #panmictic
      if(pop.type==1){
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),max(out$F[nyrs*3,]),
                out$F[nyrs,nages],out$F[nyrs*2,nages],out$F[nyrs*3,nages],out$biomass_total[1],
                out$biomass_total[nyrs],out$yield_total[nyrs],
                out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],
                out$SSB_total[1],out$SSB_total[nyrs], out$Bratio_total[nyrs])
      }
      
      
      # multi-area  
      if(pop.type==2) {
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),max(out$F[nyrs*3,]),
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
      
      
      # meta population
      if(pop.type==3) {
        temp<-c(i,max(out$F[nyrs,]),max(out$F[nyrs*2,]),max(out$F[nyrs*3,]),
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
      
      return(temp)    
      #msy_results[i,]<-as.numeric(temp)
      
  }

} #end the parallel

close(pb)
stopCluster(cl) #end the cluster for parallel processing

###########################################
#put the results in a data frame and remove sim files 
for(i in 1:nr){
  msy_results[i,]<-loop.ls[[i]]
  unlink(paste0(wd,"\\sim\\","Run",i,sep = ""),recursive = T)
}

unlink(paste0(wd,"\\sim",sep = ""),recursive = T)

################################################
#save the results as a csv
setwd(wd)
write.csv(msy_results,"TAC_allocation_stoch.csv")

} #end of whole deal


} #loop to the next folder


################################################
################################################

# quick comparision of other error level

#load in the other other file
one<-read.csv("C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SPASAM STUFF_MS1\\MS1_results\\HAKE\\PHASE 2 RUNS\\Stoch_runs\\4_pop_uMSY_survey\\TAC_allocation_stoch.csv")
head(one)


#looking at the plots
par(mfrow=c(2,1))

hist(one$yield_total)
hist(msy_results$yield_total)

hist(one$depletion_total)
hist(msy_results$depletion_total)



