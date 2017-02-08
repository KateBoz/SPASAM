######################################################################
# Script to iterate over F values by stock/region/fleet in search for FMSY
# Created by Dan Goethel
# Edited by Katelyn Bosley
# Date: 
#
######################################################################

# remove previous objects from workspace
rm(list = ls())


# set the working directory
wd<-"C:\\Users\\katelyn.bosley\\Desktop\\NOAA_POSTDOC\\SPASM\\Potentially Unuseful Code\\Spatial BRPs\\Most Recent Code"

setwd(wd)
wd<<-wd
                                                   
###############################
# Step through the MSY search
###############################

# optimize MSY by fleets within areas (=TRUE) or simply by using an F multiplier within an area (=FALSE). The latter reduces dimensionality of F permuations when number of fleets is too large to do full analysis
optimize.by.fleet<-'FALSE'             

# F values from .dat file to replace when iterating over F
F.name<-"input_F"                
F.end=.5        # Terminal F value to iterate over
interval=0.065  # Step increase in F 

#read the .dat file, set it as an object (character)
update=readLines("Spatial_BRP.dat",n=-1)

#pull population characteristics from the .dat file 
nyrs<-as.numeric(update[(grep("nyrs",update)+1)])    # number of years for the simulation
nstocks<-as.numeric(update[(grep("nstocks",update)+1)])    # number of stocks simulated
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])    # vector of the number of regions within stocks
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions)) # matrix of fleets per region per stock ####NEEDS TO BE GENERALIZED FOR UNEVEN REGIONS BY STOCK (i.e. ncol)



# LOAD LIBRARIES
#--------------------------------------------------------
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
#---------------------------------------------------------


##########################################
# Prepare for F iterations
##########################################

# Create a sequence of desired F values to iterate over

F.seq<-seq(0,F.end,interval)                                                      
if(optimize.by.fleet=='TRUE')
{
  permutation<-permutations(length(F.seq),sum(nfleets),F.seq,repeats.allowed=TRUE)  # determine all permutations of F in each stock/region/fleet
}else {
  if(optimize.by.fleet=='FALSE')
  {
    permutation<-permutations(length(F.seq),sum(nregions),F.seq,repeats.allowed=TRUE)  # determine all permutations of F in each stock/region/fleet
  }
}



# set up directories to store all the report files for each permutation
ntrials<-length(permutation[,1])
seq.length<-length(F.seq)
dir.create(paste0(wd,"/MSY Results",sep=""))
dir.create(paste0(wd,"/MSY Results/Figures",sep=""))
dir.create(paste0(wd,"/MSY Results/Report Files",sep=""))


# Function for calculating an rounded-up integer from inputs
Roundup <- function(from,to) ceiling(from/to)*to



# Function for running the iterations
MSY<-function(WD,ntrial,F.grab,perm) # params are WD, trial#, F combination and permutation number
 {
  setwd(WD) #set the WD

  dir.create(paste0(WD,"/MSY Results/Run",ntrial,sep="")) #create the results directory in the WD with run number
  
  invisible(file.copy(from=paste0(WD,"/Spatial_BRP.exe",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/Spatial_BRP.exe",sep=""))) # add the files to each run folder
  invisible(file.copy(from=paste0(WD,"/Spatial_BRP.dat",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/Spatial_BRP.dat",sep="")))
  invisible(file.copy(from=paste0(WD,"/Spatial_BRP.tpl",sep=""),to=paste0(WD,"/MSY Results/Run",ntrial,"/Spatial_BRP.tpl",sep="")))

  setwd(paste0(WD,"/MSY Results/Run",ntrial,sep="")) # now set the working direcctory as the run# file

  update=readLines("Spatial_BRP.dat",n=-1) #read in the .dat from the run# file
  
  nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
  nstocks<-as.numeric(update[(grep("nstocks",update)+1)])
  nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])
  nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nregions))],ncol=nregions))  ####NEEDS TO BE GENERALIZED FOR UNEVEN REGIONS BY STOCK (i.e. ncol)
  
  if(optimize.by.fleet=='TRUE') # set up the F iterations
  {
    F.region1<-matrix(perm[ntrial,],ncol=nfleets,nrow=sum(nregions),byrow=T)  ####NEEDS TO BE GENERALIZED FOR UNEVEN FLEETS BY REGION (i.e. ncol)
    F.region<-apply(matrix(as.character(perm[ntrial,]),ncol=nfleets,nrow=sum(nregions),byrow=T),1,paste0,sep=" ",collapse="")  ####NEEDS TO BE GENERALIZED FOR UNEVEN FLEETS BY REGION (i.e. ncol)
  }else{
    if(optimize.by.fleet=='FALSE')
    {
      F.region1<-matrix(rep(perm[ntrial,],each=nfleets[1]),ncol=nfleets,nrow=sum(nregions),byrow=T)  ####NEEDS TO BE GENERALIZED FOR UNEVEN FLEETS BY REGION (i.e. ncol)
      F.region<-apply(matrix(as.character(rep(perm[ntrial,],each=nfleets[1])),ncol=nfleets,nrow=sum(nregions),byrow=T),1,paste0,sep=" ",collapse="")  ####NEEDS TO BE GENERALIZED FOR UNEVEN FLEETS BY REGION (i.e. ncol)
      }
  }

  # set the new F vals into the .dat file
  update[(grep(F.grab,update)+1):(grep(F.grab,update)+sum(nregions))]=F.region
  writeLines(update,"Spatial_BRP.dat")
  
  ### Run ADMB with update F
  invisible(shell("Spatial_BRP -nohess",wait=T)) #show.output.on.console=FALSE))  
  
  # remove the exe from the run# folder
  invisible(file.remove(paste0(WD,"/MSY Results/Run",ntrial,"/Spatial_BRP.exe",sep="")))
  
  #copy the report file from run# and report files directory
  invisible(file.copy(from=paste0(WD,"/MSY Results/Run",ntrial,"/Spatial_BRP.rep",sep=""),to=paste0(WD,"/MSY Results/Report Files/Report",ntrial,".rep",sep="")))
  
  #read in the report file as a list
  out=readList("Spatial_BRP.rep")
  
  # list the parameter names in the .rep file
  par_names=c('overlap_switch','larval_move_switch','move_switch','F','biomass_AM','biomass_stock',
              'biomass_total','yield_fleet','yield_region','yield_stock',
              'yield_total', 'harvest_rate_region_bio','harvest_rate_stock_bio',
              'harvest_rate_total_bio','depletion_region','depletion_stock','depletion_total',
              'SSB_region','SSB_stock','SSB_total','Bratio_stock','Bratio_total','biomass_AM_overlap_region',
              'biomass_stock_overlap','biomass_natal_overlap', 'yield_region_overlap','yield_stock_overlap',
              'yield_natal_overlap','harvest_rate_region_bio_overlap','harvest_rate_stock_bio_overlap',
              'harvest_rate_natal_bio_overlap','depletion_region_overlap','depletion_stock_overlap','depletion_natal_overlap',
              'SSB_region_overlap','SSB_stock_overlap','SSB_natal_overlap','Bratio_stock_overlap','Bratio_natal_overlap')
  
  #pull out the results we are after
  result=out[par_names]
  
  #set each item as an object
  bio_region<-out$biomass_AM
  bio_stock<-out$biomass_stock
  bio_total<-out$biomass_total
  #yield_fleet<-out$yield_fleet
  yield_region<-out$yield_region
  yield_stock<-out$yield_stock
  yield_total<-out$yield_total
  u_region<-out$harvest_rate_region_bio
  u_stock<-out$harvest_rate_stock_bio
  u_total<-out$harvest_rate_total_bio
  dep_region<-out$depletion_region
  dep_stock<-out$depletion_stock
  dep_total<-out$depletion_total
  SSB_region<-out$SSB_region
  SSB_stock<-out$SSB_stock
  SSB_total<-out$SSB_total
  Bratio_stock<-out$Bratio_stock
  Bratio_total<-out$Bratio_total
  bio_region_overlap<-out$biomass_AM_overlap_region
  bio_stock_overlap<-out$biomass_stock_overlap
  bio_total_overlap<-out$biomass_natal_overlap
  yield_region_overlap<-out$yield_region_overlap
  yield_stock_overlap<-out$yield_stock_overlap
  yield_total_overlap<-out$yield_natal_overlap
  u_region_overlap<-out$harvest_rate_region_bio_overlap
  u_stock_overlap<-out$harvest_rate_stock_bio_overlap
  u_total_overlap<-out$harvest_rate_natal_bio_overlap
  dep_region_overlap<-out$depletion_region_overlap
  dep_stock_overlap<-out$depletion_stock_overlap
  dep_total_overlap<-out$depletion_natal_overlap
  SSB_region_overlap<-out$SSB_region_overlap
  SSB_stock_overlap<-out$SSB_stock_overlap
  SSB_total_overlap<-out$SSB_natal_overlap 
  Bratio_stock_overlap<-out$Bratio_stock_overlap
  Bratio_natal_overlap<-out$Bratio_natal_overlap
  
  
  
  # PUll out the results for 1 region model
  if(sum(nregions)==1)
  { 
    c(i,F.region1,bio_region[1],bio_region[nyrs],bio_stock[1],bio_stock[nyrs],bio_total[1],bio_total[nyrs],
      yield_region[1],yield_region[nyrs],yield_stock[1],yield_stock[nyrs],yield_total[1],yield_total[nyrs],
      u_region[1],u_region[nyrs],u_stock[1],u_stock[nyrs],u_total[1],u_total[nyrs],
      dep_region[1],dep_region[nyrs],dep_stock[1],dep_stock[nyrs],dep_total[1],dep_total[nyrs],
      SSB_region[1],SSB_region[nyrs],SSB_stock[1],SSB_stock[nyrs],SSB_total[1],SSB_total[nyrs],
      Bratio_stock[1],Bratio_stock[nyrs],Bratio_total[1],Bratio_total[nyrs],
      bio_region_overlap[1],bio_region_overlap[nyrs],bio_stock_overlap[1],bio_stock_overlap[nyrs],
      bio_total_overlap[1],bio_total_overlap[nyrs], yield_region_overlap[1],yield_region_overlap[nyrs],yield_stock_overlap[1],
      yield_stock_overlap[nyrs],yield_total_overlap[1],yield_total_overlap[nyrs],
      u_region_overlap[1],u_region_overlap[nyrs],u_stock_overlap[1],u_stock_overlap[nyrs],u_total_overlap[1],
      u_total_overlap[nyrs], dep_region_overlap[1],dep_region_overlap[nyrs],dep_stock_overlap[1],dep_stock_overlap[nyrs],
      dep_total_overlap[1],dep_total_overlap[nyrs], SSB_region_overlap[1],SSB_region_overlap[nyrs],SSB_stock_overlap[1],
      SSB_stock_overlap[nyrs],SSB_total_overlap[1],SSB_total_overlap[nyrs],Bratio_stock_overlap[1],
      Bratio_stock_overlap[nyrs],Bratio_natal_overlap[1],Bratio_natal_overlap[nyrs]
    )
  }
  else # pull out results when the there is one stock with more than one region - for paper 2
    {
  if(sum(nstocks)<2  & sum(nregions)>1)
   {
    c(i,rbind(F.region1[,]),rbind(bio_region[,1]),rbind(bio_region[,nyrs]),rbind(bio_stock[1]),rbind(bio_stock[nyrs]),
      bio_total[1],bio_total[nyrs], rbind(yield_region[,1]),rbind(yield_region[,nyrs]),
      rbind(yield_stock[1]),rbind(yield_stock[nyrs]),yield_total[1],yield_total[nyrs], rbind(u_region[,1]),
      rbind(u_region[,nyrs]),rbind(u_stock[1]),rbind(u_stock[nyrs]),u_total[1],u_total[nyrs],
      rbind(dep_region[,1]),rbind(dep_region[,nyrs]),rbind(dep_stock[1]),rbind(dep_stock[nyrs]),dep_total[1],
      dep_total[nyrs],rbind(SSB_region[,1]),rbind(SSB_region[,nyrs]),rbind(SSB_stock[1]),rbind(SSB_stock[nyrs]),
      SSB_total[1],SSB_total[nyrs],Bratio_stock[1],Bratio_stock[nyrs],Bratio_total[1],Bratio_total[nyrs], rbind(bio_region_overlap[,1]),
      rbind(bio_region_overlap[,nyrs]),rbind(bio_stock_overlap[1]),rbind(bio_stock_overlap[nyrs]),bio_total_overlap[1],
      bio_total_overlap[nyrs], rbind(yield_region_overlap[,1]),rbind(yield_region_overlap[,nyrs]),rbind(yield_stock_overlap[1]),
      rbind(yield_stock_overlap[nyrs]),yield_total_overlap[1],yield_total_overlap[nyrs],
      rbind(u_region_overlap[,1]),rbind(u_region_overlap[,nyrs]),rbind(u_stock_overlap[1]),rbind(u_stock_overlap[nyrs]),
      u_total_overlap[1],u_total_overlap[nyrs], rbind(dep_region_overlap[,1]),rbind(dep_region_overlap[,nyrs]),
      rbind(dep_stock_overlap[1]),rbind(dep_stock_overlap[nyrs]),dep_total_overlap[1],dep_total_overlap[nyrs], 
      rbind(SSB_region_overlap[,1]),rbind(SSB_region_overlap[,nyrs]),rbind(SSB_stock_overlap[1]),rbind(SSB_stock_overlap[nyrs]),
      SSB_total_overlap[1],SSB_total_overlap[nyrs],Bratio_stock_overlap[1],Bratio_stock_overlap[nyrs],
      Bratio_natal_overlap[1],Bratio_natal_overlap[nyrs]
    )
  }else # pull out results for more than one stock with more than one region
    {
  if(sum(nstocks)>1 & sum(nregions)>1)
    {
     matrix(c(i,as.vector(t(F.region1)),rbind(bio_region[,1]),rbind(bio_region[,nyrs]),rbind(bio_stock[,1]),
              rbind(bio_stock[,nyrs]),bio_total[1],bio_total[nyrs],rbind(yield_region[,1]),rbind(yield_region[,nyrs]),
              rbind(yield_stock[,1]),rbind(yield_stock[,nyrs]),yield_total[1],yield_total[nyrs], rbind(u_region[,1]),
              rbind(u_region[,nyrs]),rbind(u_stock[,1]),rbind(u_stock[,nyrs]),u_total[1],u_total[nyrs], 
              rbind(dep_region[,1]),rbind(dep_region[,nyrs]),rbind(dep_stock[,1]),rbind(dep_stock[,nyrs]),dep_total[1],
              dep_total[nyrs],rbind(SSB_region[,1]),rbind(SSB_region[,nyrs]),rbind(SSB_stock[,1]),rbind(SSB_stock[,nyrs]),
              SSB_total[1],SSB_total[nyrs],Bratio_stock[,1],Bratio_stock[,nyrs],Bratio_total[1],Bratio_total[nyrs],
              rbind(bio_region_overlap[,1]),rbind(bio_region_overlap[,nyrs]),rbind(bio_stock_overlap[,1]),rbind(bio_stock_overlap[,nyrs]),
              bio_total_overlap[,1],bio_total_overlap[,nyrs], rbind(yield_region_overlap[,1]),rbind(yield_region_overlap[,nyrs]),
              rbind(yield_stock_overlap[,1]),rbind(yield_stock_overlap[,nyrs]),yield_total_overlap[,1],yield_total_overlap[,nyrs],
              rbind(u_region_overlap[,1]),rbind(u_region_overlap[,nyrs]),rbind(u_stock_overlap[,1]),rbind(u_stock_overlap[,nyrs]),
              u_total_overlap[,1],u_total_overlap[,nyrs],rbind(dep_region_overlap[,1]),rbind(dep_region_overlap[,nyrs]),
              rbind(dep_stock_overlap[,1]),rbind(dep_stock_overlap[,nyrs]),dep_total_overlap[,1],dep_total_overlap[,nyrs],
              rbind(SSB_region_overlap[,1]),rbind(SSB_region_overlap[,nyrs]),rbind(SSB_stock_overlap[,1]),rbind(SSB_stock_overlap[,nyrs]),
              SSB_total_overlap[,1],SSB_total_overlap[,nyrs],Bratio_stock_overlap[,1],Bratio_stock_overlap[,nyrs],Bratio_natal_overlap[,1],
              Bratio_natal_overlap[,nyrs]
              ),nrow=1)
    }
   }
  }
}

# end of MSY_search function


######################################
# Set up and run parallel processing 
#####################################

#detect computer cores
cl <- makeSOCKcluster(detectCores()) 

# export the job to the different cores
clusterExport(cl, c("wd","ntrials","F.name","permutation"))
registerDoSNOW(cl)

# build progress bar
pb <- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Simulation Run 0 of ",ntrials,sep=""),max=100)
progress<-function(n) setWinProgressBar(pb,(n/ntrials*100),label=paste("Simulation Run", n,"of", ntrials,"Completed"))
opts<-list(progress=progress)

# set up the permutations/loops using the 'foreach' package and syntax for parallel processing and run. Collect results in t
t<- foreach(i=1:ntrials,.combine=rbind,.options.snow=opts,.packages=c('PBSmodelling','matrixStats','TeachingDemos','snowfall','parallel')
) %dopar% { MSY(wd,i,F.name,permutation) }



# Change names of the output
if(sum(nstocks)<2  & sum(nregions)>1)
  {
    colnames(t)=c("trial",rep("F",times=sum(nfleets)),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),
                  rep("bio_stock_st",times=sum(nstocks)),rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",
                  rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                  rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),"yield_tot_st","yield_tot_end",
                  rep("u_reg_st",times=sum(nregions)), rep("u_reg_end",times=sum(nregions)),rep("u_stock_st",times=sum(nstocks)), 
                  rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end",
                  rep("dep_reg_st",times=sum(nregions)), rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)), 
                  rep("dep_stock_end",times=sum(nstocks)), "dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)), 
                  rep("SSB_reg_end",times=sum(nregions)),rep("SSB_stock_st",times=sum(nstocks)), 
                  rep("SSB_stock_end",times=sum(nstocks)),"SSB_tot_st","SSB_tot_end",rep("Bratio_stock_st",times=sum(nstocks)),
                  rep("Bratio_stock_end",times=sum(nstocks)),"Bratio_tot_st","Bratio_tot_end",
                  rep("bio_reg_st_over",times=(sum(nregions))),rep("bio_reg_end_over",times=(sum(nregions))),rep("bio_stock_st_over",times=(sum(nstocks))^2),
                  rep("bio_stock_end_over",times=(sum(nstocks))^2),rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)),
                  rep("yield_reg_st_over",times=(sum(nregions))),rep("yield_reg_end_over",times=(sum(nregions))),
                  rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),rep("yield_natal_st_over",times=sum(nstocks)), 
                  rep("yield_natal_end_over",times=sum(nstocks)),rep("u_reg_st_over",times=(sum(nregions))),
                  rep("u_reg_end_over",times=(sum(nregions))),rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                  rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                  rep("dep_reg_st_over",times=(sum(nregions))), rep("dep_reg_end_over",times=(sum(nregions))),rep("dep_stock_st_over",times=(sum(nstocks))^2), 
                  rep("dep_stock_end_over",times=(sum(nstocks))^2),
                  rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),rep("SSB_reg_st_over",times=(sum(nregions))), 
                  rep("SSB_reg_end_over",times=(sum(nregions))),rep("SSB_stock_st_over",times=(sum(nstocks))^2), 
                  rep("SSB_stock_end_over",times=(sum(nstocks))^2),rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)),
                  rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),
                  rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks)))

    }    
  if(sum(nregions)>1  | sum(nstocks)>1)
   {
  colnames(t)=c("trial",rep("F",times=sum(nfleets)),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),
                rep("bio_stock_st",times=sum(nstocks)), rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",
                rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),"yield_tot_st","yield_tot_end",
                rep("u_reg_st",times=sum(nregions)), rep("u_reg_end",times=sum(nregions)),rep("u_stock_st",times=sum(nstocks)), 
                rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end", rep("dep_reg_st",times=sum(nregions)), 
                rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)), rep("dep_stock_end",times=sum(nstocks)),
                "dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)), rep("SSB_reg_end",times=sum(nregions)),
                rep("SSB_stock_st",times=sum(nstocks)), rep("SSB_stock_end",times=sum(nstocks)),"SSB_tot_st","SSB_tot_end",
                rep("Bratio_stock_st",times=sum(nstocks)), rep("Bratio_stock_end",times=sum(nstocks)),"Bratio_tot_st","Bratio_tot_end",
                rep("bio_reg_st_over",times=(sum(nregions))^2),rep("bio_reg_end_over",times=(sum(nregions))^2),rep("bio_stock_st_over",times=(sum(nstocks))^2),
                rep("bio_stock_end_over",times=(sum(nstocks))^2),rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)),
                rep("yield_reg_st_over",times=(sum(nregions))^2),rep("yield_reg_end_over",times=(sum(nregions))^2),
                rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),rep("yield_natal_st_over",times=sum(nstocks)), 
                rep("yield_natal_end_over",times=sum(nstocks)),rep("u_reg_st_over",times=(sum(nregions))^2),
                rep("u_reg_end_over",times=(sum(nregions))^2),rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                rep("dep_reg_st_over",times=(sum(nregions))^2), rep("dep_reg_end_over",times=(sum(nregions))^2),rep("dep_stock_st_over",times=(sum(nstocks))^2), 
                rep("dep_stock_end_over",times=(sum(nstocks))^2),
                rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),rep("SSB_reg_st_over",times=(sum(nregions))^2), 
                rep("SSB_reg_end_over",times=(sum(nregions))^2),rep("SSB_stock_st_over",times=(sum(nstocks))^2), 
                rep("SSB_stock_end_over",times=(sum(nstocks))^2),rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)),
                rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),
                rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks)))
  } 
#end change col names for t matrix output

#write the output as a csv in the results/figures file
write.csv(t, file=paste0(wd,"/MSY Results/Figures/Output Quantities.csv",sep=""))

#close the parallel processing
close(pb) 
stopCluster(cl)
closeAllConnections()
output3<-gc()


##################################
# MAKE PLOTS for MSY determination
##################################

# read CSV for graphing
graph<-read.csv(paste0(wd,"/MSY Results/Figures/Output Quantities.csv",sep=""))

#pull some key things for plotting
u<-graph$u_tot_end #fishing rate
MSY1<-graph$yield_tot_end #total yeild at model end
SPR<-graph$Bratio_tot_end #Bratio
bio<-graph$bio_tot_end #total biomass at the end

#MSY!
stock<-graph[which.max(graph$yield_tot_end),] # pull out the run with max yeild at the end

#create a data frame of yield and SPR
MSY.SPR<-data.frame(cbind(SPR,MSY1))
MSY.SPR<-data.table(MSY.SPR,key="SPR") #key = sort by SPR

#create a data frame of yield and SPR
MSY.u<-cbind(u,MSY1)
MSY.u<-data.table(MSY.u,key="u")

#create a data frame of biomass and yield
MSY.bio<-cbind(bio,MSY1)
MSY.bio<-data.table(MSY.bio,key="bio")

#create a data frame of biomass and fishing rate
u.bio<-cbind(bio,u)
u.bio<-data.table(u.bio,key="bio")


#create matrices for biomass, yield, SPR, u by trial and stock
bio.stock<-matrix(NA,nrow=ntrials,ncol=nstocks)
MSY.stock<-matrix(NA,nrow=ntrials,ncol=nstocks)
SPR.stock<-matrix(NA,nrow=ntrials,ncol=nstocks)
u.stock<-matrix(NA,nrow=ntrials,ncol=nstocks)

MSY.SPR.stock<-array(NA,c(ntrials,2,nstocks))
MSY.bio.stock<-array(NA,c(ntrials,2,nstocks))
MSY.u.stock<-array(NA,c(ntrials,2,nstocks))
u.bio.stock<-array(NA,c(ntrials,2,nstocks))
msy.run.stock<-matrix(1:nstocks)



#fill in the matrices from graph data.frame
for(i in 1:nstocks)
{
  bio.stock[,i]<-graph[,(11+i-1)]
  MSY.stock[,i]<-graph[,(21+i-1)]
  SPR.stock[,i]<-graph[,(57+i-1)]
  u.stock[,i]<-graph[,(31+i-1)]

  MSY.SPR.stock[,,i]<-cbind(SPR.stock[,i],MSY.stock[,i])
  MSY.u.stock[,,i]<-cbind(u.stock[,i],MSY.stock[,i])
  MSY.bio.stock[,,i]<-cbind(bio.stock[,i],MSY.stock[,i])
  u.bio.stock[,,i]<-cbind(bio.stock[,i],u.stock[,i])
  
  msy.run.stock[i]<-which.max(MSY.SPR.stock[,2,i])
}


# Create a new csv with the MSY values for F. Save to results file
MSY.values<-rbind(graph[which.max(graph$yield_tot_end),],graph[c(msy.run.stock),])
write.csv(MSY.values, file=paste0(wd,"/MSY Results/Figures/MSY Outputs.csv",sep=""))

######################
#STOCK SPECIFIC PLOTS
######################

# open a png file for results
png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (STOCK).png",sep=""))

##PLOT 1
#set plot parameters and plot
par(mfrow=c(2,2)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
par(las=1, mar=c(5,5,3,1))
plot(MSY.SPR.stock[,1,],MSY.SPR.stock[,2,], type="n", xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio by Stock',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.SPR.stock[,2,]),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR.stock[,1,]),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR.stock[,2,]),1000),Roundup(max(MSY.SPR.stock[,2,]),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR.stock[,1,]),.1), Roundup(max(MSY.SPR.stock[,1,]),.1)/5),cex.axis=.8)


#plot equilibrium SSB vs. yield
for(i in 1:nstocks)
  {
  matlines(MSY.SPR.stock[,1,i],MSY.SPR.stock[,2,i], type="p", pch=(21+i-1))
  segments(0, max(MSY.SPR.stock[,2,i]), MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
  segments(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],0, MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
  matpoints(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],max(MSY.SPR.stock[,2,i]),pch=19,col='grey70',cex=2);
  }



##PLOT 2
#set up plotting params and plot for equilibrium harvest vs. yield
plot(MSY.u.stock[,1,],MSY.u.stock[,2,], type="n", xlab = 'Equilibrium Harvest Rate', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Harvest Rate by Stock',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.u.stock[,2,]),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.u.stock[,1,]),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.u.stock[,2,]),1000),Roundup(max(MSY.u.stock[,2,]),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.u.stock[,1,]),.1), Roundup(max(MSY.u.stock[,1,]),.1)/5),cex.axis=.8)

# add data to the plot for each stock
  for(i in 1:nstocks)
  {
    matlines(MSY.u.stock[,1,i],MSY.u.stock[,2,i], type="p", pch=(21+i-1))
    segments(0, max(MSY.u.stock[,2,i]), MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
    segments(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],0, MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
    matpoints(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],max(MSY.u.stock[,2,i]),pch=19,col='grey70',cex=2);
  }




##PLOT 3
#set up plot for Equilibrium biomass vs. yeild
plot(MSY.bio.stock[,1,],MSY.bio.stock[,2,], type="n", xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass by Stock',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.bio.stock[,2,]),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio.stock[,1,i]),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.bio.stock[,2,]),1000),Roundup(max(MSY.bio.stock[,2,]),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.bio.stock[,1,]),.1), Roundup(max(MSY.bio.stock[,1,]),.1)/5),cex.axis=.8)

#add data for each plot
for(i in 1:nstocks)
{
  matlines(MSY.bio.stock[,1,i],MSY.bio.stock[,2,i], type="p", pch=(21+i-1))
  segments(0, max(MSY.bio.stock[,2,i]), MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
  segments(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],0, MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
  matpoints(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],max(MSY.bio.stock[,2,i]),pch=19,col='grey70',cex=2);
}


##PLOT 4
#set up plot for Equilbrium biomass vs. harvest rate
plot(u.bio.stock[,1,],u.bio.stock[,2,], type="n", xlab = 'Equilibrium Biomass', ylab = 'Harvest Rate', lwd=2,lty=1, main='Harvest Rate vs. Biomass by Stock',cex.main=1.15,
     ylim=c(0,Roundup(max(u.bio.stock[,2,]),.1)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(u.bio.stock[,1,i]),1000)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(u.bio.stock[,2,]),.1),Roundup(max(u.bio.stock[,2,]),.1)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(u.bio.stock[,1,]),1000), Roundup(max(u.bio.stock[,1,]),1000)/5),cex.axis=.8)

#add data
for(i in 1:nstocks)
{
  matlines(u.bio.stock[,1,i],u.bio.stock[,2,i], type="p", pch=(21+i-1))
}

dev.off()


#convexhull.xy(MSY.SPR)
#plot(hull)



#################
# MSY PLOTS
#################

png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (TOTAL).png",sep=""))
par(mfrow=c(2,2)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)

#PLOT1 
#MSY plot
par(las=1, mar=c(5,5,3,1))
plot(MSY.SPR$SPR,MSY.SPR$MSY, typ='p', xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio',cex.main=1.15,
        ylim=c(0,Roundup(max(MSY.SPR$MSY),1000)),
        cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR$SPR),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR$MSY),1000),Roundup(max(MSY.SPR$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR$SPR),.1), Roundup(max(MSY.SPR$SPR),.1)/5),cex.axis=.8)

#library(lawn)
#which.min(MSY.SPR$MSY[MSY.SPR$SPR=c(1.88:2)])
#plot(lawn_concave(MSY.SPR),add=TRUE)
#plot(convexhull.xy(MSY.SPR),add=TRUE,lwd=3)

segments(0, max(MSY.SPR$MSY), MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],0, MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],max(MSY.SPR$MSY),pch=19,col='grey70',cex=2);



#PLOT 2
# Harvest vs. yeild
par(las=1, mar=c(5,5,3,1))
plot(MSY.u$u,MSY.u$MSY, typ='p', xlab = 'Equilibrium Harvest Rate', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Harvest Rate',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.u$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.u$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.u$MSY),1000),Roundup(max(MSY.u$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.u$u),.1), Roundup(max(MSY.u$u),.1)/5),cex.axis=.8)

segments(0, max(MSY.u$MSY), MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.u$u[which.max(MSY.u$MSY)],0, MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.u$u[which.max(MSY.u$MSY)],max(MSY.u$MSY),pch=19,col='grey70',cex=2);


#PLOT 3
#Equilibrium Biomass vs. Yield

par(las=1, mar=c(5,5,3,1))
plot(MSY.bio$bio,MSY.bio$MSY, typ='p', xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.bio$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio$bio),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.bio$MSY),1000),Roundup(max(MSY.bio$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.bio$bio),1000), Roundup(max(MSY.bio$bio),1000)/5),cex.axis=.8)

segments(0, max(MSY.bio$MSY), MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.bio$bio[which.max(MSY.bio$MSY)],0, MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.bio$bio[which.max(MSY.bio$MSY)],max(MSY.bio$MSY),pch=19,col='grey70',cex=2);


#PLOT 4 
#Equilibrium Harvest vs. Biomass
par(las=1, mar=c(5,5,3,1))
plot(u.bio$u,u.bio$bio, typ='p', xlab = 'Equilibrium Harvest', ylab = 'Biomass', lwd=2,lty=1, main='Biomass vs. Harvest',cex.main=1.15,
     ylim=c(0,Roundup(max(u.bio$bio),.1)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(u.bio$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(u.bio$bio),1000),Roundup(max(u.bio$bio),1000)/10),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(u.bio$u),.1), Roundup(max(u.bio$u),.1)/10),cex.axis=.8)

dev.off() #close plotting connection




#####################
# If there are 2 stocks
#####################

if(sum(nstocks==2))
{
  dir.create(paste0(wd,"/Stock1",sep=""))
  wd1<<-paste0(wd,"/Stock1",sep="") #F:/NOAA FILES/Research/Spatial BRPs/Initial Model/F iteration'
  setwd(wd1)
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.exe",sep=""),to=paste0(wd1,"/Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.dat",sep=""),to=paste0(wd1,"/Spatial_BRP.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.tpl",sep=""),to=paste0(wd1,"/Spatial_BRP.tpl",sep="")))
  dir.create(paste0(wd1,"/MSY Results",sep=""))
  dir.create(paste0(wd1,"/MSY Results/Figures",sep=""))
  dir.create(paste0(wd1,"/MSY Results/Report Files",sep=""))
  
  permutation=cbind(F.seq,rep(0,times=length(F.seq)))
  ntrials1<-length(F.seq)
  cl <- makeSOCKcluster(detectCores())
  clusterExport(cl, c("wd1","ntrials","F.name","permutation"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Simulation Run 0 of ",ntrials1,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/ntrials1*100),label=paste("Simulation Run", n,"of", ntrials1,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:ntrials1,.combine=rbind,.options.snow=opts,.packages=c('PBSmodelling','matrixStats','TeachingDemos','snowfall','parallel')
  ) %dopar% {
    MSY(wd1,i,F.name,permutation) 
  }
  
  colnames(t)=c("trial",rep("F",times=sum(nfleets)),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),rep("bio_stock_st",times=sum(nstocks)),
                rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),"yield_tot_st","yield_tot_end",rep("u_reg_st",times=sum(nregions)),
                rep("u_reg_end",times=sum(nregions)),rep("u_stock_st",times=sum(nstocks)), rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end",
                rep("dep_reg_st",times=sum(nregions)), rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)), rep("dep_stock_end",times=sum(nstocks)),
                "dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)), rep("SSB_reg_end",times=sum(nregions)),rep("SSB_stock_st",times=sum(nstocks)), 
                rep("SSB_stock_end",times=sum(nstocks)),"SSB_tot_st","SSB_tot_end",rep("Bratio_stock_st",times=sum(nstocks)), rep("Bratio_stock_end",times=sum(nstocks)),"Bratio_tot_st","Bratio_tot_end",
                rep("bio_reg_st_over",times=(sum(nregions))^2),rep("bio_reg_end_over",times=(sum(nregions))^2),rep("bio_stock_st_over",times=(sum(nstocks))^2),
                rep("bio_stock_end_over",times=(sum(nstocks))^2),rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)),
                rep("yield_reg_st_over",times=(sum(nregions))^2),rep("yield_reg_end_over",times=(sum(nregions))^2),
                rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),rep("yield_natal_st_over",times=sum(nstocks)), 
                rep("yield_natal_end_over",times=sum(nstocks)),rep("u_reg_st_over",times=(sum(nregions))^2),
                rep("u_reg_end_over",times=(sum(nregions))^2),rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                rep("dep_reg_st_over",times=(sum(nregions))^2), rep("dep_reg_end_over",times=(sum(nregions))^2),rep("dep_stock_st_over",times=(sum(nstocks))^2), 
                rep("dep_stock_end_over",times=(sum(nstocks))^2),
                rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),rep("SSB_reg_st_over",times=(sum(nregions))^2), 
                rep("SSB_reg_end_over",times=(sum(nregions))^2),rep("SSB_stock_st_over",times=(sum(nstocks))^2), 
                rep("SSB_stock_end_over",times=(sum(nstocks))^2),rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)),
                rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks)))
  
  write.csv(t, file=paste0(wd1,"/MSY Results/Figures/Output Quantities (Stock 1 MSY).csv",sep=""))
  
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  
  dir.create(paste0(wd,"/Stock2",sep=""))
  wd2<<-paste0(wd,"/Stock2",sep="") #F:/NOAA FILES/Research/Spatial BRPs/Initial Model/F iteration'
  setwd(wd2)
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.exe",sep=""),to=paste0(wd2,"/Spatial_BRP.exe",sep="")))
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.dat",sep=""),to=paste0(wd2,"/Spatial_BRP.dat",sep="")))
  invisible(file.copy(from=paste0(wd,"/Spatial_BRP.tpl",sep=""),to=paste0(wd2,"/Spatial_BRP.tpl",sep="")))
  dir.create(paste0(wd2,"/MSY Results",sep=""))
  dir.create(paste0(wd2,"/MSY Results/Figures",sep=""))
  dir.create(paste0(wd2,"/MSY Results/Report Files",sep=""))
  
  permutation<-cbind(rep(0,times=length(F.seq)),F.seq)
  cl <- makeSOCKcluster(detectCores())
  clusterExport(cl, c("wd2","ntrials1","F.name","permutation"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste("MSY Search Progress Bar"), label=paste("Simulation Run 0 of ",ntrials1,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/ntrials1*100),label=paste("Simulation Run", n,"of", ntrials1,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:ntrials1,.combine=rbind,.options.snow=opts,.packages=c('PBSmodelling','matrixStats','TeachingDemos','snowfall','parallel')
  ) %dopar% {
    MSY(wd2,i,F.name,permutation) 
  }
  
  colnames(t)=c("trial",rep("F",times=sum(nfleets)),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),rep("bio_stock_st",times=sum(nstocks)),
                rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),"yield_tot_st","yield_tot_end",rep("u_reg_st",times=sum(nregions)),
                rep("u_reg_end",times=sum(nregions)),rep("u_stock_st",times=sum(nstocks)), rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end",
                rep("dep_reg_st",times=sum(nregions)), rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)), rep("dep_stock_end",times=sum(nstocks)),
                "dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)), rep("SSB_reg_end",times=sum(nregions)),rep("SSB_stock_st",times=sum(nstocks)), 
                rep("SSB_stock_end",times=sum(nstocks)),"SSB_tot_st","SSB_tot_end",rep("Bratio_stock_st",times=sum(nstocks)), rep("Bratio_stock_end",times=sum(nstocks)),"Bratio_tot_st","Bratio_tot_end",
                rep("bio_reg_st_over",times=(sum(nregions))^2),rep("bio_reg_end_over",times=(sum(nregions))^2),rep("bio_stock_st_over",times=(sum(nstocks))^2),
                rep("bio_stock_end_over",times=(sum(nstocks))^2),rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)),
                rep("yield_reg_st_over",times=(sum(nregions))^2),rep("yield_reg_end_over",times=(sum(nregions))^2),
                rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),rep("yield_natal_st_over",times=sum(nstocks)), 
                rep("yield_natal_end_over",times=sum(nstocks)),rep("u_reg_st_over",times=(sum(nregions))^2),
                rep("u_reg_end_over",times=(sum(nregions))^2),rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                rep("dep_reg_st_over",times=(sum(nregions))^2), rep("dep_reg_end_over",times=(sum(nregions))^2),rep("dep_stock_st_over",times=(sum(nstocks))^2), 
                rep("dep_stock_end_over",times=(sum(nstocks))^2),
                rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),rep("SSB_reg_st_over",times=(sum(nregions))^2), 
                rep("SSB_reg_end_over",times=(sum(nregions))^2),rep("SSB_stock_st_over",times=(sum(nstocks))^2), 
                rep("SSB_stock_end_over",times=(sum(nstocks))^2),rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)),
                rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks)))
  
  write.csv(t, file=paste0(wd2,"/MSY Results/Figures/Output Quantities (Stock 2 MSY).csv",sep=""))
  
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  
  
  graph1<-read.csv(paste0(wd1,"/MSY Results/Figures/Output Quantities (Stock 1 MSY).csv",sep=""))
  graph2<-read.csv(paste0(wd2,"/MSY Results/Figures/Output Quantities (Stock 2 MSY).csv",sep=""))
  
  bio.stock<-matrix(NA,nrow=ntrials1,ncol=nstocks)
  MSY.stock<-matrix(NA,nrow=ntrials1,ncol=nstocks)
  SPR.stock<-matrix(NA,nrow=ntrials1,ncol=nstocks)
  u.stock<-matrix(NA,nrow=ntrials1,ncol=nstocks)
  
  MSY.SPR.stock<-array(NA,c(ntrials1,2,nstocks))
  MSY.bio.stock<-array(NA,c(ntrials1,2,nstocks))
  MSY.u.stock<-array(NA,c(ntrials1,2,nstocks))
  u.bio.stock<-array(NA,c(ntrials1,2,nstocks))
  msy.run.stock<-matrix(1:nstocks)
  
  bio.stock[,1]<-graph1[,(11+1-1)]
  MSY.stock[,1]<-graph1[,(21+1-1)]
  SPR.stock[,1]<-graph1[,(57+1-1)]
  u.stock[,1]<-graph1[,(31+1-1)]
  
  bio.stock[,2]<-graph2[,(11+2-1)]
  MSY.stock[,2]<-graph2[,(21+2-1)]
  SPR.stock[,2]<-graph2[,(57+2-1)]
  u.stock[,2]<-graph2[,(31+2-1)]
  
  for(i in 1:nstocks)
  {
    
    MSY.SPR.stock[,,i]<-cbind(SPR.stock[,i],MSY.stock[,i])
    MSY.u.stock[,,i]<-cbind(u.stock[,i],MSY.stock[,i])
    MSY.bio.stock[,,i]<-cbind(bio.stock[,i],MSY.stock[,i])
    u.bio.stock[,,i]<-cbind(bio.stock[,i],u.stock[,i])
    
    msy.run.stock[i]<-which.max(MSY.SPR.stock[,2,i])
  }
  
  
  png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (STOCK MAX).png",sep=""))
  par(mfrow=c(2,2)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  par(las=1, mar=c(5,5,3,1))
  plot(MSY.SPR.stock[,1,],MSY.SPR.stock[,2,], type="n", xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio by Stock',cex.main=1.15,
       ylim=c(0,Roundup(max(MSY.SPR.stock[,2,]),1000)),
       cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR.stock[,1,]),.1)) ,axes=FALSE)
  axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR.stock[,2,]),1000),Roundup(max(MSY.SPR.stock[,2,]),1000)/5),cex.axis=.8)
  axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR.stock[,1,]),.1), Roundup(max(MSY.SPR.stock[,1,]),.1)/5),cex.axis=.8)
  
  for(i in 1:nstocks)
  {
    matlines(MSY.SPR.stock[,1,i],MSY.SPR.stock[,2,i], type="p", pch=(21+i-1))
    segments(0, max(MSY.SPR.stock[,2,i]), MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
    segments(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],0, MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
    matpoints(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],max(MSY.SPR.stock[,2,i]),pch=19,col='grey70',cex=2);
  }
  
  plot(MSY.u.stock[,1,],MSY.u.stock[,2,], type="n", xlab = 'Equilibrium Harvest Rate', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Harvest Rate by Stock',cex.main=1.15,
       ylim=c(0,Roundup(max(MSY.u.stock[,2,]),1000)),
       cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.u.stock[,1,]),.1)) ,axes=FALSE)
  axis( 2, pos=0,seq(0,Roundup(max(MSY.u.stock[,2,]),1000),Roundup(max(MSY.u.stock[,2,]),1000)/5),cex.axis=.8)
  axis( 1,pos=0,seq(0, Roundup(max(MSY.u.stock[,1,]),.1), Roundup(max(MSY.u.stock[,1,]),.1)/5),cex.axis=.8)
  
  for(i in 1:nstocks)
  {
    matlines(MSY.u.stock[,1,i],MSY.u.stock[,2,i], type="p", pch=(21+i-1))
    segments(0, max(MSY.u.stock[,2,i]), MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
    segments(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],0, MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
    matpoints(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],max(MSY.u.stock[,2,i]),pch=19,col='grey70',cex=2);
  }
  
  
  plot(MSY.bio.stock[,1,],MSY.bio.stock[,2,], type="n", xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass by Stock',cex.main=1.15,
       ylim=c(0,Roundup(max(MSY.bio.stock[,2,]),1000)),
       cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio.stock[,1,i]),.1)) ,axes=FALSE)
  axis( 2, pos=0,seq(0,Roundup(max(MSY.bio.stock[,2,]),1000),Roundup(max(MSY.bio.stock[,2,]),1000)/5),cex.axis=.8)
  axis( 1,pos=0,seq(0, Roundup(max(MSY.bio.stock[,1,]),.1), Roundup(max(MSY.bio.stock[,1,]),.1)/5),cex.axis=.8)
  
  for(i in 1:nstocks)
  {
    matlines(MSY.bio.stock[,1,i],MSY.bio.stock[,2,i], type="p", pch=(21+i-1))
    segments(0, max(MSY.bio.stock[,2,i]), MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
    segments(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],0, MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
    matpoints(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],max(MSY.bio.stock[,2,i]),pch=19,col='grey70',cex=2);
  }
  
  
  plot(u.bio.stock[,2,],u.bio.stock[,1,], type="n", xlab = 'Equilibrium Biomass', ylab = 'Harvest Rate', lwd=2,lty=1, main='Biomass vs. Harvest Rate',cex.main=1.15,
       ylim=c(0,Roundup(max(u.bio.stock[,1,]),1000)),
       cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(u.bio.stock[,2,i]),.1)) ,axes=FALSE)
  axis( 2, pos=0,seq(0,Roundup(max(u.bio.stock[,1,]),1000),Roundup(max(u.bio.stock[,1,]),1000)/5),cex.axis=.8)
  axis( 1,pos=0,seq(0, Roundup(max(u.bio.stock[,2,]),.1), Roundup(max(u.bio.stock[,2,]),.1)/5),cex.axis=.8)
  
  for(i in 1:nstocks)
  {
    matlines(u.bio.stock[,2,i],u.bio.stock[,1,i], type="p", pch=(21+i-1))
  }
  dev.off()


png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (Total vs. Stock).png",sep=""))
par(mfrow=c(2,2)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)

par(las=1, mar=c(5,5,3,1))
plot(MSY.SPR$SPR,MSY.SPR$MSY, typ='p', xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.SPR$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR$SPR),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR$MSY),1000),Roundup(max(MSY.SPR$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR$SPR),.1), Roundup(max(MSY.SPR$SPR),.1)/5),cex.axis=.8)

segments(0, max(MSY.SPR$MSY), MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],0, MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],max(MSY.SPR$MSY),pch=19,col='grey70',cex=2);

for(i in 1:nstocks)
{
  matlines(MSY.SPR.stock[,1,i],MSY.SPR.stock[,2,i], type="p", col="grey",lwd=2,pch=(18+i-1))
  segments(0, max(MSY.SPR.stock[,2,i]), MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
  segments(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],0, MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i], max(MSY.SPR.stock[,2,i]), lty=6, col='grey70', lwd=2)
  matpoints(MSY.SPR.stock[which.max(MSY.SPR.stock[,2,i]),1,i],max(MSY.SPR.stock[,2,i]),pch=19,col='grey70',cex=2);
}

par(las=1, mar=c(5,5,3,1))
plot(MSY.u$u,MSY.u$MSY, typ='p', xlab = 'Equilibrium Harvest Rate', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Harvest Rate',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.u$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.u$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.u$MSY),1000),Roundup(max(MSY.u$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.u$u),.1), Roundup(max(MSY.u$u),.1)/5),cex.axis=.8)

segments(0, max(MSY.u$MSY), MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.u$u[which.max(MSY.u$MSY)],0, MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.u$u[which.max(MSY.u$MSY)],max(MSY.u$MSY),pch=19,col='grey70',cex=2);

for(i in 1:nstocks)
{
  matlines(MSY.u.stock[,1,i],MSY.u.stock[,2,i], type="p", col="grey",lwd=2,pch=(18+i-1))
  segments(0, max(MSY.u.stock[,2,i]), MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
  segments(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],0, MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i], max(MSY.u.stock[,2,i]), lty=6, col='grey70', lwd=2)
  matpoints(MSY.u.stock[which.max(MSY.u.stock[,2,i]),1,i],max(MSY.u.stock[,2,i]),pch=19,col='grey70',cex=2);
}

par(las=1, mar=c(5,5,3,1))
plot(MSY.bio$bio,MSY.bio$MSY, typ='p', xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.bio$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio$bio),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.bio$MSY),1000),Roundup(max(MSY.bio$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.bio$bio),1000), Roundup(max(MSY.bio$bio),1000)/5),cex.axis=.8)

segments(0, max(MSY.bio$MSY), MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.bio$bio[which.max(MSY.bio$MSY)],0, MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.bio$bio[which.max(MSY.bio$MSY)],max(MSY.bio$MSY),pch=19,col='grey70',cex=2);

for(i in 1:nstocks)
{
  matlines(MSY.bio.stock[,1,i],MSY.bio.stock[,2,i], type="p", lwd=2,col="grey",pch=(18+i-1))
  segments(0, max(MSY.bio.stock[,2,i]), MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
  segments(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],0, MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i], max(MSY.bio.stock[,2,i]), lty=6, col='grey70', lwd=2)
  matpoints(MSY.bio.stock[which.max(MSY.bio.stock[,2,i]),1,i],max(MSY.bio.stock[,2,i]),pch=19,col='grey70',cex=2);
}


par(las=1, mar=c(5,5,3,1))
plot(u.bio$u,u.bio$bio, typ='p', xlab = 'Equilibrium Harvest', ylab = 'Biomass', lwd=2,lty=1, main='Biomass vs. Harvest Rate',cex.main=1.15,
     ylim=c(0,Roundup(max(u.bio$bio),.1)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(u.bio$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(u.bio$bio),1000),Roundup(max(u.bio$bio),1000)/10),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(u.bio$u),.1), Roundup(max(u.bio$u),.1)/10),cex.axis=.8)

for(i in 1:nstocks)
{
  matlines(u.bio.stock[,2,i],u.bio.stock[,1,i], type="p", col="grey",lwd=2,pch=(18+i-1))
}
dev.off()



png(file=paste0(wd,"/MSY Results/Figures/MSY Curve (Total With Corresponding Stock Values).png",sep=""))
par(mfrow=c(2,2)) #,  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)

par(las=1, mar=c(5,5,3,1))
plot(MSY.SPR$SPR,MSY.SPR$MSY, typ='p', xlab = 'Equilibrium SSB Ratio', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. SSB Ratio',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.SPR$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.SPR$SPR),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.SPR$MSY),1000),Roundup(max(MSY.SPR$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.SPR$SPR),.1), Roundup(max(MSY.SPR$SPR),.1)/5),cex.axis=.8)

segments(0, max(MSY.SPR$MSY), MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],0, MSY.SPR$SPR[which.max(MSY.SPR$MSY)], max(MSY.SPR$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.SPR$SPR[which.max(MSY.SPR$MSY)],max(MSY.SPR$MSY),pch=19,col='grey70',cex=2);

segments(0, stock$yield_stock_end,stock$Bratio_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
segments(stock$Bratio_stock_end, 0,stock$Bratio_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
matpoints(stock$Bratio_stock_end, stock$yield_stock_end,pch=17,col='grey60',cex=2)

segments(0, stock$yield_stock_end.1,stock$Bratio_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
segments(stock$Bratio_stock_end.1, 0,stock$Bratio_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
matpoints(stock$Bratio_stock_end.1, stock$yield_stock_end.1,pch=15,col='grey60',cex=2)


par(las=1, mar=c(5,5,3,1))
plot(MSY.u$u,MSY.u$MSY, typ='p', xlab = 'Equilibrium Harvest Rate', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Harvest Rate',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.u$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.u$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.u$MSY),1000),Roundup(max(MSY.u$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.u$u),.1), Roundup(max(MSY.u$u),.1)/5),cex.axis=.8)

segments(0, max(MSY.u$MSY), MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.u$u[which.max(MSY.u$MSY)],0, MSY.u$u[which.max(MSY.u$MSY)], max(MSY.u$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.u$u[which.max(MSY.u$MSY)],max(MSY.u$MSY),pch=19,col='grey70',cex=2);

segments(0, stock$yield_stock_end,stock$u_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
segments(stock$u_stock_end, 0,stock$u_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
matpoints(stock$u_stock_end, stock$yield_stock_end,pch=17,col='grey60',cex=2)

segments(0, stock$yield_stock_end.1,stock$u_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
segments(stock$u_stock_end.1, 0,stock$u_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
matpoints(stock$u_stock_end.1, stock$yield_stock_end.1,pch=15,col='grey60',cex=2)



legend('topright',c("Total","Stock 1", "Stock 2"), pch=c(19,17,15),col='grey70')


par(las=1, mar=c(5,5,3,1))
plot(MSY.bio$bio,MSY.bio$MSY, typ='p', xlab = 'Equilibrium Biomass', ylab = 'Yield', lwd=2,lty=1, main='MSY vs. Biomass',cex.main=1.15,
     ylim=c(0,Roundup(max(MSY.bio$MSY),1000)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(MSY.bio$bio),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(MSY.bio$MSY),1000),Roundup(max(MSY.bio$MSY),1000)/5),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(MSY.bio$bio),1000), Roundup(max(MSY.bio$bio),1000)/5),cex.axis=.8)

segments(0, max(MSY.bio$MSY), MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
segments(MSY.bio$bio[which.max(MSY.bio$MSY)],0, MSY.bio$bio[which.max(MSY.bio$MSY)], max(MSY.bio$MSY), lty=6, col='grey70', lwd=2)
matpoints(MSY.bio$bio[which.max(MSY.bio$MSY)],max(MSY.bio$MSY),pch=19,col='grey70',cex=2);

segments(0, stock$yield_stock_end,stock$bio_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
segments(stock$bio_stock_end, 0,stock$bio_stock_end,stock$yield_stock_end,lty=2,col='grey60',lwd=2)
matpoints(stock$bio_stock_end, stock$yield_stock_end,pch=17,col='grey60',cex=2)

segments(0, stock$yield_stock_end.1,stock$bio_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
segments(stock$bio_stock_end.1, 0,stock$bio_stock_end.1,stock$yield_stock_end.1,lty=2,col='grey60',lwd=2)
matpoints(stock$bio_stock_end.1, stock$yield_stock_end.1,pch=15,col='grey60',cex=2)



par(las=1, mar=c(5,5,3,1))
plot(u.bio$u,u.bio$bio, typ='p', xlab = 'Equilibrium Harvest', ylab = 'Biomass', lwd=2,lty=1, main='Biomass vs. Harvest Rate',cex.main=1.15,
     ylim=c(0,Roundup(max(u.bio$bio),.1)),
     cex.axis=1., cex.lab=1.,xlim=c(0.0, Roundup(max(u.bio$u),.1)) ,axes=FALSE)
axis( 2, pos=0,seq(0,Roundup(max(u.bio$bio),1000),Roundup(max(u.bio$bio),1000)/10),cex.axis=.8)
axis( 1,pos=0,seq(0, Roundup(max(u.bio$u),.1), Roundup(max(u.bio$u),.1)/10),cex.axis=.8)

segments(0, MSY.bio$bio[which.max(MSY.bio$MSY)], MSY.u$u[which.max(MSY.u$MSY)], MSY.bio$bio[which.max(MSY.bio$MSY)], lty=6, col='grey70', lwd=2)
segments(MSY.u$u[which.max(MSY.u$MSY)],0, MSY.u$u[which.max(MSY.u$MSY)], MSY.bio$bio[which.max(MSY.bio$MSY)], lty=6, col='grey70', lwd=2)
matpoints(MSY.u$u[which.max(MSY.u$MSY)],MSY.bio$bio[which.max(MSY.bio$MSY)],pch=19,col='grey70',cex=2);



segments(0, stock$bio_stock_end,stock$u_stock_end,stock$bio_stock_end,lty=2,col='grey60',lwd=2)
segments(stock$u_stock_end, 0,stock$u_stock_end,stock$bio_stock_end,lty=2,col='grey60',lwd=2)
matpoints(stock$u_stock_end, stock$bio_stock_end,pch=17,col='grey60',cex=2)

segments(0, stock$bio_stock_end.1,stock$u_stock_end.1,stock$bio_stock_end.1,lty=2,col='grey60',lwd=2)
segments(stock$u_stock_end.1, 0,stock$u_stock_end.1,stock$bio_stock_end.1,lty=2,col='grey60',lwd=2)
matpoints(stock$u_stock_end.1, stock$bio_stock_end.1,pch=15,col='grey60',cex=2)


dev.off()


}

