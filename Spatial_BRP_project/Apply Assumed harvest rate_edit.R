######################################################################
# Script to apply MSY to Spatial BRP
# Created by Dan Goethel
# Edited by Katelyn Bosley
# Date: 2/8/2015
######################################################################

#Dan's notes:
# Script to Apply an assumed harvest rate or TAC
# Need to use the true .dat file, turn use_TAC switch to 1 (TAC) or 2 (u), and input the corresponding value for either input (initial F has no impact)
# Need to copy MSY outputs file for both the assumed and true MSY to the current directory and rename to Assumed MSY or True MSY



# set working directory
wd<-"C:\\Users\\katelyn.bosley\\Desktop\\NOAA_POSTDOC\\SPASM\\Potentially Unuseful Code\\Spatial BRPs\\Most Recent Code\\MSY_test"
setwd(wd)


#load libraries
#--------------------------------------------------------------------------------------------------------------------------------------------------
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
#--------------------------------------------------------------------------------------------------------------------------------------------------


# create a directory for the results
dir.create(paste0(wd,"/MSY Results",sep=""))

# create a function to round up values
Roundup <- function(from,to) ceiling(from/to)*to


#run the ADMB model with .exe
  invisible(shell("Spatial_BRP -nohess",wait=T)) #show.output.on.console=FALSE))  ### Run ADMB with update F
  
#remove junk files
  invisible(file.remove(paste0(wd,"/Spatial_BRP.bar",sep="")))
  invisible(file.remove(paste0(wd,"/Spatial_BRP.log",sep="")))
  invisible(file.remove(paste0(wd,"/variance",sep="")))
  invisible(file.remove(paste0(wd,"/fmin.log",sep="")))                      
                          
# read the report file     
  out =readList("Spatial_BRP.rep")

# create the parameter names
  par_names=c('overlap_switch','larval_move_switch','move_switch','F','biomass_AM','biomass_stock','biomass_total','yield_fleet','yield_region',
              'yield_stock','yield_total', 'harvest_rate_region_bio','harvest_rate_stock_bio','harvest_rate_total_bio',
              'depletion_region','depletion_stock','depletion_total', 'SSB_region','SSB_stock','SSB_total','Bratio_stock',
              'Bratio_total', 'biomass_AM_overlap_region','biomass_stock_overlap','biomass_natal_overlap','yield_region_overlap',
              'yield_stock_overlap','yield_natal_overlap','harvest_rate_region_bio_overlap','harvest_rate_stock_bio_overlap',
              'harvest_rate_natal_bio_overlap','depletion_region_overlap','depletion_stock_overlap',
              'depletion_natal_overlap','SSB_region_overlap','SSB_stock_overlap','SSB_natal_overlap','Bratio_stock_overlap','Bratio_natal_overlap')
  
# pull out the results we are interested in
  result=out[par_names]
  
# pull out the population params
  nyrs<-out$nyrs
  nstocks<-out$nstocks
  nregions<-out$nregions
  nfleets<-out$nfleets
  nages<-out$nages
  

# assign separate objects to results components
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
  
  
# Create matrix of output values
  
# only one region
  if(sum(nregions)==1)
  { 
    F.region<-out$F[nyrs,]
   t<- matrix(c(F.region,bio_region[1],bio_region[nyrs],bio_stock[1],bio_stock[nyrs],bio_total[1],bio_total[nyrs],
      yield_region[1],yield_region[nyrs],yield_stock[1],yield_stock[nyrs],yield_total[1],yield_total[nyrs],
      u_region[1],u_region[nyrs],u_stock[1],u_stock[nyrs],u_total[1],u_total[nyrs],
      dep_region[1],dep_region[nyrs],dep_stock[1],dep_stock[nyrs],dep_total[1],dep_total[nyrs],
      SSB_region[1],SSB_region[nyrs],SSB_stock[1],SSB_stock[nyrs],SSB_total[1],SSB_total[nyrs],Bratio_stock[1],
      Bratio_stock[nyrs],Bratio_total[1],Bratio_total[nyrs],
      bio_region_overlap[1],bio_region_overlap[nyrs],bio_stock_overlap[1],bio_stock_overlap[nyrs],bio_total_overlap[1],bio_total_overlap[nyrs],
      yield_region_overlap[1],yield_region_overlap[nyrs],yield_stock_overlap[1],yield_stock_overlap[nyrs],yield_total_overlap[1],
      yield_total_overlap[nyrs],u_region_overlap[1],u_region_overlap[nyrs],u_stock_overlap[1],u_stock_overlap[nyrs],u_total_overlap[1],
      u_total_overlap[nyrs], dep_region_overlap[1],dep_region_overlap[nyrs],dep_stock_overlap[1],dep_stock_overlap[nyrs],
      dep_total_overlap[1],dep_total_overlap[nyrs], SSB_region_overlap[1],SSB_region_overlap[nyrs],SSB_stock_overlap[1],
      SSB_stock_overlap[nyrs],SSB_total_overlap[1],SSB_total_overlap[nyrs],Bratio_stock_overlap[1],Bratio_stock_overlap[nyrs],
      Bratio_natal_overlap[1],Bratio_natal_overlap[nyrs]),nrow=1)
  }else 
    {
  if(sum(nstocks)<2  & sum(nregions)>1)# if 1 stock with multiple regions
  {
    F.region<-matrix(NA,nrow=sum(nregions),ncol=nages)
    for(i in 1:sum(nregions))
    {
      F.region[i,]<-out$F[nyrs*i,] # fill in the age-specific F for each region
    }
  t<-matrix(c(as.vector(t(F.region[,])),rbind(bio_region[,1]),rbind(bio_region[,nyrs]),rbind(bio_stock[1]),rbind(bio_stock[nyrs]),
              bio_total[1],bio_total[nyrs], rbind(yield_region[,1]),rbind(yield_region[,nyrs]),rbind(yield_stock[1]),
              rbind(yield_stock[nyrs]),yield_total[1],yield_total[nyrs], rbind(u_region[,1]),rbind(u_region[,nyrs]),
              rbind(u_stock[1]),rbind(u_stock[nyrs]),u_total[1],u_total[nyrs],rbind(dep_region[,1]),rbind(dep_region[,nyrs]),
              rbind(dep_stock[1]),rbind(dep_stock[nyrs]),dep_total[1],dep_total[nyrs],rbind(SSB_region[,1]),rbind(SSB_region[,nyrs]),
              rbind(SSB_stock[1]),rbind(SSB_stock[nyrs]),SSB_total[1],SSB_total[nyrs],Bratio_stock[1],Bratio_stock[nyrs],
              Bratio_total[1],Bratio_total[nyrs],rbind(bio_region_overlap[,1]),rbind(bio_region_overlap[,nyrs]),
              rbind(bio_stock_overlap[1]),rbind(bio_stock_overlap[nyrs]),bio_total_overlap[1],bio_total_overlap[nyrs],
              rbind(yield_region_overlap[,1]),rbind(yield_region_overlap[,nyrs]),rbind(yield_stock_overlap[1]),
              rbind(yield_stock_overlap[nyrs]),yield_total_overlap[1],yield_total_overlap[nyrs], rbind(u_region_overlap[,1]),
              rbind(u_region_overlap[,nyrs]),rbind(u_stock_overlap[1]),rbind(u_stock_overlap[nyrs]),
              u_total_overlap[1],u_total_overlap[nyrs], rbind(dep_region_overlap[,1]),rbind(dep_region_overlap[,nyrs]),
              rbind(dep_stock_overlap[1]),rbind(dep_stock_overlap[nyrs]),dep_total_overlap[1],dep_total_overlap[nyrs],
              rbind(SSB_region_overlap[,1]),rbind(SSB_region_overlap[,nyrs]),rbind(SSB_stock_overlap[1]),rbind(SSB_stock_overlap[nyrs]),
              SSB_total_overlap[1],SSB_total_overlap[nyrs],Bratio_stock_overlap[1],Bratio_stock_overlap[nyrs],Bratio_natal_overlap[1],
              Bratio_natal_overlap[nyrs] 
              ),nrow=1)
  }else # if multiple stock with multiple regions
    {
  if(sum(nstocks)>1  & sum(nregions)>1)
    {
    F.region<-matrix(NA,nrow=sum(nregions),ncol=nages)
    for(i in 1:sum(nregions))
        {
         F.region[i,]<-out$F[nyrs*i,]
        }
    
     t<- matrix(c(as.vector(t(F.region[,])),rbind(bio_region[,1]),rbind(bio_region[,nyrs]),rbind(bio_stock[,1]),rbind(bio_stock[,nyrs]),
              bio_total[1],bio_total[nyrs], rbind(yield_region[,1]),rbind(yield_region[,nyrs]),rbind(yield_stock[,1]),
              rbind(yield_stock[,nyrs]),yield_total[1],yield_total[nyrs], rbind(u_region[,1]),rbind(u_region[,nyrs]),
              rbind(u_stock[,1]),rbind(u_stock[,nyrs]),u_total[1],u_total[nyrs], rbind(dep_region[,1]),rbind(dep_region[,nyrs]),
              rbind(dep_stock[,1]),rbind(dep_stock[,nyrs]),dep_total[1],dep_total[nyrs], rbind(SSB_region[,1]),
              rbind(SSB_region[,nyrs]),rbind(SSB_stock[,1]),rbind(SSB_stock[,nyrs]),SSB_total[1],SSB_total[nyrs],Bratio_stock[,1],
              Bratio_stock[,nyrs],Bratio_total[1],Bratio_total[nyrs], rbind(bio_region_overlap[,1]),rbind(bio_region_overlap[,nyrs]),
              rbind(bio_stock_overlap[,1]),rbind(bio_stock_overlap[,nyrs]),bio_total_overlap[,1],bio_total_overlap[,nyrs], 
              rbind(yield_region_overlap[,1]),rbind(yield_region_overlap[,nyrs]),rbind(yield_stock_overlap[,1]),
              rbind(yield_stock_overlap[,nyrs]),yield_total_overlap[,1],yield_total_overlap[,nyrs],rbind(u_region_overlap[,1]),
              rbind(u_region_overlap[,nyrs]),rbind(u_stock_overlap[,1]),rbind(u_stock_overlap[,nyrs]),u_total_overlap[,1],
              u_total_overlap[,nyrs], rbind(dep_region_overlap[,1]),rbind(dep_region_overlap[,nyrs]),rbind(dep_stock_overlap[,1]),
              rbind(dep_stock_overlap[,nyrs]),dep_total_overlap[,1],dep_total_overlap[,nyrs], rbind(SSB_region_overlap[,1]),
              rbind(SSB_region_overlap[,nyrs]),rbind(SSB_stock_overlap[,1]),rbind(SSB_stock_overlap[,nyrs]),SSB_total_overlap[,1],
              SSB_total_overlap[,nyrs],Bratio_stock_overlap[,1],Bratio_stock_overlap[,nyrs],Bratio_natal_overlap[,1],
              Bratio_natal_overlap[,nyrs]
              ),nrow=1)
     }
    }
   }

  
  
#########FIND OUT WHERE THESE COME FROM##########################
# load in values from the assumed MSY run 
assumed.msy.file<-read.csv(file=paste0(wd,"/MSY Assumed.csv",sep=""),nrows=1)
assumed.msy.file<-assumed.msy.file[4:length(assumed.msy.file)]

# load in values from the True MSY run
true.msy.file<-read.csv(file=paste0(wd,"/MSY True.csv",sep=""),nrows=1)
true.msy.file<-true.msy.file[4:length(true.msy.file)]

###################################################################


# fix up the col names for t

if(sum(nstocks)<2  & sum(nregions)>1) # for 1 stock with muliple regions
{
  colnames(t)=c(rep("F",times=(sum(nregions)*nages),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),
                    rep("bio_stock_st",times=sum(nstocks)), rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",
                    rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                    rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),"yield_tot_st",
                    "yield_tot_end",rep("u_reg_st",times=sum(nregions)),rep("u_reg_end",times=sum(nregions)),
                    rep("u_stock_st",times=sum(nstocks)), rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end",
                    rep("dep_reg_st",times=sum(nregions)), rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)),
                    rep("dep_stock_end",times=sum(nstocks)),"dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)),
                    rep("SSB_reg_end",times=sum(nregions)),rep("SSB_stock_st",times=sum(nstocks)), rep("SSB_stock_end",times=sum(nstocks)),
                    "SSB_tot_st","SSB_tot_end",rep("Bratio_stock_st",times=sum(nstocks)), rep("Bratio_stock_end",times=sum(nstocks)),
                    "Bratio_tot_st","Bratio_tot_end", rep("bio_reg_st_over",times=(sum(nregions))),rep("bio_reg_end_over",times=(sum(nregions))),
                    rep("bio_stock_st_over",times=(sum(nstocks))^2), rep("bio_stock_end_over",times=(sum(nstocks))^2),
                    rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)), 
                    rep("yield_reg_st_over",times=(sum(nregions))),rep("yield_reg_end_over",times=(sum(nregions))), 
                    rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),
                    rep("yield_natal_st_over",times=sum(nstocks)), rep("yield_natal_end_over",times=sum(nstocks)),
                    rep("u_reg_st_over",times=(sum(nregions))), rep("u_reg_end_over",times=(sum(nregions))),
                    rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                    rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                    rep("dep_reg_st_over",times=(sum(nregions))), rep("dep_reg_end_over",times=(sum(nregions))),
                    rep("dep_stock_st_over",times=(sum(nstocks))^2), rep("dep_stock_end_over",times=(sum(nstocks))^2),
                    rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),
                    rep("SSB_reg_st_over",times=(sum(nregions))), rep("SSB_reg_end_over",times=(sum(nregions))),
                    rep("SSB_stock_st_over",times=(sum(nstocks))^2), rep("SSB_stock_end_over",times=(sum(nstocks))^2),
                    rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)), 
                    rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),
                    rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks))))
                }    
if(sum(nregions)==1  | sum(nstocks)>1) # if more than on region or greater than one stock
{
  colnames(t)=c(rep("F",times=(sum(nregions)*nages)),rep("bio_reg_st",times=sum(nregions)),rep("bio_reg_end",times=sum(nregions)),
                    rep("bio_stock_st",times=sum(nstocks)), rep("bio_stock_end",times=sum(nstocks)),"bio_tot_st","bio_tot_end",
                    rep("yield_reg_st",times=sum(nregions)),rep("yield_reg_end",times=sum(nregions)),
                    rep("yield_stock_st",times=sum(nstocks)), rep("yield_stock_end",times=sum(nstocks)),
                    "yield_tot_st","yield_tot_end",rep("u_reg_st",times=sum(nregions)), rep("u_reg_end",times=sum(nregions)),
                    rep("u_stock_st",times=sum(nstocks)), rep("u_stock_end",times=sum(nstocks)),"u_tot_st","u_tot_end",
                    rep("dep_reg_st",times=sum(nregions)), rep("dep_reg_end",times=sum(nregions)),rep("dep_stock_st",times=sum(nstocks)),
                    rep("dep_stock_end",times=sum(nstocks)),
                    "dep_tot_st","dep_tot_end",rep("SSB_reg_st",times=sum(nregions)), rep("SSB_reg_end",times=sum(nregions)),
                    rep("SSB_stock_st",times=sum(nstocks)), rep("SSB_stock_end",times=sum(nstocks)),"SSB_tot_st","SSB_tot_end",
                    rep("Bratio_stock_st",times=sum(nstocks)), rep("Bratio_stock_end",times=sum(nstocks)),"Bratio_tot_st","Bratio_tot_end",
                    rep("bio_reg_st_over",times=(sum(nregions))^2),rep("bio_reg_end_over",times=(sum(nregions))^2),
                    rep("bio_stock_st_over",times=(sum(nstocks))^2), rep("bio_stock_end_over",times=(sum(nstocks))^2),
                    rep("bio_natal_st_over",times=sum(nstocks)),rep("bio_natal_end_over",times=sum(nstocks)),
                    rep("yield_reg_st_over",times=(sum(nregions))^2),rep("yield_reg_end_over",times=(sum(nregions))^2),
                    rep("yield_stock_st_over",times=(sum(nstocks))^2), rep("yield_stock_end_over",times=(sum(nstocks))^2),
                    rep("yield_natal_st_over",times=sum(nstocks)), rep("yield_natal_end_over",times=sum(nstocks)),
                    rep("u_reg_st_over",times=(sum(nregions))^2), rep("u_reg_end_over",times=(sum(nregions))^2),
                    rep("u_stock_st_over",times=(sum(nstocks))^2), rep("u_stock_end_over",times=(sum(nstocks))^2),
                    rep("u_natal_st_over",times=sum(nstocks)), rep("u_natal_end_over",times=sum(nstocks)),
                    rep("dep_reg_st_over",times=(sum(nregions))^2), rep("dep_reg_end_over",times=(sum(nregions))^2),
                    rep("dep_stock_st_over",times=(sum(nstocks))^2), rep("dep_stock_end_over",times=(sum(nstocks))^2),
                    rep("dep_natal_st_over",times=sum(nstocks)), rep("dep_natal_end_over",times=sum(nstocks)),
                    rep("SSB_reg_st_over",times=(sum(nregions))^2), rep("SSB_reg_end_over",times=(sum(nregions))^2),
                    rep("SSB_stock_st_over",times=(sum(nstocks))^2), rep("SSB_stock_end_over",times=(sum(nstocks))^2),
                    rep("SSB_natal_st_over",times=sum(nstocks)), rep("SSB_natal_end_over",times=sum(nstocks)),
                    rep("Bratio_stock_st_over",times=(sum(nstocks))^2), rep("Bratio_stock_end_over",times=(sum(nstocks))^2),
                    rep("Bratio_natal_st_over",times=sum(nstocks)), rep("Bratio_natal_end_over",times=sum(nstocks)))
} 


# ready for plots!
write.csv(t, file=paste0(wd,"/MSY Results/Output Quantities.csv",sep=""))
write.csv(true.msy.file, file=paste0(wd,"/MSY Results/True MSY.csv",sep=""))
write.csv(assumed.msy.file, file=paste0(wd,"/MSY Results/Assumed MSY.csv",sep=""))
