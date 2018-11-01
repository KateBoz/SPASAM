#################################################
# Wrapper for running the improved MSY search code
# Created by: Katelyn Bosley
# Date:5/20/2017
#################################################

# remove previous objects from workspace
rm(list = ls())


#load libraries
library(PBSmodelling)

#set wd
wd<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\MY DOCS\\SPASAM STUFF_MS1\\MS1_results\\HAKE\\Move_high_residency"

setwd(wd)

#set up the data save

# select the population type
# 1 - panmictic
# 2 - multiple area
# 3 - metapop
# 4 - natal homing
pop.type<-2



###################################################################################
{ #run the code

#read in the info for the population
update=readLines("GreatDana_MSY_search.dat",n=-1)

nyrs<-as.numeric(update[(grep("nyrs",update)+1)])
nstocks<-as.numeric(update[(grep("npopulations",update)+1)])
nregions<-as.numeric(update[(grep("nregions",update)+1):(grep("nregions",update)+nstocks)])

#need to adjust by number of populations
nfleets<-matrix(as.numeric(update[(grep("nfleets",update)+1):(grep("nfleets",update)+sum(nstocks))],ncol=nstocks))

####################

#panmictic
if(pop.type==1){
  #picking out only the values I want.
  msy_results<-data.frame(matrix(NA,nrow = 1,ncol=10))
  names(msy_results)<-c("Model","F","biomass_total_start","biomass_total_end",
                        "yield_total","harvest_rate_total_bio","depletion_total",
                        "SSB_total_start","SSB_total_end","Bratio_total")
}


# if using the multiple area model
if(pop.type==2){
  
  N_par_reg<-6 # number of parameters with regional values # need to fix this up..
  N_par_pop<-9 # number of parameters for stock
  
  #picking out only the values I want.
  msy_results<-data.frame(matrix(NA,nrow = 1,ncol=(N_par_reg*nregions)+N_par_pop)) # number of parameters with multiple regions+number of total pop values
  
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
  
  #picking out only the values I want.
  msy_results<-data.frame(matrix(NA,nrow = 1,ncol=((N_par_reg*nstocks)+N_par_pop))) # number of parameters with multiple regions+number of total pop values
  
  #fill in the table slowly for error checking
  names(msy_results)[1]<-c("Model")
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

#msy_results



#run the model
invisible(shell("GreatDana_MSY_search -nohess",wait=T))

#invisible(shell("Spatial_BRP",wait=T))

out=PBSmodelling::readList("GreatDana_MSY_search.rep")






#store results to a full spreadsheet

#panmictic
if(pop.type==1){
temp<-c(1,out$F_est,out$biomass_total[1],
        out$biomass_total[nyrs],out$yield_total[nyrs],
        out$harvest_rate_total_bio[nyrs],out$depletion_total[nyrs],
        out$SSB_total[1],out$SSB_total[nyrs], out$Bratio_total[nyrs])
}


# multi-area  
if(pop.type==2) {
  temp<-c(1,out$F_est,
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
  temp<-c(1,out$F_est,
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


msy_results[1,]<-temp

write.csv(msy_results,"true_msy_greatdana.csv")

#generalized output
write.csv(out$biomass_AM_age[seq(nyrs,nrow(out$biomass_AM_age),nyrs),], "biomass_at_age.csv")
write.csv(out$abundance_at_age_AM[seq(nyrs,nrow(out$abundance_at_age_AM),nyrs),], "abund_at_age.csv")

                             
} #end code


