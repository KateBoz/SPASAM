###########################################################
# The beginnings of a wrapper for the Tag-integrated OM/EM
#
# Created by Katelyn Bosley
# Date: 8/24/2017
###########################################################

##############
#getting ready
##############

#clean environment
rm(list = ls())

#load libraries
load_libraries<-function() {
 library(PBSmodelling)
 library(data.table)
}
load_libraries()

#####################################
# Run the OM
#####################################

#Set the folder for the OM
OM_folder<-"E:\\SPASAM-master\\Tag_Integrated_Model\\Operating_Model"

#set the name of the OM
OM_name<-"TIM_OM"

#run the OM and get a report file
setwd(OM_folder)


# use if need to compile
#invisible(shell(paste0("adcomp ",OM_name,sep=""),wait=T))

#run the OM
invisible(shell(paste0(OM_name," -nohess",sep=""),wait=T))

out=PBSmodelling::readList(paste0(OM_name,".rep", sep=""))

names(out)

#pulling out core specifications for formatting the data for EM
nages<-out$nages
nyrs<-out$nyrs
npops<-out$npops
nregions<-out$nregions
nfleets<-out$nfleets
nfleets_survey<-out$nfleets_survey
tag_rel_years<-out$years_of_tag_releases

#################
# REC DATA
#################
#rec_index
#write.csv(out2$rec_index_BM,"rec_index_BM.csv")

#################
#SURVEY DATA
#################
#obs_survey_fleet
OBS_survey_region_bio<-t(out$OBS_survey_region_bio)
#write.csv(OBS_survey_region_bio,"OBS_survey_bio.csv")

#setting up the obs_survey_fleet_se
#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_survey_region_bio_se<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_survey_region_bio_se[,i]<-out$sigma_survey
}
#write.csv(OBS_survey_region_bio_se,"OBS_survey_region_bio_se.csv")

#Setting up the age compositions. Needs to be (1,np,1,nreg,1,ny,1,nfs,1,na). 
temp<-grep("alt_OBS_survey_prop", names(out), value = TRUE) #pulling out the age_comps
temp2<-out[temp] 
OBS_survey_prop<-do.call("rbind", temp2)
#write.csv(OBS_survey_prop,"OBS_survey_prop.csv")

#setting up the survey_sample N
#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_survey_prop_N<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
OBS_survey_prop_N[,i]<-out$survey_sampleN
}
#write.csv(OBS_survey_prop_N,"OBS_survey_prop_N.csv")

#####################
# FISHERY CATCH DATA
#####################

#yield catch matrix
OBS_yield_region<-t(out$OBS_yield_region)
#write.csv(OBS_yield_region,"OBS_yield_fleet.csv")

#yield_se
OBS_yield_fleet_se<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_yield_fleet_se[,i]<-out$sigma_catch
}
#write.csv(OBS_yield_fleet_se,"OBS_yield_fleet_se.csv")

#fishery age compositions
temp<-grep("alt_OBS_catch_prop", names(out), value = TRUE) #pulling out the age_comps
temp2<-out[temp] 
OBS_catch_at_age_fleet_prop<-do.call("rbind", temp2)
#write.csv(OBS_catch_at_age_fleet_prop,"OBS_catch_at_age_fleet_prop.csv")

#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_catch_prop_N<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_catch_prop_N[,i]<-out$fishery_sampleN
}
#write.csv(OBS_catch_prop_N,"OBS_catch_prop_N.csv")

#####################
# TAGGING DATA
#####################

#ntags (1,np,1,nreg,1,ny_rel,1,na) #number of tags released in each year for each age
ntags<-out$ntags_matrix # this will be tricky

#not generallized
temp1<-t(ntags[1:30,])
temp2<-t(ntags[31:60,])
temp3<-t(ntags[61:90,])
ntags_matrix<-rbind(temp1,temp2,temp3)
write.csv(ntags_matrix,"ntags_matrix.csv")

#backup plan  = set number of tags manually
ntags_matrix2<-ntags_matrix
ntags_matrix2[,]<-200 #can fill in with whatever you want! Or use out$SIM_ntag
write.csv(ntags_matrix2,"ntags_matrix_ESS.csv")

######################
# ADDITIONAL PARAMS  #
######################


#pick the parameters we want to use in the estimation model in the order we want them in
select<-c(
  'nages',#number of ages
  'nyrs',#number of years
  'npops',# number of populations
  'nregions', #number of regions,
  'nfleets', #number of fleets
  'sigma_survey', #survey error
  'survey_sampleN', #survey N 
  'fishery_sampleN', #fishery N
  'rec_index_BM',
  'OBS_survey_region_bio', #observed survey biomass by fleet - only 1 fleet per region
  'OBS_yield_region', #observed landings (biomass)
  'alt_OBS_survey_prop_1', #survey age comps for pop1, region, fleet and year
  'alt_OBS_survey_prop_2', #survey age comps for pop2, region, fleet and year
  'alt_OBS_catch_prop_1', #fishery age comp for pop1, region, fleet, and year
  'alt_OBS_catch_prop_2',#fishery age comp for pop2, region, fleet, and year
  'nyears_tag_release',
  'years_of_tag_releases',
  'reporting_rate',
  'ntags_total', #total number of tags
  'ntags_matrix'
  ) #recaptures 

out2<-out[select] # pull all the data we want...

#formatting the age_comps by pop_region_fleet_year








# getting recaps from the report file
temp<-grep("recaps", names(out), value = TRUE) #pulling out the recaptures

temp2<-out[temp] #pulling all the recapture
recaps<-do.call("rbind", temp2)

out2[[length(out2)+1]] <- recaps #merge the recapture data

names(out2)[20]<-"recaps"

#all the data...

#save the file for later use
saveRDS(out2,'TM_EM_data.rds')

#load the file to get info off it or write to csv, text, etc
data<-readRDS('TM_EM_data.rds')
names(data)


#######################################################################
#data to include...
#
#Catch by area/fleet
#Age compositions and sample sizes for fishery and for survey by area/fleet
#Survey Abundance and error
#Tagging recovery data by release year, recocvery year, area
#Tag release data by year,  and area
########################################################################

########################################################################
# BRIANS COMMENTS ON THE REPORT OF 5/6-D arrays
#
#Since regions are by population, need to include in population loop
#Thus, need counters to index total number of regions.
#For example, with region, for 3 pops with 2 1 and 1 regions, region1 is region 1 of pop 1, 
#region2 is region 2 of pop 1, region3 is region 1 of pop 2, and region4 is region 1 of pop 3. 
#


#######################################################################
# Testing some params with plots
#######################################################################
years<-1:out$nyrs
na<-1:nages

#rec_index
plot(years,colSums(out2$rec_index_BM))
abline(h=15.78)

#obs_survey_fleet
plot(years,rowSums(out2$OBS_survey_region_bio))
points(years,out$true_survey_total_bio,type="l")

#OBS_survey_prop
#save as .csv for now
#aggregated
plot(na,colMeans(OBS_survey_prop[1:50,]), lwd = 2, type="l", ylim=c(0,0.5), ylab="Mean Proportion",xlab="age")
points(na,colMeans(OBS_survey_prop[51:100,]),lwd=2,  type="l", col=2)
points(na,colMeans(OBS_survey_prop[101:150,]), lwd=2, type="l", col=3)
legend("topright",legend = c("Area 1","Area 2","Area 3"), col = c(1,2,3), bty="n", lty=1, lwd = 2)

#OBS_yield_fleet
plot(years,rowSums(out2$OBS_yield_region))
points(years,out$yield_total,type="l")

#OBS_yield_prop
#save as .csv for now
#aggregated
plot(na,colMeans(OBS_catch_at_age_fleet_prop[1:50,]), lwd = 2, type="l", ylim=c(0,0.5), ylab="Mean Proportion",xlab="age")
points(na,colMeans(OBS_catch_at_age_fleet_prop[51:100,]),lwd=2,  type="l", col=2)
points(na,colMeans(OBS_catch_at_age_fleet_prop[101:150,]), lwd=2, type="l", col=3)
legend("topright",legend = c("Area 1","Area 2","Area 3"), col = c(1,2,3), bty="n", lty=1, lwd = 2)


