###########################################################
# The beginnings of a wrapper for the Tag-integrated OM/EM
#
# Created by Katelyn Bosley
# Date: 
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
OM_folder<-"F:\\SPASAM-master\\SPASAM-master\\Tag_Integrated_Model\\Operating_Model"
#set the name of the OM
OM_name<-"CODE_TO_ALTER_Spatial_BRP"

#run the OM and get a report file
setwd(OM_folder)


# use if need to compile
#invisible(shell(paste0("adcomp ",OM_name,sep=""),wait=T))

#run the OM
invisible(shell(paste0(OM_name," -nohess",sep=""),wait=T))

out=PBSmodelling::readList(paste0(OM_name,".rep", sep=""))

names(out)

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

# getting recaps from the report file
temp<-grep("recaps", names(out), value = TRUE) #pulling out the recaptures

temp2<-out[temp] #pulling all the recapture
recaps<-do.call("rbind", temp2)

out2[[length(out2)+1]] <- recaps #merrge the recapture data

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







