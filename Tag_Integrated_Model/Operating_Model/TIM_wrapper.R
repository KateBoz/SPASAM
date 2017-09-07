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
OM_folder<-"E:\\SPASAM-master\\OM_EM_tests\\OM_SABLEFISH"

#set the name of the OM
OM_name<-"TIM_OM"

#run the OM and get a report file
setwd(OM_folder)

# use if need to compile
#invisible(shell(paste0("adcomp ",OM_name,sep=""),wait=T))

#run the OM
#invisible(shell(paste0(OM_name," -nohess",sep=""),wait=T))

out=PBSmodelling::readList(paste0(OM_name,".rep", sep=""))

names(out)

#pulling out core specifications for formatting the data for EM
nages<-out$nages #ages
nyrs<-out$nyrs #simulation years
npops<-out$npops #populations
nregions<-out$nregions #regions - thankfully only 1 per pop
nfleets<-out$nfleets #nfleets fishery
nfleets_survey<-out$nfleets_survey #nfleets survey (usually 1)
tag_rel_years<-out$years_of_tag_releases # model years with tag releases
n_rel<-length(tag_rel_years) #number of releases
max_life_tags<-out$max_life_tags # number of years tally tags after release - need to output this in OM

#setting everything as a .csv for now to check dims and values

#################
# REC DATA
#################
#rec_index
write.csv(out$rec_index_BM,"rec_index_BM.csv")

#################
# Mortality 
################
#create a matrix with dimensions (1,np,1,nreg,1,ny,1,na)
write.csv(out$M,"M_matrix.csv")

#################
#SURVEY DATA
#################
#obs_survey_fleet
OBS_survey_region_bio<-t(out$OBS_survey_region_bio)
write.csv(OBS_survey_region_bio,"OBS_survey_bio.csv")

#setting up the obs_survey_fleet_se
#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_survey_region_bio_se<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_survey_region_bio_se[,i]<-out$sigma_survey
}
write.csv(OBS_survey_region_bio_se,"OBS_survey_region_bio_se.csv")

#Setting up the age compositions. Needs to be (1,np,1,nreg,1,ny,1,nfs,1,na). 
temp<-grep("alt_OBS_survey_prop", names(out), value = TRUE) #pulling out the age_comps
temp2<-out[temp] 
OBS_survey_prop<-do.call("rbind", temp2)
write.csv(OBS_survey_prop,"OBS_survey_prop.csv")

#setting up the survey_sample N
#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_survey_prop_N<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
OBS_survey_prop_N[,i]<-out$survey_sampleN
}
write.csv(OBS_survey_prop_N,"OBS_survey_prop_N.csv")

#####################
# FISHERY CATCH DATA
#####################

#yield catch matrix
OBS_yield_region<-t(out$OBS_yield_region)
write.csv(OBS_yield_region,"OBS_yield_fleet.csv")

#yield_se
OBS_yield_fleet_se<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_yield_fleet_se[,i]<-out$sigma_catch
}
write.csv(OBS_yield_fleet_se,"OBS_yield_fleet_se.csv")

#fishery age compositions
temp<-grep("alt_OBS_catch_prop", names(out), value = TRUE) #pulling out the age_comps
temp2<-out[temp] 
OBS_catch_at_age_fleet_prop<-do.call("rbind", temp2)
write.csv(OBS_catch_at_age_fleet_prop,"OBS_catch_at_age_fleet_prop.csv")

#create a matrix with dimensions (1,np,1,nreg,1,ny,1,nfs)
OBS_catch_prop_N<-matrix(NA,(npops*nregions*nfleets),nyrs)
for(i in 1:nyrs){
  OBS_catch_prop_N[,i]<-out$fishery_sampleN
}
write.csv(OBS_catch_prop_N,"OBS_catch_prop_N.csv")

#####################
# TAGGING DATA
#####################

#ntags (1,np,1,nreg,1,ny_rel,1,na) #number of tags released in each year for each age
ntags<-out$ntags_matrix # this will be tricky
#not generallized #SET UP FOR 3 POPS/AREA
#temp1<-(ntags[1:years_tally])
#temp2<-(ntags[(years_tally+1):(years_tally*2),])
#temp3<-(ntags[((years_tally*2)+1):(years_tally*3),])
#ntags_matrix<-(rbind(temp1,temp2,temp3))

write.csv(ntags,"ntags_matrix.csv")

#backup plan  = set number of tags manually
ntags_matrix2<-ntags
ntags_matrix2[,]<-200 #can fill in with whatever you want! Or use out$SIM_ntag
write.csv(ntags_matrix2,"ntags_matrix_ESS.csv")

##########################
#Getting Observed tag recaps
##########################
#ESS matrix 
temp<-grep("alt_SIM_tag_prop", names(out), value = TRUE) #pulling out the recaptures
temp2<-out[temp] #pulling all the recapture

#proportion - full matrix
temp1<-grep("OBS_tag_prop_final", names(out), value = TRUE) #pulling out the recaptures
temp1.1<-out[temp1] #pulling all the recapture
recaps<-do.call("rbind", temp1.1) # if many populations
write.csv(recaps,"tag_recaps.csv")


temp2<-grep("OBS_tag_prop_reg", names(out), value = TRUE)
length(temp2)
temp2.1<-out[temp2] #pulling all the recaptures


#move it into an array for re-arranging
index<-n_rel*nregions #length
test.array<-array(NA,dim=c(nages,max_life_tags*nregions+1,index)) #3d array 

#fill in array for further indexing
  for(i in 1:index){
    test.array[,,i]<-temp2.1[[i]]
    
  }

#pulling only the recaps
 test.array<-test.array[,1:15,]
 dim(test.array)
 
 #split the array into three separate
 
 recap1<-test.array[,1:5,] #recoveries for reg1
 recap2<-test.array[,6:10,] #recoveries in reg2
 recap3<-test.array[,11:15,] #recoveries in reg3
 
 new.array1<-array(NA,dim=c(max_life_tags,1,index,nages)) # think I have it...now filling it
 
 for(i in 1:max_life_tags){
   for(j in 1:index){
     for(k in 1:nages){
        new.array1[i,1,j,k]<-recap1[k,i,j]
 
 }}}
 
 new.array2<-array(NA,dim=c(max_life_tags,1,index,nages)) # think I have it...now filling it
 
 for(i in 1:max_life_tags){
   for(j in 1:index){
     for(k in 1:nages){
       new.array2[i,1,j,k]<-recap2[k,i,j]
       
     }}}
 
 new.array3<-array(NA,dim=c(max_life_tags,1,index,nages)) # think I have it...now filling it
 
 for(i in 1:max_life_tags){
   for(j in 1:index){
     for(k in 1:nages){
       new.array3[i,1,j,k]<-recap3[k,i,j]
       
     }}}
 
#combine arrays
tag.matrix.EM<-data.frame(reg1=as.vector(new.array1),reg2=as.vector(new.array2),reg3=as.vector(new.array3))
 
write.csv(tag.matrix.EM,"obs_tag_array.csv")



# getting true recaps from 7 d matrix from the report file
recaps_all<-grep("recaps", names(out), value = TRUE) #pulling out the recaptures
temp_recaps<-out[recaps_all[[1]]]


######################
# ADDITIONAL PARAMS  #
######################

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


