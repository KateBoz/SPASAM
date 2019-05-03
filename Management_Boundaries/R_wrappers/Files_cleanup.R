##########################################
# Quick little code to replace em and om files
# and clean run for re-do
# 
# Created by: Katelyn Bosley
# Date 5/3/2019
############################################


rm(list=(ls()))

#what is the model name
OM_Name<-"MB_OM" #name of the OM you are wanting to run
EM_Name<-"MB_EM" ###name of .dat, .tpl., .rep, etc.

#directories for model locations
OM_loc<-"C:\\Users\\katelyn.bosley\\Desktop\\MB_Runs\\Model\\Operating_Model"
EM_loc<-"C:\\Users\\katelyn.bosley\\Desktop\\MB_Runs\\Model\\Estimation_Model"

#directory with runs to be updated
direct_master<-"C:\\Users\\katelyn.bosley\\Desktop\\MB_Runs\\All_Runs"

# Create a list of files in the directory
files<-list.files(direct_master) # these folders in the master will be the individual scenarios 

#seems to be some missing .dat files so I will check
yes_dat<-rep(0,length(files)-1)

#this code chunk deletes all the old .tpl, .dat, .rep, etc files from previous runs and replaces with new files

for (i in 1:(length(files)-1)){
#OM Location

#delete old diagnostic folders
unlink(paste0(direct_master,"\\",files[i],"\\Diagnostics"),recursive = T)

#do we have a .dat?
if(file.exists(paste0(direct_master,"\\",files[i],"\\Operating_Model\\",OM_Name,".dat",sep="")))
{yes_dat[i]=1}

# remove old files and replace with new ones NOT including the .dat files
#deleting all those annoying files
OMdir<-paste0(direct_master,"\\",files[i],"\\Operating_Model")
id <- grep(paste0(OM_Name,".dat"),dir(OMdir)) #save only the .dat file
todelete <- dir(OMdir, full.names = TRUE)[-id]
unlink(todelete)

# add new file
file.copy(from=paste0(OM_loc,"\\",OM_Name,".exe",sep=""),to=paste0(direct_master,"\\",files[i],"\\Operating_Model",sep=""))



#deleting all those EM annoying files
EMdir<-paste0(direct_master,"\\",files[i],"\\Estimation_Model")
todelete <- dir(EMdir, full.names = TRUE)
unlink(todelete)

# add new file
file.copy(from=paste0(EM_loc,"\\",EM_Name,".exe",sep=""),to=paste0(direct_master,"\\",files[i],"\\Estimation_Model",sep=""))

}


#pull out the files that did not have a .dat
#list.files(direct_master)[which(yes_dat==0)]




#save this for later...
#deleting all those annoying files
#EMdir<-paste0(direct_master,"\\",files[i],"\\Estimation_Model")
#id <- c(grep(paste0(EM_Name,".dat"),dir(EMdir)),grep(paste0(EM_Name,".exe"),dir(EMdir)))
#todelete <- dir(EMdir, full.names = TRUE)[-id]
#unlink(todelete)