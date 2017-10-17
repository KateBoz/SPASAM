####################################################
# Simple initial code for visuallizing model outputs
# Created by JJD/KB
####################################################

rm(list=(ls()))


#load libraries
load_libraries<-function() {
  library(PBSmodelling)
}
load_libraries()


###inputs for running models
# Manually make changes in the OM .dat and run both OM and EM together if you want

#OM Location
OM_direct<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SPASAM-master\\Tag_Integrated_Model\\simple model\\OM"
OM_name<-"TIM_OM_TEST1" #name of the OM you are wanting to run

#EM Location
EM_direct<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SPASAM-master\\Tag_Integrated_Model\\simple model\\EM" #location of run(s)
EM_name<-"tim_hakeish" ###name of .dat, .tpl., .rep, etc.


#run the OM
setwd(OM_direct)
invisible(shell(paste0(OM_name," -nohess"),wait=T))

# move the .rep from the OM_direct and change name. 
from<-paste0(OM_direct,"\\",OM_name,".rep")
to<-paste0(EM_direct,"\\",EM_name,".dat")

file.remove(to) #remove old version if present
file.copy(from = from,  to = to)

#run the EM
setwd(EM_direct)
invisible(shell(paste0(EM_name," -nohess"),wait=T))

#########################################################
# Look at the outputs
#########################################################

#Read in model .rep
out<-readList(paste(EM_direct,paste0(EM_name,".rep"),sep="\\")) #read in .rep file


#pull info about the model
na<-out$nages
yrs<-out$nyrs
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)

#recruitment for two area
plot(years,out$recruits_BM[1,], type = "l", lwd = 2, ylab = "Recruits", ylim=c(0,max(out$recruits_BM)))
points(years,out$recruits_BM[2,], type = "l", col = "blue",lwd = 2)
##points(years,out$recruits_BM[3,], type = "l", col = "red",lwd = 2)

points(years,out$recruits_BM_TRUE[1,], type = "l", lty=2)
points(years,out$recruits_BM_TRUE[2,], type = "l", col = "blue", lty = 2)
#points(years,out$recruits_BM_TRUE[3,], type = "l", col = "red",lty = 2)

legend("topright",legend = c("Estimated","True"),lty = c(1,2), lwd = 2)

#SSB
plot(years,out$SSB_region[1,], type = "l", lwd = 2, ylab = "SSB", ylim=c(min(out$SSB_region),max(out$SSB_region)))
points(years,out$SSB_region[2,], type = "l", col = "blue",lwd = 2)

points(years,out$SSB_region_TRUE[1,], type = "l", lty=2)
points(years,out$SSB_region_TRUE[2,], type = "l", col = "blue", lty = 2)
legend("topright",legend = c("Estimated","True"),lty = c(1,2), lwd = 2)

#q
out$q_survey
out$q_survey_TRUE

#sel params
out$sel_beta1_survey
out$sel_beta1_TRUE

out$sel_beta2_survey
out$sel_beta2_TRUE


#F
F.temp<-out$F[1:yrs,]
F.mean<-colMeans(F.temp)
plot(ages,F.mean, type = "l", col = "black", ylim=c(0,1), lwd = 2)

#true
F.true<-colMeans(out$F_TRUE[1:yrs,])
points(ages,F.true, type = "b", col = "black", lwd = 1, lty = 2)


#additional plots

output=out

#function to make time series plots with truth as solid black and estimate(s) as dashed
tsplot<-function(true=NULL,est=NULL,plotname=NULL){
  plot(true,lwd=2,type='l',col="black",ylab=plotname,xlab="Year")
  lines(est,lwd=2,lty=2,col="black")
}

Fs<-matrix(NA,output$nregions,output$nyrs) #to hold fully selected Fs
TrueFs<-matrix(NA,output$nregions,output$nyrs) #to hold fully selected Fs
for(r in 1:output$nregions) {
  #Recruits
  esti<-output$recruits_BM[r,]
  true.b<-output$recruits_BM_TRUE[r,]
  tsplot(true=true.b,est=esti,plotname="Recruitment")
  #SSB
  esti<-output$SSB_region[r,]
  true.b<-output$SSB_region_TRUE[r,]
  tsplot(true=true.b,est=esti,plotname="SSB")
  #Fully selected F
  for(y in 1:output$nyrs){
  y2<-y+(output$nyrs*(r-1)) #finding the start point of each region
  Fs[r,y]<-max(output$F[y2,])
  TrueFs[r,y]<-max(output$F_TRUE[y2,])
  } #end y loop
  esti<-Fs[r,]
  true.b<-TrueFs[r,]
  tsplot(true=true.b,est=esti,plotname="Fully Selected F")
} #end r loop
