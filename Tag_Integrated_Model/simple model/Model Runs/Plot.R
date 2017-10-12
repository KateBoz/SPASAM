rm(list=(ls()))

###inputs
direct<-"C:\\Spatial_SPASM_2016\\EM_Project2\\SPASAM-master\\SPASAM-master\\Tag_Integrated_Model\\simple model\\EM" #location of run(s)
name<-"tim_hakeish" ###name of .dat, .tpl., .rep, etc.

#packages
require(PBSmodelling)

####
output <-readList(paste(direct,paste0(name,".rep"),sep="\\")) #read in .rep file

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
