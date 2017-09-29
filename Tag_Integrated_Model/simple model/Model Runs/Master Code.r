
rm(list=(ls()))
wd<-getwd()
setwd(wd)

############### INPUTS ##############################################################################################################################

nloops <-50                                                 # How many simulations to perform

#--------------------------------------------------------------------------------------------------------------------------------------------------
suppressWarnings(suppressMessages(require(PBSmodelling)))
suppressWarnings(suppressMessages(require(matrixStats)))
suppressWarnings(suppressMessages(require(TeachingDemos)))
suppressWarnings(suppressMessages(require(snowfall)))
suppressWarnings(suppressMessages(library(parallel)))
suppressWarnings(suppressMessages(library('snow')))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doSNOW)))
suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(gtools)))
suppressWarnings(suppressMessages(library(spatstat)))
suppressWarnings(suppressMessages(library(alphahull)))
suppressWarnings(suppressMessages(library(beanplot)))
#--------------------------------------------------------------------------------------------------------------------------------------------------


##########################################################################################################################################
############ Run Fmsy Simulation ###########################################################################################################

update1=readLines("SIM_TAC.dat",n=-1)


run.SIM<-function(ntrial,WD1,WD)
{

  dir.create(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  dir.create(paste0(WD,"/Assessment Results/Run",ntrial,sep=""))
  invisible(file.copy(from=paste0(WD1,"/TIM_OM.exe",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/TIM_OM.exe",sep="")))
  invisible(file.copy(from=paste0(WD1,"/TIM_OM.tpl",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/TIM_OM.tpl",sep="")))
  invisible(file.copy(from=paste0(WD1,"/TIM_OM.dat",sep=""),to=paste0(WD,"/Simulation Results/Run",ntrial,"/TIM_OM.dat",sep="")))
  
  setwd(paste0(WD,"/Simulation Results/Run",ntrial,sep=""))
  
  SIM.DAT=readLines("TIM_OM.dat",n=-1)
  SIM.DAT[(grep("myseed_yield",SIM.DAT)+1)]=ntrial
  SIM.DAT[(grep("myseed_survey",SIM.DAT)+1)]=100000+ntrial
  SIM.DAT[(grep("myseed_F",SIM.DAT)+1)]=200000+ntrial
  SIM.DAT[(grep("myseed_rec_devs",SIM.DAT)+1)]=300000+ntrial
  SIM.DAT[(grep("myseed_rec_apport",SIM.DAT)+1)]=400000+ntrial
  SIM.DAT[(grep("myseed_rec_index",SIM.DAT)+1)]=500000+ntrial
  SIM.DAT[(grep("myseed_survey_age",SIM.DAT)+1)]=600000+ntrial
  SIM.DAT[(grep("myseed_catch_age",SIM.DAT)+1)]=700000+ntrial
  SIM.DAT[(grep("myseed_tag",SIM.DAT)+1)]=800000+ntrial
  
  writeLines(SIM.DAT,"TIM_OM.dat")
  
  system(paste0("TIM_OM -nox -nohess TIM_OM.dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen
  
  invisible(shell(paste("copy TIM_OM.par TIM_OM",ntrial,".par", sep="")))
  invisible(shell(paste("copy TIM_OM.rep TIM_OM",ntrial,".rep", sep="")))

  invisible(file.copy(from=paste0(WD,"/Simulation Results/Run",ntrial,"/TIM_OM",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".dat",sep="")))

  invisible(shell(paste("del TIM_OM.par", sep="")))
  invisible(shell(paste("del TIM_OM.rep", sep="")))
  invisible(shell("del TIM_OM.exe")) 
  invisible(shell(paste("del TIM_OM.bar",sep="")))
  invisible(shell(paste("del TIM_OM.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
  
}

run.est<-function(ntrial,WD1,WD)
{

  invisible(file.copy(from=paste0(WD1,"/tim_sab.exe",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab.exe",sep="")))
  invisible(file.copy(from=paste0(WD1,"/tim_sab.tpl",sep=""),to=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab.tpl",sep="")))

  setwd(paste0(WD,"/Assessment Results/Run",ntrial,sep=""))

  system(paste0("tim_sab -nox -ind tim_sab",ntrial,".dat",sep=""),wait=TRUE,show.output.on.console=FALSE)  #-nox keeps from showing the ADMB gradient outputs on screen
  
  invisible(shell(paste("copy tim_sab.par tim_sab",ntrial,".par", sep="")))
  invisible(shell(paste("copy tim_sab.rep tim_sab",ntrial,".rep", sep="")))

  if(file.exists("tim_sab.cor")==TRUE)
  {
    print(paste0("Model Run ",ntrial," Converged"))
    invisible(shell(paste("copy tim_sab.std tim_sab",ntrial,".std", sep="")))
    invisible(shell(paste("copy tim_sab.cor tim_sab",ntrial,".cor", sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Report Files/tim_sab",ntrial,".rep",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".cor",sep=""),to=paste0(WD,"/Assessment Results/Report Files/tim_sab",ntrial,".cor",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".std",sep=""),to=paste0(WD,"/Assessment Results/Report Files/tim_sab",ntrial,".std",sep="")))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".par",sep=""),to=paste0(WD,"/Assessment Results/Report Files/tim_sab",ntrial,".par",sep="")))
    invisible(shell(paste("del tim_sab.std", sep="")))
    invisible(shell(paste("del tim_sab.cor", sep="")))  
    invisible(shell("del admodel.cov"))
  }
  if(file.exists("tim_sab.cor")==FALSE)
  {
    print(paste0("Model Run ",ntrial," Did Not Converge"))
    invisible(file.copy(from=paste0(WD,"/Assessment Results/Run",ntrial,"/tim_sab",ntrial,".rep",sep=""),to=paste0(WD,"/Assessment Results/Report Files/tim_sab",ntrial,".rep",sep="")))
  }

  
  invisible(shell(paste("del tim_sab.par", sep="")))
  invisible(shell(paste("del tim_sab.rep", sep="")))
  invisible(shell("del tim_sab.exe")) 
  invisible(shell("del admodel.dep"))
  invisible(shell("del admodel.hes")) 
  invisible(shell(paste("del tim_sab.b01",sep="")))
  invisible(shell(paste("del tim_sab.P01",sep="")))
  invisible(shell(paste("del tim_sab.R01",sep="")))
  invisible(shell(paste("del tim_sab.b02",sep="")))
  invisible(shell(paste("del tim_sab.P02",sep="")))
  invisible(shell(paste("del tim_sab.R02",sep="")))
  invisible(shell(paste("del tim_sab.b03",sep="")))
  invisible(shell(paste("del tim_sab.P03",sep="")))
  invisible(shell(paste("del tim_sab.R03",sep="")))
  invisible(shell(paste("del tim_sab.bar",sep="")))
  invisible(shell(paste("del tim_sab.eva",sep="")))
  invisible(shell(paste("del tim_sab.log",sep="")))
  invisible(shell("del fmin.log"))
  invisible(shell("del variance"))
  
}

plot.beans.combined<-function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(as.numeric(unlist(data.frame.name)),bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'), #ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Percent Bias"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  abline(h=median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),lty=2,col="red",lwd=2)
  legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),3),"%"),
                  paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
 
}   

plot.beans=function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(data.frame.name,bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),#ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Percent Bias"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  abline(h=median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),lty=2,col="red",lwd=2)
  legend("top", c(paste0("Median Bias ", signif(median(as.numeric(unlist(data.frame.name)),na.rm=TRUE),3),"%"),
                  paste0("Max Bias ", signif(max(abs(as.numeric(unlist(data.frame.name))),na.rm=TRUE),3),"%")),lty=c(4,-1),cex=.5,ncol=2,bg='white')
  
}   

plot.beans.value=function(parameter.name,data.frame.name,Converge,percent.converge,...)
{
  par(mfrow=c(1,1),  mgp=c(2.5,.6,.5), mar=c(5, 4, 3, .5) + 0.1)
  
  beanplot(data.frame.name,bw="nrd0",boxwex=1.2,what=c(0,1,0,0),col=c('grey','black','black','white'),#ylim=c(min(data.frame.name,na.rm=TRUE),max(data.frame.name,na.rm=TRUE)),
           beanlines="median",main=paste0(parameter.name," Observed and Predicted Value"),  mgp=c(2.75,.75,.5),
           sub=paste0(Converge, " runs converged out of ", nloops," (", percent.converge,"% Convergence)"),las=1,log="",
           cex.axis=.5,na.rm=TRUE,...)
  abline(h=0,lty=1,lwd=2)
  legend("top",c('True Value'),lty=c(1),pch=1,bg='white')
  
}  

result=function(WD,iteration)
{
  setwd(paste0(WD,"/Assessment Results/Report Files"))
  
  if(file.exists(paste0(WD,"/Assessment Results/Report Files/tim_sab",iteration,".cor",sep=""))==TRUE)
  {
    output =readList(paste0(WD,"/Assessment Results/Report Files/tim_sab",iteration,".rep",sep=""))
    par_names=c('SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')  #'Obj_Func','max_gradient',
    results=output[par_names]
    c(rep=iteration,
      converged=1,
      #objective=results$Obj_Func,
      #grad=results$max_gradient,
      #h_est=results$h,
      R0_est=results$R_ave,
      F_est=as.vector(results$F),
      ssb_est=as.vector(results$SSB_region),
      recruits_est=as.vector(results$recruits_BM)
    )
  }
  else
  {
    output =readList(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".rep",sep=""))
    par_names=c('SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
    results=output[par_names]
    c(rep=iteration,
      converged=NA,
      #objective=NA,
      #grad=NA,
      #h_est=NA,
      R0_est=NA,
      F_est=as.vector(rep(NA,times=length(results$fmult))),
      ssb_est=as.vector(rep(NA,times=length(results$ssb))),
      recruits_est=as.vector(rep(NA,times=length(results$recruits)))
    )
  }
}

true_value=function(WD,iteration)
{
  setwd(paste0(WD,"/Assessment Results/Report Files"))
  output =readList(paste0(WD,"/Assessment Results/Report Files/TAC_est",iteration,".rep",sep=""))
  par_names=c('SIM_R_ave','SIM_F','SIM_ssb','SIM_recruits','h','R_ave','fmult','ssb','recruits')
  results=output[par_names]
  c(rep=iteration,
    #h_true=results$SIM_h,
    R0_true=results$R_ave_TRUE,
    F_true=as.vector(results$F_TRUE),
    ssb_true=as.vector(results$SSB_region_TRUE),
    recruits_true=as.vector(results$recruits_BM_TRUE)
  )
}

#calc_bias=function(WD,iteration)
#{
#  est=cbind(sapply(1:iteration,function(i)result(WD,i)))
#  true_temp=cbind(sapply(1:iteration,function(i)true_value(WD,i)))
  
#  Converged<-(iteration-length(which(is.na(est[2,]))))
#  percent.converged<-(Converged/iteration)*100
  
#  est1<-est[-c(1:2),]
#  write.csv(est1,file=paste0(WD,"/Bias/Estimated Parameters.csv",sep=""))
  
#  setwd(paste0(WD,"/Assessment Results/Report Files"))
#  report.file=readList(paste0(WD,"/Assessment Results/Report Files/tim_sab1.rep",sep=""))
#  nyrs<-as.numeric(report.file['nyrs'])
#  nyrs_ASS<-as.numeric(report.file['nyrs'])
  
#  true=true_temp[c(2:(nyrs_ASS+3),(nyrs_ASS+3+1+(nyrs-nyrs_ASS)):(nyrs_ASS+3+1+(nyrs)-1),(nyrs_ASS+3+1+(nyrs)+(nyrs-nyrs_ASS)):(nyrs_ASS+3+1+(nyrs)-1+(nyrs))),]
#  write.csv(true_temp,file=paste0(WD,"/Bias/Simulated Parameters(ALL).csv",sep=""))
#  write.csv(true,file=paste0(WD,"/Bias/Simulated Parameters(Assess).csv",sep=""))
  
#  setwd(WD)
#  bias<-matrix(NA,nrow=length(true[,1]),ncol=iteration)
#  percent_bias<-matrix(NA,nrow=length(true[,1]),ncol=iteration)
  
#  for(i in 1:length(true[,1]))  
#  {
#    for(j in 1:iteration)
#    {
#      if(true[i,j]==0 || is.na(est1[i,j]))
#      {
#        bias[i,j]=NA
#        percent_bias[i,j]=NA
#      }
#      else{
#        bias[i,j]=est1[i,j]-true[i,j]
#        percent_bias[i,j]=bias[i,j]/true[i]*100
#      }
#    }
#  }
  
#  row.names(percent_bias)<-c('Steep','R0',rep('F',times=nyrs_ASS),rep('SSB',times=nyrs_ASS),rep('RECR',times=nyrs_ASS))
#  row.names(bias)<-c('Steep','R0',rep('F',times=nyrs_ASS),rep('SSB',times=nyrs_ASS),rep('RECR',times=nyrs_ASS))
#  write.csv(percent_bias,file=paste0(WD,"/Bias/Percent_Bias.csv",sep=""))
#  write.csv(bias,file=paste0(WD,"/Bias/Bias.csv",sep=""))
  
#  pdf(file=paste0(WD,"/Bias/Percent_Bias.pdf",sep=""))
#  plot.beans.combined('Steepness',percent_bias[1,],Converged,percent.converged,ylab='% Bias')
#  plot.beans.combined('R0',percent_bias[2,],Converged,percent.converged,ylab='% Bias')
#  plot.beans.combined('F (Combined Across Years)',percent_bias[c(3:(nyrs_ASS+2)),],Converged,percent.converged,ylab='% Bias')
#  plot.beans('F (By Year)',as.data.frame(t(percent_bias[c(3:(nyrs_ASS+2)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

#  plot.beans.value('F (By Year)',as.data.frame(t(est1[c(3:(nyrs_ASS+2)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='F')
#  matplot(1:75,true[c(3:(nyrs_ASS+2)),2],type='l',lwd=2,add=T)
#  matpoints(1:75,true[c(3:(nyrs_ASS+2)),2],pch=21,col='black',cex=0.75,bg='white')
  
#  plot.beans.combined('SSB (Combined Across Years)',percent_bias[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),],Converged,percent.converged,ylab='% Bias')
#  plot.beans('SSB (By Year)',as.data.frame(t(percent_bias[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

#  plot.beans.value('SSB (By Year)',as.data.frame(t(est1[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='SSB')
#  matplot(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),2],type='l',lwd=2,add=T)
#  matpoints(1:75,true[c((nyrs_ASS+2+1):(nyrs_ASS+2+1+(nyrs_ASS)-1)),2],pch=21,col='black',cex=0.75,bg='white')
  
#  plot.beans.combined('Recruits (Combined Across Years)',percent_bias[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),],Converged,percent.converged,ylab='% Bias')
#  plot.beans('Recruits (By Year)',as.data.frame(t(percent_bias[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='% Bias')

#  plot.beans.value('Recruits (By Year)',as.data.frame(t(est1[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),])),Converged,percent.converged,names=c(1:75),xlab='Assessment Year',ylab='Recruits')
#  matplot(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),2],type='l',lwd=2,add=T)
#  matpoints(1:75,true[c((nyrs_ASS+2+1+(nyrs_ASS)):(nyrs_ASS+2+1+(nyrs_ASS)+nyrs_ASS-1)),2],pch=21,col='black',cex=0.75,bg='white')
  
#  dev.off()
#}


Model_Outputs=function(WD,iteration)
{
 result(WD,iteration)
 true_value(WD,iteration)
# calc_bias(WD,iteration)
}




###############################################################################################################################################
###############################################################################################################################################
################ Perform Simulations and Estimation ###########################################################################################
###############################################################################################################################################
wd1<-getwd()
setwd(wd1)

  ######### CNST TAC, low F ###############################################################
  
  wd<-paste0(wd1,"/Model Runs",sep="")
  dir.create(paste0(wd,"/Simulation Results",sep=""))
  dir.create(paste0(wd,"/Bias",sep=""))
  dir.create(paste0(wd,"/Assessment Results",sep=""))
  dir.create(paste0(wd,"/Assessment Results/Report Files",sep=""))
  
  cl <- makeSOCKcluster(detectCores()-1)
  clusterExport(cl, c("wd1","wd","run.SIM","run.est"))
  registerDoSNOW(cl)
  pb <- winProgressBar(paste0("Simulation Model Progress Bar (",label,")",sep=""), label=paste("Simulation Run 0 of ",nloops,sep=""),max=100)
  progress<-function(n) setWinProgressBar(pb,(n/nloops*100),label=paste("Simulation Run", n,"of", nloops,"Completed"))
  opts<-list(progress=progress)
  
  t<- foreach(i=1:nloops,.options.snow=opts,.packages=c('PBSmodelling','beanplot','matrixStats','TeachingDemos','snowfall','parallel','tcltk2')
  ) %dopar% {
    run.SIM(i,wd1,wd) 
    run.est(i,wd1,wd) 
  }
  close(pb) 
  stopCluster(cl)
  closeAllConnections()
  output3<-gc()
  
  Model_Outputs(wd,nloops)

  
  