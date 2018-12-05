####################################################
# Model diagnostics for TIM OM/EM
# Created by: Katelyn Bosley
####################################################

#use this code when running simulations

#load libraries
load_libraries<-function() {
  library(PBSmodelling)
  library(data.table)
  library(ggplot2)
  library(reshape2)
  library(gridExtra)
  library(gplots)
  library(colorspace)
  library(RColorBrewer)
  library(dplyr)
  library(matrixStats) 
  library(PBSadmb)
  library(gtable)
  library(corrplot)
  }
load_libraries()


######################################################################
########### PLOT SETTINGS ############################################
######################################################################

##############################################
#plot set up 
##############################################
#set up the colors you want to use for TRUE and ESTIMATED
t.col="black" #true
e.col="blue"  #estimated

#set up the color that are wanted for MOVEMENT plot - used a color ramp
mycols=colorRampPalette(c("blue", "cyan","black"))

# Select Line width
line.wd=0.8

#select threshold correlation level for corrlation plot
cor.level = 0.4# can change this in the make.plots function below


##############################################
#Residual values
##############################################

#set residual switch for different values plotted
#resid.switch=2
# =1 straight residual (ESTIMATED-TRUE; not % of true)
# =2 Relative percent difference ((ESTIMATED-TRUE)/TRUE *100)


#################################################################################
########### INPUTS FOR RUNNING MODELS ###########################################
#################################################################################
  
#name OPERATING MODEL of the .exe
OM_name<-OM_name #name of the OM you are wanting to run

#name ESTIMATION MODEL of the .exe
EM_name<-EM_name ###name of .dat, .tpl., .rep, etc.


############################################################################
########### AUTOMATED...DO NOT CHANGE ######################################
############################################################################

# One touch operation...run the whole code. Comment out this and very last bracket to run sections manually


#{  

#running the OM and EM together
run.model<-function(){
  
#run the OM
setwd(OM_direct)
invisible(shell(paste0(OM_name," -nohess"),wait=T))

# move the .rep from the OM_direct and change name. 
from<-paste0(OM_direct,"\\",OM_name,".rep")
to<-paste0(EM_direct,"\\",EM_name,".dat")


#run the EM
setwd(EM_direct)

#remove exsisting cor and std files
cor.name<-paste0(EM_name,".cor")
std.name<-paste0(EM_name,".std")

if (file.exists(cor.name)) {
  file.remove(cor.name)
}

if (file.exists(std.name)) {
  file.remove(std.name)
}


#remove .dat
if (file.exists(to)) {
  file.remove(to)
}


#copy .dat over from OM
file.copy(from = from,  to = to)


#run the model
system.time( # keeping track of time for run

invisible(shell(paste0(EM_name),wait=T)))

} #end running the OM/EM together

time.elapsed<-run.model()

#use this function to run the model
#make.plots()

#########################################################
#########################################################
##### Ploting Code        ###############################
##### Look at the outputs ###############################
#########################################################
#########################################################
  
#plot and output function
make.plots<-function(direct=EM_direct){ #run diagnostics plotting
  
#Read in model .rep
out<-readList(paste(EM_direct,paste0(EM_name,".rep"),sep="\\")) #read in .rep file

#use the pbsadmb package to get the cor and std files in

cor.name<-paste0(EM_name,".cor")

if (file.exists(cor.name) && out$dummy<0) {
  cor<-readRep(EM_name, suffix=(".cor"), global=FALSE)
}


std.name<-paste0(EM_name,".std")

if (file.exists(std.name)) {
std<-readRep(EM_name, suffix=".std", global=FALSE)
}



#pull info about the model
na<-out$nages
nyrs<-out$nyrs
npops_OM<-out$npops_OM
npops<-out$npops
nreg_OM<-out$nregions_OM
nreg<-out$nregions
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)


#for running the meta pop example. Might need fixing if more complex
if(npops_OM>1){
  nreg_OM=sum(nreg_OM)}

if(npops>1){
  nreg=sum(nreg)}
   

###################################
#Set my theme for ggplot plotting
#################################

diag_theme<-
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))

######################################################################################

######################################################################################
#Likelihoods in histagram

#vector of values 
likes<-data.frame(names=c("Catch","Survey","Fishery_age","Survey_age","Recruitment","Tags"),likes=c(out$catch_like,out$survey_like,out$fish_age_like,out$survey_age_like,out$rec_like,out$tag_like))

# change the order
likes$names <- factor(likes$names, as.character(likes$names))


like.p <-ggplot(likes, aes(names, likes))+
  geom_bar(stat = "identity", col="grey60")+
  xlab(" ")+
  ylab("Values")+
  ggtitle("Likelihood Components")+
  diag_theme+
  theme(axis.text.x=element_text(angle=40,hjust=0.5,vjust=0.5),plot.margin = unit(c(0.5, 0.5, 2, 0.5),"lines"))
  

##############################
#### RECRUITMENT PLOTS #######
##############################

#rec total

#Matching - might need to refine for additional mismatches
if(nreg_OM==nreg){ 
rec.total<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Rec_Est = as.vector(t(out$recruits_BM)), Rec_True=as.vector(t(out$recruits_BM_TRUE)))}

#spatial to panmictic
if(nreg_OM>1 && nreg==1){ 
#aggregating rec true  
Rec_True=colSums(out$recruits_BM_TRUE)

rec.total<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), Rec_Est = as.vector(t(out$recruits_BM)),Rec_True=Rec_True)}
    

rec.total.plot<-melt(rec.total,id=c("Reg","Year"))
rec.total.plot$Reg<-as.factor(rec.total.plot$Reg)

rec1<-ggplot(rec.total.plot,aes(Year,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Recruitment")+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Recruitment Total")


##########################
#recruit residual plot

if(resid.switch==1){
  #calculate the resids as Relative % Diff
  rec.total$resids<-(rec.total$Rec_Est-rec.total$Rec_True)
}

if(resid.switch==2){
  #calculate the resids as Relative % Diff
  rec.total$resids<-((rec.total$Rec_Est-rec.total$Rec_True)/rec.total$Rec_True)*100
}

#preparing for plotting
rec.resids.plot<-melt(rec.total[,c(1,2,5)],id=c("Reg","Year"))
rec.resids.plot$Reg<-as.factor(rec.resids.plot$Reg)

#rec resids plot  
R.resid<-ggplot(rec.resids.plot,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("%Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Recruitment Residuals")


###################################
# Recruitment Deviations

#matching metapop/natal homing
if(npops>1 && npops_OM>1){

Rec_Dev_Est=as.vector(t(out$rec_devs))
Rec_Dev_True=as.vector(t(out$rec_devs_TRUE[,2:ncol(out$rec_devs_TRUE)]))
#Rec_Dev_True=as.vector(t(out$rec_devs_TRUE))

rec.devs<-data.frame(Year=rep(years[-1],nreg),Reg=rep(1:nreg,each=(nyrs-1)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)

#rec.devs<-data.frame(Year=rep(years,nreg),Reg=rep(1:nreg,each=(nyrs)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)


rec.devs.plot<-melt(rec.devs,id=c("Reg","Year"))
rec.devs.plot$Reg<-as.factor(rec.devs.plot$Reg)
rec_ave<-data.frame(pop=as.character(1:npops),R_ave=out$R_ave, R_ave_TRUE=out$R_ave_TRUE)


rec2<-ggplot(rec.devs.plot,aes(Year,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  theme_bw()+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Deviations")+
  #geom_hline(data=rec_ave,aes(yintercept=R_ave), col = e.col)+
  #geom_hline(data=rec_ave, aes(yintercept = R_ave_TRUE),col=t.col,lty=2)+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Recruitment Deviations")
}


#matching panmictic/metamictic or mismatch spatial to panmictic
if(npops==1 && npops_OM==1){

rec.devs<-data.frame(Year=years[-1],Rec_Dev_Est=out$rec_devs, Rec_Dev_True=out$rec_devs_TRUE[2:length(out$rec_devs_TRUE)])
  
rec.devs.plot<-melt(rec.devs,id=c("Year"))


#add the rec_ave
  
rec2<-ggplot(rec.devs.plot,aes(Year,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
    theme_bw()+
    #facet_wrap(~Reg)+
    #geom_hline(aes(yintercept=out$R_ave), col = e.col)+
    #geom_hline(aes(yintercept=out$R_ave_TRUE), col = t.col, lty=2)+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Recruitment")+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Recruitment Deviations")
  
}


#going from metamictic to metapop
if(npops>1 && npops_OM==1){
  
  Rec_Dev_Est=colMeans(out$rec_devs)
  Rec_Dev_True=out$rec_devs_TRUE[2:length(out$rec_devs_TRUE)]
  #Rec_Dev_True=as.vector(t(out$rec_devs_TRUE))
  
  rec.devs<-data.frame(Year=rep(years[-1],nreg),Reg=rep(1:nreg,each=(nyrs-1)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)
  
  #rec.devs<-data.frame(Year=rep(years,nreg),Reg=rep(1:nreg,each=(nyrs)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)
  
  
  rec.devs.plot<-melt(rec.devs,id=c("Reg","Year"))
  rec.devs.plot$Reg<-as.factor(rec.devs.plot$Reg)
  rec_ave<-data.frame(pop=as.character(1:npops),R_ave=out$R_ave, R_ave_TRUE=out$R_ave_TRUE)
  
  
  rec2<-ggplot(rec.devs.plot,aes(Year,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
    theme_bw()+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Deviations")+
    #geom_hline(data=rec_ave,aes(yintercept=R_ave), col = e.col)+
    #geom_hline(data=rec_ave, aes(yintercept = R_ave_TRUE),col=t.col,lty=2)+
    facet_wrap(~Reg)+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Recruitment Deviations")
}

#going from metamictic to metapop
if(npops==1 && npops_OM>1){
  
  Rec_Dev_Est=out$rec_devs
  Rec_Dev_t.temp<-colMeans(out$rec_devs_TRUE)
  Rec_Dev_True=Rec_Dev_t.temp[2:length(Rec_Dev_t.temp)]

  
  rec.devs<-data.frame(Year=rep(years[-1],nreg),Reg=rep(1:nreg,each=(nyrs-1)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)
  
  #rec.devs<-data.frame(Year=rep(years,nreg),Reg=rep(1:nreg,each=(nyrs)),Rec_Dev_Est=Rec_Dev_Est,Rec_Dev_True=Rec_Dev_True)
  
  
  rec.devs.plot<-melt(rec.devs,id=c("Reg","Year"))
  rec.devs.plot$Reg<-as.factor(rec.devs.plot$Reg)
  rec_ave<-data.frame(pop=as.character(1:npops),R_ave=out$R_ave, R_ave_TRUE=out$R_ave_TRUE)
  
  
  rec2<-ggplot(rec.devs.plot,aes(Year,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
    theme_bw()+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Deviations")+
    #geom_hline(data=rec_ave,aes(yintercept=R_ave), col = e.col)+
    #geom_hline(data=rec_ave, aes(yintercept = R_ave_TRUE),col=t.col,lty=2)+
    facet_wrap(~Reg)+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Recruitment Deviations")
}



######################
#Rec_dev resid plot
  

  if(resid.switch==1){
    #calculate the resids as Relative % Diff
    rec.devs$resids<-(rec.devs$Rec_Dev_Est-rec.devs$Rec_Dev_True)
  }
  

  if(resid.switch==2){
    #calculate the resids as Relative % Diff
    rec.devs$resids<-((rec.devs$Rec_Dev_Est-rec.devs$Rec_Dev_True)/rec.devs$Rec_Dev_True)*100
  }
  

if(npops==1){  
  
#preparing for plotting
  rec.dev.resids.plot<-melt(rec.devs[,c(1,4)],id="Year")
  
  #rec resids plot  
  R.dev.resid<-ggplot(rec.dev.resids.plot,aes(Year,value))+
    geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
    geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
    theme_bw()+
    scale_color_gradient2(low="red",mid="grey",high ="blue")+
    ylab("Difference (Estimated-True)")+
    diag_theme+
    theme(legend.position = "none", legend.justification = c(1,1))+
    ggtitle("Recruitment Devs Residuals")
}
  
  if(npops>1){ 
    #preparing for plotting
    rec.dev.resids.plot<-melt(rec.devs[,c(1,2,5)],id=c("Year","Reg"))
   
     #rec resids plot  
    R.dev.resid<-ggplot(rec.dev.resids.plot,aes(Year,value))+
      geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
      geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
      theme_bw()+
      scale_color_gradient2(low="red",mid="grey",high ="blue")+
      ylab("% Difference (Estimated-True)")+
      facet_wrap(~Reg)+
      diag_theme+
      theme(legend.position = "none", legend.justification = c(1,1))+
      ggtitle("Recruitment Devs Residuals")
  }
  
  
#################  
#Rec Apportionment

#matching panmictic or metamictic
if(nreg_OM==nreg){ 
rec.apport<-data.frame(Year=rep((years[2:nyrs]),nreg), Reg=rep(c(1:nreg),each=nyrs-1),Rec_Prop_Est = as.vector(t(out$Rec_Prop)), Rec_Prop_True=as.vector(t(out$Rec_Prop_TRUE)))

#rec.apport<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Rec_Prop_Est = as.vector(t(out$Rec_Prop)), Rec_Prop_True=as.vector(t(out$Rec_Prop_TRUE)))

}
  
#taking the mean for a mismatch spatial to panmictic
if(nreg_OM>1 && nreg==1){ 
  
rec.apport<-data.frame(Year=rep(years[2:nyrs],nreg_OM), Reg=rep(c(1:nreg_OM),each=nyrs-1),Rec_Prop_Est = NA, Rec_Prop_True=as.vector(t(out$Rec_Prop_TRUE))
                       
)
}



rec.apport.plot<-melt(rec.apport,id=c("Reg","Year"))
rec.apport.plot$Reg<-as.factor(rec.apport.plot$Reg)
  
  
rec.prop<-ggplot( rec.apport.plot,aes(Year,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
    facet_wrap(~Reg)+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Apportionment Proportion")+
    ylim(0,1)+
    theme_bw()+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Recruitment Apportionment")
 
  

#resids  
  if(resid.switch==1){
    #calculate the resids as Relative % Diff
    rec.apport$resids<-(rec.apport$Rec_Prop_Est-rec.apport$Rec_Prop_True)
  }
  
  if(resid.switch==2){
    #calculate the resids as Relative % Diff
    rec.apport$resids<-((rec.apport$Rec_Prop_Est-rec.apport$Rec_Prop_True)/rec.apport$Rec_Prop_True)*100
  }
  
  #preparing for plotting
  rec.apport.resids.plot<-melt(rec.apport[,c(1,2,5)],id=c("Reg","Year"))
  #rec.apport.resids.plot$Reg<-as.factor(rec.resids.plot$Reg)
  
  #rec resids plot  
  R.apport.resid<-ggplot(rec.apport.resids.plot,aes(Year,value))+
    geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
    geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
    theme_bw()+
    scale_color_gradient2(low="red",mid="grey",high ="blue")+
    ylab("% Difference (Estimated-True)")+
    facet_wrap(~Reg)+
    diag_theme+
    theme(legend.position = "none", legend.justification = c(1,1))+
    ggtitle("Recruitment Apportionment Residuals")
  

##########################
#Mortality graph
#######################

#matching population structures or metamictic to panmictic

    mort.rate<-data.frame(years=rep(1:nyrs,nreg), Reg=rep(c(1:nreg),each=nyrs))
    m.rate.temp<-data.frame(out$M)
    names(m.rate.temp)<-1:na
    m.rate<-cbind(mort.rate,m.rate.temp)
    m.rate.df<-melt(m.rate,id=c("years","Reg"))
    names(m.rate.df)<-c("years","Reg","Age","M")

  m.plot<-data.frame(m.rate.df %>%
    group_by(Reg,Age) %>%
    summarise(M.ave=mean(M)))
  
#matching
 if(nreg_OM==nreg||npops==npops_OM){ 
  m.plot$M_T<-as.vector(out$M_TRUE)
 }

#metapop to panmictic  
if(npops<npops_OM){ 
    
  M_ave_temp<-out$M_TRUE
  
  for(i in 1:na){
    M_ave_temp[,i]<-out$M_TRUE[,i]*out$abund_frac_region_OM
    m.plot$M_T<-colSums(M_ave_temp)
    }
  }
  

#calc resids
  if(resid.switch==1){
    #calculate the resids as Relative % Diff
    m.plot$resid=(m.plot$M.ave-m.plot$M_T)
  }
  
  if(resid.switch==2){
    m.plot$resid=((m.plot$M.ave-m.plot$M_T)/(m.plot$M_T))*100
  }
  
m.plot.p<-melt(m.plot,id=c("Age","Reg"))


#time and age varying M
mort.p<-ggplot(m.plot.p[!m.plot.p$variable=="resid",],aes(x=Age, y=value, group = variable))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
      facet_wrap(~Reg)+
      scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
      scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
      ylab("Mortality Rate")+
      diag_theme+
      ylim(0,1)+
      xlab("Age")+
      #theme(legend.position = c(1, 0), legend.justification = c(1,0))+
      theme(legend.position = c(1, 1), legend.justification = c(1,1))+
      ggtitle("Mortality Rate")
  

#resid plot
m.resid.p<-ggplot(m.plot.p[m.plot.p$variable=="resid",],aes(x=Age, y=value, group = variable))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Mortality Rate Residuals")


  
##########################
# Add reporting rate here
  
#matching Panmictic or matching metapop
if(npops==npops_OM && npops>1){
tag.rep.rate<-data.frame(Rel_Reg=rep(1:nreg,each=out$nyrs_release),Rel_year=rep(1:out$nyrs_release,nreg),RR_est=out$report_rate,RR_true=out$report_rate_TRUE)
}

#Matching spatial to spatial
if(npops==npops_OM && npops==1||npops_OM==nreg||npops_OM<npops){
  tag.rep.rate<-data.frame(Rel_Reg=rep(1:nreg,each=out$nyrs_release),Rel_year=rep(1:out$nyrs_release,nreg),RR_est=as.vector(out$report_rate),RR_true=as.vector(out$report_rate_TRUE))
}

#works for mismatch Spatial to panmictic or metapop OM to metamictic
if((npops_OM>npops&&nreg>nreg_OM)||(nreg_OM>nreg&&npops==1)){
  
#fill a matrix of reporting rate
RR_temp<-matrix(out$report_rate_TRUE,out$nyrs_release,npops_OM, byrow=T)
tag.rep.rate<-data.frame(Rel_Reg=rep(1:nreg,each=out$nyrs_release),Rel_year=rep(1:out$nyrs_release,nreg),RR_est=out$report_rate,RR_true=rowMeans(RR_temp))
}


rr.plot<-melt(tag.rep.rate,id=c("Rel_Reg","Rel_year"))

rr.p<-ggplot(rr.plot,aes(Rel_year, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  facet_wrap(~Rel_Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Reporting Rate")+
  diag_theme+
  ylim(0,1)+
  theme(legend.position = c(1, 0), legend.justification = c(1,0))+
  ggtitle("Reporting Rate")



#calc resids
if(resid.switch==1){
  #calculate the resids as Relative % Diff
  tag.rep.rate$resid<-(tag.rep.rate$RR_est-tag.rep.rate$RR_true)
}

if(resid.switch==2){
  #calculate the resids as Relative % Diff
  tag.rep.rate$resid<-((tag.rep.rate$RR_est-tag.rep.rate$RR_true)/tag.rep.rate$RR_true)*100
}

tag.rep.resid<-tag.rep.rate[,c(1,2,5)]
rr.resid<-melt(tag.rep.resid,id=c("Rel_Reg","Rel_year"))


#resid plot
rr.resid.p<-ggplot(rr.resid,aes(Rel_year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Rel_Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Reporting Rate Residuals")

  
  
  

#################################
#### ABUNDANCE PLOTS ############
#################################

####################
#initial abundance
  
#panmictic/metamictic matching  
if(nreg_OM==nreg && npops==1){ 
    init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),In_ab_Est = as.vector(t(out$Init_Abund)), In_ab_True=as.vector(t(out$Init_Abund_TRUE)))
}

#matching metapopulation
if(npops>1 & npops_OM>1){
  
  i.ab_temp<-data.frame(out$Init_Abund)
  i.ab_temp_T<-data.frame(out$Init_Abund_T)
  pops_temp<-rep(1:npops,each=npops)
  
  t_T<-split(i.ab_temp_T,pops_temp)
  t_E<-split(i.ab_temp,pops_temp)
  
  #build empty matrix to fill with sums
  t2_T<-matrix(NA,npops,na)
  t2_E<-matrix(NA,npops,na)
  
  for(i in 1:npops){
    t2_T[i,]<-colSums(data.frame(t_T[i]))
    t2_E[i,]<-colSums(data.frame(t_E[i]))
  }
  
  #the sums across the populations
  init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na), In_ab_Est = as.vector(t(t2_E)), In_ab_True=as.vector(t(t2_T))) 
}


#mismatch metamictic/metapop to panmictic  
if((out$npops_OM==1 && sum(out$nregions_OM)>1&& npops==1 && nreg==1)||(out$npops_OM>1 && sum(out$nregions_OM)>1 && npops==1 && nreg==1)){ 
  init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),In_ab_Est=out$Init_Abund,In_ab_True=colSums(out$Init_Abund_TRUE))
}

#mismatch metapop to metamictic
  if(npops==1&npops_OM>1&&out$nregions>1){
  
    i.ab_temp_T<-data.frame(out$Init_Abund_TRUE)
    pops_temp<-rep(1:npops_OM,each=npops_OM)
    
    t_T<-split(i.ab_temp_T,pops_temp)
    
    #build empty matrix to fill with sums
    t2_T<-matrix(NA,npops_OM,na)
    
    for(i in 1:npops_OM){
      t2_T[i,]<-colSums(data.frame(t_T[i]))
    }

#the sums across the populations
    init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na), In_ab_Est = as.vector(t(out$Init_Abund)), In_ab_True=as.vector(t(t2_T))) 
  }

#mismatch metamictic EM to metapop OM
if(npops>1&npops_OM==1&&nreg_OM>1){
  i.ab_temp<-data.frame(out$Init_Abund)
  pops_temp<-rep(1:npops,each=npops)
  t<-split(i.ab_temp,pops_temp)
  
  #build empty matrix to fill with sums
  t2<-matrix(NA,npops,na)
  
  for(i in 1:npops){
    t2[i,]<-colSums(matrix(unlist(t[i]),npops,na))
  }
  
  #the sums across the populations
  init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na), In_ab_Est = as.vector(t(t2)), In_ab_True=as.vector(t(out$Init_Abund_TRUE))) 
}

#prepare for plot
  init.abund.plot<-melt(init.abund,id=c("Reg","Age"))
  init.abund.plot$Reg<-as.factor(init.abund.plot$Reg)

#plot
  init.ab<-ggplot(init.abund.plot,aes(Age, value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
    facet_wrap(~Reg)+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Abundance")+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Initial Abundance")  
  

 #Init Abundance residuals 

if(resid.switch==1){
  init.abund$resids<-(init.abund$In_ab_Est-init.abund$In_ab_True)
}
  
if(resid.switch==2){
  init.abund$resids<-((init.abund$In_ab_Est-init.abund$In_ab_True)/init.abund$In_ab_True)*100
}
  
  
  init.abund.resid.plot<-melt(init.abund[,c(1,2,5)],id=c("Reg","Age"))
  init.abund.resid.plot$Reg<-as.factor(init.abund.resid.plot$Reg) 
  
#resid plot
  init.ab.resid<-ggplot(init.abund.resid.plot,aes(Age,value))+
    geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
    geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
    theme_bw()+
    scale_color_gradient2(low="red",mid="grey",high ="blue")+
    ylab("% Difference (Estimated-True)")+
    facet_wrap(~Reg)+
    diag_theme+
    theme(legend.position = "none", legend.justification = c(1,1))+
    ggtitle("Initial Abundance Residuals")
  

#################################
# Total Biomass
  
#matching
if(nreg_OM==nreg){ 
bio.dat<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Bio_est = as.vector(t(out$biomass_AM)), Bio_True=as.vector(t(out$biomass_AM_TRUE)))
}
  
#mismatch metapop and metamictic to panmictic
if(nreg_OM>1 && nreg==1){ 
bio.dat<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Bio_est = out$biomass_AM, Bio_True=colSums(out$biomass_AM_TRUE))
}

  
  bio.plot<-melt(bio.dat,id=c("Reg","Year"))
  bio.plot$Reg<-as.factor(bio.plot$Reg)
  
  bio.p<-ggplot(bio.plot,aes(Year,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
    facet_wrap(~Reg)+
    scale_color_manual(values = c(e.col,t.col),labels = c("Predicted","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Predicted","True"))+
    ylab("Total Biomass")+
    diag_theme+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Total Biomass")
  
  
#calculate resids
  
if(resid.switch==1){ 
    bio.dat$resid<-(bio.dat$Bio_est-bio.dat$Bio_True)
}
  
if(resid.switch==2){ 
  bio.dat$resid<-((bio.dat$Bio_est-bio.dat$Bio_True)/bio.dat$Bio_True)*100
}
  

  bio.resid.plot<-melt(bio.dat[,c(1,2,5)],id=c("Reg","Year"))
  bio.resid.plot$Reg<-as.factor(bio.resid.plot$Reg)
  
  bio.resid<-ggplot(bio.resid.plot,aes(Year,value))+
    geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
    geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
    theme_bw()+
    scale_color_gradient2(low="red",mid="grey",high ="blue")+
    ylab("% Difference (Estimated-True)")+
    facet_wrap(~Reg)+
    diag_theme+
    theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
    theme(legend.position = "none", legend.justification = c(1,1))+
    ggtitle("Biomass Residuals") 
  
  
###############################
#SSB


#matching
if(nreg_OM==nreg){ 
ssb.dat<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),SSB_est = as.vector(t(out$SSB_region)), SSB_True=as.vector(t(out$SSB_region_TRUE)))
}
  
#aggregated
if(nreg_OM>1 && nreg==1){ 
ssb.dat<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),SSB_est = out$SSB_region, SSB_True=colSums(out$SSB_region_TRUE))}
  

ssb.plot<-melt(ssb.dat,id=c("Reg","Year"))
ssb.plot$Reg<-as.factor(ssb.plot$Reg)

ssb.p<-ggplot(ssb.plot,aes(Year,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Predicted","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Predicted","True"))+
  ylab("SSB")+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("SSB")



#SSB residual plot

if(resid.switch==1){ 
ssb.dat$resid<-(ssb.dat$SSB_est-ssb.dat$SSB_True)
}

if(resid.switch==2){
  ssb.dat$resid<-((ssb.dat$SSB_est-ssb.dat$SSB_True)/ssb.dat$SSB_True)*100
}

ssb.resid.plot<-melt(ssb.dat[,c(1,2,5)],id=c("Reg","Year"))
ssb.resid.plot$Reg<-as.factor(ssb.resid.plot$Reg)

ssb.resid<-ggplot(ssb.resid.plot,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("SSB Residuals")



############################################################
########## Selectivity Parameters ##########################
############################################################

#####################
#Fishery Selectivity


if(nreg_OM==nreg){ 
f.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = as.vector(t(out$selectivity_age)), Select_T=as.vector(t(out$selectivity_age_TRUE)))

if(out$nfleets>1 && out$nregions==1){ # this is a fleets as areas model
f.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:out$nfleets),each=na),Select_Est = as.vector(out$selectivity_age), Select_T=as.vector(out$selectivity_age_TRUE))
}
}

if(nreg_OM>1 && nreg==1){ 
##if(out$nfleets>out$nfleets_OM){
##f.select<-data.frame(Age=rep(ages,nreg_OM), Reg=rep(c(1:nreg_OM),each=na), Select_Est = as.vector(out$selectivity_age), Select_T=as.vector(out$selectivity_age_TRUE))}
  

#aggregate by abund frac...the correct way do to it!
f.select.temp<-matrix(out$selectivity_age_TRUE,nreg_OM,byrow=T)
f.select.temp2<-f.select.temp
for(i in 1:na){
  f.select.temp2[,i]<-f.select.temp[,i]*out$abund_frac_region_OM
}
  
f.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = out$selectivity_age, Select_T=colSums(f.select.temp2))

 }


f.select.plot<-melt(f.select,id=c("Reg","Age"))
f.select.plot$Reg<-as.factor(f.select$Reg)

f.select.p<-ggplot(f.select.plot,aes(Age, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Selectivity")+
  diag_theme+
  theme(legend.position = c(1, 0), legend.justification = c(1,0))+
  ggtitle("Fishery Selectivity")

#fishery selectivity resid
if(resid.switch==1){
  f.select$resid<-(f.select$Select_Est-f.select$Select_T)
}

if(resid.switch==2){
f.select$resid<-((f.select$Select_Est-f.select$Select_T)/f.select$Select_T)*100
}

f.select.resid.plot<-melt(f.select[,c(1,2,5)],id=c("Reg","Age"))
f.select.resid.plot$Reg<-as.factor(f.select.resid.plot$Reg)


#resid plot
f.select.resid<-ggplot(f.select.resid.plot,aes(Age,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Fishery Selectivity Residuals")




########################
#Survey Selectivity

if(nreg_OM==nreg){ 
s.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = as.vector(t(out$survey_selectivity_age)),Select_T=as.vector(t(out$survey_selectivity_age_TRUE)))
}

if(nreg_OM>1 && nreg==1){ 
  #aggreagating
  temp1<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg_OM),each=na),Select_T=as.vector(t(out$survey_selectivity_age_TRUE)))
  temp2<-group_by(temp1,Age) %>% summarise(Select_T = mean(Select_T))
  
  s.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = out$survey_selectivity_age, Select_T=temp2$Select_T)
}


s.select.plot<-melt(s.select,id=c("Reg","Age"))
s.select.plot$Reg<-as.factor(s.select$Reg)

s.select.p<-ggplot(s.select.plot,aes(Age, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Selectivity")+
  diag_theme+
  theme(legend.position = c(1, 0), legend.justification = c(1,0))+
  ggtitle("Survey Selectivity")



# survey selectivity resids plot
if(resid.switch==1){
  s.select$resid<-(s.select$Select_Est-s.select$Select_T)
}

if(resid.switch==2){
  s.select$resid<-((s.select$Select_Est-s.select$Select_T)/s.select$Select_T)*100
}

s.select.resid.plot<-melt(s.select[,c(1,2,5)],id=c("Reg","Age"))
s.select.resid.plot$Reg<-as.factor(s.select.resid.plot$Reg)

s.select.resid<-ggplot(s.select.resid.plot,aes(Age,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Survey Selectivity Residuals")




############################################################
########## F by year #######################################
############################################################

#matching or metapop OM mismatch
if((out$npops==out$npops_OM && out$nregions==out$nregions_OM)||npops>npops_OM||npops_OM>npops){ 
    f.max<-rowMaxs(out$F)
    f.max.t<-rowMaxs(out$F_TRUE)
    F.year=data.frame(Year=rep(years,nreg),Reg=rep(c(1:nreg),each=nyrs),F_year=as.vector(f.max), F_year_T=f.max.t)
}


#mismatch spatial to panmictic
if((npops_OM>1 && nreg==1 && out$nfleets==1)||(npops==out$npops && out$nregions_OM>1&&nreg==1)){
  f.max<-rowMaxs(out$F)
  f.max.t.temp<-rowMaxs(out$F_TRUE)
  
  #aggreagating - yearly weighted average
  temp1<-matrix(f.max.t.temp,nyrs,nreg_OM,byrow=F)
  temp2<-f.max.t.temp*t(out$abund_frac_year_OM)
  F.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),F_year=f.max, F_year_T=rowSums(temp2))
}



############################
#some old stuff for mismatch
#if(nreg_OM>1 && nreg==1 && out$nfleets>1){#mismatch with FAA

#  f.max.test<-as.vector(unlist(data.frame(out$F_year)))
#  f.max.t<-as.vector(unlist(data.frame(out$F_year_TRUE)))
#  f.mat.temp<-matrix(unlist(data.frame(out$F)),nyrs*nreg_OM,na, byrow=T)
#  f.max<-rowMaxs(f.mat.temp)
#  f.mat.temp.t<-matrix(unlist(data.frame(out$F_TRUE)),nyrs*nreg_OM,na, byrow=T)
#  f.max.t<-rowMaxs( f.mat.temp.t)
#  F.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:out$nfleets),each=nyrs),F_year=f.max, F_year_T=f.max.t)
#}

#pull_F_est<-out[grep("alt_F", names(out), value = TRUE)]
#F_est<-data.frame(matrix(unlist(pull_F_est),nyrs*npops,npops,byrow=TRUE))
#f.mat.temp<-matrix(unlist(data.frame(F_est)),nyrs*npops,na, byrow=T)
#f.max<-rowMaxs(f.mat.temp)
#f.mat.temp<-matrix(unlist(data.frame(F_true)),nyrs*nreg_OM,na, byrow=T)
#f.max.t<-rowMaxs(f.mat.temp)
#F.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:out$npops),each=nyrs),F_year=as.vector(unlist(f.max)), #F_year_T=as.vector(unlist(f.max.t)))
#}
##################################

F.plot<-melt(F.year,id=c("Reg","Year"))
F.plot$Reg<-as.factor(F.plot$Reg)

F.plot.p<-ggplot(F.plot,aes(Year,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  theme_bw()+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("F")+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Fully selected F by Year")


#F resids
if(resid.switch==1){
F.year$resid<-(F.year$F_year-F.year$F_year_T)}

if(resid.switch==2){
  F.year$resid<-((F.year$F_year-F.year$F_year_T)/F.year$F_year_T)*100
}


F.resid.p<-melt(F.year[,c(1,2,5)],id=c("Reg","Year"))
F.resid.p$Reg<-as.factor(F.resid.p$Reg)

F.resid.plot<-ggplot(F.resid.p,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Fully Selected F Residuals")


############################################################
########## Movement by year ################################
############################################################


#matching Panmictic and multi-area
if((out$nregions_OM==out$nregions && npops==1&&npops_OM==1)){ 
  pull<-out[grep("T_est", names(out), value = TRUE)]
  pull.t<-out[grep("T_true", names(out), value = TRUE)]
  
  T_est<-data.frame(do.call("rbind",pull))
  T_true<-data.frame(do.call("rbind",pull.t))
  
  
if(resid.switch==1){
  T_resid<-(T_true-T_est)}
  
if(resid.switch==2){
  T_resid<-((T_true-T_est)/T_true)*100}
  
for(i in 1:nreg){
    names(T_est)[i]<-paste0("Est_",i)
    names(T_true)[i]<-paste0("True_",i)
    names(T_resid)[i]<-paste0("Resid_",i)}
    
}



#other matching metapop
if(npops==npops_OM && npops>1){
  
#combine movements
pull<-out[grep("T_est", names(out), value = TRUE)]
pull.t<-out[grep("T_true", names(out), value = TRUE)]

T_true<-data.frame(matrix(unlist(pull.t),na*nyrs*npops,npops,byrow=T)) #just setting the
T_est<-data.frame(matrix(unlist(pull),na*nyrs*npops,npops,byrow=T)) #just setting the

if(resid.switch==1){
  T_resid<-(T_est-T_true)}

if(resid.switch==2){
  T_resid<-((T_est-T_true)/T_true)*100}

  
for(i in 1:npops){
names(T_est)[i]<-paste0("Est_",i)
names(T_true)[i]<-paste0("True_",i)
names(T_resid)[i]<-paste0("Resid_",i)

}
}



# for mismatch metamictic to panmictic
if(nreg==1 && nreg_OM>1&&npops_OM==1){
  T_est<-data.frame(matrix(NA,nyrs*nreg_OM*na,nreg_OM))
  pull.t<-out[grep("T_true", names(out), value = TRUE)]
  T_true<-data.frame(do.call("rbind",pull.t))
  
  if(resid.switch==1){
    T_resid<-(T_est-T_true)}
  
  if(resid.switch==2){
    T_resid<-((T_est-T_true)/T_true)*100}
  
  for(i in 1:nreg_OM){
    names(T_est)[i]<-paste0("Est_",i)
    names(T_true)[i]<-paste0("True_",i)
    names(T_resid)[i]<-paste0("Resid_",i)}
}
  

#metamictic to metapop
if(npops_OM==1 && npops>1){ 
  
  #combine movements
  pull<-out[grep("T_est", names(out), value = TRUE)]
  T_est<-data.frame(matrix(unlist(pull),na*nyrs*npops,npops,byrow=T)) #just setting the
  
  pull.t<-out[grep("T_true", names(out), value = TRUE)]
  T_true<-data.frame(do.call("rbind",pull.t))
  
  if(resid.switch==1){
    T_resid<-(T_true-T_est)}
  
  if(resid.switch==2){
    T_resid<-((T_true-T_est)/T_true)*100}
  
  
  for(i in 1:npops){
    names(T_est)[i]<-paste0("Est_",i)
    names(T_true)[i]<-paste0("True_",i)
    names(T_resid)[i]<-paste0("Resid_",i)
    
  }}

#metapop to metamictic
if(npops_OM>1 && nreg>1 && npops==1){ 

  #combine movements
  pull.t<-out[grep("T_true", names(out), value = TRUE)]
  T_true<-data.frame(matrix(unlist(pull.t),na*nyrs*npops_OM,npops_OM,byrow=T))
  
  pull<-out[grep("T_est", names(out), value = TRUE)]
  T_est<-data.frame(do.call("rbind",pull)) #just setting the estimated movement = 1
  
  if(resid.switch==1){
    T_resid<-(T_est-T_true)}
  
  if(resid.switch==2){
    T_resid<-((T_est-T_true)/T_true)*100}
  
  
  for(i in 1:npops_OM){
    names(T_est)[i]<-paste0("Est_",i)
    names(T_true)[i]<-paste0("True_",i)
    names(T_resid)[i]<-paste0("Resid_",i)
    
  }}


#metapop to panmictic
if(npops_OM>1 && nreg==1 && npops==1){ 
  
  #combine movements
  pull<-out[grep("T_true", names(out), value = TRUE)]
  T_true<-data.frame(matrix(unlist(pull),na*nyrs*npops_OM,npops_OM,byrow=TRUE)) #just setting the
  
  T_est<-data.frame(matrix(NA,na*nyrs*npops_OM,npops_OM,byrow=TRUE)) #just setting the estimated movement = 1
  
  
  if(resid.switch==1){
    T_resid<-(T_est-T_true)}
  
  if(resid.switch==2){
    T_resid<-((T_est-T_true)/T_true)*100}
  
  
  for(i in 1:npops_OM){
    names(T_est)[i]<-paste0("Est_",i)
    names(T_true)[i]<-paste0("True_",i)
    names(T_resid)[i]<-paste0("Resid_",i)
    
  }}


#build the data frame
#matching
if(nreg_OM==nreg && npops==npops_OM){ 
T.temp<-data.frame(Year=rep(years,each=na,times=nreg), Reg=rep(c(1:nreg),each=na*nyrs),Age=rep(c(1:na),nreg*nyrs))
}

#spatial to pan
if(nreg_OM==nreg && npops==1){ 
 T.temp<-data.frame(Year=rep(years,each=na,times=nreg), Reg=rep(c(1:nreg),each=na*nyrs),Age=rep(c(1:na),nreg*nyrs))
}

#metapop to metamictic
if(nreg==1 && nreg_OM>1){
  T.temp<-data.frame(Year=rep(years,each=na,times=nreg_OM), Reg=rep(c(1:nreg_OM),each=na*nyrs),Age=rep(c(1:na),nreg_OM*nyrs))
}

#metamictic to metapop
if(npops>1 && npops_OM==1 && nreg>1){
  T.temp<-data.frame(Year=rep(years,each=na,times=nreg_OM), Reg=rep(c(1:nreg_OM),each=na*nyrs),Age=rep(c(1:na),nreg_OM*nyrs))
}


T.year<-cbind(T.temp,T_est,T_true)
T.year.resid<-cbind(T.temp,T_resid)

#primary plot
T.year.plot<-melt(T.year,id=c("Reg","Year","Age"))
T.year.plot$Reg<-as.factor(T.year.plot$Reg)

T.lines<-c(rep(1,nreg),rep(6,nreg))

#if multiple pops or multiple reg matching
if(nreg>1 || npops>1){
T.col<-rep(mycols(nreg),nreg)}


#for panmictic matching
if(nreg==1 && npops==1){
  T.col<-mycols(2)}


#if spatial to panmictic
if(nreg_OM>1 && npops==1){
  T.lines<-c(rep(1,nreg_OM),rep(6,nreg_OM))
  T.col<-rep(mycols(nreg_OM),nreg_OM)}


T.year.p<-ggplot(T.year.plot,aes(Year,value))+
  geom_line(aes(col = variable,linetype=variable),lwd=line.wd,stat = "identity")+
  theme_bw()+
  #facet_wrap(~Reg)+
  facet_grid(Reg ~ Age)+
  scale_color_manual(values = T.col)+
  #scale_linetype_manual(values=T.lines,guide=FALSE)+
  scale_linetype_manual(values=T.lines)+
  ylab("Movement Rate")+
  diag_theme+
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))+
  theme(legend.position = "right", legend.justification = c(1,1))+
  ggtitle("Yearly Movement Rate")



#T_matrix resids
T.year.resid<-melt(T.year.resid,id=c("Reg","Year","Age"))
T.year.resid$Reg<-as.factor(T.year.resid$Reg)

T.resid.plot<-ggplot(T.year.resid,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=1, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_grid(Reg ~ Age)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))+
  #theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Yearly Movement Rate Residuals")



############################################################
########## Fits to data ####################################
############################################################

##############
# Yield

##fleets as areas 
#if((out$nregions_OM==out$nregions && out$nfleets>=1)||(out$nregions_OM>out$nregions && out$nfleets>1)){ 
#Y.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:out$nfleets),each=nyrs), Estimated=as.vector(out$yield_fleet),Observed=as.vector(out$OBS_yield_fleet))
#}

#mismatch to panmictic
if(nreg_OM > 1 && sum(out$nfleets)==1){ 
  Y.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), Estimated=as.vector(out$yield_fleet),Observed=as.vector(out$OBS_yield_fleet))
}
          
#spatial to spatial           
if(out$nfleets==out$nfleets_OM && nreg>1){
Y.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), Estimated=out$yield_fleet,Observed=out$OBS_yield_fleet) 
}



Y.year.plot<-melt(Y.year, id=c("Reg","Year"))

yield.p<-ggplot(Y.year.plot,aes(Year,value,shape=variable))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  geom_point(size=2, alpha = 0.5)+
  scale_shape_manual(values=c(NA,16),labels = c("Estimated","Observed"))+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","Observed"))+
  scale_linetype_manual(values=c(1,0),labels = c("Estimated","Observed"))+
  ylab("Yield")+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Yield")


# calculate residuals

#if(resid.switch==1){
#  Y.year$resid<-(Y.year$Observed-Y.year$Estimated)}

#if(resid.switch==2){
Y.year$resid<-((Y.year$Estimated-Y.year$Observed)/Y.year$Observed)*100
#}

y.resid.p<-melt(Y.year[,c(1,2,5)],id=c("Reg","Year"))


#
Y.resid.plot<-ggplot(y.resid.p,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Yield Residuals")



################
# Survey Index

Survey.year<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), SI_Est=out$survey_fleet_bio,SI_Obs=out$OBS_survey_fleet_bio ) 

Survey.year.plot<-melt(Survey.year, id=c("Reg","Year"))

survey.p<-ggplot(Survey.year.plot,aes(Year,value,shape=variable))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=line.wd)+
  geom_point(size=2, alpha = 0.5)+
  scale_shape_manual(values=c(NA,16),labels = c("Estimated","Observed"))+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","Observed"))+
  scale_linetype_manual(values=c(1,0),labels = c("Estimated","Observed"))+
  ylab("Survey Index")+
  diag_theme+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Survey Biomass")


#Calc resids
#if(resid.switch==1){
#  Survey.year$resid<-(Survey.year$SI_Obs-Survey.year$SI_Est)}

#if(resid.switch==2){
Survey.year$resid<-((Survey.year$SI_Est-Survey.year$SI_Obs)/Survey.year$SI_Obs)*100
#}

surv.resid.p<-melt(Survey.year[,c(1,2,5)],id=c("Reg","Year"))


#
Surv.resid.plot<-ggplot(surv.resid.p,aes(Year,value))+
  geom_hline(aes(yintercept=0), col = "grey20", lty = 2)+
  geom_point(aes(color=value),size=2, alpha = 0.9, pch=16)+
  theme_bw()+
  scale_color_gradient2(low="red",mid="grey",high ="blue")+
  ylab("% Difference (Estimated-True)")+
  facet_wrap(~Reg)+
  diag_theme+
  theme(legend.position = "none", legend.justification = c(1,1))+
  ggtitle("Survey Index Residuals")




###################################################
# Age Compositions
##################################################

####################
#survey age comps

survey.comps.resid<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs))


if(npops==1){
survey.prop.resid<-data.frame((out$OBS_survey_prop-out$EST_survey_prop))
survey.comps.resid<-cbind(survey.comps.resid,survey.prop.resid)
survey.long<-melt(survey.comps.resid,id.vars=c("Year","Reg"))
}


if(npops>1){
  
#combine movements age_comps
  pull.surv.obs<-out[grep("OBS_survey_prop", names(out), value = TRUE)]
  pull.surv.est<-out[grep("EST_survey_age_prop", names(out), value = TRUE)]
  
  surv_obs<-data.frame(do.call("rbind",pull.surv.obs))
  surv_est<-data.frame(do.call("rbind",pull.surv.est))
  
  survey.prop.resid<-(surv_obs-surv_est)
  survey.comps.resid<-cbind(survey.comps.resid,survey.prop.resid)
  survey.long<-melt(survey.comps.resid,id.vars=c("Year","Reg"))
  }


survey.comp.plot<-
  ggplot(survey.long, aes(x = as.numeric(variable), y = Year)) + 
  geom_raster(aes(fill=value)) + 
  #scale_fill_gradient2(low="red",mid="grey99",high ="blue",limits=c(-100, 100))+
  scale_fill_gradient2(low="red",mid="grey99",high ="blue")+
  labs(fill = "% Dif")+
  scale_y_continuous(trans = "reverse")+
  labs(x="Age", y="Year", title="Survey Age Comp Residuals") +
  facet_grid(Reg~.)+
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11))+
  #theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))




######################
#fishery age comps

fishery.comps.resid<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs))

#rewrite if fleets as areas
if(out$nfleets_OM==1 && out$nfleets>1 && out$nregions==1){ 
fishery.comps.resid<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:out$nfleets),each=nyrs))
}


fishery.comps.resid<-data.frame(Year=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs))


if(npops==1){

fishery.prop.resid<-data.frame(out$EST_catch_age_fleet_prop-out$OBS_catch_prop)
#    fishery.prop.resid<-data.frame((out$OBS_catch_prop-out$EST_catch_age_fleet_prop)/out$OBS_catch_prop)
fishery.comps.resid<-cbind(fishery.comps.resid,fishery.prop.resid)
fishery.long<-melt(fishery.comps.resid,id.vars=c("Year","Reg"))
}



if(npops>1){
  #combine movements age_comps
  pull.catch.obs<-out[grep("OBS_catch_prop", names(out), value = TRUE)]
  pull.catch.est<-out[grep("EST_catch_age_fleet_prop", names(out), value = TRUE)]
  
  catch_obs<-data.frame(do.call("rbind",pull.catch.obs))
  catch_est<-data.frame(do.call("rbind",pull.catch.est))
  
  fishery_resid<-(catch_est-catch_obs)
  
  fishery.comps.resid<-cbind(fishery.comps.resid,fishery_resid)
  fishery.long<-melt(fishery.comps.resid,id.vars=c("Year","Reg"))
}


fishery.comp.plot<-
  ggplot(fishery.long, aes(x = as.numeric(variable), y = Year)) + 
  geom_raster(aes(fill=value)) + 
  #scale_fill_gradient2(low="red",mid="grey99",high ="blue",limits=c(-100, 100))+
  labs(fill = "% Dif")+
  scale_fill_gradient2(low="red",mid="grey99",high ="blue")+
  scale_y_continuous(trans = "reverse")+
  labs(x="Age", y="Year", title="Fishery Age Comp Residuals") +
  facet_grid(Reg~.)+
  theme_bw() + 
  theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11))+
  #theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))



############################
# Fits to Tag data

#building the matrix
tag.prop.resid<-data.frame(Rel_Reg=rep(1:nreg,each=out$nyrs_release*na),Rel_year=rep(1:out$nyrs_release,each = na),Rel_age=rep(1:na,(out$nyrs_release*nreg)))

#create column names 
colnamesindex1<-rep(1:nreg,out$max_life_tags)/10
colnamesindex2<-rep(1:out$max_life_tags,each=nreg)
colnamesindex<-as.character(c((colnamesindex1+colnamesindex2),"NoCapture"))


#not there yet
if(npops==1){
  tags_obs<-out$OBS_tag_prop_final
  #tags_est<-out$TRUE_tag_prop_final
  tags_est<-out$EST_tag_prop_final
  tag_resid<-(tags_est-tags_obs)
  tag_resid<-data.frame(tag_resid)
  names(tag_resid)<-colnamesindex
   
tag_resid$Rel_Reg<-tag.prop.resid[,1]
tag_resid$Rel_year<-tag.prop.resid[,2]
tag_resid$Rel_age<-tag.prop.resid[,3]
tags.long<-melt(tag_resid,id.vars=c("Rel_Reg","Rel_year","Rel_age"))

}

#write.csv(tag_resid,"tags.csv")

if(npops>1){
  #combine movements age_comps
  pull.tags.obs<-out[grep("OBS_tag_prop", names(out), value = TRUE)]
  pull.tags.est<-out[grep("EST_tag_prop", names(out), value = TRUE)]
  
  tags_obs<-do.call(rbind, pull.tags.obs)
  tags_est<-do.call(rbind, pull.tags.est)

  
  tag_resid<-(tags_est-tags_obs)
  tag_resid<-data.frame(tag_resid)
  names(tag_resid)<-colnamesindex

  tag_resid$Rel_Reg<-tag.prop.resid[,1]
  tag_resid$Rel_year<-tag.prop.resid[,2]
  tag_resid$Rel_age<-tag.prop.resid[,3]
  tags.long<-melt(tag_resid,id.vars=c("Rel_Reg","Rel_year","Rel_age"))
}


# subset the data by release region

split.by.reg<-split(tags.long,tags.long$Rel_Reg)


#not sure what I am doing here but whatever
tags.plot<-function(j) {
  
  ggplot(split.by.reg[[j]], aes(x = Rel_age, y=as.factor(variable))) + 
  geom_raster(aes(fill=value)) + 
  #scale_fill_gradient2(low="red",mid="grey99",high ="blue",limits=c(-10, 10))+
  scale_fill_gradient2(low="red",mid="grey99",high ="blue")+
  labs(fill = "% Dif")+
  #scale_y_continuous(trans = "reverse")+
  labs(x="Age of Release", y="Time-at-Liberty . Recap Region", title="Release Event") +
  #facet_grid(Rel_Reg ~ Rel_year)+
  facet_grid(Rel_Reg~Rel_year)+
  theme_bw() + 
  theme(axis.text.x=element_text(size=6, angle=0, vjust=0.3),
        #axis.text.x=element_blank(),
        axis.text.y=element_text(size=9),
        plot.title=element_text(size=11),
        #legend.title = element_blank()),
    strip.text.x = element_text(size = 10, colour = "black", face="bold"),
    panel.spacing = unit(0.05, "lines"))

}

#tags.plot(3)

##############################################
# Generate Table of Error Values from OM .dat
##############################################

#for OM
OM_dat<-readLines(paste0(OM_direct,"\\",OM_name,".dat"))

#create a table of values

OM_error_table<-data.frame(Parameter=c("Sigma_Rec","Sigma_Rec_Apport","Sigma_F","Rec_Index_sigma","Survey_Sigma","Catch_Sigma","SIM_N_Catch","SIM_N_Survey"))

#set up matrix to fill in the values
temp<-matrix(NA,dim(OM_error_table)[1],nreg_OM)
names(temp)<-1:nreg_OM
OM_error_table<-cbind(OM_error_table,temp)

###building table

if(npops_OM==1){
rec_om<-grep("_sigma_recruit",OM_dat, fixed = T)+1
OM_error_table[1,2]<-OM_dat[rec_om]
  
rec_app_om<-grep("_sigma_rec_prop",OM_dat, fixed = T)+1
OM_error_table[2,2]<-OM_dat[rec_app_om]
}


if(npops_OM>1){
rec_om<-grep("_sigma_recruit",OM_dat, fixed = T)+1
OM_error_table[1,2:(1+nreg_OM)]<-OM_dat[rec_om:(rec_om+(nreg_OM-1))]

rec_app_om<-grep("_sigma_rec_prop",OM_dat, fixed = T)+1
OM_error_table[2,2:(1+nreg_OM)]<-OM_dat[rec_app_om:(rec_app_om+(nreg_OM-1))]
}


F_om<-grep("_sigma_F",OM_dat, fixed = T)+1
OM_error_table[3,2:(1+nreg_OM)]<-OM_dat[F_om:(F_om+(nreg_OM-1))]

rec_ind_om<-grep("_rec_index_sigma",OM_dat, fixed = T)+1
OM_error_table[4,2:(1+nreg_OM)]<-OM_dat[rec_ind_om:(rec_ind_om+(nreg_OM-1))]

surv_ind_om<-grep("_sigma_survey",OM_dat, fixed = T)+1
OM_error_table[5,2:(1+nreg_OM)]<-OM_dat[surv_ind_om:(surv_ind_om+(nreg_OM-1))]

catch_om<-grep("_sigma_catch",OM_dat, fixed = T)+1
OM_error_table[6,2:(1+nreg_OM)]<-OM_dat[catch_om:(catch_om+(nreg_OM-1))]

ncatch_om<-grep("_SIM_ncatch",OM_dat, fixed = T)+1
OM_error_table[7,2:(1+nreg_OM)]<-OM_dat[ncatch_om:(ncatch_om+(nreg_OM-1))]

nsurvey_om<-grep("_SIM_nsurvey",OM_dat, fixed = T)+1
OM_error_table[8,2:(1+nreg_OM)]<-OM_dat[nsurvey_om:(nsurvey_om+(nreg_OM-1))]


OM_table<-tableGrob(OM_error_table)
h <- grobHeight(OM_table)
w <- grobWidth(OM_table)
title <- textGrob("OM Error Parameters", y=unit(0.5,"npc") + 1.0*h, 
                  vjust=0, gp=gpar(fontsize=11))

gt <- gTree(children=gList(OM_table, title))


##################################
##EM Error Parameters  
##################################

# probably don't need this if the errors are going to be matching


############################################
# CONVERGENCE etc
###########################################

#did model converge?
if(file.exists(paste0(EM_name,".cor"))==TRUE)
     {cor.text<-paste("Model Converged")}

if(file.exists(paste0(EM_name,".cor"))==FALSE)
{cor.text<-paste("Model did NOT Converge")}

#max gradient

grad<-readLines(paste0(EM_name,".par"))[1]
pos = regexpr('Max', grad)
grad.val<-substr(grad,pos,nchar(grad))

tgrob.mod <- textGrob(paste(cor.text," ",grad.val,sep = "\n"),just = "centre")

####################################################
##Table of True and Estimated values for parameters
####################################################

{
#Estimated parameter values
# not set up for domed selectivity...can add that later
EM_est_table<-data.frame(Parameter=c("q","R_ave","beta 1 fishery","beta 2 fishery","beta 3 fishery","beta 4 fishery","beta 1 survey","beta 2 survey","beta 3 survey","beta 4 survey"))

#set up matrix to fill in the values

if((out$nfleets_OM==out$nfleets && npops_OM==npops && out$nfleets>1)){ #FAA matching
  temp_est<-data.frame(matrix(NA,dim(EM_est_table)[1],out$nfleets))
  names(temp_est)<-as.character(1:out$nfleets)}

if(out$nfleets_OM==1 && out$nfleets>1){ #mismatch with FAA
temp_est<-data.frame(matrix(NA,dim(EM_est_table)[1],out$nfleets))
names(temp_est)<-as.character(1:out$nfleets)}

if((nreg==nreg_OM) && (out$nfleets_OM!=out$nfleets)){ #mismatch with metamictic
temp_est<-data.frame(matrix(NA,dim(EM_est_table)[1], nreg_OM))
names(temp_est)<-as.character(1:nreg)}

if((nreg==nreg_OM) && (out$nfleets_OM=out$nfleets) && (nreg>1)){ #matching
  temp_est<-data.frame(matrix(NA,dim(EM_est_table)[1], nreg_OM))
  names(temp_est)<-as.character(1:nreg_OM)}

if((nreg<nreg_OM) && (out$nfleets_OM==out$nfleets) && out$nfleets==1){ #spatial to panmictic
  temp_est<-data.frame(matrix(NA,dim(EM_est_table)[1],nreg))
  names(temp_est)<-as.character(1:nreg)}

EM_est_table<-cbind(EM_est_table,temp_est)
#params
EM_est_table[1,2:ncol(EM_est_table)]<-round(out$q_survey,2)

if(npops==1){
  EM_est_table[2,2]<-round(out$R_ave,2)
}

if(npops>1){
  EM_est_table[2,2:ncol(EM_est_table)]<-round(out$R_ave,2)
}




# adding in selectivity values based on specification
if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==0)
{
  EM_est_table[3,2:ncol(EM_est_table)]<-NA
  EM_est_table[4,2:ncol(EM_est_table)]<-NA
  EM_est_table[5,2:ncol(EM_est_table)]<-NA
  EM_est_table[6,2:ncol(EM_est_table)]<-NA
}

if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==0)
{
  EM_est_table[7,2:ncol(EM_est_table)]<-NA
  EM_est_table[8,2:ncol(EM_est_table)]<-NA
  EM_est_table[9,2:ncol(EM_est_table)]<-NA
  EM_est_table[10,2:ncol(EM_est_table)]<-NA
}

if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==1)
{
  EM_est_table[3,2:ncol(EM_est_table)]<-round(out$sel_beta1,2)
  EM_est_table[4,2:ncol(EM_est_table)]<-round(out$sel_beta2,2)
  EM_est_table[5,2:ncol(EM_est_table)]<-NA
  EM_est_table[6,2:ncol(EM_est_table)]<-NA
}

if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==1)
{
EM_est_table[7,2:ncol(EM_est_table)]<-round(out$sel_beta1_survey,2)
EM_est_table[8,2:ncol(EM_est_table)]<-round(out$sel_beta2_survey,2)
EM_est_table[9,2:ncol(EM_est_table)]<-NA
EM_est_table[10,2:ncol(EM_est_table)]<-NA
}

if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==2)
{
  EM_est_table[3,2:ncol(EM_est_table)]<-round(out$sel_beta1,2)
  EM_est_table[4,2:ncol(EM_est_table)]<-round(out$sel_beta2,2)
  EM_est_table[5,2:ncol(EM_est_table)]<-round(out$sel_beta3,2)
  EM_est_table[6,2:ncol(EM_est_table)]<-round(out$sel_beta4,2)
}


if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==2)
{
EM_est_table[7,2:ncol(EM_est_table)]<-round(out$sel_beta1_survey,2)
EM_est_table[8,2:ncol(EM_est_table)]<-round(out$sel_beta2_survey,2)
EM_est_table[9,2:ncol(EM_est_table)]<-round(out$sel_beta3_survey,2)
EM_est_table[10,2:ncol(EM_est_table)]<-round(out$sel_beta4_survey,2)
}

est_params<-tableGrob(EM_est_table,theme = ttheme_default(base_size = 11))
h <- grobHeight(est_params)
w <- grobWidth(est_params)
title_est <- textGrob("Estimated Parameter Values", y=unit(0.5,"npc") + 1.0*h, 
                  vjust=0, gp=gpar(fontsize=12))

gt_est <- gTree(children=gList(est_params, title_est))


##################################
#True Values
#################################

OM_true_table<-data.frame(Parameter=c("q","R_ave","beta 1 fishery","beta 2 fishery","beta 3 fishery","beta 4 fishery","beta 1 survey","beta 2 survey","beta 3 survey","beta 4 survey"))

#set up matrix to fill in the values
if((OM_dat[grep("fleets_as_areas_switch",OM_dat, fixed = T)+3]==1 && nreg_OM==1)){
temp_om<-data.frame(matrix(NA,dim(OM_true_table)[1],out$nfleets_OM))
names(temp_om)<-1:out$nfleets_OM
}

if(OM_dat[grep("fleets_as_areas_switch",OM_dat, fixed = T)+3]==1 && nreg_OM>1){
  temp_om<-data.frame(matrix(NA,dim(OM_true_table)[1],nreg_OM))
  names(temp_om)<-1:nreg_OM
}

if(OM_dat[grep("fleets_as_areas_switch",OM_dat, fixed = T)+3]==0){
  temp_om<-data.frame(matrix(NA,dim(OM_true_table)[1],nreg_OM))
  names(temp_om)<-1:nreg_OM
}

OM_true_table<-cbind(OM_true_table,temp_om)
OM_true_table[1,2:ncol(OM_true_table)]<-round(out$q_survey_TRUE,2)

if(npops_OM==1){
OM_true_table[2,2]<-round(out$R_ave_TRUE,2)
}

if(npops_OM>1){
OM_true_table[2,2:ncol(OM_true_table)]<-round(out$R_ave_TRUE,2)
}

# adding in selectivity values based on specification
if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==0)
{
  OM_true_table[3,2:ncol(OM_true_table)]<-NA
  OM_true_table[4,2:ncol(OM_true_table)]<-NA
  OM_true_table[5,2:ncol(OM_true_table)]<-NA
  OM_true_table[6,2:ncol(OM_true_table)]<-NA
}

if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==0)
{
  OM_true_table[7,2:ncol(OM_true_table)]<-NA
  OM_true_table[8,2:ncol(OM_true_table)]<-NA
  OM_true_table[9,2:ncol(OM_true_table)]<-NA
  OM_true_table[10,2:ncol(OM_true_table)]<-NA
}

if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==1)
{
  OM_true_table[3,2:ncol(OM_true_table)]<-round(out$sel_beta1_TRUE,2)
  OM_true_table[4,2:ncol(OM_true_table)]<-round(out$sel_beta2_TRUE,2)
  OM_true_table[5,2:ncol(OM_true_table)]<-NA
  OM_true_table[6,2:ncol(OM_true_table)]<-NA
}

if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==1)
{
  OM_true_table[7,2:ncol(OM_true_table)]<-round(out$sel_beta1_survey_TRUE,2)
  OM_true_table[8,2:ncol(OM_true_table)]<-round(out$sel_beta2_survey_TRUE,2)
  OM_true_table[9,2:ncol(OM_true_table)]<-NA
  OM_true_table[10,2:ncol(OM_true_table)]<-NA
}

if(OM_dat[grep("select_switch_EM",OM_dat, fixed = T)+4]==2)
{
  OM_true_table[3,2:ncol(OM_true_table)]<-round(out$sel_beta1_TRUE,2)
  OM_true_table[4,2:ncol(OM_true_table)]<-round(out$sel_beta2_TRUE,2)
  OM_true_table[5,2:ncol(OM_true_table)]<-round(out$sel_beta3_TRUE,2)
  OM_true_table[6,2:ncol(OM_true_table)]<-round(out$sel_beta4_TRUE,2)
}


if(OM_dat[grep("select_switch_survey_EM",OM_dat, fixed = T)+4]==2)
{
  OM_true_table[7,2:ncol(OM_true_table)]<-round(out$sel_beta1_survey_TRUE,2)
  OM_true_table[8,2:ncol(OM_true_table)]<-round(out$sel_beta2_survey_TRUE,2)
  OM_true_table[9,2:ncol(OM_true_table)]<-round(out$sel_beta3_survey_TRUE,2)
  OM_true_table[10,2:ncol(OM_true_table)]<-round(out$sel_beta4_survey_TRUE,2)
}

true_params<-tableGrob(OM_true_table,theme = ttheme_default(base_size = 11))
h <- grobHeight(true_params)
w <- grobWidth(true_params)
title_true <- textGrob("True Parameter Values", y=unit(0.5,"npc") + 1.0*h, 
                      vjust=0, gp=gpar(fontsize=12))


gt_true <- gTree(children=gList(true_params, title_true))

}


##########################################
# SAVE THE OUTPUTS 
##########################################

#provide some information on the run

if(out$npops_OM>1 & out$npops>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Metapopulation",sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Spatial Heterogeneity",sep ="\n")
}
if(out$npops_OM==1 && out$nregions_OM==1 && out$npops==1 && out$nfleets==1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","Population Structure: Panmictic", sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM==1 && out$nfleets>1 && out$nfleets>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Matching","OM Population Structure: Fleets-as-Areas", "EM Population Structure: Fleets-as-Areas", sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions==1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Panmictic", sep = "\n")
}

if(out$npops_OM>1 && out$npops==1 && out$nregions==1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Metapopulation", "EM Population Structure: Panmictic", sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM>1 && out$npops==1 && out$nregions==1 && out$nfleets>1 ){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Fleets-as-Areas", sep = "\n")
}

if(out$npops_OM==1 && out$nregions_OM>1 && out$npops>1 && out$nregions==1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Spatial Heterogeneity", "EM Population Structure: Metapopulation", sep = "\n")
}

if(out$npops_OM>1 && sum(out$nregions_OM)>1 && out$npops==1 && out$nregions>1){
  text1<-paste("MODEL STRUCTURE","Match/Mismatch: Mismatch","OM Population Structure: Metapopulation", "EM Population Structure: Spatial Heterogeneity", sep = "\n")
}



#is this a diagnostics run?
diagnostic<-OM_dat[(grep("diagnostics_switch",OM_dat, fixed = T)+3)]

if(diagnostic==0)
{text2<-"Diagnostic Run: NO. Uses OBS values as data inputs"}

if(diagnostic==1)
{text2<-"Diagnostic Run: YES. Uses TRUE values as data inputs"}

if(resid.switch==1)
{text3<-"\nResidual plots represent STRAIGHT difference; (ESTIMATED-TRUE)"}

if(resid.switch==2)
{text3<-"Residual plots represent RELATIVE difference; ((ESTIMATED-TRUE)/TRUE)*100"}


text4<-paste0("\nRun Time of Estimation: ",round((as.vector(time.elapsed)[3]/60),5)," minutes")


text.all<-paste(text1,text2,text3,text4,sep = "\n")


# Create a text grob
tgrob <- textGrob(text.all,just = "centre")


################################################
###############################################
#print these plots to pdf in the EM folder

setwd(EM_direct)

#generate pdf with plots
pdf("Model_Diagnostics.pdf",paper='letter') 

grid.arrange(ncol=1,
             top=textGrob("TIM Diagnostics", gp=gpar(fontsize=18,font=3)),
            tgrob,
            gt)

# for first layout
ly<-rbind(c(1,2),c(3,4))
grid.arrange(layout_matrix=ly,
             tgrob.mod,like.p,gt_est,gt_true,heights = c(60,60))

grid.arrange(ncol = 1,
  #top="Recruitment",
  rec1, R.resid)

grid.arrange(ncol = 1,
             #top="Recruitment",
             rec2, R.dev.resid)

grid.arrange(ncol = 1,
             #top="Recruitment ",
             rec.prop, R.apport.resid)

grid.arrange(ncol = 1,
             #top="Recruitment ",
             mort.p, m.resid.p)

grid.arrange(ncol = 1,
             #top="Reporting Rate ",
             rr.p, rr.resid.p)

grid.arrange(ncol = 1,
  #top="Initial Abundance",
  init.ab,init.ab.resid
)

grid.arrange(ncol = 1,
             #top="Abundance"
             bio.p,bio.resid)

grid.arrange(ncol = 1,
    #top="Abundance"
    ssb.p,ssb.resid)

#selectivity plots
grid.arrange(ncol = 1,
      #top="Selectivity",
      s.select.p, s.select.resid)

grid.arrange(ncol = 1,
      #top="Selectivity",
      f.select.p, f.select.resid)

#other stuff
grid.arrange(ncol = 1,
    top = "Fishing Mortality",
    F.plot.p, F.resid.plot)

#other stuff
grid.arrange(ncol = 1,
    #top = "Movement",
    T.year.p)

grid.arrange(ncol = 1,
             #top = "Movement",
            T.resid.plot)

#other stuff
grid.arrange(ncol = 1,
  top ="Fits to Data",
   yield.p,Y.resid.plot)

grid.arrange(ncol = 1,
             top ="Fits to Data continued",
             survey.p,Surv.resid.plot)

grid.arrange(ncol = 1,
    top ="Fits to Data continued",
    survey.comp.plot, fishery.comp.plot)

dev.off()

################################################
#Plot corrlation Matrix
################################################

#if cor file exists
if (file.exists(cor.name)) {

cor.mat<-as.matrix(cor)

#removing diag
cor.mat[cor.mat==1]<-0
#removing low correlations
cor.mat[abs(cor.mat)<cor.level]<-0

rm<-which(colSums(cor.mat)==0)

#remove the parameter combinations with weak correlations.
cor2<-as.matrix(cor[-rm,-rm])


#create a corplot
pdf("Correlation_Matrix.pdf", paper = "a4", width = 8, height = 11)

if(nrow(cor2)==0){
  corrplot(cor.mat, type = "upper", order = "original", 
                            tl.col = "black", tl.srt = 45,tl.cex = 0.5,diag = F)
  }

if(nrow(cor2)>0){
corrplot(cor2, type = "upper", order = "original", 
         tl.col = "black", tl.srt = 45,tl.cex = 0.5,diag = F)
}

#save as PDF
dev.off()
}


#####################
# Plot Tag residuals
######################

pdf("Tag_Residuals.pdf",paper='a4r',width=11,height=8) 
for(k in 1:nreg){
  print(tags.plot(k))}

dev.off()

#############################################
# Build table of Standard Errors
#############################################

#take the values off the log scale
#std[,3:4]<-exp(std[,3:4])

if (file.exists(std.name)) {
  
std_params<-tableGrob(std,theme = ttheme_minimal(base_size = 10),rows = NULL)

std_params<-gtable_add_grob(std_params,grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                            t = 2, b = nrow(std_params), l = 1, r = ncol(std_params)) 
std_params <-gtable_add_grob(std_params,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 1, l = 1, r = ncol(std_params))
                             
#h.st <- grobHeight(std_params)
#title_std <- textGrob("Standard Errors (log scale)", y=unit(0.5,"npc") + 1.25*h.st, 
#                       vjust=0, gp=gpar(fontsize=12))
#gt_std <- gTree(children=gList(std_params, title_std))


#print over many pages
fullheight <- convertHeight(sum(std_params$heights), "cm", valueOnly = TRUE)
margin <- unit(0.61,"in")
margin_cm <- convertHeight(margin, "cm", valueOnly = TRUE)
a4height <- 29.7 - margin_cm
nrows <- nrow(std_params)
npages <- ceiling(fullheight / a4height)

heights <- convertHeight(std_params$heights, "cm", valueOnly = TRUE) 
rows <- cut(cumsum(heights), include.lowest = FALSE,
            breaks = c(0, cumsum(rep(a4height, npages))))
groups <- split(seq_len(nrows), rows)
gl <- lapply(groups, function(id) std_params[id,])


pdf("Standard_Error_table.pdf", paper = "a4", width = 0, height = 0)
ml<-marrangeGrob(grobs=gl, ncol=1, nrow=1, top="Standard Errors (log scale)")

grid.draw(ml)
dev.off()
}


####################################
#save all the outputs to a text file
#####################################

sink(file = "output.txt",
     append = F, type = "output")

#model structure
print("$nages")
print(out$nages)
print("nyrs")
print(out$nyrs)
print("npops")
print(out$npops)
print("nregions")
print(out$nregions)

#rec params
print("$R_ave")
print(out$R_ave)

print("$R_ave_TRUE")
print(out$R_ave_TRUE)

print("$R_apport")
print(out$Rec_Prop)

print("$R_apport_True")
print(out$Rec_Prop_TRUE)

print("$steep")
print(out$steep)

print("$steep_TRUE")
print(out$steep_TRUE)

print("$q_survey")
print(out$q_survey)

print("$q_survey_TRUE")
print(out$q_survey_TRUE)

print("$sel_beta1_survey")
print(out$sel_beta1_survey)

print("$sel_beta1_survey_TRUE")
print(out$sel_beta1_survey_TRUE)

print("$sel_beta2_survey")
print(out$sel_beta2_survey)

print("$sel_beta2_survey")
print(out$sel_beta2_survey_TRUE)

print("$sel_beta1")
print(out$sel_beta1)

print("$sel_beta1_TRUE")
print(out$sel_beta1_TRUE)

print("$sel_beta2")
print(out$sel_beta2)

print("$sel_beta2_TRUE")
print(out$sel_beta2_TRUE)

print("$Rec_BM")
print(out$recruits_BM)

print("$Rec_BM_TRUE")
print(out$recruits_BM_TRUE)

print("$Rec_devs")
print(out$rec_devs)

print("$Rec_devs_TRUE")
print(out$rec_devs_TRUE)

print("$init_abund")
print(out$Init_Abund)

print("$init_abund_TRUE")
print(out$Init_Abund_TRUE)

print("$selectivity_age")
print(out$selectivity_age)

print("$selectivity_age_TRUE")
print(out$selectivity_age_TRUE)

print("$survey_selectivity_age")
print(out$survey_selectivity_age)

print("$survey_selectivity_age")
print(out$survey_selectivity_age_TRUE)

print("$report_rate")
print(out$report_rate)

if(npops==1){
print("$T_year")
print(out$T_year)

print("$T_year_TRUE")
print(out$T_year_TRUE)
}

if(npops>1){
  print("$T_year")
  print(T_est)
  print("$T_true")
  print(T_true)
}

print("$yield_fleet")
print(out$yield_fleet)

print("$OBS_yield_fleet")
print(out$OBS_yield_fleet)

print("$survey_fleet_bio")
print(out$survey_fleet_bio)

print("$OBS_survey_fleet_bio")
print(out$OBS_survey_fleet_bio)

sink()

}

make.plots()

#end running the whole code


  

