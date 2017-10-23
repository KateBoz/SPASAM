####################################################
# Simple initial code for visuallizing model outputs
# Created by JJD/KB
####################################################

#remove junk from workspace
{rm(list=(ls()))

#load libraries
load_libraries<-function() {
  library(PBSmodelling)
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(gplots)
  library(colorspace)
  library(RColorBrewer)
}
load_libraries()
}


#################################################################################
########### INPUTS FOR RUNNING MODELS ###########################################
#################################################################################

# Manually make changes in the OM .dat and run both OM and EM together if you want

  
######### USER INPUTS...NEED TO CHANGE EACH RUN ##################################  
multiple_reg<-1 #running model with multiple regions, 1==yes, 0==no
multiple_pop<-0 #running model with multiple populations, 1==yes, 0==no

#DO NOT RUN MODELS WITH MUTLIPLE REGIONS AND POPULATIONS, NOT EQUIPPED TO ESTIMATE MOVEMENT AMONG REGIONS AND POPULATIONS...YOU TRY CODING THAT ESTIMATED ARRAY

#OM Location
OM_direct<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SPASAM-master (1)\\SPASAM-master\\Tag_Integrated_Model\\Operating_Model"
OM_name<-"TIM_OM" #name of the OM you are wanting to run

#EM Location
EM_direct<-"C:\\Users\\katelyn.bosley.NMFS\\Desktop\\SPASAM-master (1)\\SPASAM-master\\Tag_Integrated_Model\\Estimation_Model" #location of run(s)
EM_name<-"TIM_EM" ###name of .dat, .tpl., .rep, etc.
########################################################################################################

##############################################
#plot color set up if what something different
##############################################
#set up the colors you want to use for TRUE and Predicted
t.col="black"
e.col="blue"

#set up the color that are wanted for movement plot - used a color ramp
mycols=colorRampPalette(c("black", "cyan","blue"))


########################################################################################################
########### AUTOMATED...DO NOT CHANGE ##################################################################
########################################################################################################
{ #run this section of code

if(multiple_reg==1)
   {  
    EM_direct<-paste0(EM_direct,'\\Mult Regs',sep="")
   }
if(multiple_pop==1)
{  
  EM_direct<-paste0(EM_direct,'\\Mult Pops',sep="")
}
if(multiple_reg==0 & multiple_pop==0)
{  
  EM_direct<-paste0(EM_direct,'\\panmictic',sep="")
}
if(multiple_reg==1 & multiple_pop==1)
{  
  print("Can't Estimate Movement Among Multiple Populations and Regions Simultaneously")
  stop()
}

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
invisible(shell(paste0(EM_name),wait=T))

}#end running the OM/EM together



#########################################################
##### Look at the outputs ###############################
#########################################################

make.plots<-function(direct=EM_direct){ #run diagnostics plotting
  
#Read in model .rep
out<-readList(paste(EM_direct,paste0(EM_name,".rep"),sep="\\")) #read in .rep file


#pull info about the model
na<-out$nages
nyrs<-out$nyrs
npops<-out$npops
nreg<-out$nregions
years<-seq(1:out$nyrs)
ages<-seq(1:out$nages)


#to use the current code with meta pop example
if(npops>1){
  nreg=npops}

##############################
#### RECRUITMENT PLOTS #######
##############################

#rec total
rec.total<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Rec_Est = as.vector(t(out$recruits_BM)), Rec_True=as.vector(t(out$recruits_BM_TRUE)))

rec.total.plot<-melt(rec.total,id=c("Reg","Years"))
rec.total.plot$Reg<-as.factor(rec.total.plot$Reg)


rec1<-ggplot(rec.total.plot,aes(Years,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Recruitment")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Recruitment Total")

  
#rec devs
if(npops>1){
rec.devs<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),Rec_Est = as.vector(t(out$rec_devs)*out$R_ave), Rec_True=as.vector(t(out$rec_devs_TRUE[2:length(out$rec_devs_TRUE)]*out$R_ave_TRUE)))

rec.devs.plot<-melt(rec.devs,id=c("Reg","Years"))
rec.devs.plot$Reg<-as.factor(rec.devs.plot$Reg)

rec2<-ggplot(rec.total.plot,aes(Years,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  theme_bw()+
  facet_wrap(~Reg)+
  scale_color_manual(values = c("black","grey60"),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  geom_hline(yintercept = out$R_ave, col = e.col)+
  geom_hline(yintercept = out$R_ave_TRUE,col=t.col,lty=2)+
  ylab("Recruitment")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Recruitment Deviations")
}

if(npops==1){
  rec.devs<-data.frame(Years=years[-1],Rec_Est=out$rec_devs*out$R_ave, Rec_True=out$rec_devs_TRUE[2:length(out$rec_devs_TRUE)]*out$R_ave_TRUE)
  
rec.devs.plot<-melt(rec.devs,id=c("Years"))
  
rec2<-ggplot(rec.devs.plot,aes(Years,value))+
    geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
    theme_bw()+
    #facet_wrap(~Reg)+
    scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
    scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
    ylab("Recruitment")+
    geom_hline(yintercept = out$R_ave, col = e.col)+
    geom_hline(yintercept = out$R_ave_TRUE,col=t.col,lty=2)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          legend.background = element_rect(fill="transparent"),
          panel.border = element_rect(colour = "black"))+
    theme(legend.title = element_blank())+
    theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
    theme(legend.position = c(1, 1), legend.justification = c(1,1))+
    ggtitle("Recruitment Deviations")
  
}


#################################
#### ABUNDANCE PLOTS ############
#################################

#SSB
ssb.dat<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),SSB_est = as.vector(t(out$SSB_region)), SSB_True=as.vector(t(out$SSB_region_TRUE)))

ssb.plot<-melt(ssb.dat,id=c("Reg","Years"))
ssb.plot$Reg<-as.factor(ssb.plot$Reg)

ssb.p<-ggplot(ssb.plot,aes(Years,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("SSB")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("SSB")



#init abundance
init.abund<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),In_ab_Est = as.vector(t(out$Init_Abund)), In_ab_True=as.vector(t(out$Init_Abund_TRUE)))


init.abund.plot<-melt(init.abund,id=c("Reg","Age"))
init.abund.plot$Reg<-as.factor(init.abund.plot$Reg)

init.ab<-ggplot(init.abund.plot,aes(Age, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Abundance")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Initial Abundance")


############################################################
########## Selectivity Parameters ##########################
############################################################

#Fishery Selectivity
f.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = as.vector(t(out$selectivity_age)), Select_T=as.vector(t(out$selectivity_age_TRUE)))


f.select.plot<-melt(f.select,id=c("Reg","Age"))
f.select.plot$Reg<-as.factor(f.select$Reg)

f.select.p<-ggplot(f.select.plot,aes(Age, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Selectivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 0), legend.justification = c(1,0))+
  ggtitle("Fishery Selectivity")


#Survey Selectivity
s.select<-data.frame(Age=rep(ages,nreg), Reg=rep(c(1:nreg),each=na),Select_Est = as.vector(t(out$survey_selectivity_age)),Select_T=as.vector(t(out$survey_selectivity_age_TRUE)))


s.select.plot<-melt(s.select,id=c("Reg","Age"))
s.select.plot$Reg<-as.factor(s.select$Reg)

s.select.p<-ggplot(s.select.plot,aes(Age, value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("Selectivity")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 0), legend.justification = c(1,0))+
  ggtitle("Survey Selectivity")



############################################################
########## F by year #######################################
############################################################


#might need to change the dims if working with multi pop

F.year<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs),F_year=out$F_year, F_year_T=out$F_year_TRUE)

F.plot<-melt(F.year,id=c("Reg","Years"))
F.plot$Reg<-as.factor(F.plot$Reg)

F.plot.p<-ggplot(F.plot,aes(Years,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  theme_bw()+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","True"))+
  scale_linetype_manual(values=c(1,2),labels = c("Estimated","True"))+
  ylab("F")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Fully selected F by Year")


############################################################
########## Movement by year ################################
############################################################

#not really generalized yet..the best I can do for now. Will allow for yearly movement
T_est<-data.frame(out$T_year)
T_true<-data.frame(out$T_year_TRUE)

for(i in 1:nreg){
names(T_est)[i]<-paste0("Est_",i)
names(T_true)[i]<-paste0("True_",i)}

T.year<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs)) 
T.year<-cbind(T.year,T_est,T_true)

T.year.plot<-melt(T.year, id=c("Reg","Years"))

T.year.plot$Reg<-as.factor(T.year.plot$Reg)
T.lines<-c(rep(1,nreg),rep(2,nreg))
T.col<-rep(mycols(nreg),nreg)

T.year.p<-ggplot(T.year.plot,aes(Years,value))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  theme_bw()+
  facet_wrap(~Reg)+
  scale_color_manual(values = T.col)+
  #scale_linetype_manual(values=T.lines,guide=FALSE)+
  scale_linetype_manual(values=T.lines)+
  ylab("Movement Rate")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = "right", legend.justification = c(1,1))+
  ggtitle("Yearly Movement Rate")

############################################################
########## Fits to data ####################################
############################################################

# Yields

Y.year<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), Estimated=out$yield_fleet,Observed=out$OBS_yield_fleet ) 

Y.year.plot<-melt(Y.year, id=c("Reg","Years"))

yield.p<-ggplot(Y.year.plot,aes(Years,value,shape=variable))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  geom_point(size=2, alpha = 0.5)+
  scale_shape_manual(values=c(NA,16),labels = c("Estimated","Observed"))+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","Observed"))+
  scale_linetype_manual(values=c(1,0),labels = c("Estimated","Observed"))+
  ylab("Yield")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Yield")


# Survey Index

Survey.year<-data.frame(Years=rep(years,nreg), Reg=rep(c(1:nreg),each=nyrs), SI_Est=out$survey_fleet_bio,SI_True=out$OBS_survey_fleet_bio ) 

Survey.year.plot<-melt(Survey.year, id=c("Reg","Years"))

survey.p<-ggplot(Survey.year.plot,aes(Years,value,shape=variable))+
  geom_line(aes(col = variable,linetype=variable), stat = "identity", lwd=1)+
  geom_point(size=2, alpha = 0.5)+
  scale_shape_manual(values=c(NA,16),labels = c("Estimated","Observed"))+
  facet_wrap(~Reg)+
  scale_color_manual(values = c(e.col,t.col),labels = c("Estimated","Observed"))+
  scale_linetype_manual(values=c(1,0),labels = c("Estimated","Observed"))+
  ylab("Survey Index")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.background = element_rect(fill="transparent"),
        panel.border = element_rect(colour = "black"))+
  theme(legend.title = element_blank())+
  theme(strip.text.x = element_text(size = 10, colour = "black", face="bold"))+
  theme(legend.position = c(1, 1), legend.justification = c(1,1))+
  ggtitle("Survey Biomass")




# age comps coming soon...


#print these plots to pdf in the EM folder
setwd(direct)

#generate pdf with plots
pdf("Model_Diagnostics.pdf",width=6,height=4,paper='special') 
print(rec1)
print(rec2)
print(init.ab)
print(ssb.p)
print(s.select.p)
print(f.select.p)
print(F.plot.p)
print(T.year.p)
print(yield.p)
print(survey.p)

dev.off()


#save the outputs to a text file
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

print("$selectivity_age")
print(out$selectivity_age)

print("$selectivity_age_TRUE")
print(out$selectivity_age_TRUE)

print("$survey_selectivity_age")
print(out$survey_selectivity_age)

print("$survey_selectivity_age")
print(out$survey_selectivity_age_TRUE)

print("$T_year")
print(out$T_year)

print("$T_year_TRUE")
print(out$T_year_TRUE)

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

###################################
#JJD plots
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

