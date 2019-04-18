
rm(list=(ls()))



####### Where to store figures ############
Figure_dir<-"F:\\NOAA FILES\\Outreach\\Journal Submissions\\2019\\CAPAM SI\\Tagging manuscript\\manuscript\\Figures\\figures"


#load required libraries

  library(ggplot2)
  library(gridExtra)
  library(grid)

setwd(Figure_dir)
sel1<-c(0.05, 0.10, 0.25, 0.68, 0.95, 0.99, 1.00, 1.00,0.01, 0.05, 0.37, 0.84, 0.97, 1.00, 1.00, 1.00)
wt1<-c(0.1, 0.75, 1.5, 3.0, 5.0, 7.0, 7.5, 7.7,0.5, 0.65, 2.0, 4.0, 5.5, 7.2, 7.9, 8.3)
ages<-rep(seq(1:8),times=2)
sel3<-as.data.frame(cbind(ages,sel1))
sel3$pop<-rep(c("pop1","pop2"),each=8)
colnames(sel3)<-c("ages","sel","pop")

wt3<-as.data.frame(cbind(ages,wt1))
wt3$pop<-rep(c("pop1","pop2"),each=8)
colnames(wt3)<-c("ages","wt","pop")

my_theme.fig4<-
  (#theme_bw()+
    theme(
      plot.margin = unit(c(0, .15, 0, .15), "cm"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(colour="black",size=6,angle=0,hjust=0.45,vjust=1,face="plain"),
      axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0.3,face="plain"),
      axis.title.x = element_text(size=6),
      axis.title.y = element_text(size=6),
      #strip.text.y = element_text(angle = 180),
      plot.title=element_text(hjust=0.5,vjust=0,size=8,margin=margin(2,2,3,2)),
      plot.subtitle=element_text(hjust=0.5,size=8),
      strip.text.x = element_text(margin = margin(0,0,.15,0, "cm"),
                                  size = 5, color = "black"),
      strip.text.y = element_text(
        size = 5, color = "black", face = "italic"),
      strip.background = element_rect(
        fill=NA),
      legend.position = c(0.5, 0.5),
      legend.text=element_text(size=5),
      legend.title=element_text(size=5),
      legend.key.size = unit(.15,"cm"),
      panel.background = element_blank(), #element_rect(colour = "black", size=.75, fill=NA),
      #axis.line.x = element_line(),
      #axis.line.y = element_line(),
      legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid',size=.25),
      #legend.title=element_text(),
      legend.margin=margin(c(1,2,1,2)),
      panel.border = element_rect(colour = "black", fill=NA, size=.75)
    ))

wt.plot<-ggplot(wt3, aes(x=as.factor(ages), y=wt)) +
  geom_line(data=subset(wt3,pop %in% 'pop1'), aes(x=ages,y=wt,color="Population 1"),lty=1) + 
  geom_line(data=subset(wt3,pop %in% 'pop2'), aes(x=ages,y=wt,color="Population 2"),lty=3) + 
  scale_x_continuous(breaks=seq(1:8))+ #,labels=c("1","5","10","15","20","25","30"))+
  ylab('Weight (kg)')+
  xlab('Age')+
  #ggtitle("Weight-at-Age")+
  ylim(0,9)+
  my_theme.fig4+ 
  labs(color="Population")+
  scale_colour_manual(values = c("black","black"))+
  theme(legend.position=c(.5,.9),legend.title=element_blank(),legend.margin=margin(c(4,4,4,4)),legend.text=element_text(size=8),
        legend.key.size = unit(.3,"cm"),legend.spacing.x = unit(.2, 'cm'))+
  guides (colour = guide_legend (ncol=2,override.aes = 
                                   list(linetype = c(1,3))))
  
  sel.plot<-ggplot(sel3, aes(x=as.factor(ages), y=sel)) +
    geom_line(data=subset(sel3,pop %in% 'pop1'), aes(x=ages,y=sel),lty=1) + 
    geom_line(data=subset(sel3,pop %in% 'pop2'), aes(x=ages,y=sel),lty=3) + 
    scale_x_continuous(breaks=seq(1:8))+ #,labels=c("1","5","10","15","20","25","30"))+
    ylab('Selectivity')+
    xlab('Age')+
    #ggtitle("Weight-at-Age")+
    ylim(0,1)+
    my_theme.fig4+ 
   # labs(color="Population")+
    scale_colour_manual(values = c("black","black"))
   # theme(legend.position=c(.5,.9),legend.title=element_blank(),legend.margin=margin(c(4,4,4,4)),legend.text=element_text(size=8),
    #      legend.key.size = unit(.3,"cm"),legend.spacing.x = unit(.2, 'cm'))+
    #guides (colour = guide_legend (ncol=2,override.aes = 
                                 #    list(linetype = c(1,3))))
  

  
  tiff("Fig 1_Operating Model Weight and Selectivity.tif",width=190,height=100, unit='mm',res=500)
par(mar=c(6,6,6,6))

grid.arrange(ncol = 2,top=textGrob("Operating Model Weight and Selectivity",gp=gpar(fontsize=15,font=3)), 
             #bottom=textGrob("Year",gp=gpar(fontsize=12,font=1)), 
             arrangeGrob(wt.plot,top=textGrob("Weight-at-Age",gp=gpar(fontsize=9,font=2),hjust=0.22),nrow=1), 
             arrangeGrob(sel.plot,top=textGrob("Selectivity-at-Age",gp=gpar(fontsize=9,font=2),hjust=0.37),nrow=1)
)
dev.off()
