library(readr)
setwd("C:/Users/Dana/Downloads/")
df<-read.csv("MSY_results.csv")
head(df)

df2<-df[,c(3,4,5,11)]
df2[,1]<-jitter(as.numeric(df2[,1]))
df2[,2]<-jitter(as.numeric(df2[,2]))
df2[,3]<-jitter(as.numeric(df2[,3]))

df2$mean<-rowMeans(df[,c(3,4,5)])
df2$YieldCut<-cut(df2$yield_total,breaks=10)
levels(df2$YieldCut)<-as.character(seq(0.1,1,by=0.1))
df3<-df2[which(df2$yield_total>0.9*max(df2$yield_total)),]
df3$YieldCut<-cut(df3$yield_total,5)
ggplot(df3, aes(factor(YieldCut), F.1)) + 
  geom_violin(aes(y = F.1, colour = "1"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.2, colour = "2"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.3, colour = "3"), fill = NA, size=0.5)+
  ylab("Fishing mortality")+theme_bw()+scale_color_manual(c("Area"),values=c("red","blue","dark green"))+
  xlab("Proportion of maximum yield")+ggtitle("Combinations of areal fishing mortality")

dm<-melt(df2,id=c("yield_total","mean","YieldCut"))
names(dm)<-c("Yield","MeanF","YieldCut","Area","F")         
ggplot(dm)+geom_point(aes(x=F,y=Yield,color=Area,shape=Area),
        size=0.6,alpha=0.4)+theme_bw()+
      geom_line(aes(x=MeanF,y=Yield))+stat_smooth()
cut(dm$Yield,breaks=10)
p <- ggplot(dm, aes(factor(YieldCut), F)) + 
  geom_violin(aes(colour = ), fill = NA, size=1) + 
  geom_violin(aes(y = coco2, colour = "2"), fill = NA, size=2) +
  geom_violin(aes(y = coco3, colour = "3"), fill = NA, size=2) +
  geom_violin(aes(y = coco4, colour = "4"), fill = NA, size=2)

ggplot(df2, aes(factor(YieldCut), F.1)) + 
  geom_violin(aes(y = F.1, colour = "1"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.2, colour = "2"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.3, colour = "3"), fill = NA, size=0.5)+
  ylab("Fishing mortality")+theme_bw()+scale_color_manual(c("Area"),values=c("red","blue","dark green"))+
  xlab("Proportion of maximum yield")+ggtitle("Combinations of areal fishing mortality")

ggsave("violin.png",dpi=300)

setwd("C:/Users/Dana/Downloads/")
df<-read.csv("MSY_results_menhaden_nomove.csv")
head(df)

df2<-df[,c(3,4,5,11)]
df2[,1]<-jitter(as.numeric(df2[,1]))
df2[,2]<-jitter(as.numeric(df2[,2]))
df2[,3]<-jitter(as.numeric(df2[,3]))

df2$mean<-rowMeans(df[,c(3,4,5)])
df2$YieldCut<-cut(df2$yield_total,breaks=10)
levels(df2$YieldCut)<-as.character(seq(0.1,1,by=0.1))

dm<-melt(df2,id=c("yield_total","mean","YieldCut"))
names(dm)<-c("Yield","MeanF","YieldCut","Area","F")         
ggplot(dm)+geom_point(aes(x=F,y=Yield,color=Area,shape=Area),
                      size=0.6,alpha=0.4)+theme_bw()+
  geom_line(aes(x=MeanF,y=Yield))+stat_smooth()
cut(dm$Yield,breaks=10)
p <- ggplot(dm, aes(factor(YieldCut), F)) + 
  geom_violin(aes(colour = ), fill = NA, size=1) + 
  geom_violin(aes(y = coco2, colour = "2"), fill = NA, size=2) +
  geom_violin(aes(y = coco3, colour = "3"), fill = NA, size=2) +
  geom_violin(aes(y = coco4, colour = "4"), fill = NA, size=2)

ggplot(df2, aes(factor(YieldCut), F.1)) + 
  geom_violin(aes(y = F.1, colour = "1"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.2, colour = "2"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.3, colour = "3"), fill = NA, size=0.5)+
  ylab("Fishing mortality")+theme_bw()+scale_color_manual(c("Area"),values=c("red","blue","dark green"))+
  xlab("Proportion of maximum yield")+ggtitle("Combinations of areal fishing mortality")

ggsave("violin_men.png",dpi=300)

df3<-df2[which(df2$yield_total>0.9*max(df2$yield_total)),]
df3$YieldCut<-cut(df3$yield_total,5)
ggplot(df3, aes(factor(YieldCut), F.1)) + 
  geom_violin(aes(y = F.1, colour = "1"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.2, colour = "2"), fill = NA, size=0.5) +
  geom_violin(aes(y = F.3, colour = "3"), fill = NA, size=0.5)+
  ylab("Fishing mortality")+theme_bw()+scale_color_manual(c("Area"),values=c("red","blue","dark green"))+
  xlab("Proportion of maximum yield")+ggtitle("Combinations of areal fishing mortality")

