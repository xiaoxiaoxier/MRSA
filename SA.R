rm(list=ls())#1234512345
library(plyr)
library(ROCR)
library(pROC)
library(sva)
library(ggplot2)
library(MASS)
library(RColorBrewer)
library(keras)
library(VennDiagram)
library(svglite)
SA=read.csv("SA.csv")
colcha=c('Intensity','mz')
SA[,colcha]=lapply(SA[,colcha],as.character)
data=na.omit(subset(SA,select=c(Lab_number,mz,Intensity,GENDER,AGE,date_Day,Penicillin,Oxacillin,Clindamycin,Erythromycin,
                                Fusidic_acid,SXT)))
######pie chart for SA with the ratio of R/S#########
par(mfrow=c(6,1),mar = c(0,0,0,0))
table(SA$Penicillin)
df=data.frame(value=c(25151,1701),group=c("R","S"))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

table(SA$Oxacillin)
df=data.frame(value=c(13860,12992),group=c("R","S"))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

table(SA$Erythromycin)
df=data.frame(value=c(15068,11784),group=c("R","S"))
#df=arrange(df,desc(value))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

table(SA$Clindamycin)
df=data.frame(value=c(11638,15214),group=c("R","S"))
#df=arrange(df,desc(value))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

table(SA$Fusidic_acid)
df=data.frame(value=c(2233,24619),group=c("R","S"))
#df=arrange(df,desc(value))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

table(SA$SXT)
df=data.frame(value=c(3794,23058),group=c("R","S"))
#df=arrange(df,desc(value))
df$color=c(brewer.pal(8,"Blues")[7],brewer.pal(8,"Blues")[2])
df =df[c(2:nrow(df),1),]
labs = paste0(df$group,"\n(",round(df$value/sum(df$value)*100,2),"%)")
pie(df$value,labels = labs,init.angle=90,col=df$color,border="black")

######bar chart for resistance group#########
par(mfrow=c(2,1),mar = c(4,4,2,2))
SA_resistance <- c(25151,15068,13860,11683,3794,2233)

anti=c("Penicillin","Erythromycin","Oxacillin","Clindamycin","Sulfamethoxazole-Trimethoprim","Fusidic acid")
# Plot the bar chart.
barplot(SA_resistance,names.arg=anti,xlab="Antibiotics",ylab="The number of resistant samples",
        col=brewer.pal(8,"Blues")[7],ylim=c(0,25000),space=0.5)
######bar chart for resistance to 1/2/3/4/5/6 antibiotics#########
data_1=subset(SA,select=c(Penicillin,Oxacillin,Clindamycin,Erythromycin,Fusidic_acid,SXT))
length(which(rowSums(data_1=="R")==1))#7766
length(which(rowSums(data_1=="R")==2))#4279
length(which(rowSums(data_1=="R")==3))#3078
length(which(rowSums(data_1=="R")==4))#6352
length(which(rowSums(data_1=="R")==5))#2940
length(which(rowSums(data_1=="R")==6))#1013
length(which(rowSums(data_1=="R")==0))#1424

data_1_resistance <- c(1424,7766,4279,3078,6352,2940,1013)

anti_num=c("0","1","2","3","4","5","6")
# Plot the bar chart.
barplot(data_1_resistance,names.arg=anti_num,xlab="MDR",ylab="The number of resistant samples",
        col=brewer.pal(8,"Blues")[7],ylim=c(0,8000),space=0.5)

#--------------plot Venn Diagram--------------------
Penicillin_R=as.character(data$Lab_number[which(data$Penicillin=="R")])
Erythromycin_R=as.character(data$Lab_number[which(data$Erythromycin=="R")])
Oxacillin_R=as.character(data$Lab_number[which(data$Oxacillin=="R")])
Clindamycin_R=as.character(data$Lab_number[which(data$Clindamycin=="R")])
SXT_R=as.character(data$Lab_number[which(data$SXT=="R")])
Fusidicacid_R=as.character(data$Lab_number[which(data$Fusidic_acid =="R")])
myplot=venn.diagram(list(Penicillin=Penicillin_R,Oxacillin=Oxacillin_R),fill=c(brewer.pal(7,"Set2")[1:2]),resolution = 600,filename = NULL)
myplot=venn.diagram(list(Erythromycin=Erythromycin_R,Clindamycin=Clindamycin_R),fill=c(brewer.pal(7,"Set2")[3:4]),resolution = 600,filename = NULL)
myplot=venn.diagram(list(Erythromycin=Erythromycin_R,Clindamycin=Clindamycin_R,Fusidic_acid=Fusidicacid_R),fill=c(brewer.pal(7,"Set2")[3:5]),resolution = 600,filename = NULL)
myplot=venn.diagram(list(Oxacillin=Oxacillin_R,Clindamycin=Clindamycin_R,SulfamethoxazoleTrimethoprim=SXT_R),fill=c(brewer.pal(7,"Set2")[2],brewer.pal(7,"Set2")[4],brewer.pal(7,"Set2")[6]),resolution = 600,filename = NULL)
ggsave(myplot, file="my_plot.svg", device = "svg")

#-------------plot year informationwith R or S ---------------------
data$date_Day=as.factor(substr(data$date_Day,1,4))
ggplot(data,aes(x=factor(date_Day)))+geom_bar(stat ="count",width=0.6,fill="steelblue",colour="black")+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("the number of samples")+xlab("sample collection year")+ylim(0, 6000)

DF<-ddply(data,.(date_Day),transform,P=length(which(Penicillin=="R"))/length(Penicillin),OX=length(which(Oxacillin=="R"))/length(Oxacillin),
          E=length(which(Erythromycin=="R"))/length(Erythromycin),CC=length(which(Clindamycin=="R"))/length(Clindamycin),
          FA=length(which(Fusidic_acid=="R"))/length(Fusidic_acid),SXT1=length(which(SXT=="R"))/length(SXT)) 
ggplot(DF,aes(x=date_Day,y=P,color=date_Day))+geom_line(size=2)

write.csv(subset(DF,select=c(date_Day,P,E,CC,OX,FA,SXT1)),"figure_year.csv")

ratio=read.csv("figure_year.csv")
ggplot(ratio,aes(x=factor(Year),y=ratio_R,color=Antibiotic,group=Antibiotic))+geom_line(size=1)+geom_point(size=1.5)+xlab("sample collection year")+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  ylab("the ratio of resistance samples")+ylim(0, 1)+scale_color_manual(values=c(brewer.pal(7, "Set1")[1:5],brewer.pal(7, "Set1")[7]))
#-------------------- Multivariable logistic regression analysis
mylogit <- glm(formula = Erythromycin ~ date_Day + GENDER + AGE,
              data=data, family=binomial)
summary(mylogit)                                                                

lreg.or <-exp(cbind(OR = coef(mylogit), confint(mylogit)))
round(lreg.or, digits=4)

#-------------------------------------------------------------
data_OX_R=subset(data,Oxacillin=="R")
data_OX_S=subset(data,Oxacillin=="S")
test_mz_R=data_OX_R$mz
test_in_R=data_OX_R$Intensity
test_mz_S=data_OX_S$mz
test_in_S=data_OX_S$Intensity
myFun=function(x){strsplit(x,split=';')}
mz_R=lapply(test_mz_R, myFun)
Intens_R=lapply(test_in_R, myFun)
mz_S=lapply(test_mz_S, myFun)
Intens_S=lapply(test_in_S, myFun)

mz_num_R=ceiling(as.numeric(unlist(mz_R)))
mz_num_S=ceiling(as.numeric(unlist(mz_S)))
hist(mz_num_R,col="red",border = "red",xlim = c(2000,20000),breaks = 18000,ylim = c(0,5000))
hist(mz_num_S,col="blue",border = "blue",xlim = c(2000,20000),breaks = 18000,ylim = c(0,5000))
freq_R=table(mz_num_R)
freq_S=table(mz_num_S)

Intens_num_R=as.numeric(unlist(Intens_R))
mz_num_R=ceiling(as.numeric(unlist(mz_R)))
mz_num_R=as.numeric(unlist(mz_R))


Intens_num_S=as.numeric(unlist(Intens_S))
mz_num_S=ceiling(as.numeric(unlist(mz_S)))
mz_num_S=as.numeric(unlist(mz_S))
par(mfrow=c(2,1))
plot(mz_num_S[which(mz_num_S<=(2000+i)&mz_num_S>(1999+i))],log(Intens_num_S[which(mz_num_S<=(2000+i)&mz_num_S>(1999+i))]),pch=19,cex=0.1,ylim = c(2,12))
plot(mz_num_R[which(mz_num_R<=(2000+i)&mz_num_R>(1999+i))],log(Intens_num_R[which(mz_num_R<=(2000+i)&mz_num_R>(1999+i))]),pch=19,cex=0.1,ylim = c(2,12))

p=vector()
location_diff=vector()
for (i in 1:17269) 
  {
  if(length(which(mz_num_S<=(2000+i)&mz_num_S>(1999+i)))<100||length(which(mz_num_R<=(2000+i)&mz_num_R>(1999+i)))<100)
  {p[i]=1}
  
  else
  {
  x=log(Intens_num_S[which(mz_num_S<=(2000+i)&mz_num_S>(1999+i))])
  y=log(Intens_num_R[which(mz_num_R<(2000+i)&mz_num_R>(1999+i))])
  kk=wilcox.test(x, y,alternative = c("two.sided"),paired = FALSE,conf.int = TRUE, conf.level = 0.95)
  p[i]=kk$p.value
  location_diff[i]=kk$estimate
  print(i)
  }
}

hist(p)
which(p<0.0000000005)
# [1]  208  209  210  211  212  226  227  229  230  265  386  387  388  389  390  396  397  398  399  400
# [21]  412  413  414  415  416  417  418  429  430  431  432  433  434  435  436  437  438  439  450  451
# [41]  452  453  454  455  456  491  546  547  548  643  762  765  851  852  878  879  880  881  917  918
# [61]  919  979  980  981 1006 1007 1008 1011 1029 1030 1031 1036 1037 1038 1039 1040 1045 1046 1047 1210
# [81] 1211 1212 1218 1219 1220 1221 1222 1265 1277 1278 1295 1296 1297 1298 1299 1305 1306 1307 1452 1453
# [101] 1454 1455 1489 1490 1640 1641 1763 1764 1765 1766 1769 1770 2095 2096 2097 2486 2487 2488 2498 2499
# [121] 2568 2807 2808 2809 3285 3286 3438 3541 3777 3778 3779 3780 4420 4421 4422 4499 4500 4501 4521 4522
# [141] 4523 4524 4525 4526 4527 4528 4549 4550 4551 4552 4569 4570 4571 4572 4591 4592 4593 4594 4595 4596
# [161] 4597 4598 4632 4633 4634 4635 4886 4887 7613

boxplot(log(Intens_num_R[which(mz_num_R<=2212&mz_num_R>2207)]),outline=FALSE,ylim = c(3,9))
boxplot(log(Intens_num_S[which(mz_num_S<=2212&mz_num_S>2207)]),outline=FALSE,ylim = c(3,9))
ks.test(x,alternative = c("two.sided"))

hist(log(Intens_num_R[which(mz_num_R<=2212&mz_num_R>2207)]))
hist(log(Intens_num_S[which(mz_num_S<=2212&mz_num_S>2207)]))

x=log(Intens_num_R[which(mz_num_R<=2212&mz_num_R>2207)])#4070
y=log(Intens_num_S[which(mz_num_S<=2212&mz_num_S>2207)])#2571
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2207~2212)")

x=log(Intens_num_R[which(mz_num_R<=2230&mz_num_R>2225)])#2704
y=log(Intens_num_S[which(mz_num_S<=2230&mz_num_S>2225)])#2630
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2225~2230)")


x=log(Intens_num_R[which(mz_num_R<=2390&mz_num_R>2385)])#5479
y=log(Intens_num_S[which(mz_num_S<=2390&mz_num_S>2385)])#3502
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2385~2390)")

x=log(Intens_num_R[which(mz_num_R<=2400&mz_num_R>2395)])#4880
y=log(Intens_num_S[which(mz_num_S<=2400&mz_num_S>2395)])#3164
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2395~2400)")

x=log(Intens_num_R[which(mz_num_R<=2418&mz_num_R>2411)])#5725
y=log(Intens_num_S[which(mz_num_S<=2418&mz_num_S>2411)])#2949
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2411~2418)")

x=log(Intens_num_R[which(mz_num_R<=2439&mz_num_R>2428)])#10300
y=log(Intens_num_S[which(mz_num_S<=2439&mz_num_S>2428)])#8227
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2428~2439)")


x=log(Intens_num_R[which(mz_num_R<=2456&mz_num_R>2449)])#7329
y=log(Intens_num_S[which(mz_num_S<=2456&mz_num_S>2449)])#5182
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2449~2456)")

x=log(Intens_num_R[which(mz_num_R<=2548&mz_num_R>2545)])#5697
y=log(Intens_num_S[which(mz_num_S<=2548&mz_num_S>2545)])#4420
Oxacillin=c(rep("R",length(x)),rep("S",length(y)))
Intensity_OX=c(x,y)
simulate_gamma<- data.frame(Oxacillin=Oxacillin, Intensity=Intensity_OX)
ggplot(simulate_gamma,aes(x = Intensity, colour=Oxacillin)) + scale_colour_manual(values=c('red','black'))+
  geom_histogram(aes(y=..density.., fill=Oxacillin), alpha=0.6,position="identity",colour="NA",binwidth=0.1) +
  stat_density(geom = "line",position = "identity",size=1) + scale_fill_manual(values=c('red','black'))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  xlab("log(Intensity)   m/z(2545~2548)")
