#Staphylococcus aureus positive data

###########################################################
rm(list=ls())

#5186 samples ----remove all 58 replicates
kaohsiung=read.csv("kaohsiung_filter.csv")
kaohsiung=kaohsiung[,-1]
col2cha=c('Intensity','mz')
kaohsiung[,col2cha]=lapply(kaohsiung[,col2cha],as.character)
#Consider SXT(Trimethoprim/Sulfamethoxazole) for Staphylococcus aureus 
#---------------------------------------------------------------------
data_SXT=na.omit(subset(kaohsiung,select=c(Lab_number,mz,Intensity,Gender,AGE,SXT)))
#data_RSA=data_SXT[which(data_SXT$SXT=="R"),]##454 samples
#data_SSA=data_SXT[which(data_SXT$SXT=="S"),]##4674 samples
#Calculate the frequence of the peaks of every bins between 2000-20000 with 1 interval
#iron 1:2000-2001     peptide iron:19999-20000
#----------------------------------------------
mz=as.numeric(strsplit(data_SXT$mz,split=';')[[1]])

intensity=as.numeric(strsplit(data_SXT$Intensity,split=';')[[1]])
n=0
label=ceiling(mz[which(intensity>n)]-2000)
vec=rep(0, times=18000)
vec[label]=1
for (i in 2:nrow(data_SXT))
{
  vec1=vec
  z=as.numeric(strsplit(data_SXT$mz,split=';')[[i]])
  
  y=as.numeric(strsplit(data_SXT$Intensity,split=';')[[i]])
  
  if(sum(which(y>n))>0)
  {label=ceiling(z[which(y>n)]-2000)
  vec=rep(0, times=18000)
  vec[label]=1
  vec=rbind(vec1,vec)}
  else{vec=rep(0, times=18000)
  vec=rbind(vec1,vec)
  }
  print(i)
}
colnames(vec)=paste0("iron",2001:20000)
rownames(vec)=data_SXT$Lab_number

peak_kaohsiung=vec
#######peak_kaohsiung was made and saved in peak_kaohsiung.csv
#######same step for peak_Linkou.csv
#-----------------input set from linkou--------
Linkou=read.csv("Linkou_filter.csv")##21130 samples
Linkou=Linkou[,-1]
col2cha=c('Intensity','mz')
Linkou[,col2cha]=lapply(Linkou[,col2cha],as.character)

SA_mz=as.numeric(strsplit(Linkou$mz,split=';')[[1]])

SA_intensity=as.numeric(strsplit(Linkou$Intensity,split=';')[[1]])
n=0
label=ceiling(SA_mz[which(SA_intensity>n)]-2000)
vec=rep(0, times=18000)
vec[label]=1

for (i in 2:nrow(Linkou))
{
  vec1=vec
  z=as.numeric(strsplit(Linkou$mz,split=';')[[i]])
  
  y=as.numeric(strsplit(Linkou$Intensity,split=';')[[i]])
  
  if(sum(which(y>n))>0)
  {label=ceiling(z[which(y>n)]-2000)
  vec=rep(0, times=18000)
  vec[label]=1
  vec=rbind(vec1,vec)}
  else{vec=rep(0, times=18000)
  vec=rbind(vec1,vec)
  }
  print(i)
}

colnames(vec)=paste0("iron",2001:20000)
rownames(vec)=Linkou$Lab_number

peak_Linkou=vec
#######peak_Linkou was made and saved in peak_Linkou.csv
#peak_kaohsiung(Linkou) and kaohsiung(Linkou) have the same row name.
peak_kaohsiung_R=peak_kaohsiung[which(kaohsiung$SXT=="R"),]##454 samples
peak_kaohsiung_S=peak_kaohsiung[which(kaohsiung$SXT=="S"),]##4674 samples
peak_whole=t(rbind(peak_kaohsiung_R,peak_kaohsiung_S))
peak_whole_filter=peak_whole[which(rowSums(peak_whole)>0),]




#---of splitting the data into two parts (notnecessarily of equal size), with one part used to divine
#the fitted model and the other part reserved for validation




group=c(rep(1,nrow(peak_kaohsiung_R)),rep(0,nrow(peak_kaohsiung_S)))
age=kaohsiung$AGE[match(colnames(peak_whole_filter),kaohsiung$Lab_number)]
whole=data.frame(t(peak_whole_filter),age,group=as.factor(group))
p=vector()
for (i in 1:nrow(peak_whole_filter)) {
  p[i]=chisq.test(whole[,i],group)$p.value
}
pvalue=data.frame(p)
rownames(pvalue)=rownames(peak_whole_filter)
DE_iron=rownames(peak_whole_filter)[which(p<0.00000000000005)]
rownames(pvalue)=rownames(peak_whole_filter)
DE_peak=peak_whole_filter[match(DE_iron,rownames(peak_whole_filter)),]
whole_logistic=data.frame(t(DE_peak),age,group=as.factor(group))

anes1<- glm(formula = as.formula(paste0("group~age+",paste0("`",DE_iron,"`", collapse = "+"))), 
            family = binomial(logit), data = whole_logistic)
summary(anes1)
anes1.AIC = stepAIC(anes1,direction="backward")
summary(anes1.AIC)
pre=predict(anes1,whole_logistic,type='response')
modelroc=roc(whole_logistic$group,pre)
plot(modelroc,print.auc=TRUE,print.thres=TRUE)

pre=predict(anes1.AIC,whole_logistic,type='response',main="training set")
modelroc=roc(whole_logistic$group,pre)
plot(modelroc,print.auc=TRUE,print.thres=TRUE,main="training set")

ttt=peak_Linkou[,DE_iron]
age_sample=Linkou$AGE[match(rownames(peak_Linkou),Linkou$Lab_number)]
group_sample=Linkou$SXT[match(rownames(peak_Linkou),Linkou$Lab_number)]
Linkou_predict=data.frame(ttt,age=as.numeric(age_sample),group=as.factor(group_sample))


pre=predict(anes1,Linkou_predict,type='response')
modelroc=roc(Linkou_predict$group,pre)
plot(modelroc,print.auc=TRUE,print.thres=TRUE)

pre=predict(anes1.AIC,Linkou_predict,type='response')
modelroc=roc(Linkou_predict$group,pre)
plot(modelroc,print.auc=TRUE,print.thres=TRUE,main="test set")
#----------------
DE_iron_AIC=colnames(anes1.AIC$model)[3:length(colnames(anes1.AIC$model))]
peak_RSA_AIC=peak_RSA[,DE_iron_AIC]
peak_SSA_sample_AIC=peak_SSA_sample[,DE_iron_AIC]
colSums(peak_RSA_AIC)
colSums(peak_SSA_sample_AIC)
iron=c(colnames(peak_RSA_AIC),colnames(peak_SSA_sample_AIC))
peak=c(colSums(peak_RSA_AIC),colSums(peak_SSA_sample_AIC))
group=c(rep("R",ncol(peak_RSA_AIC)),rep("S",ncol(peak_SSA_sample_AIC)))

peak_number=data.frame(iron=iron,peak=peak,group=group)
ggplot(peak_number, aes(x=iron, y=peak,fill=group)) +
  geom_bar(stat='identity',position="dodge",width = 0.5)+scale_fill_manual(values=c('red','blue'))+
  labs(title = "Kaohsiung Branch")+xlab("iron")+ylab("the number of peaks")+
  theme(axis.text.x  = element_text(angle=30, vjust=0.5))+
  theme(panel.grid.major =element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
