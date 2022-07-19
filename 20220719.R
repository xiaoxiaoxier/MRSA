###Step 1 preprocessing Linkou data and Kaohsiung data separately###
setwd("drugresistance")

library(dplyr)
SA=read.csv("SA_Linkou.csv")
#Linkou dim 26852*23
SA$date_Day=as.factor(substr(SA$date_Day,1,4))

colcha=c('Intensity','mz')
SA[,colcha]=lapply(SA[,colcha],as.character)
colfac=c('Oxacillin','Clindamycin','Erythromycin','SXT')
SA[,colfac]=lapply(SA[,colfac],as.factor)

data_intensity=SA$Intensity
myFun=function(x){strsplit(x,split=';')}
Intens=lapply(data_intensity, myFun)
intenx_index=lapply(Intens, function(x){length(unlist(x))>0&length(unlist(x))<600})
data=SA[unlist(intenx_index),]


mz=lapply(data$mz, myFun)
mz_quzheng=lapply(mz, function(x){ceiling(as.numeric(unlist(x)))})
original_matrix=matrix(0, nrow(data),18000)
intensity=lapply(data$Intensity, myFun)
intensity_log=lapply(intensity, function(x){as.numeric(unlist(x))})
for (i in 1: nrow(data)) {
  original_matrix[i,(mz_quzheng[[i]]-2000)]=intensity_log[[i]]
}
#########################################

SA_label=data[,colfac]

#SA_label <- SA_label %>%
#mutate(Penicillin = ifelse(Penicillin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Erythromycin = ifelse(Erythromycin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Oxacillin = ifelse(Oxacillin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Clindamycin = ifelse(Clindamycin == "S",0,1))
SA_label <- SA_label %>%
  mutate(SXT = ifelse(SXT == "S",0,1))

matrix_linkou = original_matrix
label_linkou = SA_label

SA=read.csv("SA_Kaohsiung.csv")
#Kaohsiung dim 4955*27
SA$date_Day=as.factor(substr(SA$date_Day,1,4))

colcha=c('Intensity','mz')
SA[,colcha]=lapply(SA[,colcha],as.character)
colfac=c('Oxacillin','Clindamycin','Erythromycin','SXT')
SA[,colfac]=lapply(SA[,colfac],as.factor)

data_intensity=SA$Intensity
myFun=function(x){strsplit(x,split=';')}
Intens=lapply(data_intensity, myFun)
intenx_index=lapply(Intens, function(x){length(unlist(x))>0&length(unlist(x))<600})
data=SA[unlist(intenx_index),]


mz=lapply(data$mz, myFun)
mz_quzheng=lapply(mz, function(x){ceiling(as.numeric(unlist(x)))})
original_matrix=matrix(0, nrow(data),18000)
intensity=lapply(data$Intensity, myFun)
intensity_log=lapply(intensity, function(x){as.numeric(unlist(x))})
for (i in 1: nrow(data)) {
  original_matrix[i,(mz_quzheng[[i]]-2000)]=intensity_log[[i]]
}
#########################################

SA_label=data[,colfac]

#SA_label <- SA_label %>%
#mutate(Penicillin = ifelse(Penicillin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Erythromycin = ifelse(Erythromycin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Oxacillin = ifelse(Oxacillin == "S",0,1))
SA_label <- SA_label %>%
  mutate(Clindamycin = ifelse(Clindamycin == "S",0,1))
SA_label <- SA_label %>%
  mutate(SXT = ifelse(SXT == "S",0,1))
#SA_label <- SA_label %>%
#mutate(Fusidic_acid = ifelse(Fusidic_acid == "S",0,1))
matrix_kaohsiung = original_matrix
label_kaohsiung = SA_label
###step 2:binning with bin size 20##########

pseudo_ion_matrix <- function(original_matrix,bin,mz_range){
  uu=matrix(0, nrow(original_matrix),mz_range/bin)
  for (i in 1:(mz_range/bin)) {
    a=rowSums(original_matrix[,(bin*(i-1)+1):(bin*i)]);
    b=apply(original_matrix[,(bin*(i-1)+1):(bin*i)],1,function(x) length(which(x!=0)))
    c=apply(original_matrix[,(bin*(i-1)+1):(bin*i)],1,function(x) sqrt(sum(x^2)))
    uu[,i]=a/b+c
    
    print(i)
  }
  colnames(uu)= paste0("pseudo-ion", seq(1,(mz_range/bin)))
  uu[is.na(uu)]<-0
  return(uu)
}

ion_matrix_linkou=pseudo_ion_matrix(original_matrix = matrix_linkou ,bin = 20,mz_range = 18000)
ion_matrix_kaohsiung = pseudo_ion_matrix(original_matrix = matrix_kaohsiung ,bin = 20,mz_range = 18000)

#######Step 3 XGBoost classifier######
library(xgboost)
library(pROC)
set.seed(1)
params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.1, gamma=0.1, max_depth=10, min_child_weight=1, subsample=1, colsample_bytree=1,alpha=0.01,lambda=0.01)
bstSparse <- xgboost(params = params,data = ion_matrix_linkou, label = label_linkou$Oxacillin,nrounds = 120,eval_metric="auc",set.seed=123)
pred_OXA <- predict(bstSparse, ion_matrix_kaohsiung)
roc_OXA <- roc(response = label_kaohsiung$Oxacillin,predictor=pred_OXA,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_OXA
pred_linkou_OXA = predict(bstSparse, ion_matrix_linkou)
roc_linkou_OXA <- roc(response = label_linkou$Oxacillin,predictor=pred_linkou_OXA,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_linkou_OXA

params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.05, gamma=0, max_depth=8, min_child_weight=1, subsample=1, colsample_bytree=1)
bstSparse <- xgboost(params = params,data = ion_matrix_linkou, label = label_linkou$Clindamycin,nrounds = 100,eval_metric="auc",set.seed=123)
pred_CLI = predict(bstSparse,ion_matrix_kaohsiung)
roc_CLI <- roc(response = label_kaohsiung$Clindamycin,predictor=pred_CLI,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_CLI
pred_linkou_CLI = predict(bstSparse, ion_matrix_linkou)
roc_linkou_CLI <- roc(response = label_linkou$Clindamycin,predictor=pred_linkou_CLI,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_linkou_CLI

params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.1, gamma=0, max_depth=8, min_child_weight=1, subsample=1, colsample_bytree=1)
bstSparse <- xgboost(params = params,data = ion_matrix_linkou, label = label_linkou$Erythromycin,nrounds = 120,eval_metric="auc",set.seed=123)
pred_ERY = predict(bstSparse,ion_matrix_kaohsiung)
roc_ERY <- roc(response = label_kaohsiung$Erythromycin,predictor=pred_ERY,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_ERY
pred_linkou_ERY = predict(bstSparse, ion_matrix_linkou)
roc_linkou_ERY <- roc(response = label_linkou$Erythromycin,predictor=pred_linkou_ERY,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_linkou_ERY

params <- list(booster = "gbtree", objective = "binary:logistic", eta=0.05, gamma=0, max_depth=6, min_child_weight=1, subsample=1, colsample_bytree=1,alpha=0)
bstSparse <- xgboost(params = params,data = ion_matrix_linkou, label = label_linkou$SXT,nrounds = 50,eval_metric="auc",set.seed=123)
pred_SXT = predict(bstSparse,ion_matrix_kaohsiung)
roc_SXT <- roc(response = label_kaohsiung$SXT,predictor=pred_SXT,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_SXT
pred_linkou_SXT = predict(bstSparse, ion_matrix_linkou)
roc_linkou_SXT <- roc(response = label_linkou$SXT,predictor=pred_linkou_SXT,ci=TRUE, plot=TRUE,auc = TRUE,aur = TRUE)
roc_linkou_SXT





XGBoost_predict = data.frame(OXA = pred_OXA,CLI = pred_CLI, ERY=pred_ERY,SXT = pred_SXT,label_kaohsiung,multidrug=rowSums(label_kaohsiung))
XGBoost_predict_linkou = data.frame(OXA = pred_linkou_OXA,CLI = pred_linkou_CLI, ERY=pred_linkou_ERY,SXT = pred_linkou_SXT,label_linkou,multidrug=rowSums(label_linkou))


#####linkou#############
library(tidyr)
library(ggplot2)
library(ggcorrplot)

density_per_linkou = gather(XGBoost_predict_linkou[,1:4], key = "varname",value= "value",1:4)
ggplot(density_per_linkou)+theme_bw()+geom_density(aes(value),fill = "red",alpha=0.5)+facet_wrap(.~varname,scale = "free")

XGB_pred_cor_linkou = as.data.frame(cor(subset(XGBoost_predict_linkou,select=c("OXA","CLI","ERY","SXT","multidrug"))))

ggcorrplot(XGB_pred_cor_linkou,method = "square", lab = TRUE)+theme(axis.text.x = element_text(size = 10),axis.text.y = element_text(size = 10))

data_lm_linkou = subset(XGBoost_predict_linkou,select=c("OXA","CLI","ERY","SXT","multidrug"))
lm_linkou = lm(multidrug~.,data = data_lm_linkou)
summary(lm_linkou)


#par(mfrow = c(1,2))
parameter = lm_linkou$coefficients
score_linkou = parameter[1]+parameter[2]*XGBoost_predict_linkou$OXA+parameter[3]*XGBoost_predict_linkou$CLI+parameter[4]*XGBoost_predict_linkou$ERY+parameter[5]*XGBoost_predict_linkou$SXT
norm_score_linkou = (score_linkou-min(score_linkou))/max((score_linkou-min(score_linkou)))


###score plot after normalization
data_boxplot_linkou_norm = data.frame(multidrug=XGBoost_predict_linkou$multidrug,score= norm_score_linkou)
#ggplot(data_boxplot_linkou_norm,aes(x = multidrug, y = score,group = multidrug))+geom_boxplot()

# Draw the boxplot. Note result is also stored in a object called boundaries
boundaries <- boxplot(data_boxplot_linkou_norm$score ~ data_boxplot_linkou_norm$multidrug , col="#69b3a2" , ylim=c(0,1.1),outline = FALSE,main = "Linkou",
                      xlab = "the number of drugs that each sample is resistant to",ylab = "multi-drug resistance risk evaluation score")
# Add sample size on top
nbGroup <- nlevels(as.factor(data_boxplot_linkou_norm$multidrug))
text( 
  x=c(1:nbGroup), 
  y=boundaries$stats[nrow(boundaries$stats),] + 0.07, 
  paste("n = ",table(data_boxplot_linkou_norm$multidrug),sep="")  )

library(ggpubr)
multidurg_linkouplot=ggviolin(data_boxplot_linkou_norm ,x = "multidrug",y = "score", fill = "multidrug",palette = "npg",add = "boxplot",add.params =list(fill="white",width=0.05),width = 1.5 )

#####kaohsiung########

score_kaohsiung = parameter[1]+parameter[2]*XGBoost_predict_kaohsiung$OXA+parameter[3]*XGBoost_predict_kaohsiung$CLI+parameter[4]*XGBoost_predict_kaohsiung$ERY+parameter[5]*XGBoost_predict_kaohsiung$SXT
norm_score_kaohsiung = (score_kaohsiung-min(score_kaohsiung))/max((score_kaohsiung-min(score_kaohsiung)))

###score plot after normalization###
data_boxplot_kaohsiung_norm = data.frame(multidrug=XGBoost_predict_kaohsiung$multidrug,score= norm_score_kaohsiung)
#ggplot(data_boxplot_kaohsiung_norm,aes(x = multidrug, y = score,group = multidrug))+geom_boxplot()

# Draw the boxplot. Note result is also stored in a object called boundaries
boundaries <- boxplot(data_boxplot_kaohsiung_norm$score ~ data_boxplot_kaohsiung_norm$multidrug , col="#69b3a2" , ylim=c(0,1.1),outline = FALSE,main = "Kaohsiung",
                      xlab = "the number of drugs that each sample is resistant to",ylab = "multi-drug resistance risk evaluation score")
# Add sample size on top
nbGroup <- nlevels(as.factor(data_boxplo_kaohsiung_normt$multidrug))
text( 
  x=c(1:nbGroup), 
  y=boundaries$stats[nrow(boundaries$stats),] + 0.07, 
  paste("n = ",table(data_boxplot_kaohsiung_norm$multidrug),sep="")  )

multidurg_kaohsiungplot=ggviolin(data_boxplot_kaohsiung_norm ,x = "multidrug",y = "score", fill = "multidrug",palette = "npg",add = "boxplot",add.params =list(fill="white",width=0.05),width = 1.5,ylim =c(0,1) )

####linkou_norm_score_matrix###
OXA_score = data.frame(drug = XGBoost_predict_linkou$Oxacillin,score = norm_score_linkou)
CLI_score = data.frame(drug = XGBoost_predict_linkou$Clindamycin,score = norm_score_linkou)
ERY_score = data.frame(drug = XGBoost_predict_linkou$Erythromycin,score = norm_score_linkou)
SXT_score = data.frame(drug = XGBoost_predict_linkou$SXT.1,score = norm_score_linkou)

OXA_R_score = subset(OXA_score,drug=="1");OXA_R_score$drug = "OXA"
CLI_R_score = subset(CLI_score,drug=="1");CLI_R_score$drug = "CLI"
ERY_R_score = subset(ERY_score,drug=="1");ERY_R_score$drug = "ERY"
SXT_R_score = subset(SXT_score,drug=="1");SXT_R_score$drug = "SXT"

total_linkou = rbind.data.frame(OXA_R_score,CLI_R_score,ERY_R_score,SXT_R_score)
ggplot(total_linkou,aes(score,fill= drug,color = drug))+geom_density(alpha = 0.1)
ggboxplot(total_linkou ,x = "drug",y = "score", fill = "drug",palette = "npg",outlier.shape = NA )
singledrug_linkouplot=ggviolin(total_linkou ,x = "drug",y = "score", fill = "drug",palette = "lancet",add = "boxplot",add.params = list(fill="white",width=0.02),width = 1,ylim =c(0,1))


####kaohsiung_norm_score_matrix###
OXA_score = data.frame(drug = XGBoost_predict_kaohsiung$Oxacillin,score = norm_score_kaohsiung)
CLI_score = data.frame(drug = XGBoost_predict_kaohsiung$Clindamycin,score = norm_score_kaohsiung)
ERY_score = data.frame(drug = XGBoost_predict_kaohsiung$Erythromycin,score = norm_score_kaohsiung)
SXT_score = data.frame(drug = XGBoost_predict_kaohsiung$SXT.1,score = norm_score_kaohsiung)

OXA_R_score = subset(OXA_score,drug=="1");OXA_R_score$drug = "OXA"
CLI_R_score = subset(CLI_score,drug=="1");CLI_R_score$drug = "CLI"
ERY_R_score = subset(ERY_score,drug=="1");ERY_R_score$drug = "ERY"
SXT_R_score = subset(SXT_score,drug=="1");SXT_R_score$drug = "SXT"

total_kaohsiung = rbind.data.frame(OXA_R_score,CLI_R_score,ERY_R_score,SXT_R_score)
ggplot(total_kaohsiung,aes(score,fill= drug,color = drug))+geom_density(alpha = 0.1)
ggboxplot(total_kaohsiung ,x = "drug",y = "score", fill = "drug",palette = "npg",outlier.shape = NA )
singledrug_kaohsiungplot=ggviolin(total_kaohsiung ,x = "drug",y = "score", fill = "drug",palette = "lancet",add = "boxplot",add.params = list(fill="white",width=0.02),width = 1,ylim = c(0,1))


library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 2)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(multidurg_linkouplot, vp = vplayout(1, 1))
print(multidurg_kaohsiungplot, vp = vplayout(1, 2))
print(singledrug_linkouplot, vp = vplayout(2, 1))
print(singledrug_kaohsiungplot, vp = vplayout(2, 2))


