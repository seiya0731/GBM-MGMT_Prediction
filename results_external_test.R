library(parallel)
library(glmnet)
library(rms)
library(Hmisc)
library(lattice)
library(Formula)
library(ggplot2)
library(foreign)
library(psych)
library(pROC)
library(sampling)
library(ggpubr)
library(Matrix)
library(e1071)
library(mRMRe)
library(caret)
library(ResourceSelection)
library(rmda)
library(regplot)
library(car)
library(data.table)
library(randomForest)
###############################################################
setwd("C:/Users/pc01/Desktop/PresentWork/GBM/MGMT/revision/data")
getwd()
################################################################
logit_sig <- function(var_y,varlist,data){
  in_formula <- as.formula(paste(var_y,"~",varlist)) 
  p <- glm(in_formula,family=binomial(link=logit),data=data)
  coeff <- summary(p)$coefficients
  beta <- coeff[,1]
  LCI <- coeff[,1] - coeff[,2]*1.96 
  UCI <- coeff[,1] + coeff[,2]*1.96 
  OR <- exp(beta)
  OR_LCI <- exp(LCI)
  OR_UCI <- exp(UCI)
  p_value <- coeff[,4]
  name <- var
  data_var <- data.frame(OR,OR_LCI,OR_UCI,p_value)
  data_var <- data_var[-1,]
  return(data_var)
}

logit_sig_p <- function(varlist, label, data){
  in_formula <- as.formula(paste(label, "~", varlist)) 
  p <- glm(in_formula,family=binomial(link=logit), data=data)
  return(summary(p)$coefficients[2,4])
}

ttest_sig_p <- function(varlist, label, data){
  in_formula <- as.formula(paste(varlist, "~", label)) 
  t <- t.test(in_formula, data=data)
  return(t$p.value)
}

mrmr_calc <- function(data, y, feture_num){
  mrmr_feature<-data
  mrmr_feature$y <-y
  target_indices = which(names(mrmr_feature)=='y')
  for (m in which(sapply(mrmr_feature, class)!="numeric")){
    mrmr_feature[,m]=as.numeric(mrmr_feature[,m])
  }
  Data <- mRMR.data(data = data.frame(mrmr_feature))
  mrmr=mRMR.ensemble(data = Data, 
                     target_indices = target_indices, 
                     feature_count = feture_num, 
                     solution_count = 1)
  index_mrmr=mrmr@filters[[as.character(mrmr@target_indices)]]
  return(index_mrmr)
}

lasso_calc <- function(data_train, data_test, y_train, y_test){
  X_train <- as.matrix(data_train)
  Y_train <- y_train
  X_test  <- as.matrix(data_test)
  Y_test  <- y_test
  
  cv.fit <- cv.glmnet(X_train,Y_train,alpha=1,family='binomial', type.measure = "class")
  fit<-glmnet(X_train,Y_train,alpha=1,family='binomial')
  
  Coefficients <- coef(fit, s = cv.fit$lambda.min)
  Active.Index <- which(Coefficients != 0)
  Active.Weighted <- Coefficients[Active.Index]
  Active.Feature <-row.names(Coefficients)[Active.Index]
  index_lasso <- Active.Index[-1]-1
  return(index_lasso)
}

roc_metrics_calc <- function(score){
  vName <- colnames(score)
  roc_sum <- data.frame(matrix(nrow = 0, ncol = 5))
  for (i in 2:ncol(score)) {
    ROC <- roc(score[,1], as.numeric(score[, i]))
    auc95 <- paste0(round(ci(ROC)[2],3), 
                    ' (', round(ci(ROC)[1],3), 
                    '-', round(ci(ROC)[3],3), ')')
    uniroc<-data.frame('model' = vName[i],
                       'AUC(95%CI)' = auc95)
    youden <- coords(ROC, "best", 
                     ret=c("accuracy", "sensitivity", "specificity"))
    
    uniroc <- cbind(uniroc, round(youden, 3))
    roc_sum <- rbind(roc_sum, uniroc)
  }
  return(roc_sum)
}
#########################################################################
#########################################################################
#########################################################################
data <- as.data.frame(fread("data_dev.csv"))
data <- data[, -1]
data$MGMT <- factor(data$MGMT)
data$Gender <- factor(data$Gender)
data$Deep.White.Matter.Invasion <- factor(data$Deep.White.Matter.Invasion)
#########################################################################
data_external <- as.data.frame(fread("data_test.csv"))
data_external <- data_external[, -1]
data_external$MGMT <- factor(data_external$MGMT)
data_external$Gender <- factor(data_external$Gender)
data_external$Deep.White.Matter.Invasion <- factor(data_external$Deep.White.Matter.Invasion)
#########################################################################
#########################################################################
p = 0.7
N1 <- round(length(which(data$MGMT==0))*p)
N2 <- round(length(which(data$MGMT==1))*p) 
cmodels <- c("t1c_", "t1_", "t2_", "t2f_", "_1", "_2","_3", "t")
modals <- c("t1c_", "t1_", "t2_", "t2f_")
labels <- c("1", "2","3")
core <- makeCluster(4)
#########################################################################
#########################################################################
#########################################################################
g=7400
set.seed(g)
sub_str<-sampling::strata(data,stratanames=("MGMT"),size=c(N1,N2),method="srswor")
data_train<-data[sub_str$ID_unit, ]
data_test<-data[-sub_str$ID_unit, ]
score_train <- data.frame(matrix(nrow = 100, ncol = 0))
score_test  <- data.frame(matrix(nrow = 43, ncol = 0))
score_train$MGMT <- data_train$MGMT
score_test$MGMT  <- data_test$MGMT
auc_sum <- data.frame(matrix(nrow = 0, ncol = 70))
##########################step.1 Variance###############################
##########################step.2 LR#####################################
vName<-colnames(data_train[,-(1:5)])
print(paste0(g, ": LR==Started==============="))
pValue <- parSapply(core, vName, logit_sig_p, label="MGMT", data=data_train)
Index_logit<-which(pValue<0.05)
#############################################################################
for (s in c(501:10000)) {
  set.seed(s)
  sample_indices <- sample(1:nrow(data_external), 43, replace = FALSE)
  data_ext <- data_external[sample_indices, ]
  score_ext  <- data.frame(matrix(nrow = nrow(data_ext), ncol = 0))
  score_ext$MGMT  <- data_ext$MGMT
  auc  <- data.frame(matrix(nrow = 1, ncol = 0))
  auc["s"] <- s
  # cFeatures <- c()
  for (i in 1:length(modals)) {
    for (j in 1:length(labels)) {
      # i<- 3
      # j<- 3
      keyw <- paste0(modals[i], labels[j])
      cName <- paste0("score_", keyw)
      cNameTrain <- paste0(keyw, "_train")
      cNameTest <- paste0(keyw, "_test")
      cNameExt <- paste0(keyw, "_ext")
      #############################################################################
      Index_keyw <- grep(keyw, vName[Index_logit])
      Index_Rad <- Index_logit[Index_keyw]
      if (length(Index_keyw) < 5) {
        score_train[cName] <- NA
        score_test[cName]  <- NA
        auc[cNameTrain]    <- 0
        auc[cNameTest]     <- 0
        next
      }
      #######################step.3 MRMR#####################################
      if (length(Index_Rad) > 25) {
        Index_mrmr <- mrmr_calc(data_train[, Index_Rad+5], data_train$MGMT, 20)
        if(max(Index_mrmr) < length(Index_keyw))
          Index_Rad <- Index_Rad[Index_mrmr]
      }
      #######################step.4 LASSO####################################
      # set.seed(g)
      # Index_lasso <- lasso_calc(data_train[, Index_Rad + 21], 
      #                           data_test[, Index_Rad + 21], 
      #                           data_train$MGMT, 
      #                           data_test$MGMT)
      # Index_Rad <- Index_Rad[Index_lasso]
      X_train <- as.matrix(data_train[, Index_Rad + 5])
      Y_train <- data_train$MGMT
      X_test  <- as.matrix(data_test[, Index_Rad + 5])
      Y_test  <- data_test$MGMT
      X_ext  <- as.matrix(data_ext[, Index_Rad + 5])
      Y_ext  <- data_test$MGMT
      ######################################################################
      set.seed(g)
      cv.fit <- cv.glmnet(X_train, Y_train, alpha=1, family='binomial', type.measure = "class")
      # plot(cv.fit)
      # abline(v=log(c(cv.fit$lambda.min, cv.fit$lambda.1se)),
      #        col =c("red","black"),
      #        lty=c(2,2))
      ####################################################################
      fit<-glmnet(X_train, Y_train,alpha=1, family='binomial')
      # plot(fit, xvar = "lambda", label = TRUE)
      # abline(v=log(c(cv.fit$lambda.min, cv.fit$lambda.1se)),
      #        col =c("red","black"),
      #        lty =c(2,2))
      ####################################################################
      # Coefficients <- coef(fit, s = cv.fit$lambda.min)
      # Active.Index <- which(Coefficients != 0)
      # Active.Weighted <- Coefficients[Active.Index]
      # Active.Feature<-row.names(Coefficients)[Active.Index]
      # if (length(Active.Index)< 1) {
      #   score_train[cName] <- NA
      #   score_test[cName]  <- NA
      #   auc[cNameTrain]    <- 0
      #   auc[cNameTest]     <- 0
      #   next
      # }
      # cFeatures <- c(cFeatures, Active.Feature[-1])
      ##################################################################################
      #########################output#####################################
      
      score_train[cName] <-predict(fit,type="response", newx=X_train,
                                   s=cv.fit$lambda.min)
      score_test[cName]  <-predict(fit,type="response", newx=X_test,
                                   s=cv.fit$lambda.min)
      score_ext[cName]  <-predict(fit,type="response", newx=X_ext,
                                   s=cv.fit$lambda.min)
      auc[cNameTrain]  <- auc(roc(score_train$MGMT, as.numeric(score_train[[cName]])))
      auc[cNameTest]   <- auc(roc(score_test$MGMT, as.numeric(score_test[[cName]])))
      auc[cNameExt]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext[[cName]])))
      ####################################################################
      print(paste0(s, ": ", keyw, ": Finished=============================================="))
    }#end for label(j)
  }#end for model(i)
  ##################################################################################
  # Index_cFeat <- match(cFeatures, vName)
  ##################################################################################
  # for (k in 1:length(cmodels)) {
  #   keyw <- cmodels[k]
  #   cName <- paste0("score_", keyw)
  #   cNameTrain <- paste0(keyw, "_train")
  #   cNameTest <- paste0(keyw, "_test")
  #   cNameExt <- paste0(keyw, "_ext")
  #   
  #   Index_keyw <- grep(keyw, vName[Index_cFeat])
  #   Index_Rad <- Index_cFeat[Index_keyw]
  #   #######################step.3 MRMR#####################################
  #   # if (length(Index_Rad) > 25) {
  #   #   Index_mrmr <- mrmr_calc(data_train[, Index_Rad+21], data_train$MGMT, 20)
  #   #   Index_Rad <- Index_Rad[Index_mrmr]
  #   # }
  #   ####################################################################
  #   X_train <- as.matrix(data_train[, Index_Rad + 5])
  #   Y_train <- data_train$MGMT
  #   X_test  <- as.matrix(data_test[, Index_Rad + 5])
  #   Y_test  <- data_test$MGMT
  #   X_ext  <- as.matrix(data_ext[, Index_Rad + 5])
  #   Y_ext  <- data_test$MGMT
  #   ######################################################################
  #   set.seed(g)
  #   cv.fit <- cv.glmnet(X_train, Y_train, alpha=1, family='binomial', type.measure = "class")
  #   # plot(cv.fit)
  #   # abline(v=log(c(cv.fit$lambda.min, cv.fit$lambda.1se)),
  #   #        col =c("red","black"),
  #   #        lty=c(2,2))
  #   ####################################################################
  #   fit<-glmnet(X_train, Y_train,alpha=1, family='binomial')
  #   # plot(fit, xvar = "lambda", label = TRUE)
  #   # abline(v=log(c(cv.fit$lambda.min, cv.fit$lambda.1se)),
  #   #        col =c("red","black"),
  #   #        lty =c(2,2))
  #   ####################################################################
  #   Coefficients <- coef(fit, s = cv.fit$lambda.min)
  #   Active.Index <- which(Coefficients != 0)
  #   Active.Weighted <- Coefficients[Active.Index]
  #   Active.Feature<-row.names(Coefficients)[Active.Index]
  #   if (length(Active.Index)< 1) {
  #     score_train[cName] <- NA
  #     score_test[cName]  <- NA
  #     auc[cNameTrain]    <- 0
  #     auc[cNameTest]     <- 0
  #     next
  #   }
  #   #########################output#####################################
  #   score_train[cName] <-predict(fit,type="response", newx=X_train,
  #                                s=cv.fit$lambda.min)
  #   score_test[cName]  <-predict(fit,type="response", newx=X_test,
  #                                s=cv.fit$lambda.min)
  #   score_ext[cName]  <-predict(fit,type="response", newx=X_ext,
  #                                s=cv.fit$lambda.min)
  #   auc[cNameTrain]  <- auc(roc(score_train$MGMT, as.numeric(score_train[[cName]])))
  #   auc[cNameTest]   <- auc(roc(score_test$MGMT, as.numeric(score_test[[cName]])))
  #   auc[cNameExt]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext[[cName]])))
  #   ####################################################################
  #   print(paste0(g, ": ", keyw, ": Finished=============================================="))
  # }#end for k
  ##################################################################################
  sName<- colnames(score_train)
  ##################################################################################
  #t1c_LR
  Index_t1c <- grep("t1c_1|t1c_2|t1c_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_t1c]))==0)
  FML_t1c <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_t1c[Index_nna]],collapse = '+')))
  model_t1c <- glm(FML_t1c, family=binomial(link=logit),data=score_train)
  score_train["score_t1c_LR"] <-predict(model_t1c, type="response", newdata = score_train)
  score_test["score_t1c_LR"]  <-predict(model_t1c, type="response", newdata = score_test)
  score_ext["score_t1c_LR"]  <-predict(model_t1c, type="response", newdata = score_ext)
  auc["score_t1c_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_t1c_LR)))
  auc["score_t1c_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_t1c_LR)))
  auc["score_t1c_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_t1c_LR)))
  print(paste0(s, ": t1c_LR Finished=============================================="))
  ##################################################################################
  #t1_LR
  Index_t1 <- grep("t1_1|t1_2|t1_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_t1]))==0)
  FML_t1 <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_t1[Index_nna]],collapse = '+')))
  model_t1 <- glm(FML_t1, family=binomial(link=logit),data=score_train)
  score_train["score_t1_LR"] <-predict(model_t1, type="response", newdata = score_train)
  score_test["score_t1_LR"]  <-predict(model_t1, type="response", newdata = score_test)
  score_ext["score_t1_LR"]  <-predict(model_t1, type="response", newdata = score_ext)
  auc["score_t1_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_t1_LR)))
  auc["score_t1_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_t1_LR)))
  auc["score_t1_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_t1_LR)))
  print(paste0(s, ": t1_LR Finished=============================================="))
  ##################################################################################
  #t2_LR
  Index_t2 <- grep("t2_1|t2_2|t2_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_t2]))==0)
  FML_t2 <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_t2[Index_nna]],collapse = '+')))
  model_t2 <- glm(FML_t2, family=binomial(link=logit),data=score_train)
  score_train["score_t2_LR"] <-predict(model_t2, type="response", newdata = score_train)
  score_test["score_t2_LR"]  <-predict(model_t2, type="response", newdata = score_test)
  score_ext["score_t2_LR"]  <-predict(model_t2, type="response", newdata = score_ext)
  auc["score_t2_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_t2_LR)))
  auc["score_t2_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_t2_LR)))
  auc["score_t2_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_t2_LR)))
  print(paste0(s, ": t2_LR Finished=============================================="))
  ##################################################################################
  #t2f_LR
  Index_t2f <- grep("t2f_1|t2f_2|t2f_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_t2f]))==0)
  FML_t2f <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_t2f[Index_nna]],collapse = '+')))
  model_t2f <- glm(FML_t2f, family=binomial(link=logit),data=score_train)
  score_train["score_t2f_LR"] <-predict(model_t2f, type="response", newdata = score_train)
  score_test["score_t2f_LR"]  <-predict(model_t2f, type="response", newdata = score_test)
  score_ext["score_t2f_LR"]  <-predict(model_t2f, type="response", newdata = score_ext)
  auc["score_t2f_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_t2f_LR)))
  auc["score_t2f_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_t2f_LR)))
  auc["score_t2f_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_t2f_LR)))
  print(paste0(s, ": t2f_LR Finished=============================================="))
  ##################################################################################
  #ROI1_LR
  Index_ROI1 <- grep("t1c_1|t1_1|t2_1|t2f_1", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_ROI1]))==0)
  FML_ROI1 <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_ROI1[Index_nna]],collapse = '+')))
  model_ROI1 <- glm(FML_ROI1, family=binomial(link=logit),data=score_train)
  score_train["score_ROI1_LR"] <-predict(model_ROI1, type="response", newdata = score_train)
  score_test["score_ROI1_LR"]  <-predict(model_ROI1, type="response", newdata = score_test)
  score_ext["score_ROI1_LR"]  <-predict(model_ROI1, type="response", newdata = score_ext)
  auc["score_ROI1_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_ROI1_LR)))
  auc["score_ROI1_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_ROI1_LR)))
  auc["score_ROI1_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_ROI1_LR)))
  print(paste0(s, ": ROI1_LR Finished=============================================="))
  ##################################################################################
  #ROI2_LR
  Index_ROI2 <- grep("t1c_2|t1_2|t2_2|t2f_2", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_ROI2]))==0)
  FML_ROI2 <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_ROI2[Index_nna]],collapse = '+')))
  model_ROI2 <- glm(FML_ROI2, family=binomial(link=logit),data=score_train)
  score_train["score_ROI2_LR"] <-predict(model_ROI2, type="response", newdata = score_train)
  score_test["score_ROI2_LR"]  <-predict(model_ROI2, type="response", newdata = score_test)
  score_ext["score_ROI2_LR"]  <-predict(model_ROI2, type="response", newdata = score_ext)
  auc["score_ROI2_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_ROI2_LR)))
  auc["score_ROI2_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_ROI2_LR)))
  auc["score_ROI2_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_ROI2_LR)))
  print(paste0(s, ": ROI2_LR Finished=============================================="))
  ##################################################################################
  #ROI3_LR
  Index_ROI3 <- grep("t1c_3|t1_3|t2_3|t2f_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_ROI3]))==0)
  FML_ROI3 <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_ROI3[Index_nna]],collapse = '+')))
  model_ROI3 <- glm(FML_ROI3, family=binomial(link=logit),data=score_train)
  score_train["score_ROI3_LR"] <-predict(model_ROI3, type="response", newdata = score_train)
  score_test["score_ROI3_LR"]  <-predict(model_ROI3, type="response", newdata = score_test)
  score_ext["score_ROI3_LR"]  <-predict(model_ROI3, type="response", newdata = score_ext)
  auc["score_ROI3_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_ROI3_LR)))
  auc["score_ROI3_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_ROI3_LR)))
  auc["score_ROI3_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_ROI3_LR)))
  print(paste0(s, ": ROI3_LR Finished=============================================="))
  ##################################################################################
  #All_LR
  Index_All <- grep("t1c_1|t1_1|t2_1|t2f_1|t1c_2|t1_2|t2_2|t2f_2|t1c_3|t1_3|t2_3|t2f_3", sName)
  Index_nna <- which(colSums(is.na(score_train[,Index_All]))==0)
  FML_All <- as.formula(paste0("MGMT ~ ", paste0(sName[Index_All[Index_nna]],collapse = '+')))
  model_All <- glm(FML_All, family=binomial(link=logit),data=score_train)
  score_train["score_all_LR"] <-predict(model_All, type="response", newdata = score_train)
  score_test["score_all_LR"]  <-predict(model_All, type="response", newdata = score_test)
  score_ext["score_all_LR"]  <-predict(model_All, type="response", newdata = score_ext)
  auc["score_all_LR_train"]  <- auc(roc(score_train$MGMT, as.numeric(score_train$score_all_LR)))
  auc["score_all_LR_test"]   <- auc(roc(score_test$MGMT, as.numeric(score_test$score_all_LR)))
  auc["score_all_LR_ext"]   <- auc(roc(score_ext$MGMT, as.numeric(score_ext$score_all_LR)))
  print(paste0(s, ": All_LR Finished=============================================="))
  ##################################################################################
  model_All.step <- step(model_All,direction="both")
  score_train["score_all_step"] <-predict(model_All.step, type="response", 
                                          newdata = score_train)
  score_test["score_all_step"]  <-predict(model_All.step, type="response", 
                                          newdata = score_test)
  score_ext["score_all_step"]  <-predict(model_All.step, type="response", 
                                          newdata = score_ext)
  auc["score_all_step_train"]  <- auc(roc(score_train$MGMT, 
                                          as.numeric(score_train$score_all_step)))
  auc["score_all_step_test"]   <- auc(roc(score_test$MGMT, 
                                          as.numeric(score_test$score_all_step)))
  auc["score_all_step_ext"]   <- auc(roc(score_ext$MGMT, 
                                          as.numeric(score_ext$score_all_step)))
  print(paste0(s, ": All.step Finished=============================================="))
  ##################################################################################
  # Index_Clinical <- c(2,3,5)
  # clName<-colnames(data_train[,Index_Clinical])
  # print(paste0(g, ": LR_clinical==Started==============="))
  # pValue <- parSapply(core, clName, logit_sig_p, label="MGMT", data=data_train)
  # index<-which(pValue<0.05)
  
  model_clinical <- glm(MGMT~Gender+Deep.White.Matter.Invasion,
                   family=binomial(link=logit),data=data_train)
  score_train["Clinical"] <-predict(model_clinical, type="response", newdata = data_train)
  score_test["Clinical"]  <-predict(model_clinical, type="response", newdata = data_test)
  score_ext["Clinical"]  <-predict(model_clinical, type="response", newdata = data_ext)
  auc["Clinical_train"]  <- auc(roc(score_train$MGMT, 
                                          as.numeric(score_train$Clinical)))
  auc["Clinical_test"]   <- auc(roc(score_test$MGMT, 
                                          as.numeric(score_test$Clinical)))
  auc["Clinical_ext"]   <- auc(roc(score_ext$MGMT, 
                                    as.numeric(score_ext$Clinical)))
  print(paste0(s, ": Clinical Finished=============================================="))
  ##################################################################################
  data_train$score <-score_train$score_all_LR
  data_test$score <-score_test$score_all_LR
  data_ext$score <-score_ext$score_all_LR
  
  model_comb <- glm(MGMT~Gender+Deep.White.Matter.Invasion+score,
                 family=binomial(link=logit),data=data_train)
  score_train["Combined"] <-predict(model_comb, type="response", newdata = data_train)
  score_test["Combined"]  <-predict(model_comb, type="response", newdata = data_test)
  score_ext["Combined"]  <-predict(model_comb, type="response", newdata = data_ext)
  auc["Combined_train"]  <- auc(roc(score_train$MGMT, 
                                    as.numeric(score_train$Combined)))
  auc["Combined_test"]   <- auc(roc(score_test$MGMT, 
                                    as.numeric(score_test$Combined)))
  auc["Combined_ext"]   <- auc(roc(score_ext$MGMT, 
                                    as.numeric(score_ext$Combined)))
  print(paste0(s, ": Combined Finished=============================================="))
  ##################################################################################
  auc_sum <- rbind(auc_sum,auc)
}#end for s
write.csv(auc_sum, file = "auc_sum.csv")
  ##################################################################################
write.csv(score_train, file = "score_train.csv")
write.csv(score_test, file = "score_test.csv")
write.csv(score_ext, file = "score_ext.csv")
roc_metrics_train <- roc_metrics_calc(score_train)
roc_metrics_test <- roc_metrics_calc(score_test)
roc_metrics <- data.frame(matrix(nrow = 0, ncol = 5))
nRoc <- nrow(roc_metrics_train)
for (i in 1:nRoc) {
  roc_metrics <- rbind(roc_metrics, roc_metrics_train[i,])
  roc_metrics <- rbind(roc_metrics, roc_metrics_test[i,])
}
write.csv(roc_metrics_train, file = "roc_metrics_train.csv")
write.csv(roc_metrics_test, file = "roc_metrics_test.csv")
write.csv(roc_metrics, file = "roc_metrics.csv")
##################################################################################
clinical <- data[, 2:21]
clinical_train <- data_train[, 2:21]
clinical_test <- data_test[, 2:21]
index0 <- which(clinical_train$MGMT==0)
clinical_train_0 <- clinical_train[index0, ]
index1 <- which(clinical_train$MGMT==1)
clinical_train_1 <- clinical_train[index1, ]

index0 <- which(clinical_test$MGMT==0)
clinical_test_0 <- clinical_test[index0, ]
index1 <- which(clinical_test$MGMT==1)
clinical_test_1 <- clinical_test[index1, ]

length(which(clinical_train_0$Gender==0))

X<-c(14,9,14,6)
dim(X)<-c(2,2)

chisq.test(X, correct = TRUE)
fisher.test(X)




