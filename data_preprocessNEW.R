library(data.table)
################################################################################
rad_minmax <- function(Rad){
  Rad_max <- as.double(apply(Rad,2, max))
  Rad_min <- as.double(apply(Rad,2, min))
  for (i in 1:length(Rad)) {
    d <- as.double(Rad[,i])
    d <- (d-Rad_min[i])/(Rad_max[i]-Rad_min[i])
    Rad[,i]<- d
  }
  return(Rad)
}

rad_zscoring1 <- function(Rad){
  Rad_mean <- as.double(apply(Rad,2, mean))
  Rad_sd <- as.double(apply(Rad,2, sd))
  for (i in 1:length(Rad)) {
    d <- Rad[,i]
    d <- (d-Rad_mean[i])/Rad_sd[i]
    Rad[,i]<- d
  }
  return(Rad)
}

rad_zscoring2 <- function(Rad, Rad_mean, Rad_sd){
  Rad_mean <- as.double(Rad_mean)
  Rad_sd <- as.double(Rad_sd)
  for (i in 1:length(Rad)) {
    d <- Rad[,i]
    d <- (d-Rad_mean[i])/Rad_sd[i]
    Rad[,i]<- d
  }
  return(Rad)
}

rad_extreme <- function(Rad){
  for (i in 1:length(Rad)) {
    d <- as.double(Rad[,i])
    d[which(d>boxplot.stats(d)$stats[5])]<- boxplot.stats(d)$stats[5]
    d[which(d<boxplot.stats(d)$stats[1])]<- boxplot.stats(d)$stats[1]
    Rad[,i]<- d
  }
  return(Rad)
}

rad_varianceThreshold <- function(Rad, T){
  Rad_mean <- as.double(apply(Rad,2, mean))
  for (i in 1:length(Rad)) {
    d <- Rad[,i]
    d[i]<- d[i]/Rad_mean[i]
    Rad[,i]<- d
  }
  
  Rad_sd <- as.double(apply(Rad,2, sd))
  Rad <- Rad[, which(Rad_sd>T)]
  return(Rad)
}
##############################set_path#########################
setwd("C:/Users/pc01/Desktop/PresentWork/GBM/MGMT/revision/data")
getwd()
data_dev <- as.data.frame(fread("data_dev.csv"))
data_test <- as.data.frame(fread("data_test.csv"))
#####################################################################
Rad1 <- data_dev[,6:14681]
Rad2 <- data_test[,6:14681]
Rad1 <- rad_extreme(Rad1)
# Rad2 <- rad_extreme(Rad2)
#Rad_mean <- apply(Rad1,2,mean)
Rad_sd <- apply(Rad1,2,sd)
index <- which(Rad_sd <= 1e-4)
Rad1 <- Rad1[, -index]
Rad2 <- Rad2[, -index]
Rad_mean <- apply(Rad1,2,mean)
Rad_sd <- apply(Rad1,2,sd)
Rad1 <- rad_zscoring1(Rad1)
Rad2 <- rad_zscoring2(Rad2, Rad_mean, Rad_sd)
###################################################################
data_dev <- data_dev[,c(1:5)]
data_dev <-cbind(data_dev, Rad1)
data_dev <- na.omit(data_dev)
data_test <- data_test[,c(1:5)]
data_test <-cbind(data_test, Rad2)
data_test <- na.omit(data_test)
#################################saving#################################
write.csv(data_dev, file = "data_Dev.csv")
write.csv(data_test, file = "data_Test.csv")