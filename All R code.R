##All code used for analyses for my thesis
#Mary Anne Schoenhardt, August 9th, 2023
#All files are copied from original folders into MAS_Files for R

setwd('C:/Users/Ryan Danby/Desktop/Mary Anne_Lab Computer/MAS_Files for R')
set.seed(28)

#####Packages needed#####
library(ggplot2)
library(randomForest)
library(dplyr)
library(corrplot)
library(pdp)
library(raster)
library(rgdal)
library(caTools)
library(tiff)
library(plyr)
library(epiDisplay)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(factoextra)
library(vegan)
library(rpart)
library(terra)
library(spatialEco)
library(EnvStats)
library(data.table)
library(RColorBrewer)
library(effectsize)
library(MASS)
library(mvnormtest)
library(gridExtra)
library(imager)
library(cowplot)
library(permimp)


######Random forest regression######
#See RamNurs.R file in 'Random forest regression' for original code and 
#Each file contains random points from composing samples from 10% of the pixels on each range
#Rows correspond to an individual point and columns to predictor variables and KD value associated with each point
donj_points_ramnurs <- read.csv(
  'Random forest regression/DonjekPoints_RamNurs.csv', header = TRUE)
TD_points_ramnurs <- read.csv(
  'Random forest regression/TDPoints_RamNurs.csv', header = TRUE)
vulc_points_ramnurs <- read.csv(
  'Random forest regression/VulcanPoints_RamNurs.csv', header = TRUE)
auri_points_ramnurs <- read.csv(
  'Random forest regression/AuriolPoints_RamNurs.csv', header = TRUE)

#Remove unneeded columns and reorder to alphabetical 
donj_points_ramnurs <- donj_points_ramnurs[,-(1:2)] 
donj_points_ramnurs <- donj_points_ramnurs[,order(colnames(donj_points_ramnurs))]
TD_points_ramnurs <- TD_points_ramnurs[,-(1:2)]
TD_points_ramnurs <- TD_points_ramnurs[,order(colnames(TD_points_ramnurs))]
vulc_points_ramnurs <- vulc_points_ramnurs[,-(1:2)]
vulc_points_ramnurs <- vulc_points_ramnurs[,order(colnames(vulc_points_ramnurs))]
auri_points_ramnurs <- auri_points_ramnurs[,-(1:2)]
auri_points_ramnurs <- auri_points_ramnurs[,order(colnames(auri_points_ramnurs))]

#Pairwise pearson correlation and plots
donj.cor.test <- cor(donj_points.test, method = "pearson", use = "complete.obs")
View(donj.cor.test)
corrplot(donj.cor.test)
TD.cor.ramnurs <- cor(TD_points_ramnurs, method = "pearson", use = "complete.obs")
View(TD.cor.ramnurs)
corrplot(TD.cor.ramnurs)
vulc.cor.ramnurs <- cor(vulc_points_ramnurs, method = "pearson", use = "complete.obs")
View(vulc.cor.ramnurs)
corrplot(vulc.cor.ramnurs)
auri.cor.ramnurs <- cor(auri_points_ramnurs, method = "pearson", use = "complete.obs")
View(auri.cor.ramnurs)
corrplot(auri.cor.ramnurs)

#Cleaning for random forest regression
#KD for all individuals is in sheep/m2. Convert to sheep/km2 and add to df
donj_points_ramnurs <- na.omit(donj_points_ramnurs)
TD_points_ramnurs <- na.omit(TD_points_ramnurs)
vulc_points_ramnurs <- na.omit(vulc_points_ramnurs)
auri_points_ramnurs <- na.omit(auri_points_ramnurs)

donj_kd_km2 <- donj_points_ramnurs$KD_Donj_250_sr910*1000000 
td_kd_km2 <- TD_points_ramnurs$KD_TD_250_sr670*1000000
vulc_kd_km2 <- vulc_points_ramnurs$KD_Vulc_250_sr1280*1000000
auri_kd_km2 <- auri_points_ramnurs$KD_Auri_250_sr1240*1000000

donj_points_ramnurs <- cbind(donj_points_ramnurs, donj_kd_km2)
TD_points_ramnurs <- cbind(TD_points_ramnurs, td_kd_km2)
vulc_points_ramnurs <- cbind(vulc_points_ramnurs, vulc_kd_km2)
auri_points_ramnurs <- cbind(auri_points_ramnurs, auri_kd_km2)

donj_points_ramnurs <- donj_points_ramnurs[,-8] #Remove sheep/m2
TD_points_ramnurs <- TD_points_ramnurs[,-6]
vulc_points_ramnurs <- vulc_points_ramnurs[,-6]
auri_points_ramnurs <- auri_points_ramnurs[,-8]

names(donj_points_ramnurs)[names(donj_points_ramnurs) == 'Altitude'] <- 'Elevation'
names(TD_points_ramnurs)[names(TD_points_ramnurs) == 'Altitude'] <- 'Elevation'
names(vulc_points_ramnurs)[names(vulc_points_ramnurs) == 'Altitude'] <- 'Elevation'
names(auri_points_ramnurs)[names(auri_points_ramnurs) == 'Altitude'] <- 'Elevation'

#Create random noise columns to be added to df
noise_donj<-sample(100, size = nrow(donj_points_ramnurs), replace = TRUE)
donj_points_ramnurs <- cbind(donj_points_ramnurs, noise_donj)
noise_td<-sample(100, size = nrow(TD_points_ramnurs), replace = TRUE)
TD_points_ramnurs <- cbind(TD_points_ramnurs, noise_td)
noise_vulc<-sample(100, size = nrow(vulc_points_ramnurs), replace = TRUE)
vulc_points_ramnurs <- cbind(vulc_points_ramnurs, noise_vulc)
noise_auri <- sample(100, size = nrow(auri_points_ramnurs), replace = TRUE)
auri_points_ramnurs <- cbind(auri_points_ramnurs, noise_auri)

#Create test and training sets
donj_points.test <- donj_points_ramnurs[1:2200,]
TD_points.test <- TD_points_ramnurs[1:1900,]
vulc_points.test <- vulc_points_ramnurs[1:4100,]
auri_points.test <- auri_points_ramnurs[1:3700,]

donj_points_ramnurs <- donj_points_ramnurs[-(1:2200),]
TD_points_ramnurs <- TD_points_ramnurs[-(1:1900),]
vulc_points_ramnurs <- vulc_points_ramnurs[-(1:4100),]
auri_points_ramnurs <- auri_points_ramnurs[-(1:3700),]

donj_points_ramnurs <- donj_points_ramnurs[,order(colnames(donj_points_ramnurs))]
TD_points_ramnurs <- TD_points_ramnurs[,order(colnames(TD_points_ramnurs))]
vulc_points_ramnurs <- vulc_points_ramnurs[,order(colnames(vulc_points_ramnurs))]
auri_points_ramnurs <- auri_points_ramnurs[,order(colnames(auri_points_ramnurs))]

donj_points.test <- donj_points.test[,order(colnames(donj_points.test))]
TD_points.test <- TD_points.test[,order(colnames(TD_points.test))]
vulc_points.test <- vulc_points.test[,order(colnames(vulc_points.test))]
auri_points.test <- auri_points.test[,order(colnames(auri_points.test))]

#Run RF models. Original models are saved in folder with .csv's
#Donj
donj.rf.keepbag <- randomForest(donj_kd_km2~Dist2_IceSnow+Eastness+Elevation+noise_donj+
                                  Northness+Slope+Solar+VRM, data = donj_points_ramnurs, 
                                xtest = donj_points.test[,-c(1,3:5,8,11:12)], ytest = donj_points.test[,3],
                                na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                keep.forest = T, keep.inbag = T)
save(donj.rf.keepbag, file = "Random forest regression/donj.rf.keepbag.RData")

donj.rf.nurs.keepbag <- randomForest(donj_kd_km2_nurs ~ Dist2_IceSnow+Eastness+Elevation+noise_donj+Northness+
                                       Slope+Solar+VRM, data = donj_points_ramnurs, 
                                     xtest = donj_points.test[,-c(1,3:5,8,11:12)], ytest = donj_points.test[,4],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(donj.rf.nurs.keepbag, file = "Random forest regression/donj.rf.nurs.keepbag.RData")

donj.rf.rams.keepbag <- randomForest(donj_kd_km2_rams ~ Dist2_IceSnow+Eastness+Elevation+noise_donj+Northness+
                                       Slope+Solar+VRM, data = donj_points_ramnurs, 
                                     xtest = donj_points.test[,-c(1,3:5,8,11:12)], ytest = donj_points.test[,5],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(donj.rf.rams.keepbag, file = "Random forest regression/donj.rf.rams.keepbag.RData")


#TD
TD.rf.keepbag <- randomForest(td_kd_km2 ~ Dist2_IceSnow+Eastness+Elevation+noise_td+
                                Northness+Slope+Solar+VRM, data = TD_points_ramnurs, 
                              xtest = TD_points.test[,-c(1,5,8:9,12:14)], ytest = TD_points.test[,12],
                              na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                              keep.forest = T, keep.inbag = T)
save(TD.rf.keepbag, file = "Random forest regression/TD.rf.keepbag.RData")

TD.rf.nurs.keepbag <- randomForest(td_kd_km2_nurs ~ Dist2_IceSnow+Eastness+Elevation+noise_td+Northness+
                                     Slope+Solar+VRM, data = TD_points_ramnurs, 
                                   xtest = TD_points.test[,-c(1,5,8:9,12:14)], ytest = TD_points.test[,13],
                                   na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                   keep.forest = T, keep.inbag = T)
save(TD.rf.nurs.keepbag, file = "Random forest regression/TD.rf.nurs.keepbag.RData")

TD.rf.rams.keepbag <- randomForest(td_kd_km2_rams ~ Dist2_IceSnow+Eastness+Elevation+noise_td+Northness+
                                     Slope+Solar+VRM,
                                   data = TD_points_ramnurs, xtest = TD_points.test[,-c(1,5,8:9,12:14)], ytest = TD_points.test[,14],
                                   na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                   keep.forest = T, keep.inbag = T)
save(TD.rf.rams.keepbag, file = "Random forest regression/TD.rf.rams.keepbag.RData")

#Vulc
vulc.rf.keepbag <- randomForest(vulc_kd_km2 ~ Dist2_IceSnow+Eastness+Elevation+noise_vulc+
                                  Northness+Slope+Solar+VRM, data = vulc_points_ramnurs, 
                                xtest = vulc_points.test[,-c(1,5,8:9,13:15)], ytest = vulc_points.test[,13],
                                na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                keep.forest = T, keep.inbag = T)
save(vulc.rf.keepbag, file = "Random forest regression/vulc.rf.keepbag.RData")

vulc.rf.nurs.keepbag <- randomForest(vulc_kd_km2_nurs ~ Dist2_IceSnow+Eastness+Elevation+noise_vulc+Northness+
                                       Slope+Solar+VRM, data = vulc_points_ramnurs, 
                                     xtest = vulc_points.test[,-c(1,5,8:9,13:15)], ytest = vulc_points.test[,14],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(vulc.rf.nurs.keepbag, file = "Random forest regression/vulc.rf.nurs.keepbag.RData")

vulc.rf.rams.keepbag <- randomForest(vulc_kd_km2_rams ~ Dist2_IceSnow+Eastness+Elevation+noise_vulc+Northness+
                                       Slope+Solar+VRM, data = vulc_points_ramnurs, 
                                     xtest = vulc_points.test[,-c(1,5,8:9,13:15)], ytest = vulc_points.test[,15],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(vulc.rf.rams.keepbag, file = "Random forest regression/vulc.rf.rams.keepbag.RData")

#Auri
auri.rf.keepbag <- randomForest(auri_kd_km2 ~ Dist2_IceSnow+Eastness+Elevation+noise_auri+
                                  Northness+Slope+Solar+VRM, data = auri_points_ramnurs, 
                                xtest = auri_points.test[,-c(1:4,8,11:12)], ytest = auri_points.test[,2],
                                na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                keep.forest = T, keep.inbag = T)
save(auri.rf.keepbag, file = "Random forest regression/auri.rf.keepbag.RData")

auri.rf.nurs.keepbag <- randomForest(auri_kd_km2_nurs ~ Dist2_IceSnow+Eastness+Elevation+noise_auri+Northness+
                                       Slope+Solar+VRM, data = auri_points_ramnurs, 
                                     xtest = auri_points.test[,-c(1:4,8,11:12)], ytest = auri_points.test[,3],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(auri.rf.nurs.keepbag, file = "Random forest regression/auri.rf.nurs.keepbag.RData")

auri.rf.rams.keepbag <- randomForest(auri_kd_km2_rams ~ Dist2_IceSnow+Eastness+Elevation+noise_auri+Northness+
                                       Slope+Solar+VRM, data = auri_points_ramnurs, 
                                     xtest = auri_points.test[,-c(1:4,8,11:12)], ytest = auri_points.test[,4],
                                     na.action = na.omit, mtry = 6, ntree = 1000, importance = T, 
                                     keep.forest = T, keep.inbag = T)
save(auri.rf.rams.keepbag, file = "Random forest regression/auri.rf.rams.keepbag.RData")

#Conditional permutation importance for each model
#Original values saved in 'Random forest regression/Conditional Permutation Importance'
permimp.donj <- permimp(donj.rf.keepbag)
capture.output(permimp.donj, file = "Permutation Importance Donjek.csv")
permimp.donj.nurs <- permimp(donj.rf.nurs.keepbag)
capture.output(permimp.donj.nurs, file = "Permutation Importance Donjek Nurs.csv")
permimp.donj.rams <- permimp(donj.rf.rams.keepbag)
capture.output(permimp.donj.rams, file = "Permutation Importance Donjek Rams.csv")

permimp.TD <- permimp(TD.rf.keepbag)
capture.output(permimp.TD, file = "Permutation Importance TD.csv")
permimp.TD.nurs <- permimp(TD.rf.nurs.keepbag)
capture.output(permimp.TD.nurs, file = "Permutation Importance TD Nurs.csv")
permimp.TD.rams <- permimp(TD.rf.rams.keepbag)
capture.output(permimp.TD.rams, file = "Permutation Importance TD Rams.csv")

permimp.vulc <- permimp(vulc.rf.keepbag)
capture.output(permimp.vulc, file = "Permutation Importance Vulcan.csv")
permimp.vulc.nurs <- permimp(vulc.rf.nurs.keepbag)
capture.output(permimp.vulc.nurs, file = "Permutation Importance Vulcan Nurs.csv")
permimp.vulc.rams <- permimp(vulc.rf.rams.keepbag)
capture.output(permimp.vulc.rams, file = "Permutation Importance Vulcan Rams.csv")

permimp.auri <- permimp(auri.rf.keepbag)
capture.output(permimp.auri, file = "Permutation Importance Auriol.csv")
permimp.auri.nurs <- permimp(auri.rf.nurs.keepbag)
capture.output(permimp.auri.nurs, file = "Permutation Importance Auriol Nurs.csv")
permimp.auri.rams <- permimp(auri.rf.rams.keepbag)
capture.output(permimp.auri.rams, file = "Permutation Importance Auriol Rams.csv")

#Partial dependence plots. Takes a VERY long time to run
par(mfrow=c(1,4))

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Elevation, main = 'Donjek', cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Elevation, main = 'TD', cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Elevation, main = 'Vulcan', cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Elevation, main = 'Auriol', cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.keepbag, donj_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.keepbag, TD_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.keepbag, vulc_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.keepbag, auri_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)

#For rams
partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Elevation, main = 'Donjek', cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Elevation, main = 'TD', cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Elevation, main = 'Vulcan', cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Elevation, main = 'Auriol', cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.rams.keepbag, donj_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.rams.keepbag, TD_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.rams.keepbag, vulc_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.rams.keepbag, auri_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)

#For nursery groups
partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Elevation, main = 'Donjek', cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Elevation, main = 'TD', cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Elevation, main = 'Vulcan', cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Elevation, main = 'Auriol', cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Dist2_IceSnow, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Eastness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Northness, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Slope, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, Solar, cex.axis = 2, cex.main = 3)

partialPlot(donj.rf.nurs.keepbag, donj_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(TD.rf.nurs.keepbag, TD_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(vulc.rf.nurs.keepbag, vulc_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)
partialPlot(auri.rf.nurs.keepbag, auri_points_ramnurs, VRM, cex.axis = 2, cex.main = 3)

#####Analysing NDVI images#####
#NDVI function with standardization
#Files from R have plot photo as bands 1-3 and standard photo bands as 4-6
vifun_norm <- function(img, stan) {
  b1 <- img[[1]]
  b2 <- img[[2]]
  stan1 <- cellStats(stan[[1]], stat='mean')
  stan2 <- cellStats(stan[[2]], stat='mean')
  vi_norm <- ((b1/stan1)-(b2/stan2))/((b1/stan1)+(b2/stan2))
  return(vi_norm)
}

#Build dataframe for average ndvi values
col <- c("Site", "Plot", "NDVI")
MeanNDVI <- data.frame(matrix(nrow = 0, ncol = length(col)))
colnames(MeanNDVI) <- col

#Masked photos (using ROI in ENVI) are all saved in folders corresponding to plot level naming
#Need to run the following lines of code for all individual plots and for each site
#Four folders contain NDVI photos - three from TD (Sheep Bullion, Congdon Creek, Sheep Mountain) and one from Donjek
plot.stack <- read.ENVI("Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/plot_mask")
plot.brick <- brick(plot.stack)
plot.bricksubset <- clamp(plot.brick, lower=0.0001, upper=Inf, useValues=FALSE)

stan.stack <- read.ENVI("Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/stan_mask")
stan.brick <- brick(stan.stack)
stan.bricksubset <- clamp(stan.brick, lower=0.0001, upper=Inf, useValues=FALSE)

ndvi <- vifun_norm(plot.bricksubset, stan.bricksubset) #calculate NDVI
#writeRaster(ndvi, "Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/ndvi.tiff", format="GTiff", overwrite = TRUE)

#Add mean ndvi to table
x <- c('TD_1sd47', '6', cellStats(ndvi, stat = 'mean'))
MeanNDVI <- rbind(x, MeanNDVI)
#Mean NDVI has already been calculated for all plots and is saved in 'Field data'

######Vegetation data#####
#Cover file is composed of all data collected except presence of forage species (Hoefs, 1976)
#Each row contains data from an individual plot
#NDVI (from above) has already been calculated and added to the cover file
cover <- read.csv(
  "Field data/Field Data_heights and cover.csv", 
  header = TRUE, na.strings = "")
colnames(cover)[colnames(cover) == "?..Range"] <- "Range"
#Species data contains presence/absence data for all species listed by Hoefs (1976)
#Each row contains data from an individual plot, with 1 indicating presence and 0 absence
species <- read.csv(
  "Field data/Field Data_species presence.csv", 
  header = TRUE, na.strings = "")

#Organize cover data
cover.highuse <- subset(cover, SD>=4)
cover.lowuse <- subset(cover, SD<=3) #All low and medium use sites are grouped together as low use

cover.site <- cover %>% #Group by site ID
  group_by(Site.ID) %>% 
  dplyr::summarize(across(Range:Slope, ~first(.)), across(Tree.Cover:LichMoss.Height, ~mean(.)), across(NDVI, ~mean(.)))
Use <- ifelse(cover.site$SD>=4, 'high','low')
cover.site <- cbind(cover.site, Use)
cover.site$Use <- as.factor(cover.site$Use)

#View statistical distribution (boxplots) for high vs low use sites
cover.site.short <- cover.site[,-(1:7)] #Remove columns of site level characteristics
cover.site.short <- cover.site.short[,-(8:9)] #Remove tree and tall shrub cover (very rarely observed)
cover.site.percent <- cover.site.short[,-(8:12)] #Remove species heights so only percent cover remains
colnames(cover.site.percent) <- c("Trees", "Tall Shrub", "Dwarf Shrub", "Graminoids", "Forbs",
                                  "Lichen and Moss", "Bare Ground", "Use")
cover.site.height <- cover.site.short[,-(1:7)]
cover.site.height <- cover.site.height[,-5] #Remove NDVI
colnames(cover.site.height) <- c("Dwarf Shrub", "Graminoids", "Forbs", "Lichen and Moss", "Use")
ndvi <- cover.site.short[,(12:13)]

cover.long <- cover.site.percent %>%
  pivot_longer(-Use, names_to = "variables", values_to = "Value") #Convert to long type
height.long <- cover.site.height %>%
  pivot_longer(-Use, names_to = "variables", values_to = "Value")
ndvi.long <- ndvi %>%
  pivot_longer(-Use, names_to = "variables", values_to = "Value")

stat.test.cover <- cover.long %>% #Run t-tests for all variables
  group_by(variables) %>%
  t_test(Value ~ Use) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.height <- height.long %>%
  group_by(variables) %>%
  t_test(Value ~ Use) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
stat.test.ndvi <- ndvi.long %>%
  group_by(variables) %>%
  t_test(Value ~ Use) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

coverplot <- ggboxplot( #Percent cover boxplots
  cover.long, x = "Use", y = "Value",
  fill = "Use", palette = c("#E60000", "#F9F200"), legend = "none", main = "Percent Cover",
  ggtheme = theme_pubr(border = TRUE)
) + facet_wrap(~variables)
stat.test.cover <- stat.test.cover %>% add_xy_position(x = "Use") #Add statistical differences
coverplot + 
  stat_pvalue_manual(stat.test.cover, label = "p.adj.signif") +
  theme(plot.title = element_text(hjust = 0.5))

heightplot <- ggboxplot( #Height boxplots
  height.long, x = "Use", y = "Value",
  fill = "Use", palette = c("#E60000", "#F9F200"), legend = "none", main = "Maximum Height (cm)",
  ggtheme = theme_pubr(border = TRUE)
) + facet_wrap(~variables)
stat.test.height <- stat.test.height %>% add_xy_position(x = "Use")
heightplot + 
  stat_pvalue_manual(stat.test.height, label = "p.adj.signif") + 
  theme(plot.title = element_text(hjust = 0.5))

ndviplot <- ggboxplot( #NDVI boxplot
  ndvi.long, x = "Use", y = "Value",
  fill = "Use", palette = c("#E60000", "#F9F200"), legend = "none", 
  main = "Normalized Difference Vegetation Index (NDVI)",
  ggtheme = theme_pubr(border = TRUE)
) + facet_wrap(~variables)
# Add statistical test p-values
stat.test.ndvi <- stat.test.ndvi %>% add_xy_position(x = "Use")
ndviplot + 
  stat_pvalue_manual(stat.test.ndvi, label = "p.adj.signif") + 
  theme(plot.title = element_text(hjust = 0.5))


#####PCA and OLS with veg data#####
#File contains cover data and species presence for all species observed at 5+ sites
TransectCenters <- read.csv(
  "Field data/Transect_centers_wKD.csv", header = T)
TransectCenters[is.na(TransectCenters)] <- 0
PCAdata <- TransectCenters[,-(1:6)] #Remove site level characteristics
PCAdata <- PCAdata[,-(29:30)] #Individual KD columns
PCAdata <- PCAdata[,-13] #Sheep use

#Transform percent variables to decimals
PCAdata$Dwarf_Cover <- PCAdata$Dwarf_Cover/100
PCAdata$Gram_Cover <- PCAdata$Gram_Cover/100
PCAdata$Forb_Cover <- PCAdata$Forb_Cover/100
PCAdata$LichMoss_Cover <- PCAdata$LichMoss_Cover/100
PCAdata$Bare_Cover <- PCAdata$Bare_Cover/100
PCAdata$Dwarf_Height <- PCAdata$Dwarf_Height/max(PCAdata$Dwarf_Height)
PCAdata$Gram_Height <- PCAdata$Gram_Height/max(PCAdata$Gram_Height)
PCAdata$Forb_Height <- PCAdata$Forb_Height/max(PCAdata$Forb_Height)

#Run PCA
VegPCA <- prcomp(PCAdata[,-28], center = T, scale. = T)
summary(VegPCA)
#Add first three PC's to veg CSV
PCAcol <- VegPCA$x[,(1:3)]
PCAdata <- cbind(PCAdata, PCAcol) 

#Linear model with PC's
kd.pc <- lm(All.KD~PC1+PC2+PC3, data = PCAdata)
plot(kd.pc)
summary(kd.pc)
anova(kd.pc)


#####Pixelwise Thiel-Sen analysis#####
#List all .tif files in the given directory. Files include KD distributions for each year of survey data, as well as blank rasters for all unsurveyed years
#Need to repeat for ram and nurs files in each range folder, as well as for each range
allfiles <- list.files("Pixelwise TheilSen/Donjek", 
                       pattern = ".tif$", full.names = T)
#Initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
Donj.spatrast <- rast(empty)

#Creates a raster stack of all file names stored in allfiles
for (i in 1:length(allfiles)){
  temp_data <- rast(allfiles[i]) #each file will be read in
  Donj.spatrast <- c(Donj.spatrast, temp_data) #for each iteration, bind the new data to the building dataset
}
#Check that nlyr=45 for output raster (46 for Vulcan and TD)
Donj.spatrast

#Calculate Thiel-Sen estimator for each pixel
Donj.TS <- raster.kendall(Donj.spatrast, p.value = TRUE)
writeRaster(Donj.TS[[1]], "Pixelwise TheilSen/Donjek/Donj_TS_slope.tif") #Slope
writeRaster(Donj.TS[[2]], "Pixelwise TheilSen/Donj_TS_Pval.tif") #P-value

#Calculate mean and SD for yearly values
Donj_mean <- app(Donj.spatrast, fun = mean)
Donj_SD <- app(Donj.spatrast, fun = sd)
writeRaster(Donj_mean, "Pixelwise TheilSen/Donjek/Donj_mean.tif")
writeRaster(Donj_SD, "Pixelwise TheilSen/Donjek/Donj_sd.tif")

#####Characteristics of regions of change#####
#Each file contains points from 10% of pixels in each of the positive, negative, and no change classes for each range and group
allfiles <- list.files(path = "Change characteristics", pattern = "*.csv", full.names = TRUE) 
df <- lapply(allfiles, read.csv)

#Group together all files (positive, negative, and no change) by range and group
auri.all <- rbind(df[[1]],df[[4]],df[[7]])
auri.nurs <- rbind(df[[2]],df[[5]],df[[8]])
auri.rams <- rbind(df[[3]],df[[6]],df[[9]])
donj.all <- rbind(df[[10]],df[[13]],df[[16]])
donj.nurs <- rbind(df[[11]],df[[14]],df[[17]])
donj.rams <- rbind(df[[12]],df[[15]],df[[18]])
TD.all <- rbind(df[[19]],df[[22]],df[[25]])
TD.nurs <- rbind(df[[20]],df[[23]],df[[26]])
TD.rams <- rbind(df[[21]],df[[24]],df[[27]])
vulc.all <- rbind(df[[28]],df[[31]],df[[34]])
vulc.nurs <- rbind(df[[29]],df[[32]],df[[35]])
vulc.rams <- rbind(df[[30]],df[[33]],df[[36]])

#Organize files
file.list <- list(auri.all, auri.nurs, auri.rams, donj.all, donj.nurs, donj.rams,
                  TD.all, TD.nurs, TD.rams, vulc.all, vulc.nurs, vulc.rams) #Create list with all files
file.list <- lapply(file.list, as.data.table) #Convert to data table to for use with melt function
file.list <- lapply(file.list, setcolorder, neworder = 
                      c("Change", "CID", "DEM", "Dist2_IceSnow", "Eastness", "KD", "Northness",
                        "OID_", "Range", "Slope", "Solar", "VRM")) #Set file order
#file.list <- lapply(file.list, setnames, c("Change", "CID", "Elevation (m.a.s.l)", "Distance to ice (m)", "Eastness", "KD", "Northness",
#                                           "OID_", "Range", "Slope (%)", "Solar (w/m2)", "Ruggedness (VRM)")) #Give all files the same names

#Create testing and training sets for all files in file.list
#Runs random forest classifier and outputs conditional permutation importance
#Remove marked lines if comparing to regions of no change (as opposed to just increase vs decrease)
rf_fun <- function(df) {
  sample.neg <- sample(1:(nrow(df[df$Change=="Negative"])), 
                       size = (nrow(df[df$Change=="Negative"])/2), replace = F)
  sample.mid <- sample(nrow(df[df$Change=="Negative"]):(nrow(df)-nrow(df[df$Change=="Positive"])), 
                       size = (nrow(df[df$Change=="None"])/2), replace = F) #Remove if comparing to no change regions
  sample.pos <- sample((nrow(df)-nrow(df[df$Change=="Positive"])):nrow(df), 
                       size = (nrow(df[df$Change=="Positive"]))/2, replace = F)
  range.name <- df[1,"Range"]
  test <- c(sample.neg, sample.mid, sample.pos)
  test.set <- df[test,-c(2,6,8:9)]
  test.set <- na.omit(test.set)
  test.set <- as.data.frame(test.set)
  test.set <- test.set[!(test.set$Change=="None"),] #Remove if comparing to no change
  test.set$Change <- as.factor(test.set$Change)
  train.set <- df[-test,-c(2,6,8:9)]
  train.set <- na.omit(train.set)
  train.set <- as.data.frame(train.set)
  train.set <- train.set[!(train.set$Change=="None"),] #Remove if comparing to no change
  train.set$Change <- as.factor(train.set$Change)
  changeclass.rf <- randomForest(Change~., data = train.set,
                                 xtest = test.set[,-1], ytest = test.set[,1],
                                 ntree = 500, mtry = 3,
                                 importance = T, keep.forest = T, keep.inbag = T)
  permimp.rf <- (permimp(changeclass.rf, threshold = 0.8))
  plot(permimp.rf)
  return(permimp.rf)
}

lapply(file.list, rf_fun)

#Create boxplots across all ranges and groups for topographic characteristics by change group
long_fun <- function(x) { #Removes unneeded columns
  x <- x[,-("OID_")]
  x <- x[,-("CID")]
  x <- x[,-("KD")]
  x <- x[,-("Range")]
}
list.long <- lapply(file.list, long_fun)
list.long <- lapply(list.long, melt, id = "Change") #Convert to long type grouping based on change category

graph.fun <- function(l){
  ggplot(l, aes(x = variable, y = value, fill = Change, color = Change)) +
    geom_boxplot(alpha = 0.75) +
    scale_color_manual(values = c("#3f8ebf", "#5d5f45", "#e36007")) +
    scale_fill_manual(values = c("#74add1", "#8a8d67", "#fa9857")) +
    facet_wrap(~variable, scale="free") +
    labs(y = "Topographic Variable Value", x = "") +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_blank(),
          panel.grid.major.y = element_line("white"),
          panel.grid.minor.y = element_line("white"),
          panel.grid.major.x = element_line("white"),
          panel.spacing = unit(2, "lines"),
          legend.key.size = unit(30, "pt"),
          legend.position = "bottom")
}
lapply(list.long, graph.fun)

