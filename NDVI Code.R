library(raster)
library(rgdal)
library(caTools)
library(tiff)

#setwd("C:/YukonPhotos")

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


#Build dataframe for mean ndvi values from each plot to be added to
col <- c("Site", "Plot", "NDVI")
MeanNDVI <- data.frame(matrix(nrow = 0, ncol = length(col)))
colnames(MeanNDVI) <- col

--------------
#Import file and format
plot.stack <- read.ENVI("E:/Mary Anne_Lab Computer/Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/plot_mask")
plot.brick <- brick(plot.stack) #Converting the stach to a brick allows for faster processing
plot.bricksubset <- clamp(plot.brick, lower=0.0001, upper=Inf, useValues=FALSE) #Ensures any values of 0 are removed. These correspond to areas outside the ROI

stan.stack <- read.ENVI("E:/Mary Anne_Lab Computer/Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/stan_mask")
stan.brick <- brick(stan.stack)
stan.bricksubset <- clamp(stan.brick, lower=0.0001, upper=Inf, useValues=FALSE)


#Calculate NDVI
ndvi <- vifun_norm(plot.bricksubset, stan.bricksubset)
writeRaster(ndvi, "E:/Mary Anne_Lab Computer/Yukon!/Tachal Dhal Photos/Sheep Bullion/NDVI Photos/TD_1sd47_6/ndvi.tiff", format="GTiff", overwrite = TRUE)

#Add mean ndvi to table
x <- c('TD_1sd47', '6', cellStats(ndvi, stat = 'mean')) #Change site and plot labels to correspond with the photos analyzed
MeanNDVI <- rbind(x, MeanNDVI) #Adds the mean NDVI for a given plot to the ongoing dataframe

#Run once NDVI has been calculated for all photos
write.csv(MeanNDVI, file = "E:/Mary Anne_Lab Computer/Yukon!/Tachal Dhal Photos/MeanNDVI_Bullion.csv")








