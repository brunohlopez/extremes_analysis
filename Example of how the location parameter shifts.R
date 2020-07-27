setwd("C:/Users/Luno Bropez/Downloads/Research_Desktop") #sets the working directory to my desktop

library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal) #Loads all neccesary packages
library(lattice)
library(RColorBrewer)
library(rasterVis)
library(evd)

nc_file <- nc_open('nc_data_psd.nc', write = T) #loads ncdf
time <- ncvar_get(nc_file ,"time")
sst <- ncvar_get(nc_file, "tmp")


long_point <- (80) * 4
lat_point <- 23 * 4
the_point <- sst[long_point, lat_point,]

max_point <- which.max(sst)
max_location <- arrayInd(max_point, c(1440,720)) #gets the longitude and latitude of the highest loctioon parameter

the_max <- sst[361,399,]


mean_extremes <- matrix(data = NA, nrow = 1, ncol = length(the_max))


for (i in 1:length(the_max)) {
  mean_extremes[i] <- mean(sst[,,i], na.rm = T)
}

max_data <- data.frame(t(rbind(time, the_max, mean_extremes)))
names(max_data) <- c('time', 'the_max', 'mean_extremes')

x <- seq(-4,6,by = .1)
y <- dnorm(x,mean = -0.06, sd = 1)
z <- dnorm(x, mean = 2.28, sd = 1)
plot(x,y,type = 'l', xlab = 'SST (degrees Celsius)', ylab = 'Density', col = 'blue')
plot(x,z,type = 'l', xlab = 'SST (degrees Celsius)', ylab = 'Density', col = 'red')
lines(x,y,type = 'l', xlab = 'SST (degrees Celsius)', col = 'blue', ylab = 'Density')
legend(-4, 0.3, legend = c('1982', '2018'), col = c('blue','red'),
       title = "Yearly Maximum", lty = 1)