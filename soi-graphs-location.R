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

nc_file <- nc_open('nc_data_psd.nc', write = T) #loads ncdf

lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
time <- ncvar_get(nc_file ,"time")
sst <-  ncvar_get(nc_file, "tmp")

lon <- lon - 180 #converts the longitude to a 0-180 degrees scale

#to make the arrays
nlon <- length(lon)
nlat <- length(lat)
nt <- length(time)
trend <- 1:nt #So that we are able to fit the trend
grid = expand.grid(lon = lon, lat = lat)

#===================load the data from the ao analysis script

load("soi_effect.RData")
load("lr_p_value_soi.RData")
sig_or_not_soi <- matrix(data=NA,nrow=nlon,ncol=nlat)
load("soi_trend.RData")
load('sig_soi_and_trend.RData')

#makes vectors of all of the soi vectors to make it easier to analyze
soi_effect_vector <- as.vector(soi_effect)
lr_p_value_soi_vector <- as.vector(lr_p_value_soi)
lr_p_value_soi_trend1_vector <- as.vector(lr_p_value_soi_trend1)
the_soi_effect_vector <- as.vector(the_soi_effect)
trend_effect_vector <- as.vector(trend_effect)

#gets which points are significant
ind_soi <- which(lr_p_value_soi < 0.05)
ind_soi_not <- which(lr_p_value_soi > 0.05)
sig_or_not_soi[ind_soi] <- 1
sig_or_not_soi[ind_soi_not] <- 0
sig_vector <- as.vector(sig_or_not_soi)

#the total amount of significantr points
yes <- which(sig_vector == 1)
non <- which(sig_vector == 0 )
total_significance <- length(yes) / (length(yes) + length(non))

soi_effect_only_sig <- as.vector(soi_effect)
#dataframe of all data combined
soi_data <- data.frame(grid, soi_effect_vector, lr_p_value_soi_vector, lr_p_value_soi_trend1_vector, the_soi_effect_vector, trend_effect_vector, sig_vector, soi_effect_only_sig )

#=================================plots the significance======================================
ggplot(data = soi_data, aes(x = lon, y = lat,fill = sig_vector )) +
  geom_raster(interpolate = T, show.legend = F) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("             Significance of SOI as a Covariate") +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "black", na.value = "gainsboro", guide = "legend")

#soi effect 

ind_soi_remove <- which(soi_effect_vector > 1.0)
soi_effect_vector[ind_soi_remove] <- NA
cuts_soi_trend <- seq(from = -1, to = 1, by = 0.2)

#plot of the significance of the SOI
levelplot(soi_effect_vector ~ lon * lat, data = grid, cuts = 9, at = cuts_soi_trend, 
          xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.8),
          useRaster = T, col.regions= rev(brewer.pal(10, "RdBu")))

#significant points in the SOI location
zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 6))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

cuts_signficant_soi <- seq(from = -1, to=  1, by = .2)

soi_effect_only_sig[ind_soi_not] <- -1

#SOI with the Significant Points
soi_sig <- levelplot(soi_effect_only_sig ~ lon * lat, data = grid, cuts = 8, at = cuts_signficant_soi , par.settings = mytheme,
                     xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.7),
                     useRaster = T)


#soi with trend and soi both as covariates

load(file = 'Test1SOI.RData')

ind_which_bad_trend <- which(trend_effect < -1.5)
ind_which_bad_soi <- which(the_soi_effect > 2)
trend_effect[ind_which_bad_trend] <- NA
the_soi_effect[ind_which_bad_soi] <- NA

cuts_trend_Soi <- seq(from = -0.25, to = 0.25, by = 0.05)
cuts_soi_trend <- seq(from = -1, to = 1, by = 0.2)


#Which of these points were deemed as stastically significant
soi_trend_sig <- which(lr_p_value_soi_trend1 < 0.05)
soi_trend_sig_not <- which(lr_p_value_soi_trend1 > 0.05)
sig_or_not_soi_trend1[soi_trend_sig] <- 1
sig_or_not_soi_trend1[soi_trend_sig_not] <- 0

sig_or_not_soi_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)

#custom color palette
zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 5))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

cuts_soi_trend <- seq(from = -.25, to = .25, by = 0.05)
trend_effect[soi_trend_sig_not] <- -0.25

#plot fo the trend + the ao
trend_soi <- levelplot(trend_effect ~ lon * lat, data = grid, cuts = 9, at = cuts_soi_trend , par.settings = mytheme,
                       xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.7),
                       useRaster = T)

grid.arrange(soi_sig, trend_soi, nrow = 2)
