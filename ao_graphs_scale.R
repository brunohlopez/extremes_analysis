setwd("C:/Users/Luno Bropez/Downloads/Research_Desktop") #sets the working directory to my desktop
#Loads all neccesary packages
library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal) 
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
sig_or_not_ao <- matrix(data=NA,nrow=nlon,ncol=nlat)


load('ap_effect_scale.RData')
load('lr_p_value_ao_scale.RData')
load('ao_trend_scale_parameter.Rdata')

#gets the significance of the ao on the trend
ind_ao <- which(lr_p_value_ao < 0.05)
ind_ao_not <- which(lr_p_value_ao > 0.05)
sig_or_not_ao[ind_ao] <- 1
sig_or_not_ao[ind_ao_not] <- 0
sig_vector <- as.vector(sig_or_not_ao)
ao_effect_only_sig <- as.vector(ao_effect)




#turn everything into vecotrs to turn into a dataframe
ao_effect_vector <- as.vector(ao_effect)
lr_p_value_ao_vector <- as.vector(lr_p_value_ao)
lr_p_value_ao_trend1_vector <- as.vector(lr_p_value_ao_trend1)
the_ao_effect_vector <- as.vector(the_ao_effect)
trend_effect_vector <- as.vector(trend_effect)

ao_data <- data.frame(grid, ao_effect_vector, lr_p_value_ao_vector, lr_p_value_ao_trend1_vector,
                      the_ao_effect_vector, trend_effect_vector, sig_vector, ao_effect_only_sig)

non <- which(sig_vector == 0) #how many points are significant
yes <- which(sig_vector == 1) #vice versa
significant_points <- length(yes) / (length(yes) + length(non))

#In total about 16.78 grid cells are significant  

#=================================plots the significance======================================
#I would expect the significance to be moret orwards the north
ggplot(data = ao_data, aes(x = lon, y = lat,fill = sig_vector )) +
  geom_raster(interpolate = T, show.legend = F) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("             Significance of AO as a Covariate in the Scale Parameter") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "black", na.value = "gainsboro", guide = "legend")

#cut to see the AO effect better
cuts = seq(from = -1, to=  1, by = 0.2)

levelplot(ao_data$ao_effect_vector ~ lon * lat, data = grid, cuts = 8, at = cuts, 
          xlab = "Longitude", ylab = "Latitude", main = paste(" ", sep = " "),
          useRaster = T, col.regions= rev(brewer.pal(10, "RdBu")), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.8))

zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 5))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

cuts_signficant_ao <- seq(from = -1, to=  1, by = .2)
ao_effect_only_sig[ind_ao_not] <- -1

plot_ao_sig <- levelplot(ao_effect_only_sig ~ lon * lat, data = grid, cuts = 8, at = cuts_signficant_ao , par.settings = mytheme,
                         xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = 'Degrees Celsius', vjust = -0.75, cex = 0.75),
                         useRaster = T)

#This is just the regular trend from before
ind_trend_ao <- which(trend_effect < -1.0)
trend_effect[ind_trend_ao] <- NA
cuts_trend_ao <- seq(from = -.25, to = .25, by = 0.05)

levelplot(trend_effect ~ lon * lat, data = grid, cuts =9, at = cuts_trend_ao,
          xlab = "Longitude", ylab = "Latitude", main = paste("Trend with AO & Trend as Covariates", sep = " "),
          useRaster = T, col.regions= rev(brewer.pal(10, "RdBu")))

sig_or_not_ao_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)

sig_points <- which(lr_p_value_ao_trend1 < 0.05)
not_sig_points <- which(lr_p_value_ao_trend1 > 0.05)
sig_or_not_ao_trend1[sig_points] <- 1
sig_or_not_ao_trend1[not_sig_points] <- 0

zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 5))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

trend_effect[not_sig_points] <- -0.25 #If you want to set back to original value reload data line 94
cuts_trend_ao <- seq(from = -.25, to = .25, by = 0.05)

#Shows onlt the significant trend points when ao and trend
trend_ao_plot <- levelplot(trend_effect ~ lon * lat, data = grid, cuts = 9, at = cuts_trend_ao , par.settings = mytheme,
                           xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.75),
                           useRaster = T)

grid.arrange(plot_ao_sig, trend_ao_plot, nrow = 2)



