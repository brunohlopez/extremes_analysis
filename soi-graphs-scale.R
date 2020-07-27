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
sig_or_not_soi <- matrix(data=NA,nrow=nlon,ncol=nlat)

load('soi_effect_scale.RData')
load('lr_p_value_soi_scale.RData')
load('soi_trend_scale_parameter.Rdata')

#gets the significance of the ao on the trend
ind_soi <- which(lr_p_value_soi < 0.05)
ind_soi_not <- which(lr_p_value_soi > 0.05)
sig_or_not_soi[ind_soi] <- 1
sig_or_not_soi[ind_soi_not] <- 0
sig_vector <- as.vector(sig_or_not_soi)
soi_effect_only_sig <- as.vector(soi_effect)

#turn everything into vecotrs to turn into a dataframe
soi_effect_vector <- as.vector(soi_effect)
lr_p_value_soi_vector <- as.vector(lr_p_value_soi)
lr_p_value_soi_trend1_vector <- as.vector(lr_p_value_soi_trend1)
the_soi_effect_vector <- as.vector(the_soi_effect)
trend_effect_vector_soi <- as.vector(trend_effect_soi)


soi_data <- data.frame(grid, soi_effect_vector, lr_p_value_soi_vector, lr_p_value_soi_trend1_vector,
                      the_soi_effect_vector, trend_effect_vector_soi, sig_vector, soi_effect_only_sig)

non <- which(sig_vector == 0)
yes <- which(sig_vector == 1)
significant_points <- length(yes) / (length(yes) + length(non)) #only about 15% percent of points are significant

ggplot(data = soi_data, aes(x = lon, y = lat,fill = sig_vector )) +
  geom_raster(interpolate = T, show.legend = F) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("             Significance of SOI as a Covariate in the Scale Parameter") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "black", na.value = "gainsboro", guide = "legend")

cuts = seq(from = -1, to=  1, by = 0.2)

levelplot(soi_data$soi_effect_vector ~ lon * lat, data = grid, cuts = 8, at = cuts, 
          xlab = "Longitude", ylab = "Latitude", main = paste(" ", sep = " "),
          useRaster = T, col.regions= rev(brewer.pal(10, "RdBu")), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.8))


zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 5))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

cuts_signficant_soi <- seq(from = -1, to=  1, by = .2)
soi_effect_only_sig[ind_soi_not] <- -1

plot_soi_sig <- levelplot(soi_effect_only_sig ~ lon * lat, data = grid, cuts = 8, at = cuts_signficant_soi , par.settings = mytheme,
                         xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = 'Degrees Celsius', vjust = -0.75, cex = 0.75),
                         useRaster = T)


#This is just the regular trend from before
ind_trend_soi <- which(trend_effect_soi < -1.0)
trend_effect[ind_trend_ao] <- NA
cuts_trend_soi <- seq(from = -.25, to = .25, by = 0.05)

levelplot(trend_effect_soi ~ lon * lat, data = grid, cuts =9, at = cuts_trend_soi,
          xlab = "Longitude", ylab = "Latitude", main = paste("Trend with soi & Trend as Covariates", sep = " "),
          useRaster = T, col.regions= rev(brewer.pal(10, "RdBu")))


#Which points are significant for the trend and the SOI effect
sig_or_not_soi_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)

sig_points <- which(lr_p_value_soi_trend1 < 0.05)
not_sig_points <- which(lr_p_value_soi_trend1 > 0.05)
sig_or_not_soi_trend1[sig_points] <- 1
sig_or_not_soi_trend1[not_sig_points] <- 0

zeroCol <- "#050505" #the color black
red <- brewer.pal('Reds', n = 7)
blue <- rev(brewer.pal("Blues", n = 5))
mytheme <- rasterTheme(region = c(zeroCol, blue, red))
mytheme$panel.background$col <- 'gray'

trend_effect_soi[not_sig_points] <- -0.25 #If you want to set back to original value reload data line 94
cuts_trend_soi <- seq(from = -.25, to = .25, by = 0.05)

#Shows onlt the significant trend points when ao and trend
trend_soi_plot <- levelplot(trend_effect_soi ~ lon * lat, data = grid, cuts = 9, at = cuts_trend_soi , par.settings = mytheme,
                           xlab = "Longitude", ylab = "Latitude", main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.75),
                           useRaster = T)

grid.arrange(plot_soi_sig, trend_soi_plot, nrow = 2)



