library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal)
library(RColorBrewer)
library(rasterVis)
library(gridExtra)
library(grid)


nc_file <- nc_open('nc_data_psd.nc', write = T) #loads netcdf

#extracts the variables
lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
time <- ncvar_get(nc_file ,"time")
sst <-  ncvar_get(nc_file, "tmp")
lonlat <- as.matrix(expand.grid(lon, lat)) #makes the lat and lon into columns
lonlatdata <- data.frame(lonlat) #then convers them into a dataframe so we can add info


lon <- lon - 180 #makes the latitude on a -180-180 scale

#the dimensions of the arrays
nlon <- length(lon)
nlat <- length(lat)
nt <- length(time)
lintr <- 1:nt #So that we are able to fit the trend


loc1 <- matrix(data=NA,nrow=nlon,ncol=nlat) #location parameter
sca1 <- matrix(data=NA,nrow=nlon,ncol=nlat)  #scale parameter
sha1 <- matrix(data=NA,nrow=nlon,ncol=nlat)  #shape parameter
tr_est <- matrix(data=NA,nrow=nlon,ncol=nlat)  #trend
tr_est_1 <- matrix(data=NA,nrow=nlon,ncol=nlat)
lr_p_value <- matrix(data=NA,nrow=nlon,ncol=nlat) #p-value of the trend
sig_or_not <- matrix(data=NA,nrow=nlon,ncol=nlat) #trend signifiance

load('org_parameters.RData')


for (i in 1:nlon) {
  for(j in 1:nlat) {
    
    x = sst[i,j,]
    
    if ( sum(is.na(x))  == 0 && is.na(loc1[i,j]) == F && is.na(sca1[i,j]) == F && is.na(sha1[i,j]) == F) {
      if (sd(x) >0){
        
        outlt = fevd(x,units = "deg C", scale.fun = ~lintr) #fits the linear trend on scale
        tr_est[i,j] = outlt$results$par[2] #scale parameter
        tr_est_1[i,j] = outlt$results$par[3] #mu 1 (influence on the scale from the Scale)
        
      }
    } 
  }
}

save(tr_est, tr_est_1, file = 'Scale_analysis.RData')
write.csv(tr_est_1, 'scale_analysis.csv')

lonlat <- as.matrix(expand.grid(lon, lat)) #makes the lat and lon into columns
lonlatdata <- data.frame(lonlat) #then convers them into a dataframe so we can add info
names(lonlatdata) <- c("lon", "lat") #changes the name

#mkaes them into a 1 column vector each so I can put into an excel file and then plot

tr_est_scale <- as.vector(tr_est_1) #the trend of the scale is made into a vector

trend_data <- data.frame(cbind(lonlatdata, tr_est_scale))
names(trend_data) <- c('Lon', 'Lat', 'Trend_scale')

trend_ind <- which(trend_data$Trend_scale > 3 ) #the values that are bigger than this
#are ones that also had an issue when doing the location parameter

trend_data$Trend_scale[trend_ind] <- NA #sets the big values to NA

#how to have the data in both of the formats
save(trend_data, file = 'scale_analysis_trend') #the trend of the scale analysis
write.csv(trend_data, 'scale_analysis_trend.csv')

#comparing the originakl model tro the model with the trend in the scale parameter
for (c in 1:nlon) {
  for (d in 1:nlat) {
    
    x_lr = sst[c,d,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc1[c,d]) == F && is.na(sca1[c,d]) == F && is.na(sha1[c,d]) == F) {
      if (sd(x_lr) >0){
        
        out_lin = fevd(x_lr, units = "deg C")
        outlt = fevd(x_lr,units = "deg C",scale.fun = ~lintr)
        lr_test_method = lr.test(out_lin,outlt) #using the likelihod ratio test compare the two models
        lr_p_value[c,d] = lr_test_method$p.value #use the p-value to determine signficance
        
      }
    }
  } 
}

save(lr_p_value, file = 'p-value-scale.RData')
write.csv(lr_p_value, 'p-value-scale.Rdata')

ind_p = which(lr_p_value < 0.05)
ind_not_p = which(lr_p_value >= 0.05)
sig_or_not[ind_p] = 1
sig_or_not[ind_not_p] = 0

save(sig_or_not, file = 'significant_scale_parameter.RData')
image(lon, lat,sig_or_not)

length(ind_p) / (length(ind_p) + length(ind_not_p)) #the amount of significant points
#The scale parameter has around 20% of points that are signifcant

#Plots for the signifcant scale parameter locations


load('significant_scale_parameter.RData')
load('Scale_analysis.RData')

cuts_trend <- seq(from = -0.25, to = 0.25, by = 0.05) #How many brteaks to show in the graphs
grid <- expand.grid(lon = lon, lat = lat) #makes every possible combination of lon and lat

#the trend of the scale parameter kind of all over the place but not as much as the location
scale_no_sig <- levelplot(tr_est_1 ~ lon * lat, data=grid, at = cuts_trend,  cuts= 10, pretty=T, xlab = "Longitude", ylab = "Latitude",
          col.regions=(rev(brewer.pal(10, "RdBu"))), main = paste("", sep = " "), ylab.right = list(label = "Degrees Celsius", vjust = -0.75, cex = 0.8))

sig_ind_keep <- which(sig_or_not == 0)
tr_est_1[sig_ind_keep] <- -0.25 #makes all the non-significant points black dots

#For making a plot of the significance and of the trend
zeroCol <-"#050505" #the color black
reds <- brewer.pal('Reds', n = 7) #to use for the red colors
blues <- rev(brewer.pal('Blues', n = 5)) #to use for blues
mytheme <- rasterTheme(region = c(zeroCol,blues,reds)) #makes a custom pallete
mytheme$panel.background$col = 'gray' #changes the NA color

#plot of the trend of the scale parameter and the trend 

scale_sig <- levelplot(tr_est_1 ~ lon * lat, data=grid, cuts= 10, at = cuts_trend ,par.settings = mytheme, xlab = "Longitude",
          ylab = "Latitude", useRaster = T, ylab.right = list(label = 'Degrees Celsius', vjust = -0.75, cex = 0.8))
  
grid.arrange(scale_no_sig, scale_sig)
