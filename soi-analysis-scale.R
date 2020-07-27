library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal) #Loads all neccesary packages

#functions===============================================

#finds the mean of the index(SOI or AO) from 1982-2018

dec_feb_mean <- function(text_file, net_cdf) {
  indices <- read.table(text_file)
  names(indices) <- c("Year", "January", "February", "March", "April", "May", "June", "July",
                      "August", "September", "October", "November", "December")
  
  indices <- filter(indices, Year >= 1981 & Year < 2020)
  dec <- indices$December
  jan <- indices$January
  feb <- indices$February
  year <- indices$Year
  the_year <- 1982:2018
  mean_pna_year <- apply(indices[,-1], 1, mean) #Gets the mean of the indices
  mean_pna_year <- as.double(mean_pna_year) #  converts that value to a double
  
  nc_file <- nc_open(net_cdf)
  sst <- ncvar_get(nc_file, "tmp")
  yearly_sst_mean <- matrix(data = NA, nrow = length(the_year), ncol = 1)
  
  for (i in 1:length(the_year)) {
    yearly_sst_mean[i] <- mean(sst[,,i], na.rm = T)
  }
  
  dec_feb <- data.frame(cbind(year,jan,feb,dec))
  dec_feb[1,2:3] <- NA
  dec_feb[39,4] <-  NA
  
  dec_1981_2017 <- dec_feb$dec[1:37] #from the years 1981-2017
  jan_1982_2018 <- dec_feb$jan[2:38] #from the years 1982-2018
  feb_1982_2018 <- dec_feb$feb[2:38] #from the years 1982-2018
  
  three_month_anom <- data.frame(dec_1981_2017,jan_1982_2018,feb_1982_2018)
  three_month_ind <- apply(three_month_anom,1,mean) #getting the means of the three months dec-feb
  
  return(three_month_ind) }

nc_file <- nc_open('nc_data_psd.nc', write = T) #loads ncdf

#extracts the parameters from the NETCDF file
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

load('org_parameters.RData')

soi <- dec_feb_mean("SOI.txt", "nc_data_psd.nc") #this is the ao values means

soi_effect <- matrix(data=NA,nrow=nlon,ncol=nlat)
lr_p_value_soi <- matrix(data=NA,nrow=nlon,ncol=nlat) 
sig_or_not_soi <- matrix(data=NA,nrow=nlon,ncol=nlat)

#influence of the SOI on the scale parameter

for (c in 1:nlon) {
  for (d in 1:nlat) {
    
    x_lr = sst[c,d,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc1[c,d]) == F && is.na(sca1[c,d]) == F && is.na(sha1[c,d]) == F) {
      if (sd(x_lr) >0){
        
        out_lin_soi = fevd(x_lr, units = "deg C")
        outlt_soi = fevd(x_lr,units = "deg C",scale.fun = ~soi) #looks at the scale parameter
        soi_effect[c,d] = outlt_soi$results$par[3] #third results is the soi influence
        lr_test_method = lr.test(out_lin_soi,outlt_soi)
        lr_p_value_soi[c,d] = lr_test_method$p.value
        
      }
      
    }
  } 
}

save(soi_effect, file = "soi_effect_scale.RData ")
save(lr_p_value_soi, file = 'lr_p_value_soi_scale.RData')


#========================================trend vs trend + SOI================================

lr_p_value_soi_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)
sig_or_not_soi_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)
the_soi_effect <- matrix(data = NA, nrow = nlon, ncol = nlat)
trend_effect_soi <- matrix(data = NA, nrow = nlon, ncol = nlat)

load('org_parameters.RData')


for (e in 1:nlon) {
  for (f in 1:nlat) {
    
    x_lr = sst[e,f,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc1[e,f]) == F && is.na(sca1[e,f]) == F && is.na(sha1[e,f]) == F) {
      if (sd(x_lr) >0){
        
        out_lin_soi = fevd(x_lr, units = "deg C", scale.fun = ~soi + trend)
        outlt_soi = fevd(x_lr,units = "deg C", scale.fun = ~trend)
        lr_test_method = lr.test(out_lin_soi,outlt_soi)
        lr_p_value_soi_trend1[e,f] = lr_test_method$p.value
        the_soi_effect[e,f] = out_lin_soi$results$par[3] #the affect of the ao
        trend_effect_soi[e,f] = out_lin_soi$results$par[4] #the affect of the location
        
        
      }
      
    }
  } 
}

ind_soi_trend = which(lr_p_value_soi_trend1 < 0.05)
ind_soi_trend_not = which(lr_p_value_soi_trend1 > 0.05)
sig_or_not_soi_trend1[ind_soi_trend] = 1
sig_or_not_soi_trend1[ind_soi_trend_not] = 0


save(lr_p_value_soi_trend1,sig_or_not_soi_trend1 ,the_soi_effect, trend_effect_soi, file = "soi_trend_scale_parameter.RData")

