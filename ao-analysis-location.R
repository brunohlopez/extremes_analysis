library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal) #Loads all neccesary packages

#======================functions===================================================

dec_feb_index <- function(text_file, net_cdf, x_lab, y_lab, a_title) {
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
  
  print(ggplot(three_month_anom, aes(three_month_ind, yearly_sst_mean)) +
          geom_point(col = "green", size = 4) +
          xlab(x_lab) +
          ylab(y_lab) +
          ggtitle(a_title) +
          geom_label_repel(aes(label = the_year)) +
          theme_bw())
  
  corr_dec_feb <- cor.test(three_month_ind, yearly_sst_mean, method = "s")
  return(corr_dec_feb) }


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
#Creates all ot the required maticies
loc1 <- matrix(data=NA,nrow=nlon,ncol=nlat) 
sca1 <- matrix(data=NA,nrow=nlon,ncol=nlat) 
sha1 <- matrix(data=NA,nrow=nlon,ncol=nlat) 
tr_est <- matrix(data=NA,nrow=nlon,ncol=nlat)
lr_p_value <- matrix(data=NA,nrow=nlon,ncol=nlat) 
sig_or_not <- matrix(data=NA,nrow=nlon,ncol=nlat)


ao <- dec_feb_mean("AO_index.txt", "nc_data_psd.nc") #this is the ao values means


#=======================================COmparing AO to original Data===========================================
#loop 1
ao_effect <- matrix(data=NA,nrow=nlon,ncol=nlat)
lr_p_value_ao <- matrix(data=NA,nrow=nlon,ncol=nlat) 
sig_or_not_ao <- matrix(data=NA,nrow=nlon,ncol=nlat)


for (c in 1:nlon) {
  for (d in 1:nlat) {
    
    x_lr = sst[c,d,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc_ao[c,d]) == F && is.na(sca_ao[c,d]) == F && is.na(sha_ao[c,d]) == F) {
      if (sd(x_lr) >0){
        
        out_lin_ao = fevd(x_lr, units = "deg C")
        outlt_ao = fevd(x_lr,units = "deg C",location.fun = ~ao)
        ao_effect[c,d] = outlt_ao$results$par[2]
        lr_test_method = lr.test(out_lin_ao,outlt_ao)
        lr_p_value_ao[c,d] = lr_test_method$p.value
        
      }
      
    }
  } 
}

#This is where I learned RData files are a lot easier to work with then csv
save(ao_effect, file = "ao_effect.RData")
save(lr_p_value_ao, file = "lr_p_value_ao.RData")

ind_ao = which(lr_p_value_ao < 0.05)
ind_not_ao = which(lr_p_value_ao > 0.05)
sig_or_not_ao[ind_ao] = 1
sig_or_not_ao[ind_not_ao] = 0

#===================================== trend vs trend + a0==================================================
#loop 2

lr_p_value_ao_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)
sig_or_not_ao_trend1 <- matrix(data=NA,nrow=nlon,ncol=nlat)
the_ao_effect <- matrix(data = NA, nrow = nlon, ncol = nlat)
trend_effect <- matrix(data = NA, nrow = nlon, ncol = nlat)

load('org_parameters.RData')

for (e in 1:nlon) {
  for (f in 1:nlat) {
    
    x_lr = sst[e,f,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc1[e,f]) == F && is.na(sca1[e,f]) == F && is.na(sha1[e,f]) == F) {
      if (sd(x_lr) >0){
        
        out_lin_ao = fevd(x_lr, units = "deg C", location.fun = ~ao + trend)
        outlt_ao = fevd(x_lr,units = "deg C",location.fun = ~trend)
        lr_test_method = lr.test(out_lin_ao,outlt_ao)
        lr_p_value_ao_trend1[e,f] = lr_test_method$p.value
        the_ao_effect[e,f] = out_lin_ao$results$par[2]
        trend_effect[e,f] = out_lin_ao$results$par[3]
        
        
      }
      
    }
  } 
}

save(lr_p_value_ao_trend1, the_ao_effect, trend_effect, file = "ao_trend1.RData")

ind_ao_trend = which(lr_p_value_ao_trend1 < 0.05)
ind_ao_trend_not = which(lr_p_value_ao_trend1 > 0.05)
sig_or_not_ao_trend1[ind_ao_trend] = 1
sig_or_not_ao_trend1[ind_ao_trend_not] = 0