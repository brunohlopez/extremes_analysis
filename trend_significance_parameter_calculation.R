library(ncdf4)
library(ggplot2)
library(extRemes)
library(distillery)
library(Lmoments)
library(dplyr)
library(raster)
library(rgdal)


nc_file <- nc_open('nc_data_psd.nc', write = T) #loads netcdf

#extracts the variables
lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
time <- ncvar_get(nc_file ,"time")
sst <-  ncvar_get(nc_file, "tmp")

lon <- lon - 180 #makes the latitude on a -180-180 scale

#the dimensions of the arrays
nlon <- length(lon)
nlat <- length(lat)
nt <- length(time)

loc1 <- matrix(data=NA,nrow=nlon,ncol=nlat) #location parameter
sca1 <- matrix(data=NA,nrow=nlon,ncol=nlat)  #scale parameter
sha1 <- matrix(data=NA,nrow=nlon,ncol=nlat)  #shape parameter
tr_est <- matrix(data=NA,nrow=nlon,ncol=nlat)  #trend
lr_p_value <- matrix(data=NA,nrow=nlon,ncol=nlat) #p-value of the trend
sig_or_not <- matrix(data=NA,nrow=nlon,ncol=nlat) #trend signifiance

#Finds the location,scale and shape using maximum likelihood estimation

for (i in 1:nlon) {
  for(j in 1:nlat) {
    
    x = sst[i,j,]
    y = abs(max(x) - min(x))
    
    if (sum(is.na(x))  == 0 & ( y != 0) ) {
      out <- fevd(x,  units = "deg C")
      loc1[i,j] <- out$results$par[1]
      sca1[i,j] <- out$results$par[2]
      sha1[i,j] <- out$results$par[3]
    } 
    
    else {
      #this isn't neccesary I just like looking at it to not be confused
      loc1[i,j] <- NA
      sca1[i,j] <- NA
      sha1[i,j] <- NA
    }
  }
}


#========================================================
#Take out the weird values (probably didn't converge")
#Maybe a problem from fitting using the MLE estimation

ind = which(loc1 < -10) # the regions with large parameters are locally consistent
loc1[ind]= NA
sca1[ind]= NA
sha1[ind]= NA

ind1=which(sca1 > 100) 
loc1[ind1]=NA
sca1[ind1]=NA
sha1[ind1]=NA

ind2=which(sha1 > 100) 
loc1[ind2]=NA
sca1[ind2]=NA
sha1[ind2]=NA

#=============================================================
#For loop to get all the values with the trend fitted

for (i in 1:nlon) {
  for(j in 1:nlat) {
    
    x = sst[i,j,]
    
    if ( sum(is.na(x))  == 0 && is.na(loc1[i,j]) == F && is.na(sca1[i,j]) == F && is.na(sha1[i,j]) == F) {
      if (sd(x) >0){
        
        outlt = fevd(x,units = "deg C",location.fun = ~lintr) #fits the linear trend
        tr_est[i,j] = outlt$results$par[2]
        
      }
    } 
  }
}


lonlat <- as.matrix(expand.grid(lon, lat)) #makes the lat and lon into columns
lonlatdata <- data.frame(lonlat) #then convers them into a dataframe so we can add info
names(lonlatdata) <- c("lon", "lat") #changes the name

#mkaes them into a 1 column vector each so I can put into an excel file and then plot
locvec <- as.vector(loc1)
scavec <- as.vector(sca1)
shavec <- as.vector(sha1)
trest <- as.vector(tr_est)


#a new dataframe of the values up above
trend_data <- data.frame(cbind(lonlatdata, locvec, scavec, shavec, trest))

names(trend_data) <- c("Lon", "Lat", "Location", "Scale", "Shape", "trend")

#There were like 4 weird trend values, so I used this to take them out
tre_ind <- which(trend_data$trend < -10)
trend_data$Location[tre_ind] <- NA
trend_data$Scale[tre_ind] <- NA
trend_data$Shape[tre_ind] <- NA
trend_data$trend[tre_ind] <- NA


write.csv(trend_data, "Trend_Data.csv")

#Comparing the original model to the model with the trend in the location parameter

for (c in 1:nlon) {
  for (d in 1:nlat) {
    
    x_lr = sst[c,d,]
    
    if ( sum(is.na(x_lr))  == 0 && is.na(loc1[c,d]) == F && is.na(sca1[c,d]) == F && is.na(sha1[c,d]) == F) {
      if (sd(x_lr) >0){
        
        out_lin = fevd(x_lr, units = "deg C")
        outlt = fevd(x_lr,units = "deg C",location.fun = ~trend)
        lr_test_method = lr.test(out_lin,outlt)
        lr_p_value[c,d] = lr_test_method$p.value
        
      }
    }
  } 
}


ind_p = which(lr_p_value < 0.05) #If the P-value is less than 0.05 than the test is significant (non-stationary)
ind_not_p = which(lr_p_value >= 0.05) #if greater than it is not, doesn't change (stationary)
sig_or_not[ind_p] = 1 #points that are signficanr
sig_or_not[ind_not_p] = 0

#copies of the longitudes and latitudes to be made for 
#creating the netcdf files
lat2 <- lat 
lon2 <- lon
nlat <- length(lat)
nlon <- length(lon)

#these will be used to make the netcdf files
loc_array <- array(loc1, dim = c(nlon,nlat))
sca_array <- array(sca1, dim = c(nlon,nlat))
sha_array <- array(sha1, dim = c(nlon,nlat))
sig_array <- array(sig_or_not, dim = c(nlon,nlat))
p_array <- array(lr_p_value, dim = c(nlon,nlat))

#creation of the netcdf file with the different parameters
#these were location,scale,shape,trend_singificance and trend values
parameter_nc <- "parameter.nc"

londim <- ncdim_def("lon", "degrees east", as.double(lon))
latdim <- ncdim_def("lat", "degrees north", as.double(lat))

fill_value <- NA

loc_name <- "Location vector of Yearly max"
loc_def <- ncvar_def("Location", "unitless", list(londim, latdim), fill_value, loc_name, prec = "single")

sca_name <- "Scale vector of Yearly max"
sca_def <- ncvar_def("Scale", "unitless", list(londim, latdim), fill_value, sca_name, prec = "single")

sha_name <- "Shape vector of Yearly max"
sha_def <- ncvar_def("Shape", "unitless", list(londim, latdim), fill_value, sha_name, prec = "single")

p_name <- "P values of having a Trend Vs No Trend"
p_def <- ncvar_def("P-value", "Unitless", list(londim,latdim), fill_value, p_name, prec = "single")

sig_name <- "Significance of Trend Value"
sig_def <- ncvar_def("Significant", "1 = Significant, 0 + Not", list(londim,latdim), fill_value, sig_name, prec = "single")


nc_out <- nc_create(parameter_nc, list(loc_def, sca_def, sha_def, p_def, sig_def))

ncvar_put(nc_out, loc_def, loc_array)
ncvar_put(nc_out, sca_def, sca_array)
ncvar_put(nc_out, sha_def,sha_array)
ncvar_put(nc_out,p_def, p_array)
ncvar_put(nc_out,sig_def, sig_array)

title1 <- "Parameters From 1982-2018"
ncatt_put(nc_out, 0, "title", title1)
