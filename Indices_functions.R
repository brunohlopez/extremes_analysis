library(ggplot2)
library(ggrepel)
library(dplyr)
library(ncdf4)

"""
The function yearly_mean_index is going to take in a text file that has the recordings of that specific index 
Example (PNO) and it is going to take the yearly average of those indices and then it is going to plot it against 
yearly sea surface temperature to see how well it fits, returns the correlation test of how well those were fitted
against each other as well as the p-value
"""

yearly_mean_index <- function(text_file, net_cdf, x_lab, y_lab, a_title) {
  indices <- read.table(text_file)
  names(indices) <- c("Year", "January", "February", "March", "April", "May", "June", "July",
                      "August", "September", "October", "November", "December")
  
  indices <- filter(indices, Year >= 1982 & Year < 2019)
  year <- indices$Year
  mean_pna_year <- apply(indices[,-1], 1, mean) #Gets the mean of the indices
  mean_pna_year <- as.double(mean_pna_year) #  converts that value to a double
  
  nc_file <- nc_open(net_cdf)
  sst <- ncvar_get(nc_file, "tmp")
  yearly_sst_mean <- matrix(data = NA, nrow = length(year), ncol = 1)
  
  for (i in 1:length(year)) {
    yearly_sst_mean[i] <- mean(sst[,,i], na.rm = T)
  }
  
  sst_ind_mean <- data.frame(year, mean_pna_year, yearly_sst_mean)
  
  print(ggplot(sst_ind_mean, aes(mean_pna_year, yearly_sst_mean)) +
          geom_point(col = "red", size = 4) +
          xlab(x_lab) +
          ylab(y_lab) +
          ggtitle(a_title) +
          geom_label_repel(aes(label = year)) +
          theme_bw())
  
  corr_data <- cor.test(mean_pna_year,yearly_sst_mean, method = "s") #which method shoud I use?
  
  return(corr_data) }

"""
This function gets takes everything like the text file above but this only gets the average of december, january,
and february of the indices instead of the entire year. and then fits those against the yearly average. returns the
correlation (Spearmans correlation)
"""
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

"""
 This is the same as the function above except it takes three consecutive months of Jan-March
 """

three_months_same_year <- function(text_file,net_cdf,x_lab, y_lab, a_title) {
  
  indices <- read.table(text_file)
  names(indices) <- c("Year", "January", "February", "March", "April", "May", "June", "July",
                      "August", "September", "October", "November", "December")
  indices <- filter(indices, Year > 1981 & Year < 2019)
  month_1 <- indices$January
  month_2 <- indices$February
  month_3 <- indices$March
  year <- indices$Year
  mean_pna_year <- apply(indices[,-1], 1, mean) #Gets the mean of the indices
  mean_pna_year <- as.double(mean_pna_year) #  converts that value to a double
  
  nc_file <- nc_open(net_cdf)
  sst <- ncvar_get(nc_file, "tmp")
  yearly_sst_mean <- matrix(data = NA, nrow = length(year), ncol = 1)
  
  for (i in 1:length(year)) {
    yearly_sst_mean[i] <- mean(sst[,,i], na.rm = T)
  }
  print(length(year))
  print(length(month_1))
  print(length(month_2))
  print(length(month_3))
  three_months <- data.frame(year,month_1,month_2,month_3)
  three_months_no_year <- data.frame(month_1,month_2,month_3)
  three_ind <- apply(three_months_no_year,1,mean)
  final_frame <- data.frame(year,three_ind,yearly_sst_mean)
  
  print(ggplot(final_frame, aes(three_ind,yearly_sst_mean)) +
          geom_point(col = "orange", size = 4) +
          xlab(x_lab) +
          ylab(y_lab) +
          ggtitle(a_title) +
          geom_label_repel(aes(label = year)) +
          theme_bw())
  print(length(three_ind))
  print(length(yearly_sst_mean))
  corr_three <- cor.test(three_ind, yearly_sst_mean, method = "s")
  return(corr_three)
}

sst_org <- function(net_cdf) {
  
  nc_file <- nc_open(net_cdf, write = T) #loads ncdf
  
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
  
  ind = which(loc1 < -10) # the regions with large parameters are locally consistent
  # length(ind)/(1440*720) # small portion of the overall data
  # set them to NAs for now
  loc1[ind]= NA
  sca1[ind]= NA
  sha1[ind]= NA
  
  ind1=which(sca1 > 100) 
  # length(ind)/(1440*720) # small portion of the overall data
  loc1[ind1]=NA
  sca1[ind1]=NA
  sha1[ind1]=NA
  
  ind2=which(sha1 > 100) 
  loc1[ind2]=NA
  sca1[ind2]=NA
  sha1[ind2]=NA
  
  write.csv(loc1, "location_org.csv")
  write.csv(sha1, "shape_org.csv")
  write.csv(sca1, "scale_org.csv")
}

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


