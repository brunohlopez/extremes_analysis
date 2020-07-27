library(ncdf4)
library(lattice)


#this functions gets the max temperature and returns it as a vector
#after changing the negative infinite values into NA
max_temp_point <- function(maxfilename, year) {
  maxyear <- read.csv(maxfilename)[2:721] #The first column is not needed
  
  tmp_array <- as.matrix(maxyear)
  tmp_array[tmp_array <= -Inf] <- NA
  
  tmp_vec <- as.vector(tmp_array)
  names(tmp_vec) <- c(year)
  
  return(tmp_vec)
}

#gets all of the max values and converts them into a array and compares them against all of the longitudes
maxyear82 <- max_temp_point("max82.csv", "1982")
maxyear83 <- max_temp_point("max83.csv", "1983")
maxyear84 <- max_temp_point("max84.csv", "1984")
maxyear85 <- max_temp_point("max85.csv", "1985")
maxyear86 <- max_temp_point("max86.csv", "1986")
maxyear87 <- max_temp_point("max87.csv", "1987")
maxyear88 <- max_temp_point("max88.csv", "1988")
maxyear89 <- max_temp_point("max89.csv", "1989")
maxyear90 <- max_temp_point("max90.csv", "1990")
maxyear91 <- max_temp_point("max91.csv", "1991")
maxyear92 <- max_temp_point("max92.csv", "1992")
maxyear93 <- max_temp_point("max93.csv", "1993")
maxyear94 <- max_temp_point("max94.csv", "1994")
maxyear95 <- max_temp_point("max95.csv", "1995")
maxyear96 <- max_temp_point("max96.csv", "1996")
maxyear97 <- max_temp_point("max97.csv", "1997")
maxyear98 <- max_temp_point("max98.csv", "1998")
maxyear99 <- max_temp_point("max99.csv", "1999")
maxyear00 <- max_temp_point("max00.csv", "2000")
maxyear01 <- max_temp_point("max01.csv", "2001")
maxyear02 <- max_temp_point("max02.csv", "2002")
maxyear03 <- max_temp_point("max03.csv", "2003")
maxyear04 <- max_temp_point("max04.csv", "2004")
maxyear05 <- max_temp_point("max05.csv", "2005")
maxyear06 <- max_temp_point("max06.csv", "2006")
maxyear07 <- max_temp_point("max07.csv", "2007")
maxyear08 <- max_temp_point("max08.csv", "2008")
maxyear09 <- max_temp_point("max09.csv", "2009")
maxyear10 <- max_temp_point("max10.csv", "2010")
maxyear11 <- max_temp_point("max11.csv", "2011")
maxyear12 <- max_temp_point("max12.csv", "2012")
maxyear13 <- max_temp_point("max13.csv", "2013")
maxyear14 <- max_temp_point("max14.csv", "2014")
maxyear15 <- max_temp_point("max15.csv", "2015")
maxyear16 <- max_temp_point("max16.csv", "2016")
maxyear17 <- max_temp_point("max17.csv", "2017")
maxyear18 <- max_temp_point("max18.csv", "2018")


#puts all of the max years according to their longitude and latitude, compared against the other ones
#This is what I will be analyzing
all_data <- data.frame(cbind(lonlatdata, maxyear82, maxyear83, maxyear84, maxyear85, maxyear86, maxyear87, maxyear88,
                             maxyear89, maxyear90, maxyear91, maxyear92, maxyear93, maxyear94, maxyear95, maxyear96,
                             maxyear97, maxyear98, maxyear99, maxyear00, maxyear01, maxyear02, maxyear03, maxyear04,
                             maxyear05, maxyear06, maxyear07, maxyear08, maxyear09, maxyear10, maxyear11, maxyear12,
                             maxyear13, maxyear14, maxyear15, maxyear16, maxyear17, maxyear18))


#This is the data to analyze without themissing values
all_data_no_na <- na.omit(all_data)

#gets the min and the max for all of the combos of longitudes and latitudes
#basicall the max of the max and the min of the max for al years
all_data$max <- apply(all_data[3:39], 1, max) #max temperature
all_data$min <- apply(all_data[3:39], 1, min) #minimum temperatures

head(na.omit(all_data)) #so we can view some of the data from the console

#Does the same thing as above
all_data_no_na <- na.omit(all_data)


#--------------------------------------------------------
#Makes the graphs of the Yearly Maximums

lat2 <- lat
lon2 <- lon
time2 <- time
nlat <- length(lat) #the length of lat's used to creat the netcdf file
nlon <- length(lon)
nt <- length(time)

grid <- expand.grid(lon = lon, lat = lat) #makes every combination of longitudes and latitudes
cutpts <- c(-5, 5, 10, 15, 20, 25, 30, 35, 40) #cutpoints of temepratures for making graphs


#turns it back into an array
tmp_mat_2 <- as.matrix(all_data[3:39])
dim(tmp_mat_2) #confrims the dimensions

#makes a new dataframe that is gonna be made into the netcdf file
tmp_ncdf <- array(tmp_mat_2, dim = c(nlon, nlat, nt)) #Max temperatures for each year
dim(tmp_ncdf)

max_ncdf <- array(all_data$max, dim = c(nlon, nlat)) #maximum temperatures of every year
dim(max_ncdf)

min_ncdf <- array(all_data$min, dim = c(nlon, nlat)) #minimum temperatures of every year
dim(min_ncdf)

#Plots all of the different max temperatures for each grid cell around the world
for (i in range(1:length(time))) {
  
  levelplot(tmp_ncdf[,,i] ~ lon * lat, data=grid, at=cutpts, cuts= 10, pretty=T, 
            col.regions=(rev(brewer.pal(10, "RdBu"))), main = paste("Max Temperature of", time[i], sep = " "))
  
  
}

#----------------------------------------------------------
#Makes the NETCDF File
nc_data_psd <- "nc_data_psd.nc" #the name of the netcdf file

#the dimensions for the longitude, latitude and time
londim <- ncdim_def("lon", "degrees east", as.double(lon))
latdim <- ncdim_def("lat", "degrees north", as.double(lat))
timedim <- ncdim_def("time", "years", as.double(time))

fill_value <- NA #The fill value for the netcdf file

#Inserting the data into the NetCDF file
max_temp <- "Max Temperatures of Each Year"
tmp_def <- ncvar_def("tmp", "deg C", list(londim, latdim, timedim), fill_value, max_temp, prec = "single")

max_all <- "Max Temperature of All Years"
max_def <- ncvar_def("Max Temp", "deg C", list(londim, latdim), fill_value,max_all, prec = "single")

min_all <- "The Minimum Maximum of All Years"
minmax_all <- ncvar_def("Min Max Temp", "deg C", list(londim, latdim), fill_value, min_all, prec = "single")

#creating the netcdf file
nc_out <- nc_create(nc_data_psd, list(tmp_def, max_def, minmax_all))

#putting the variables into the netcdf file
ncvar_put(nc_out, tmp_def, tmp_ncdf)
ncvar_put(nc_out, max_def, max_ncdf)
ncvar_put(nc_out, minmax_all, min_ncdf)

#including the attributes
ncatt_put(nc_out, "lon", "axis", "X")
ncatt_put(nc_out, "lat", "axis", "Y")
ncatt_put(nc_out, "time", "axis", "T")

title1 <- "Max Temps From 1982-2018"
ncatt_put(nc_out, 0, "title", title1)

nc_out #the final product

