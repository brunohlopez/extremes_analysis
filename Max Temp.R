library(ncdf4)

'''
Finds the max temperature of every grid cell for each year from 1982-2018
'''


#Running this function can crash computer if a lot of programs are open
#Or there is not enough memory on CPU (Learned hard way)
get_max_temp <- function(filename) {
  
  data <- nc_open(filename)
  
  lon <- ncvar_get(data, varid = "lon")
  lat <- ncvar_get(data, varid = "lat")
  time <- ncvar_get(data, varid = "time")
  
  sst <- ncvar_get(data, varid = "sst")
  
  #Gets the max for each grid cell for that year
  yearly_max <- apply(sst, 1:2, na.rm = T) 
  
  return(yearly_max)
  
}

#I run all of these seperatley so my computer does not crash
#If running on a 8-16gb memory pc allow 10-30 minutes for each function to run/

#For the year 1982
max82 <- get_max_temp('sst.day.mean.1982.nc')
write.csv(max82, file = "max82.csv")

#For the year 1983
max83 <- get_max_temp('sst.day.mean.1983.nc')
write.csv(max83, file = "max83.csv")

#For the year 1984
max84 <- get_max_temp('sst.day.mean.1984.nc')
write.csv(max84, file = "max84.csv")

#For the year 1985
max85 <- get_max_temp('sst.day.mean.1985.nc')
write.csv(max85, file = "max85.csv")

#For the year 1986
max86 <- get_max_temp('sst.day.mean.1986.nc')
write.csv(max86, file = "max86.csv")

#For the year 1987
max87 <- get_max_temp('sst.day.mean.1987.nc')
write.csv(max87, file = "max87.csv")

#For the year 1988
max88 <- get_max_temp('sst.day.mean.1988.nc')
write.csv(max88, file = "max88.csv")

#For the year 1989
max89 <- get_max_temp('sst.day.mean.1989.nc')
write.csv(max89, file = "max89.csv")

#For the year 1990
max90 <- get_max_temp('sst.day.mean.1990.nc')
write.csv(max90, file = "max90.csv")

#For the year 1991
max91 <- get_max_temp('sst.day.mean.1991.nc')
write.csv(max91, file = "max91.csv")

#For the year 1992
max92 <- get_max_temp('sst.day.mean.1992.nc')
write.csv(max92, file = "max92.csv")

#For the year 1993
max93 <- get_max_temp('sst.day.mean.1993.nc')
write.csv(max93, file = "max93.csv")

#For the year 1994
max94 <- get_max_temp('sst.day.mean.1994.nc')
write.csv(max94, file = "max94.csv")

#For the year 1995
max95 <- get_max_temp('sst.day.mean.1995.nc')
write.csv(max95, file = "max95.csv")

#For the year 1996
max96 <- get_max_temp('sst.day.mean.1996.nc')
write.csv(max96, file = "max96.csv")

#For the year 1997
max97 <- get_max_temp('sst.day.mean.1997.nc')
write.csv(max97, file = "max97.csv")

#For the year 1998
max98 <- get_max_temp('sst.day.mean.1998.nc')
write.csv(max98, file = "max98.csv")

#For the year 1999
max99 <- get_max_temp('sst.day.mean.1999.nc')
write.csv(max99, file = "max99.csv")

#For the year 2000
max00 <- get_max_temp('sst.day.mean.2000.nc')
write.csv(max00, file = "max00.csv")

#For the year 2001
max01 <- get_max_temp('sst.day.mean.2001.nc')
write.csv(max01, file = "max01.csv")

#For the year 2002
max02 <- get_max_temp('sst.day.mean.2002.nc')
write.csv(max02, file = "max02.csv")

#For the year 2003
max03 <- get_max_temp('sst.day.mean.2003.nc')
write.csv(max03, file = "max03.csv")

#For the year 2004
max04 <- get_max_temp('sst.day.mean.2004.nc')
write.csv(max04, file = "max04.csv")

#For the year 2005
max05 <- get_max_temp('sst.day.mean.2005.nc')
write.csv(max05, file = "max05.csv")

#For the year 2006
max06 <- get_max_temp('sst.day.mean.2006.nc')
write.csv(max06, file = "max06.csv")

#For the year 2007
max07 <- get_max_temp('sst.day.mean.2007.nc')
write.csv(max07, file = "max07.csv")

#For the year 2008
max08 <- get_max_temp('sst.day.mean.2008.nc')
write.csv(max08, file = "max08.csv")

#For the year 2009
max09 <- get_max_temp('sst.day.mean.2009.nc')
write.csv(max09, file = "max09.csv")

#For the year 2010
max10 <- get_max_temp('sst.day.mean.2010.nc')
write.csv(max10, file = "max10.csv")

#For the year 2011
max11 <- get_max_temp('sst.day.mean.2011.nc')
write.csv(max11, file = "max11.csv")

#For the year 2012
max12 <- get_max_temp('sst.day.mean.2012.nc')
write.csv(max12, file = "max12.csv")

#For the year 2013
max13 <- get_max_temp('sst.day.mean.2013.nc')
write.csv(max13, file = "max13.csv")

#For the year 2014
max14 <- get_max_temp('sst.day.mean.2014.nc')
write.csv(max14, file = "max14.csv")

#For the year 2015
max15 <- get_max_temp('sst.day.mean.2015.nc')
write.csv(max15, file = "max15.csv")

#For the year 2016
max16 <- get_max_temp('sst.day.mean.2016.nc')
write.csv(max16, file = "max16.csv")

#For the year 2017
max17 <- get_max_temp('sst.day.mean.2017.nc')
write.csv(max17, file = "max17.csv")

#For the year 2018
max18 <- get_max_temp('sst.day.mean.2018.nc')
write.csv(max18, file = "max18.csv")
