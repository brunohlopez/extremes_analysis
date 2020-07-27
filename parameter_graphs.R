library(ncdf4)
library(raster)
library(rgdal) #Loads all neccesary packages
library(lattice)
library(RColorBrewer)
library(ggplot2)
library(rasterVis)
library(gridExtra)
library(grid)

nc_file1 <- nc_open('parameter.nc', write = T) #loads ncdf

lon <- ncvar_get(nc_file1, "lon")
lat <- ncvar_get(nc_file1, "lat")
significane <- ncvar_get(nc_file1, "Significant" )
p_val <- ncvar_get(nc_file1, "P-value")

lonlat <- as.matrix(expand.grid(lon, lat)) #makes the lat and lon into columns
lonlatdata <- data.frame(lonlat) #then convers them into a dataframe so we can add info
names(lonlatdata) <- c("lon", "lat") #changes the name
significance_vec <- as.vector(significane) #vector of o or 1
p_val_vec <- as.vector(p_val) #vector of the p-value not as important

Sig_Data <- data.frame(cbind(lonlatdata, significance_vec, p_val_vec))

names(Sig_Data) <- c("Lon", "Lat", "Significant", "P_Value")


#plot of the significance of the location parameter trend
ggplot(data = Sig_Data, aes(x = Lon, y = Lat,fill = significance_vec )) +
  geom_raster(interpolate = T, show.legend = F) +
  xlab("Longitude") +
  ylab("Latitude") +
  ggtitle("             Trend Significance of Global Ocean Extremes from 1982-2018") +
  theme_update(plot.title = element_text(hjust = 0.5)) +
  scale_fill_gradient(low = "white", high = "black", na.value = "gainsboro", guide = "legend")  +
  theme(panel.grid.major =  element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) 

#loads the trend data csv
trend_data <- read.csv("Trend_Data.csv")

#gets the location,scale,shape and trend from the previously created data frame
location = trend_data$Location
shape = trend_data$Shape
scale = trend_data$Scale
trend = trend_data$trend

nc_file <- nc_open('nc_data_psd.nc') #loads ncdf

lon <- ncvar_get(nc_file, "lon")
lat <- ncvar_get(nc_file, "lat")
lon <- lon - 180

trend_and_sig <- cbind(trend_data, significance_vec) #Binds the significance with the trend data

grid <- expand.grid(lon = lon, lat = lat) #makes every possible combination of lon and lat

#makes cuts for the graphs so values above a certain value will be treated as the max or min cutoff point

#These specific values are used so that the graph is not skewed
cuts_shape <- seq(from = -1, to = 1, by = 0.25)
cuts_scale <- seq(from = -0.1, to = 2, by =  0.3)
cuts_trend <- seq(from = -0.25, to = 0.25, by = 0.05)

#Makes the final plots, of location, scale, shape, trend
mytheme <- rasterTheme(panel.background = "Gainsboro")

plot_location <- levelplot(location ~ lon * lat, data=grid, cuts= 10, pretty=T, xlab = 'Longitude', ylab = 'Latitude', interpolate = T,
                           col.regions=(rev(brewer.pal(10, "RdBu"))), main = paste("", sep = " "))

plot_shape <- levelplot(shape ~ lon * lat, data=grid, at = cuts_shape, cuts= 10, pretty=T, xlab = 'Longitude', ylab = 'Latitude',
                        col.regions=(rev(brewer.pal(10, "YlGnBu"))), main = paste("", sep = " "))

plot_scale <- levelplot(scale ~ lon * lat, data=grid, at = cuts_scale, cuts= 5, pretty=T, xlab = 'Longitude', ylab = 'Latitude',
                        col.regions=(rev(brewer.pal(10, "Spectral"))), main = paste("", sep = " "))

#arranges them into one plot
grid.arrange(plot_location, plot_shape,plot_scale, nrow = 3)

#cuts for making graphs of the different subsections of
#central america and the tasman sea
the_at <- seq(from = 25, to = 35, by = 1.11)

#plot off the coast of Asia
plot_location_asia <- levelplot(location ~ lon * lat, data=grid, cuts= 9, pretty=T, xlab = 'Longitude', ylab = 'Latitude', interpolate = T,
                                main = paste("", sep = " "), at = the_at,  col.regions=(brewer.pal(9, "YlOrRd")),
                                xlim = c(-180,40), ylim= c(-30,40), ylab.right = list(label ="Degrees Celsius", vjust = -0.75, cex = 0.8))

#plot of Central America
plot_location_central <- levelplot(location ~ lon * lat, data=grid, cuts= 9, pretty=T, xlab = 'Longitude', ylab = 'Latitude', interpolate = T,
                                   main = paste("", sep = " "), at = the_at,  col.regions=(brewer.pal(9, "YlOrRd")),
                                   xlim = c(65,110), ylim= c(0,40), ylab.right = list(label = "Degrees Celsius", vjust =- 0.75, cex = 0.8)) 

#plots together
grid.arrange(plot_location_asia, plot_location_central)


#The trend without taking in significanceinto account
Tasman_sea <- levelplot(trend ~ lon * lat, data=grid, at = cuts_trend,  cuts= 10, pretty=T, xlab = "Longitude", ylab = "Latitude",
                        col.regions=(rev(brewer.pal(10, "RdBu"))), main = paste("", sep = " "), ylab.right = list(label = "Increase in SST extremes Per Year", vjust = -0.75, cex = 0.8),
                        xlim = c(-80,0), ylim = c(-50,0))

Canada_sea <- levelplot(trend ~ lon * lat, data=grid, at = cuts_trend,  cuts= 10, pretty=T, xlab = "Longitude", ylab = "Latitude",
                        col.regions=(rev(brewer.pal(10, "RdBu"))), main = paste("", sep = " "), ylab.right = list(label = "Increase in SST extremes Per Year", vjust = -0.75, cex = 0.8),
                        xlim = c(-180,0), ylim = c(30,90))

#plots of the trend located at Canada and Tasman Sea
#Where the trend is prevalent
grid.arrange(Canada_sea, Tasman_sea, nrow = 2)


