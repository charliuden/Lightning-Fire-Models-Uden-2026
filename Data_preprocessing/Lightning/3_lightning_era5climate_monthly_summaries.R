#Charlotte Uden
#August 12 2024

library(ggplot2)
library(dplyr)
library(terra)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

# Load configuration files
# - location of model drivers, rds files (containing models), tables of parameter estimates and performance metrics, and figures. 

read_properties <- function(file_path) {
  lines <- readLines(file_path)
  lines <- lines[grepl("=", lines)]  # only keep lines with key-value pairs
  key_vals <- strsplit(lines, "=")
  props <- setNames(
    trimws(sapply(key_vals, `[`, 2)),
    trimws(sapply(key_vals, `[`, 1))
  )
  return(props)
}

#config <- read_properties("/path/to/config/file/my_config.properties")  # path to your config file

# Extract values from config
drivers_root <- config[["drivers_root"]]
rds_root <- config[["rds_root"]]
table_root <- config[["table_root"]]
predictions_root <- config[["predictions_root"]]
figure_root <- config[["figure_root"]]
data_root <- config[["data_root"]]
netcdf_root <- config[["netcdf_root"]]

#climate data file location
climate_file_path <- file.path(data_root, "ERA5_climate_and_fuel_2002_2018.csv")

#comment out to apply code to each period (2002-2011 versus 2012-2018)
#---2002-2011---
input_file_path <- file.path(data_root, "lightning_counts_2002_2011.csv")
lightning_count_file_path <- file.path(data_root, "lightning_monthlySummaries_2002_2011.csv")
output_file_path <- file.path(drivers_root, "lightning_prediction_drivers_2002_2011.csv")
#---2012-2018---
input_file_path <- file.path(data_root, "lightning_counts_2012_2018.csv")
lightning_count_file_path <- file.path(data_root, "lightning_monthlySummaries_2012_2018.csv")
output_file_path <- file.path(drivers_root, "lightning_prediction_drivers_2012_2018.csv")

#read in lightning data
df <- read.csv(input_file_path)
df <- df[,2:5]
df$date <- as.Date(df$date)
head(df)
str(df)


#------------------------------------------------------------------------------------------------------
#-----1. get gridcell area -----
#------------------------------------------------------------------------------------------------------

# Get USA polygon
usa <- ne_states(country = "United States of America", returnclass = "sf")
alaska <- usa[usa$name == "Alaska", ]
alaska <- st_transform(alaska, crs = 4326)
map <- ggplot() + geom_sf(data = alaska, fill = "lightgray", color="white") +
  coord_sf(xlim = c(-180, -120), ylim = c(47, 80))
map

#convert points to raster and set crs

#get grid
points <- filter(df, date=="2001-06-01")
map + geom_point(data=points, aes(x=lon, y=lat, col=strikes, alpha=0.01)) +
  scale_color_gradient(low = "white", high = "blue") 

#add unique id to each point
points$id <- seq(1:nrow(points))
#convert to raster 
r <- rast(points[,c("lon", "lat", "strikes")])
plot(r)
crs(r) <- "+init=epsg:4326"
#now make a raster that includes the id column
r <- rast(points[,c("lon", "lat", "id")])
plot(r)
crs(r) <- "+init=epsg:4326"
#convert to polygons
p = as.polygons(r, aggregate=FALSE) #if aggregate=false, each raster cell is a polygon 
plot(p)    
#get area of each polygon
area <- expanse(p, unit="km", transform=T) #	If transform=TRUE, planar CRS are transformed to lon/lat for accuracy. can also use "m" instead of "km"
p$area <- area

plot(p, "area")
plot(p, "id")

d <- data.frame(values(p))
head(d)

#merge grid and area values by id
d <- merge(points[,c("lat", "lon", "id")], d, by="id")

#merge area values with data by lat lon
df <- left_join(df, d, by=c("lat", "lon"))

#map to check it worked
points <- filter(df, date=="2006-06-01")
map + geom_point(data=points, aes(x=lon, y=lat, alpha = 0.01, col=strikes)) + scale_color_gradient(low = "white", high = "blue") 
map + geom_point(data=points, aes(x=lon, y=lat, col=area)) + scale_color_gradient(low = "white", high = "blue") 


#------------------------------------------------------------------------------------------------------------------------------
# calculate monthly (JJA) and summer lightning strike rate (km-2 month-1) from daily strike count and grid cell area -----
#------------------------------------------------------------------------------------------------------------------------------

#df$date <- as.Date(df$date)
df$month <- as.numeric(format(df$date,'%m'))
df$year <- as.numeric(format(df$date,'%Y'))
#df$day <- as.numeric(format(df$date,'%d'))
#df$week <- week(df$date)

head(df)

cols <- c("lon", "lat", "year", "month", "area")

#sum strike count for each summer month 
summary_strike <- df %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(strikes_monthly_sum = sum(strikes), .groups = 'drop')

#calculate strike rate for each month (km-2 month-1)
summary_strike <- mutate(summary_strike, monthly_strike_rate = strikes_monthly_sum/area)

hist(summary_strike$monthly_strike_rate)

#for figure 6:
write.csv(summary_strike, lightning_count_file_path)

#caluclate summer average stirke rate

cols <- c("lon", "lat", "year")

#calculate mean strike rate across all (summer) months
summary_strike <- summary_strike %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(summer_monthly_strike_rate = mean(monthly_strike_rate), .groups = 'drop')

hist(summary_strike$summer_monthly_strike_rate)


#------------------------------------------------------------------------------------------------------
# ----- calculate climate monthly (june, july, august) averages -----
#------------------------------------------------------------------------------------------------------

#-----------------
#get climate data
#-----------------

climate <- read.csv(climate_file_path)
climate$date <- as.Date(climate$date)
climate$month <- as.numeric(format(climate$date,'%m'))
climate$year <- as.numeric(format(climate$date,'%Y'))


#------------------------------------------------------------------------------------------------------
#calculate climate summer averages -----
#------------------------------------------------------------------------------------------------------

cols <- c("lon", "lat", "year")

#calculate total precipitation across all summer months
summary_precip <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(precip = mean(precipitation), .groups = 'drop')

#temperature
summary_temp<- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(temp = mean(temperature), .groups = 'drop')

#wind
summary_wind <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(wind = mean(wind), .groups = 'drop')

#msdwswrf
summary_swr <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(swr = mean(swr), .groups = 'drop')

#sp
summary_sp <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(sp = mean(sp), .groups = 'drop')

#rh
summary_rh <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(rh = mean(relative_humidity), .groups = 'drop')

#cape
summary_cape <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(cape = mean(cape), .groups = 'drop')

#cxp
summary_cxp <- climate %>% 
  group_by(across(all_of(cols))) %>% 
  summarize(cxp = mean(cxp), .groups = 'drop')


#merge summaries
summary_climate <- merge(summary_temp, summary_precip, by=cols)
summary_climate <- merge(summary_climate, summary_wind, by=cols)
summary_climate <- merge(summary_climate, summary_swr, by=cols)
summary_climate <- merge(summary_climate, summary_sp, by=cols)
summary_climate <- merge(summary_climate, summary_rh, by=cols)
summary_climate <- merge(summary_climate, summary_cape, by=cols)
summary_climate <- merge(summary_climate, summary_cxp, by=cols)
head(summary_climate)

#---clip to alaska polygon:---

#make summary df an sf object
summary_climate_sf <- st_as_sf(summary_climate, coords = c("lon", "lat"), crs = 4326)

# Filter to only points within Alaska polygon
summary_climate_sf <- summary_climate_sf[st_within(summary_climate_sf, alaska, sparse = FALSE), ]

ggplot() +
  geom_sf(data=filter(summary_climate_sf, year==2001), aes(color=precip, alpha=0.01)) +
  theme_minimal() +
  guides(alpha = FALSE) +
  scale_color_gradient(low = "white", high = "blue") +
  geom_sf(data = alaska, fill = NA, color = "black") +
  coord_sf(xlim = c(-180, -120), ylim = c(47, 80)) 

#back to dataframe
summary_climate_trimmed <- summary_climate_sf %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()


df <- merge(summary_strike, summary_climate_trimmed, by=c("year", "lon", "lat"))

#write to csv file
write.csv(df, output_file_path)
