#------------------------------
#RASTER METHOD
#------------------------------
#July 26th 2024

#Script to aggregate Vaisala lightning strike data to PalEON data

# Lightning data are downloaded from:
#   https://fire.ak.blm.gov/predsvcs/maps.php

#2.1.1 Lightning 
# Lightning strike data, including location and time, were obtained from the Alaska Lightning 
# Detection Network (ALDN), maintained by the Alaska Fire Service (AFS) (U.S. Department of 
# the Interior, Bureau of Land Management (BLM), 2024, 2025). We used two datasets, the second 
# reflecting a major upgrade, to the detection network: an earlier dataset (2002-2011) and a newer 
# dataset (2012-2018). From 2002-2011, lightning was detected using an impact-based system 
# that reported cloud-to-ground lightning primarily as flashes. During this period, changes 
# in sensor technology, processing software, and network configuration resulted in variable 
# detection efficiency and relatively coarse spatial accuracy, particularly in remote regions 
# of Alaska. As a result, this dataset is characterized by lower positional accuracy. 
# Beginning in 2012, the ALDN transitioned to a time-of-arrival (TOA)-based system. 
# The TOA system records individual strokes and includes cloud-to-ground. 
# Due to differences in stroke versus flash reporting, expanded sensor coverage, and improved 
# detection range, the AFS reports detecting ~2.25 times more lightning events relative to the 
# earlier systems; there are 1,882,857 strikes in the earlier dataset and 1,784,663 strikes 
# in the later dataset. Due to these differences in detection methods and accuracy, 
# we train and test our models on each time period separately. 

# RASTER METHOD: turn the point grid into a raster and find number of lightning strikes (points) that fall into 
#each raster cell each day. 

library(raster)
library(rdgal)
library(dplyr)
library(sp)
library(ggplot2)
library(patchwork)

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

#output directory - comment out based on 2002-2011 vs 2012-2018 
data_file_path <- file.path(data_root,"lightning_counts_2002_2011.csv")
data_file_path <- file.path(data_root,"lightning_counts_2012_2018.csv")

#----------------------
#lightning data for Alaska 2012 - 2018
#----------------------

df <- read.table(file.path(data_root,"Historical_Lightning_2012_2024_TOA_WGS84.txt"), sep=",")
colnames(df) <- df[1,] #apply column names, which are in the first row
df <- df[-1,] #remove first row
df <- df %>% filter(
  STROKETYPE=="GROUND_STROKE") 
df <- df[c("UTCDATETIME", "LATITUDE", "LONGITUDE")] #only need time and location of each strike
df$UTCDATETIME <- as.POSIXct(df$UTCDATETIME, format = "%m/%d/%Y %H:%M:%S", tx="UTC")
df$UTCDATETIME <- as.Date(df$UTCDATETIME)
df$LATITUDE <- as.numeric(df$LATITUDE)
df$LONGITUDE <- as.numeric(df$LONGITUDE)

#only need june, july, august
nrow(df)
df$month <- as.integer(format(df$UTCDATETIME, "%m"))
df <- filter(df, month %in% c(6, 7, 8))
df <- df[c("UTCDATETIME", "LATITUDE", "LONGITUDE")] 
nrow(df)

#----------------------
#lightning data for Alaska 2001 - 2012
#----------------------

df <- read.table(file.path(data_root,"Historical_Lightning_1986_2012_ImpactSystem_NAD83.txt"), sep=",")

#tidy data
colnames(df) <- df[1,] #apply column names, which are in the first row
df <- df[-1,] #remove first row
df <- df[c("STRIKETIME", "LAT", "LON")] #only need time and location of each strike
df$STRIKETIME <- as.POSIXct(df$STRIKETIME, format = "%Y/%m/%d %H:%M", tx="UTC")
df$STRIKETIME <- as.Date(df$STRIKETIME) #only want day month year for the for loop
df$LAT <- as.numeric(df$LAT)
df$LON <- as.numeric(df$LON)
head(df)

#only need 2001 - 2012
df$year <- as.integer(format(df$STRIKETIME, "%Y"))
df <- filter(df, year<2012)
df <- filter(df, year>=2001)
# Filter for June, July, August
df$month <- as.integer(format(df$STRIKETIME, "%m"))
df <- filter(df, month %in% c(6, 7, 8))
df <- df[c("STRIKETIME", "LAT", "LON")] #only need these columns
head(df)
str(df)

#----------------------
#get era5 grid 
#----------------------
grid <- read.csv(file.path(data_root, "era5_grid.csv"))
head(grid)
grid <- grid[,2:3]
colnames(grid) <- c("lon", "lat")
plot(grid$lon, grid$lat)

#assign index value to grid cells 
grid <- cbind(grid, index=seq(1:nrow(grid)))
head(grid)
#use paleon grid to make spatial object, then raster
r <- grid
#coordinates(r) <- ~ lon + lat 
#r <- raster(r, ncol=26, nrow=14)
r <- rasterFromXYZ(r)
#values(r) <- runif(ncell(r))
crs(r) <- CRS("+proj=longlat +datum=WGS84")

#view point grid over raster
p <- grid
coordinates(p) <- ~ lon + lat
plot(r)
points(p)

#get lightning data
data <- df
colnames(data) <- c("date", "lat", "lon")
str(data)

#function to count points that fall into each raster cell
pointcount = function(r, pts){
  # make a raster of zeroes like the input
  r2 = r
  r2[] = 0
  # get the cell index for each point and make a table:
  counts = table(cellFromXY(r,pts))
  # fill in the raster with the counts from the cell index:
  r2[as.numeric(names(counts))] = counts
  return(r2)
}

#empty dataframe to hold output
lightningStrikes <- data.frame(matrix(nrow=0, ncol=4))
colnames(lightningStrikes) <- c("x", "y", "layer", "date")

#loop through days - select years that apply to the lightning dataset...
dates <- seq(as.Date("2012-06-01"), as.Date("2018-12-31"), by="days")
#dates <- seq(as.Date("2002-06-01"), as.Date("2011-08-31"), by="days")

# Keep only dates in June, July, or August
dates <- dates[format(dates, "%m") %in% c("06", "07", "08")]

for(day in dates){
  
  pts <- filter(data, date==day)
  
  if (nrow(pts) > 0) {
    
    #convert day's data to spatial object
    coordinates(pts) <- ~ lon + lat
    crs(pts) <- CRS("+proj=longlat +datum=WGS84")
    
    #use point count function from above
    #r2 <- pointcount(r, pts)
    #--
    r2 = r
    r2[] = 0
    # get the cell index for each point and make a table:
    counts = table(cellFromXY(r,pts))
    # fill in the raster with the counts from the cell index:
    r2[as.numeric(names(counts))] = counts
    #--
    
  } else if (nrow(pts) == 0){
    
    #if there was no lightning that day, make a raster of 0s
    r2 <- r
    r2[] <- 0
    
  }
  
  df <- as.data.frame(r2, xy=TRUE)
  df <- cbind(df, date=rep(day, nrow(df)))
  lightningStrikes <- rbind(lightningStrikes, df)
}

#convert date column from numeric to date, rename columns 
lightningStrikes$date <- as.Date(lightningStrikes$date, origin="1970-01-01")
colnames(lightningStrikes) <- c("lon", "lat", "strikes", "date")
head(lightningStrikes)
hist(lightningStrikes$strikes)

#----------
#map to check that it worked:

library(rnaturalearth)
library(rnaturalearthdata)

# Load Alaska polygon
usa <- ne_states(country = "United States of America", returnclass = "sf")
alaska <- usa[usa$name == "Alaska", ]
alaska <- st_transform(alaska, crs = 4326)  # Ensure it's in WGS84

lightningStrikes_sf <- st_as_sf(filter(lightningStrikes, date=="2012-06-01"), coords = c("lon", "lat"), crs = 4326)
#lightningStrikes_sf <- st_as_sf(filter(lightningStrikes, date=="2002-07-01"), coords = c("lon", "lat"), crs = 4326)

# Filter to only points within Alaska polygon
lightningStrikes_sf <- lightningStrikes_sf[st_within(lightningStrikes_sf, alaska, sparse = FALSE), ]

ggplot() +
  geom_sf(data=lightningStrikes_sf, aes(color=strikes, alpha=0.01)) +
  labs(title = "Strike counts on a single day") +
  theme_minimal() +
  guides(alpha = FALSE) +
  scale_color_gradient(low = "white", high = "blue") +
  geom_sf(data = alaska, fill = NA, color = "black") +
  coord_sf(xlim = c(-180, -120), ylim = c(47, 80)) 

#---------------
#filter only points in lightingStrikes that are in the Alaska polygon
# Filter to only points within Alaska polygon
lightningStrikes_sf <- st_as_sf(lightningStrikes, coords = c("lon", "lat"), crs = 4326)
lightningStrikes_sf <- lightningStrikes_sf[st_within(lightningStrikes_sf, alaska, sparse = FALSE), ]

ggplot() +
  geom_sf(data=filter(lightningStrikes_sf, date=="2012-06-01"), aes(color=strikes, alpha=0.01)) +
  labs(title = "Strike counts on a single day") +
  theme_minimal() +
  guides(alpha = FALSE) +
  scale_color_gradient(low = "white", high = "blue") +
  geom_sf(data = alaska, fill = NA, color = "black") +
  coord_sf(xlim = c(-180, -120), ylim = c(47, 80)) 

#back to dataframe:
lightningStrikes_trimmed <- lightningStrikes_sf %>%
  mutate(lon = st_coordinates(.)[, 1],
         lat = st_coordinates(.)[, 2]) %>%
  st_drop_geometry()

#---------------
#write to csv 
write.csv(lightningStrikes_trimmed, data_file_path)



