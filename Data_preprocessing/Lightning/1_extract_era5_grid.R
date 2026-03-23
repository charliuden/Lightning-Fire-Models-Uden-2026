#Script to extract gridlist from ERA5 climate data
#June 24th 2025

#This gridlist will be used to aggregate lightning points to the same spatial grid as the 
#ERA5 climate data used to drive the lightning prediction model. 

library(ncdf4)
library(ggplot2)
library(dplyr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

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

####################
# Load Alaska polygon
####################
usa <- ne_states(country = "United States of America", returnclass = "sf")
alaska <- usa[usa$name == "Alaska", ]
alaska <- st_transform(alaska, crs = 4326)  # Ensure it's in WGS84

####################
# Open NetCDF
####################
nc_file <- file.path(netcdf_root, "file_name.nc") #replace this with your netcdf file name 

# Read lat/lon
lat <- ncvar_get(nc, "latitude")   # 89, descending
lon <- ncvar_get(nc, "longitude")  # 169

# Get first time slice of temperature (optional for masking)
t2m <- ncvar_get(nc, "t2m", start = c(1, 1, 1), count = c(-1, -1, 1))

# Close file
nc_close(nc)

# Flip latitude and t2m if needed
lat <- rev(lat)
t2m <- t2m[, ncol(t2m):1]  # adjust t2m to match flipped lat

# Create grid of centroid points
lon_grid <- rep(lon, each = length(lat))
lat_grid <- rep(lat, times = length(lon))

df_points <- data.frame(lon = lon_grid, lat = lat_grid)
sf_points <- st_as_sf(df_points, coords = c("lon", "lat"), crs = 4326)

# Filter to only points within Alaska polygon
sf_land_points <- sf_points[st_within(sf_points, alaska, sparse = FALSE), ]

# Plot
ggplot() +
  geom_sf(data = alaska, fill = NA, color = "black") +
  geom_sf(data = sf_land_points, size = 0.5, color = "blue") +
  labs(title = "ERA5 Grid Centroids Within Alaska") +
  theme_minimal() +
  coord_sf(xlim = c(-180, -120), ylim = c(47, 80))

#Extract lon/lat again
land_coords <- st_coordinates(sf_land_points)
head(land_coords)

#write to csv
write.csv(land_coords, file.path(data_root, "era5_grid.csv")
