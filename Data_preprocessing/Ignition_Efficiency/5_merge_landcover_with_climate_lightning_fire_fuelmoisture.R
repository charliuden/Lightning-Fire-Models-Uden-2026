library(ggplot2)
library(patchwork)
library(dplyr)
library(lubridate)


#Load configuration files
#- location of model drivers, rds files (containing models), tables of parameter estimates and performance metrics, and figures. 

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

# Use the values to build full paths. Comment out the period not being processed:
#--- 2002_2011 ---
climate_fuel_path <- file.path(data_root, "binary_lightning_fire_climate_matched_2002_2011.csv")
data_out_path <- file.path(drivers_root, "ignition_efficiency_drivers_2002_2011.csv")
#--- 2012_2018 ---
climate_fuel_path <- file.path(data_root, "binary_lightning_fire_climate_matched_2012_2018.csv")
data_out_path <- file.path(drivers_root, "ignition_efficiency_drivers_2012_2018.csv")

#------------------------------------------------
#First, peek at land cover data for each year:
#------------------------------------------------
file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2005.csv")
df <- read.csv(file_path)

str(df)

p1 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=C)) + ggtitle("2005 conifer")

p2 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=B)) + ggtitle("2005 broadleaf")

p3 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=U)) + ggtitle("2005 open land")

p1 + p2 + p3

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2010.csv")
df <- read.csv(file_path)

p4 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=C)) + ggtitle("2010 conifer")

p5 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=B)) + ggtitle("2010 broadleaf")

p6 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=U)) + ggtitle("2010 open land")

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2015.csv")
df <- read.csv(file_path)

p7 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=C)) + ggtitle("2015 conifer")

p8 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=B)) + ggtitle("2015 broadleaf")

p9 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=U)) + ggtitle("2015 open land")

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2020.csv")
df <- read.csv(file_path)

p10 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=C)) + ggtitle("2020 conifer")

p11 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=B)) + ggtitle("2020 broadleaf")

p12 <- ggplot()+geom_point(data=df, aes(x=lon, y=lat, color=U)) + ggtitle("2020 open land")


(p1 | p2 | p3) / (p4 | p5 | p6) / (p7 | p8 | p9) / (p10 | p10 | p12) 

#------------------------------------------------
#join all years
#------------------------------------------------

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2020.csv")
df_2005 <- read.csv(file_path)
year <- rep(2005, nrow(df_2005))
df_2005 <- cbind(df_2005, year)

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2010.csv")
df_2010 <- read.csv(file_path)
year <- rep(2010, nrow(df_2010))
df_2010 <- cbind(df_2010, year)

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2015.csv")
df_2015 <- read.csv(file_path)
year <- rep(2015, nrow(df_2015))
df_2015 <- cbind(df_2015, year)

file_path <- file.path(data_root, "lc_on_era5_grid/landcover_ERA5_2020.csv")
df_2020 <- read.csv(file_path)
year <- rep(2020, nrow(df_2020))
df_2020 <- cbind(df_2020, year)

df <- rbind(df_2005, df_2010, df_2015, df_2020)


ggplot() + geom_boxplot(data=df, aes(x=as.factor(year), y=U))

ggplot() + geom_boxplot(data=df, aes(x=as.factor(year), y=C))

ggplot() + geom_boxplot(data=df, aes(x=as.factor(year), y=B))

#write.csv(df, file.path(data_root, "landcover_ERA5_2005-2020.csv")

#------------------------------------------------
#read in binary lightning caused fire, climate, and fuel moisture
#------------------------------------------------
df_2 <- read.csv(climate_fuel_path)

head(df_2)
#note that lat.1, lon.1, and date.1 are the locations of the locations of the era5 cliamte data. 
#the other lat, lon, date columns are for the lightning locations. 

#we want to match the locations of land cover data with lat.1, lon.1, date.1, 
#since this is the ERA5 grid, the same grid as the landcover data.

# df_2: lightning rows with ERA5 grid coords in lat.1/lon.1 and date.1
# df: landcover rows with lon/lat/C/B/U/year for years 2005/2010/2015/2020

#------------------------------------------------
#merge binary lightning outcome, landcover, climate, and fuel moisture by year and location
#------------------------------------------------
landcover_years <- c(2005, 2010, 2015, 2020)

# Helper: map any year -> nearest in landcover_years
nearest_landcover_year <- function(y, available = landcover_years) {
  # If tie (e.g., 2012.5), this chooses the earlier year; adjust if you prefer later.
  available[ max.col(-abs(outer(y, available, "-")), ties.method = "first") ]
}

# Make sure types are consistent
df_lc <- df %>%
  mutate(
    year = as.integer(year),
    lat  = as.numeric(lat),
    lon  = as.numeric(lon),
    C = as.numeric(C),
    B = as.numeric(B),
    U = as.numeric(U)
  )

df_merged <- df_2 %>%
  mutate(
    date.1 = as.Date(date.1),                 # if it's already Date, harmless
    era5_year = year(date.1),
    lc_year = nearest_landcover_year(era5_year),
    lat_key = as.numeric(lat.1),
    lon_key = as.numeric(lon.1)
  ) %>%
  left_join(
    df_lc %>%
      transmute(lon_key = lon, lat_key = lat, lc_year = year, C, B, U),
    by = c("lon_key", "lat_key", "lc_year")
  )

# Quick diagnostics
cat("Rows in df_2:", nrow(df_2), "\n")
cat("Rows after merge:", nrow(df_merged), "\n")
cat("Missing landcover (C is NA):", sum(is.na(df_merged$C)), "\n")

head(df_merged)

write.csv(df_merged, data_out_path)







