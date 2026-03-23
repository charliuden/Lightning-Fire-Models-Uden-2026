import pandas as pd
from scipy.spatial import cKDTree
from datetime import datetime
from tqdm import tqdm
from pyproj import Transformer

#---------------------------------------------------------------------
#function to read config file
#---------------------------------------------------------------------
def read_properties(file_path):
    props = {}
    
    with open(file_path, "r") as f:
        for line in f:
            line = line.strip()
            
            # Skip empty lines and comments
            if not line or "=" not in line:
                continue
            
            key, value = line.split("=", 1)  # split only on first "="
            props[key.strip()] = value.strip()
    
    return props


# Read config
config = read_properties("/path/to/config/file/my_config.properties")

# Extract values
drivers_root = config.get("drivers_root")
rds_root = config.get("rds_root")
table_root = config.get("table_root")
predictions_root = config.get("predictions_root")
figure_root = config.get("figure_root")
data_root = config.get("data_root")
netcdf_root = config.get("netcdf_root")

#build file paths
from pathlib import Path

#--- climate data---
climate_data_path = Path(data_root) / "ERA5_climate_and_fuel_2002_2018.csv"

#--- lightning data--- comment out depending on period
# ---2002-2011 ---
data_in_path = Path(data_root) / "binary_fire_lightning_matches_2002_2011.csv"
start_date = pd.to_datetime("2002-01-01")
end_date = pd.to_datetime("2011-12-31")
data_out_path = Path(data_root) / "inary_lightning_fire_climate_matched_2002_2011.csv"

# ---2012-1018 ---
#data_in_path = Path(data_root) / "binary_fire_lightning_matches_2012_2018.csv"
#start_date = pd.to_datetime("2012-01-01")
#end_date = pd.to_datetime("2018-12-31")
#data_out_path = Path(data_root) / "inary_lightning_fire_climate_matched_2012_2018.csv"

# --- Load Data ---
# Load lightning data
#sampled negative strikes, 1:1.5
lightning = pd.read_csv(data_in_path, parse_dates=["date"])
print(lightning.columns)
print(lightning.head)
lightning["date"] = pd.to_datetime(lightning["date"], format="%Y-%m-%d", errors="coerce")

# Load cliamte data
climate = pd.read_csv(climate_data_path, parse_dates=["date"])
print(climate.columns)
climate["date"] = pd.to_datetime(climate["date"], format="%Y-%m-%d", errors="coerce")

print(lightning.head)
print(climate.head)

# --- Filter dataset to start and end dates ---
lightning = lightning[(lightning["date"] >= start_date) & (lightning["date"] <= end_date)]
climate = climate[(climate["date"] >= start_date) & (climate["date"] <= end_date)]

# --- Group by date ---
climate_groups = climate.groupby("date")
lightning_groups = lightning.groupby("date")

matched_rows = []

print("Matching lightning to nearest climate grid point by date...")

for date, strikes_on_day in tqdm(lightning_groups):
    if date not in climate_groups.groups:
        continue

    daily_climate = climate_groups.get_group(date).reset_index(drop=True)
    tree = cKDTree(daily_climate[["lat", "lon"]].values)

    distances, indices = tree.query(strikes_on_day[["lat", "lon"]].values, k=1)
    matched_climate = daily_climate.loc[indices].reset_index(drop=True)
    
    merged = pd.concat([strikes_on_day.reset_index(drop=True), matched_climate.reset_index(drop=True)], axis=1)
    matched_rows.append(merged)

# --- Combine all matched data ---
final_df = pd.concat(matched_rows, ignore_index=True)

# --- Select and rename relevant columns ---
final_df = final_df[[
    "lat", "lon", "date", "caused_fire", "temperature", "dewpoint", "precipitation",
    "wind", "cape", "cxp", "relative_humidity", "FFMC", "DMC", "DC"
]]

print(final_df.head)
#final_df["date"] = pd.to_datetime(final_df["date"], format="%Y-%m-%d", errors="coerce")
#final_df["year"] = final_df["date"].dt.year

# --- Save as CSV ---
final_df.to_csv(data_out_path, index=False)
print(f"✅ Matched and saved {len(final_df)} rows to 'lightning_climate_matched.csv'")
