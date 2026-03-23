import pandas as pd
import numpy as np
import fastparquet

#script to combine drivers for P_fire - we need all dirvers for r_strike and P_ignited in one dataframe

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
config = read_properties("/raid/cuden/config_files/lightning_fire_config.properties")

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

met_path = Path(data_root) / "Alaska_ERA5_climate_and_fuel_2002_2018.csv"
veg_path = Path(data_root) / "landcover_ERA5grid_2005-2020.csv"
met_r_strike_path = Path(drivers_root) / "lightning_prediction_drivers_2002_2011.csv"

#output paths
daily_out_path = Path(drivers_root) / "fire_probability_daily_drivers_2002_2011.parquet"
seasonal_out_path = Path(drivers_root) / "fire_probability_seasonal_drivers_2002_2011.parquet")

#---------------------------------------------------------------------

df_met = pd.read_csv(met_path)
df_veg = pd.read_csv(veg_path)
df_met_r_strike = pd.read_csv(met_r_strike_path)

#---------------------------------------------------------------------

def attach_veg_to_met(df_met, df_veg,
                      year_min=2002,
                      year_max=2011,
                      lat_decimals=2,
                      lon_decimals=2):
    """
    Merge daily met/fire drivers with vegetation data by (lat, lon)
    and nearest vegetation year, restricted to [year_min, year_max].
    """

    df_met = df_met.copy()
    df_veg = df_veg.copy()

    # Drop the index-like column if present
    df_veg = df_veg.drop(columns=['Unnamed: 0'], errors='ignore')

    # --- Ensure types are sensible ---
    # 1. Ensure date is datetime and extract "year"
    df_met['date'] = pd.to_datetime(df_met['date'])
    df_met['year'] = df_met['date'].dt.year.astype(int)

    # If veg year is float or string, coerce to int as well
    df_veg['year'] = df_veg['year'].astype(int)

    # 2. Keep only years 2002–2011 (or whatever range you specify)
    mask = (df_met['year'] >= year_min) & (df_met['year'] <= year_max)
    df_met = df_met.loc[mask].copy()

    # 3. Round lat/lon (helps avoid float precision mismatches)
    for df in (df_met, df_veg):
        df['lat'] = df['lat'].round(lat_decimals)
        df['lon'] = df['lon'].round(lon_decimals)

    # 4. Map each met year to the closest vegetation year
    veg_years = np.array(sorted(df_veg['year'].unique()))
    unique_met_years = np.sort(df_met['year'].unique())

    # Quick print so you can see what’s going on
    print("Veg years:", veg_years)
    print("Met years (after filter):", unique_met_years)

    year_to_veg_year = {
        y: veg_years[np.abs(veg_years - y).argmin()]
        for y in unique_met_years
    }
    print("Year mapping (met -> veg):", year_to_veg_year)

    df_met['veg_year'] = df_met['year'].map(year_to_veg_year)

    # 5. Rename veg 'year' to 'veg_year' and merge
    df_veg = df_veg.rename(columns={'year': 'veg_year'})

    df_merged = df_met.merge(
        df_veg,
        how='left',
        on=['lat', 'lon', 'veg_year'],
        suffixes=('', '_veg')
    )

    # --- Debug: how many rows actually matched? ---
    if 'C' in df_veg.columns:
        n_non_null = df_merged['C'].notna().sum()
        print(f"Non-null C after merge: {n_non_null} / {len(df_merged)} rows")
    else:
        print("Warning: 'C' not found in df_veg columns:", df_veg.columns)

    return df_merged


#fucntion to add unique key to lat/lon/year
def add_cell_key(df,
                 year_col="year",
                 lat_col="lat",
                 lon_col="lon",
                 key_name="cell_id",
                 lat_decimals=2,
                 lon_decimals=2):
    """
    Add a stable key like '2005_65.25_-150.00' based on (year, lat, lon).

    Rounds lat/lon to avoid floating point mismatch and formats them
    with fixed decimals so the string keys line up between dataframes.
    """
    df = df.copy()

    # round to ensure consistency
    df[lat_col] = df[lat_col].round(lat_decimals)
    df[lon_col] = df[lon_col].round(lon_decimals)

    df[key_name] = (
        df[year_col].astype(int).astype(str) + "_" +
        df[lat_col].map(lambda x: f"{x:.{lat_decimals}f}") + "_" +
        df[lon_col].map(lambda x: f"{x:.{lon_decimals}f}")
    )

    return df

# --- Amerge vegetation df and cliamte/fuel moisture df ---
df_merged = attach_veg_to_met(df_met, df_veg, 
                              year_min=2002,
                              year_max=2011,
                              lat_decimals=2,
                              lon_decimals=2)

# --- Add keys to the daily dataframe (df_merged) ---
df_merged = add_cell_key(df_merged,
                         year_col="year",
                         lat_col="lat",
                         lon_col="lon",
                         key_name="cell_id",
                         lat_decimals=2,
                         lon_decimals=2)

# --- Add keys to the summer-summaries dataframe (df_met_r_strike) ---

# Drop index-like column if present
df_met_r_strike = df_met_r_strike.drop(columns=["Unnamed: 0"], errors="ignore")

# select columns we need
cols = ["year","lon","lat", "temp", "precip", "wind", "swr", "sp", "rh"]
df_met_r_strike = df_met_r_strike[cols]

# Make sure lat/lon and year are in the expected form
df_met_r_strike["year"] = df_met_r_strike["year"].astype(int)

df_met_r_strike = add_cell_key(df_met_r_strike,
                               year_col="year",
                               lat_col="lat",
                               lon_col="lon",
                               key_name="cell_id",
                               lat_decimals=2,
                               lon_decimals=2)
# After the merge (df_merged already created):

n_total = len(df_merged)
n_with_veg = df_merged['C'].notna().sum()

print("Total rows:", n_total)
print("Rows with veg data (C not null):", n_with_veg)
print("Fraction with veg data:", n_with_veg / n_total)

met_cells = df_merged[['lat', 'lon']].drop_duplicates()
veg_cells = df_veg[['lat', 'lon']].drop_duplicates()

print("Unique met cells:", len(met_cells))
print("Unique veg cells:", len(veg_cells))

# How many met cells are in veg grid at all?
merged_cells = met_cells.merge(veg_cells, on=['lat', 'lon'], how='inner')
print("Unique overlapping cells:", len(merged_cells))
print("Overlap fraction (met):", len(merged_cells) / len(met_cells))

print(df_merged.head)
print(df_met_r_strike.head)
print(df_merged.columns)
print(df_met_r_strike.columns)

#export
df_merged.to_parquet(daily_out_path, 
                     index=False, 
                     engine="fastparquet"
                     )

df_met_r_strike.to_parquet(
    seasonal_out_path,
    index=False,
    engine="fastparquet"
)

