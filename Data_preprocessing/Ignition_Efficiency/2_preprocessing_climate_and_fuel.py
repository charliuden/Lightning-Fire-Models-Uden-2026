#Charlotte Uden
#script to extract era5 climate variables from netcdf files. 
#then calculate daily mean relative humidity, fuel moisture codes
#and daily values of mean temp, precip, shortwave radiation, surface pressure, and wind. 

import numpy as np
import pandas as pd

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
netcdf_path = Path(netcdf_root)
data_path = Path(data_root) / "ERA5_climate_and_fuel_2002_2018.csv"

#---------------------------------------------------------------------
#read in data
#---------------------------------------------------------------------
def read_and_merge_climate_csvs(base_filenames, start_year=2002, end_year=2018, path=""):
    """
    Reads and merges CSVs for multiple variables across multiple years,
    renaming columns to concise names.

    Parameters:
        base_filenames (list of str): Base names like ['2m_temperature', ...].
        start_year (int): Start year.
        end_year (int): End year.
        path (str): Directory path for CSVs.

    Returns:
        pd.DataFrame: Merged DataFrame with renamed climate columns.
    """
    # Mapping from raw variable names to short column names
    rename_map = {
        "2m_temperature": "temperature",
        "2m_dewpoint_temperature": "dewpoint",
        "mean_total_precipitation_rate": "precipitation", #total precip
        "instantaneous_10m_wind_gust": "wind",
        "mean_surface_downward_short_wave_radiation_flux": "swr",
        "surface_pressure": "sp",
        "convective_available_potential_energy" : "cape"
    }

    all_years_data = []

    for year in range(start_year, end_year + 1):
        yearly_data = None

        for var in base_filenames:
            filename = f"{path}{var}_{year}_JJA_dailygrid.csv"

            try:
                df = pd.read_csv(filename)

                # Drop 'number' if it exists
                df = df.drop(columns=['number'], errors='ignore')

                # Identify the variable column (exclude date/lat/lon)
                value_col = [col for col in df.columns if col not in ["date", "latitude", "longitude"]][0]

                # Rename that column to the simple name
                simple_name = rename_map.get(var, var)  # fallback: use original var if no match
                df = df.rename(columns={
                    "latitude": "lat",
                    "longitude": "lon",
                    value_col: simple_name
                })

                if yearly_data is None:
                    yearly_data = df
                else:
                    yearly_data = pd.merge(yearly_data, df, on=["date", "lat", "lon"], how="inner")

            except FileNotFoundError:
                print(f"Warning: Missing file for {var} in {year} — skipping.")

        if yearly_data is not None:
            yearly_data["year"] = year
            all_years_data.append(yearly_data)

    final_df = pd.concat(all_years_data, ignore_index=True)
    return final_df

variables = [
    "2m_temperature",
    "2m_dewpoint_temperature",
    "mean_total_precipitation_rate", #total precip
    "instantaneous_10m_wind_gust",
    "mean_surface_downward_short_wave_radiation_flux",
    "surface_pressure",
    "convective_available_potential_energy"
]

df = read_and_merge_climate_csvs(variables, path=file_path)

#convert date to datetime structure
df['date'] = pd.to_datetime(df['date'])

print(df.head())

#---------------------------------------------------------------------
#convert total precip from meters to mm
#---------------------------------------------------------------------
#df["precipitation"] = (
    #df.precipitation * 1000
#)

#print(df.head())

#---------------------------------------------------------------------
#calculate CAPE x precip
#---------------------------------------------------------------------
df["cxp"] = (
    df.precipitation * df.cape
)

print(df.head())

#---------------------------------------------------------------------
##convert temp and dew point temp in kelvin to celcius
#---------------------------------------------------------------------
df["dewpoint"] = (
    df.dewpoint - 273.15
)

df["temperature"] = (
    df.temperature - 273.15
)

#---------------------------------------------------------------------
#calculate Relative humidity (%) from dew point temperature and temperature
#---------------------------------------------------------------------
beta = 17.625
gamma = 243.04

df["relative_humidity"] = (
    100*(np.exp((beta * df.dewpoint) / (gamma + df.dewpoint)) / np.exp((beta * df.temperature)/(gamma + df.temperature)))
)

#limit rh (%) to values between 0 adn 100
df["relative_humidity"] = df["relative_humidity"].clip(0, 100)

print(df.head())


#---------------------------------------------------------------------
#compute fine fuel moisture code
#https://wikifire.wsl.ch/tiki-index91f7.html?page=Fine+fuel+moisture+code
#---------------------------------------------------------------------

def compute_ffmc(df, initial_ffmc=85.0):
    df = df.sort_values(['lat', 'lon', 'date']).copy()

    # Ensure FFMC column exists for shifting (will be overwritten anyway)
    if 'FFMC' not in df.columns:
        df['FFMC'] = np.nan

    # Get previous day's FFMC per location, defaulting to initial_ffmc for first day
    df['prev_ffmc'] = df.groupby(['lat', 'lon'])['FFMC'].shift(1)
    df['prev_ffmc'] = df['prev_ffmc'].fillna(initial_ffmc)

    # Step 1: Convert FFMC to moisture content
    mo = 147.2 * (101.0 - df['prev_ffmc']) / (59.5 + df['prev_ffmc'])

    # Step 2: Effective rainfall
    rf = df['precipitation']
    re = np.where(rf > 0.5, rf - 0.5, 0.0)

    # Step 3: Moisture content after rain (with safe masking)
    mo_r = mo.copy()
    mask = re > 0
    safe_re = np.where(mask, re, 1.0)  # avoid div by zero in the unused branch

    mo_r[mask] = mo[mask] + 42.5 * safe_re[mask] * np.exp(-100.0 / (251.0 - mo[mask])) * (1 - np.exp(-6.93 / safe_re[mask]))

    # Step 4: Drying (evaporation and wind)
    ed = 0.942 * (df['relative_humidity'] ** 0.679) + (11.0 * np.exp((df['temperature'] - 32.8) / 17.3)) + 0.18 * (21.1 - df['temperature']) * (1 - (df['relative_humidity'] / 100.0) ** 0.5)
    ko = 0.424 * (1.0 - ((100.0 - df['relative_humidity']) / 100.0) ** 1.7) + (0.0694 * df['wind'] ** 0.5) * (1.0 - np.exp(-0.115 * mo_r))
    k = ko * 0.581 * np.exp(0.0365 * df['temperature'])

    # Step 5: Moisture content after drying
    m = mo_r + (ed - mo_r) * 10 ** (-k)

    # Step 6: Convert back to FFMC
    ffmc_new = 59.5 * (250.0 - m) / (147.2 + m)

    # Clip FFMC to [0, 101]
    df['FFMC'] = np.clip(ffmc_new, 0, 101)

    # Clean up intermediate columns
    df = df.drop(columns=['prev_ffmc'])

    return df

df = compute_ffmc(df, initial_ffmc=85.0)

print(df.head())

print("Fine Fuel Moisture Code Calculated")

#---------------------------------------------------------------------
#compute duff moisture code
#https://wikifire.wsl.ch/tiki-index9436.html?page=Duff+moisture+code
#---------------------------------------------------------------------
def compute_dmc(df, initial_dmc=6.0):
    df = df.sort_values(['lat', 'lon', 'date']).copy()

    # Ensure DMC column exists
    if 'DMC' not in df.columns:
        df['DMC'] = np.nan

    # Get previous day's DMC
    df['prev_dmc'] = df.groupby(['lat', 'lon'])['DMC'].shift(1)
    df['prev_dmc'] = df['prev_dmc'].fillna(initial_dmc)

    # Ensure month column
    df['month'] = pd.to_datetime(df['date']).dt.month

    # Monthly day length factor L_f (for Northern Hemisphere, from standard FWI tables)
    L_f_table = {
        1: 6.5, 2: 7.5, 3: 9.0, 4: 12.8, 5: 14.8, 6: 16.0,
        7: 16.0, 8: 14.0, 9: 12.0, 10: 10.5, 11: 9.0, 12: 7.0
    }
    df['L_f'] = df['month'].map(L_f_table)

    # Effective rainfall
    rf = df['precipitation']
    re = np.where(rf > 1.5, rf - 1.5, 0.0)

    # Moisture content from previous DMC
    m_prev = 20.0 + np.exp(5.6348 - df['prev_dmc'] / 43.43)

    # Updated moisture after rain
    b = np.where(df['prev_dmc'] <= 33.0,
                 100.0 / (0.5 + 0.3 * df['prev_dmc']),
                 6.2 * np.log(df['prev_dmc']) - 17.2)

    m_rain = np.where(re > 0,
                      m_prev + (1000.0 * re) / (48.77 + b * re),
                      m_prev)

    # Convert moisture back to DMC
    dmc_rain = 43.43 * (5.6348 - np.log(m_rain - 20.0))

    # Drying rate k
    k = (1.894 * (df['temperature'] + 1.1) *
         (100.0 - df['relative_humidity']) *
         df['L_f'] * 1e-6)

    # Final DMC: apply drying or use post-rain value
    dmc_new = np.where(re > 0,
                       dmc_rain,
                       df['prev_dmc'] + k)

    # Clip to valid range
    df['DMC'] = np.clip(dmc_new, 0, None)

    # Drop intermediate columns
    df = df.drop(columns=['prev_dmc', 'month', 'L_f'])

    return df

df = compute_dmc(df, initial_dmc=6.0)

print(df.head())

print("Duff Moisture Code Calculated")

#---------------------------------------------------------------------
#compute drought code
#https://wikifire.wsl.ch/tiki-indexd5c6.html?page=Drought+code
#---------------------------------------------------------------------
def compute_dc(df, initial_dc=15.0):
    df = df.sort_values(['lat', 'lon', 'date']).copy()

    if 'DC' not in df.columns:
        df['DC'] = np.nan

    df['prev_dc'] = df.groupby(['lat', 'lon'])['DC'].shift(1)
    df['prev_dc'] = df['prev_dc'].fillna(initial_dc)

    df['month'] = pd.to_datetime(df['date']).dt.month

    # Monthly day length adjustment factor L_f
    L_f_table = {
        1: 6.5, 2: 7.5, 3: 9.0, 4: 12.8, 5: 14.8, 6: 16.0,
        7: 16.0, 8: 14.0, 9: 12.0, 10: 10.5, 11: 9.0, 12: 7.0
    }
    df['L_f'] = df['month'].map(L_f_table)

    rf = df['precipitation']
    re = np.where(rf > 2.8, rf - 2.8, 0.0)

    Q = 800 * np.exp(-df['prev_dc'] / 400.0)
    Q_r = Q + 3.937 * re
    dc_rain = 400 * np.log(800.0 / Q_r)

    # Drying rate V
    temp = df['temperature']
    V = np.where(temp > -2.8,
                 0.36 * (temp + 2.8) + df['L_f'],
                 0.0)

    # Final DC
    dc_new = np.where(re > 0, dc_rain, df['prev_dc'] + V)
    df['DC'] = np.clip(dc_new, 0, None)

    # Drop helper columns
    df = df.drop(columns=['prev_dc', 'month', 'L_f'])

    return df

df = compute_dc(df, initial_dc=15.0)

print("Drought Code Calculated")

print(df.head())

df.to_csv(data_path, index=False)
