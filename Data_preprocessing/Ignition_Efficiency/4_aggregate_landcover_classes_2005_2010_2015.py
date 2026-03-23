#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pathlib import Path
import subprocess
import sqlite3
import tempfile
import os

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
landcover_root = config.get("landcover_root")

#build file paths
from pathlib import Path

# -----------------------------
# User settings
# -----------------------------
ERA5_GRID_CSV = Path(data_root) / "era5_grid.csv"

LANDCOVER_FILES = {
    2005: Path(landcover_root) / "alaska_landcover_2005_250m.csv",
    2010: Path(landcover_root) / "alaska_landcover_2010_250m.csv",
    2015: Path(landcover_root) / "alaska_landcover_2015_30m.gpkg",  # huge
    # 2020 missing -> omit
}

GPKG_LAYER_2015 = "landcover2015"
OUTDIR = Path(data_root) / "lc_on_era5_grid"
LAT_BAND_DEG = 2.0

# -----------------------------
# NALCMS mapping
# -----------------------------
CONIFER_CLASSES = {1, 2}
BROADLEAF_CLASSES = {3, 4, 5}
MIXED_CLASS = 6
OPEN_CLASSES = {7, 8, 9, 10, 11, 12, 13, 14, 19}

MASK_CLASSES = {0, 15, 16, 17, 18}  # your request + background + water

# -----------------------------
# Helpers
# -----------------------------
def to_360(lon_series: pd.Series) -> pd.Series:
    lon = lon_series.astype(float).to_numpy()
    lon = np.where(lon < 0, lon + 360.0, lon)
    lon = np.where(lon >= 360.0, lon - 360.0, lon)
    return pd.Series(lon, index=lon_series.index)

def infer_resolution(values: np.ndarray) -> float:
    u = np.unique(np.round(values.astype(float), 10))
    u.sort()
    diffs = np.diff(u)
    diffs = diffs[diffs > 1e-10]
    return float(np.median(diffs))

def assign_to_grid_centers(coord: np.ndarray, centers: np.ndarray, res: float) -> np.ndarray:
    c0 = float(np.min(centers))
    idx = np.round((coord - c0) / res).astype(int)
    return c0 + idx * res

def lat_band(lat: pd.Series, band_deg: float) -> pd.Series:
    return (np.floor(lat.astype(float) / band_deg) * band_deg).astype(float)

def compute_pure_ratio_from_counts(n_con: float, n_brd: float) -> float:
    denom = n_con + n_brd
    return np.nan if denom <= 0 else (n_con / denom)

# -----------------------------
# CSV processing (small enough)
# -----------------------------
def process_csv_year(year: int, lc_csv: str, era5_grid: pd.DataFrame) -> pd.DataFrame:
    lc = pd.read_csv(lc_csv).rename(columns={"landcover": "lc"})
    lc = lc[["lat", "lon", "lc"]].copy()
    lc["lat"] = lc["lat"].astype(float)
    lc["lon"] = lc["lon"].astype(float)
    lc["lc"] = lc["lc"].astype(int)

    lc["lon360"] = to_360(lc["lon"])
    lc = lc[~lc["lc"].isin(MASK_CLASSES)].copy()

    era5_lons = era5_grid["lon360"].values
    era5_lats = era5_grid["lat"].values
    res_lon = infer_resolution(era5_lons)
    res_lat = infer_resolution(era5_lats)

    lc["era5_lon"] = assign_to_grid_centers(lc["lon360"].values, era5_lons, res_lon)
    lc["era5_lat"] = assign_to_grid_centers(lc["lat"].values, era5_lats, res_lat)

    counts = (
        lc.groupby(["era5_lon", "era5_lat"])["lc"]
          .value_counts()
          .unstack(fill_value=0)
          .reset_index()
          .rename(columns={"era5_lon": "lon360", "era5_lat": "lat"})
    )
    return counts

# -----------------------------
# GPKG processing (SQL aggregation)
# -----------------------------
def aggregate_gpkg_to_counts(gpkg_path: str, layer: str, res_lon: float, res_lat: float) -> pd.DataFrame:
    """
    Uses ogr2ogr to create a temporary sqlite with an aggregated table:
      lon360_center, lat_center, landcover, n
    then pivots to wide counts in pandas.
    """
    gpkg_path = str(gpkg_path)

    # Temporary sqlite file
    tmpdir = tempfile.mkdtemp(prefix="lc_agg_")
    sqlite_path = os.path.join(tmpdir, "agg.sqlite")

    # Expression to map lon to 0..360 BEFORE binning
    # - lon in [-180,180]
    # - lon360 = CASE WHEN lon < 0 THEN lon+360 ELSE lon END
    # Then bin to nearest ERA5 center:
    # center = round((lon360 - lon0)/res)*res + lon0
    # But we don't know lon0 in SQL; simplest: bin to res grid using:
    #   center = round(lon360 / res) * res
    # This matches ERA5 centers if your ERA5 grid is on exact 0.25 increments from 0.
    #
    # For ERA5 0.25°, centers are ... 0.00, 0.25, 0.50 ... which matches round(lon360/res)*res.
    #
    # Same for lat.

    sql = f"""
    SELECT
      ROUND((CASE WHEN ST_X(geom) < 0 THEN ST_X(geom)+360.0 ELSE ST_X(geom) END) / {res_lon}) * {res_lon} AS lon360,
      ROUND(ST_Y(geom) / {res_lat}) * {res_lat} AS lat,
      landcover AS lc,
      COUNT(*) AS n
    FROM "{layer}"
    WHERE landcover NOT IN ({",".join(str(x) for x in sorted(MASK_CLASSES))})
    GROUP BY lon360, lat, lc
    """

    # Run ogr2ogr to execute the SQL and write results into sqlite table "agg"
    cmd = [
        "ogr2ogr",
        "-f", "SQLite", sqlite_path,
        gpkg_path,
        "-sql", sql,
        "-nln", "agg",
        "-dialect", "SQLITE"
    ]
    print("Running:", " ".join(cmd[:6]) + " ... (sql omitted)")
    subprocess.run(cmd, check=True)

    # Read back the small aggregated table
    con = sqlite3.connect(sqlite_path)
    agg = pd.read_sql_query("SELECT lon360, lat, lc, n FROM agg", con)
    con.close()

    # Pivot to wide counts per ERA5 cell
    counts = (
        agg.pivot_table(index=["lon360", "lat"], columns="lc", values="n", fill_value=0, aggfunc="sum")
           .reset_index()
    )
    return counts

# -----------------------------
# Convert class counts -> C,B,U with mixed splitting
# -----------------------------
def counts_to_CBU(counts_wide: pd.DataFrame, era5_grid: pd.DataFrame) -> pd.DataFrame:
    # Ensure all expected columns exist
    for cls in range(0, 20):
        if cls not in counts_wide.columns:
            counts_wide[cls] = 0

    # Pure counts
    counts_wide["n_con"] = counts_wide[1] + counts_wide[2]
    counts_wide["n_brd"] = counts_wide[3] + counts_wide[4] + counts_wide[5]
    counts_wide["n_mix"] = counts_wide[6]
    counts_wide["n_open"] = sum(counts_wide[c] for c in OPEN_CLASSES)

    # Within-cell ratio
    counts_wide["cell_ratio"] = counts_wide.apply(
        lambda r: compute_pure_ratio_from_counts(r["n_con"], r["n_brd"]),
        axis=1
    )

    # Latitude band ratio (regional)
    counts_wide["lat_band"] = lat_band(counts_wide["lat"], LAT_BAND_DEG)
    band = (
        counts_wide.groupby("lat_band", as_index=False)[["n_con","n_brd"]]
        .sum()
    )
    band["band_ratio"] = band.apply(lambda r: compute_pure_ratio_from_counts(r["n_con"], r["n_brd"]), axis=1)
    counts_wide = counts_wide.merge(band[["lat_band","band_ratio"]], on="lat_band", how="left")

    # Global ratio fallback
    global_ratio = compute_pure_ratio_from_counts(
        counts_wide["n_con"].sum(),
        counts_wide["n_brd"].sum()
    )

    mix_ratio = counts_wide["cell_ratio"].copy()
    mix_ratio = mix_ratio.fillna(counts_wide["band_ratio"])
    mix_ratio = mix_ratio.fillna(global_ratio)
    mix_ratio = mix_ratio.fillna(0.5)

    # Allocate mixed
    counts_wide["C_eff"] = counts_wide["n_con"] + counts_wide["n_mix"] * mix_ratio
    counts_wide["B_eff"] = counts_wide["n_brd"] + counts_wide["n_mix"] * (1.0 - mix_ratio)
    counts_wide["U_eff"] = counts_wide["n_open"]

    denom = counts_wide["C_eff"] + counts_wide["B_eff"] + counts_wide["U_eff"]
    counts_wide["C"] = np.where(denom > 0, counts_wide["C_eff"] / denom, np.nan)
    counts_wide["B"] = np.where(denom > 0, counts_wide["B_eff"] / denom, np.nan)
    counts_wide["U"] = np.where(denom > 0, counts_wide["U_eff"] / denom, np.nan)

    out = counts_wide[["lon360", "lat", "C", "B", "U"]].copy()

    # Ensure full ERA5 grid coverage
    out = era5_grid.merge(out, on=["lon360", "lat"], how="left")

    # Output lon in -180..180
    out["lon"] = np.where(out["lon360"] > 180.0, out["lon360"] - 360.0, out["lon360"])
    out = out[["lon", "lat", "C", "B", "U"]].sort_values(["lat","lon"]).reset_index(drop=True)
    return out

# -----------------------------
# Main
# -----------------------------
def main():
    outdir = Path(OUTDIR)
    outdir.mkdir(parents=True, exist_ok=True)

    # ERA5 grid
    era5 = pd.read_csv(ERA5_GRID_CSV)
    if "X" in era5.columns and "Y" in era5.columns:
        era5 = era5.rename(columns={"X":"lon", "Y":"lat"})
    era5["lon"] = era5["lon"].astype(float)
    era5["lat"] = era5["lat"].astype(float)
    era5["lon360"] = to_360(era5["lon"])
    era5 = era5[["lon360","lat"]].drop_duplicates().reset_index(drop=True)

    # infer ERA5 resolution from grid itself
    res_lon = infer_resolution(era5["lon360"].values)
    res_lat = infer_resolution(era5["lat"].values)
    print(f"Inferred ERA5 resolution: lon={res_lon}°, lat={res_lat}°")

    for year, path in LANDCOVER_FILES.items():
        if path is None or str(path).strip() == "" or not Path(path).exists():
            print(f"Skipping {year}: no file found")
            continue

        suffix = Path(path).suffix.lower()
        print(f"Processing {year}: {path}")

        if suffix == ".csv":
            counts = process_csv_year(year, path, era5)
        elif suffix == ".gpkg":
            counts = aggregate_gpkg_to_counts(path, GPKG_LAYER_2015, res_lon=res_lon, res_lat=res_lat)
        else:
            raise ValueError(f"Unsupported file type: {suffix}")

        out = counts_to_CBU(counts, era5)
        out_csv = outdir / f"landcover_ERA5_{year}.csv"
        out.to_csv(out_csv, index=False)
        print(f"  wrote {out_csv}")

if __name__ == "__main__":
    main()
