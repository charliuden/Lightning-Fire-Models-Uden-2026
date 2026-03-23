#!/usr/bin/env python3
"""
Aggregate NALCMS 2020 categorical landcover GeoTIFF to ERA5 grid and output
C (conifer), B (broadleaf), U (open/unforested) proportions per ERA5 cell.

Designed for huge rasters: streaming, window-based, no full raster load.

Inputs:
  - era5_grid.csv (columns X,Y or lon,lat; centers)
  - alaska_landcover_2020_30m.tif (categorical uint8)
Output:
  - landcover_ERA5_2020.csv with columns lon, lat, C, B, U

Dependencies (conda-forge recommended):
  conda install -c conda-forge rasterio geopandas shapely pyproj rtree
"""

from __future__ import annotations

import os
import sys
import time
from collections import defaultdict
from pathlib import Path
import numpy as np
import pandas as pd

import rasterio
from rasterio.features import rasterize
from shapely.geometry import box
import geopandas as gpd

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
# USER SETTINGS
# -----------------------------
ERA5_GRID_CSV = Path(data_root) / "era5_grid.csv"
TIF_2020 = Path(landcover_root) / "alaska_landcover_2020_30m.tif"
OUT_CSV = Path(data_root) / "lc_on_era5_grid" / "landcover_ERA5_2020.csv"

# Mixed-forest regional splitting:
LAT_BAND_DEG = 2.0  # 2-degree bands is a good default

# Monitoring / performance knobs:
PRINT_EVERY_WINDOWS = 200  # print progress every N windows
MAX_WINDOWS = None         # set to an int for a short test run (e.g., 50)

# If you want to test quickly on a subset of the raster, set a pixel window here.
# Example: TEST_WINDOW = (row_off, col_off, height, width)
TEST_WINDOW = None  # e.g., (100000, 100000, 5000, 5000)


# -----------------------------
# NALCMS class mapping
# -----------------------------
CONIFER_CLASSES = {1, 2}
BROADLEAF_CLASSES = {3, 4, 5}
MIXED_CLASS = 6
OPEN_CLASSES = {7, 8, 9, 10, 11, 12, 13, 14, 19}

# Mask (your request) + also water/background
MASK_CLASSES = {0, 15, 16, 17, 18}
# Raster NoData for your file (from metadata)
RASTER_NODATA_FALLBACK = 127


# -----------------------------
# Helpers
# -----------------------------
def to_360(lon: np.ndarray) -> np.ndarray:
    lon = lon.astype(float).copy()
    lon[lon < 0] += 360.0
    lon[lon >= 360.0] -= 360.0
    return lon

def infer_resolution(values: np.ndarray) -> float:
    u = np.unique(np.round(values.astype(float), 10))
    u.sort()
    if len(u) < 2:
        raise ValueError("Cannot infer resolution from grid: not enough unique coords.")
    diffs = np.diff(u)
    diffs = diffs[diffs > 1e-10]
    if len(diffs) == 0:
        raise ValueError("Cannot infer resolution from grid: diffs all ~0.")
    return float(np.median(diffs))

def lat_band(lat: np.ndarray, band_deg: float) -> np.ndarray:
    return (np.floor(lat.astype(float) / band_deg) * band_deg).astype(float)

def compute_ratio(n_con: float, n_brd: float) -> float:
    denom = n_con + n_brd
    return np.nan if denom <= 0 else (n_con / denom)

def pairs_to_counts(cid: np.ndarray, cls: np.ndarray) -> dict[int, int]:
    """
    Efficient joint histogram. Hash pair as cid*256 + cls (cls is uint8).
    Returns dict: hashed_pair -> count
    """
    pairs = cid.astype(np.int64) * 256 + cls.astype(np.int64)
    vals, cnts = np.unique(pairs, return_counts=True)
    return {int(v): int(c) for v, c in zip(vals, cnts)}


def build_era5_cells(era5: pd.DataFrame, res_lon: float, res_lat: float) -> gpd.GeoDataFrame:
    """
    Build ERA5 grid-cell polygons from grid centers.
    We build in lon=-180..180 coordinates for dateline sanity, then project to raster CRS later.
    """
    half_lon = res_lon / 2.0
    half_lat = res_lat / 2.0

    # Ensure lon360 + lon(-180..180)
    lon360 = era5["lon360"].to_numpy()
    lon = np.where(lon360 > 180.0, lon360 - 360.0, lon360)
    lat = era5["lat"].to_numpy()

    polys = [box(x - half_lon, y - half_lat, x + half_lon, y + half_lat) for x, y in zip(lon, lat)]
    gdf = gpd.GeoDataFrame(
        {"cell_id": np.arange(len(era5), dtype=np.int64), "lon360": lon360, "lat": lat},
        geometry=polys,
        crs="EPSG:4326",
    )
    return gdf


# -----------------------------
# Core aggregation
# -----------------------------
def aggregate_tif_to_class_counts(
    tif_path: str,
    era5_grid: pd.DataFrame,
    res_lon: float,
    res_lat: float,
) -> pd.DataFrame:
    """
    Stream through raster windows; rasterize ERA5 cell IDs per window; accumulate counts.
    Returns a wide DataFrame with columns lon360, lat, and class-code columns containing counts.
    """
    tif_path = str(tif_path)

    era5_cells_ll = build_era5_cells(era5_grid, res_lon=res_lon, res_lat=res_lat)

    with rasterio.open(tif_path) as src:
        if src.crs is None:
            raise ValueError("Raster CRS is missing; cannot align to ERA5 grid.")
        raster_crs = src.crs

        # Project ERA5 polygons into raster CRS (one-time cost)
        era5_cells = era5_cells_ll.to_crs(raster_crs)

        # Spatial index for quick window intersection
        sindex = era5_cells.sindex

        nodata = src.nodata
        if nodata is None:
            nodata = RASTER_NODATA_FALLBACK

        # Accumulator: hashed_pair -> count
        # where hashed_pair = cell_id*256 + class
        joint = defaultdict(int)

        # Choose window iterator
        if TEST_WINDOW is not None:
            row_off, col_off, height, width = TEST_WINDOW
            windows = [(((0, 0)), rasterio.windows.Window(col_off, row_off, width, height))]
        else:
            windows = list(src.block_windows(1))  # may be many; still ok because it's just metadata

        t0 = time.time()
        for i, (ji, window) in enumerate(windows, start=1):
            if MAX_WINDOWS is not None and i > MAX_WINDOWS:
                break

            # Read window
            arr = src.read(1, window=window)

            # Valid landcover mask (not nodata and not masked classes)
            valid = (arr != nodata)
            if not np.any(valid):
                continue

            # Also remove mask classes from valid (so they don't count toward denom)
            # (arr is uint8, MASK_CLASSES ints)
            if MASK_CLASSES:
                for m in MASK_CLASSES:
                    valid &= (arr != m)
                if not np.any(valid):
                    continue

            # Window bounds in raster CRS
            w_bounds = rasterio.windows.bounds(window, src.transform)
            w_box = box(*w_bounds)

            # Find intersecting ERA5 polygons
            cand_idx = list(sindex.intersection(w_box.bounds))
            if len(cand_idx) == 0:
                continue
            sub = era5_cells.iloc[cand_idx]
            sub = sub[sub.intersects(w_box)]
            if sub.empty:
                continue

            # Rasterize cell_id into this window
            shapes = [(geom, int(cid)) for geom, cid in zip(sub.geometry, sub["cell_id"].to_numpy())]
            cell_id_raster = rasterize(
                shapes=shapes,
                out_shape=(window.height, window.width),
                transform=rasterio.windows.transform(window, src.transform),
                fill=-1,
                dtype="int32",
            )

            inside = cell_id_raster >= 0
            use = valid & inside
            if not np.any(use):
                continue

            cls = arr[use].astype(np.uint8)
            cid = cell_id_raster[use].astype(np.int64)

            local = pairs_to_counts(cid, cls)
            for k, v in local.items():
                joint[k] += v

            if (i % PRINT_EVERY_WINDOWS) == 0:
                elapsed = time.time() - t0
                print(f"[{i}/{len(windows)}] windows processed | elapsed {elapsed/60:.1f} min | joint keys {len(joint):,}")

        # Unpack joint counts into rows (cell_id, class, n)
        rows = []
        for key, n in joint.items():
            cell_id = key // 256
            lc_class = key % 256
            rows.append((cell_id, lc_class, n))

    if len(rows) == 0:
        # No overlap or all masked/nodata
        return pd.DataFrame(columns=["lon360", "lat"])

    df = pd.DataFrame(rows, columns=["cell_id", "lc", "n"])

    # Pivot to wide counts per cell_id
    wide = (
        df.pivot_table(index="cell_id", columns="lc", values="n", fill_value=0, aggfunc="sum")
          .reset_index()
    )

    # Attach lon360/lat
    wide = wide.merge(
        era5_cells_ll[["cell_id", "lon360", "lat"]],
        on="cell_id",
        how="left"
    ).drop(columns=["cell_id"])

    # Make sure lon360/lat are first
    class_cols = [c for c in wide.columns if c not in ("lon360", "lat")]
    wide = wide[["lon360", "lat"] + class_cols]
    return wide


def counts_to_CBU(counts_wide: pd.DataFrame, era5_grid: pd.DataFrame) -> pd.DataFrame:
    """
    Convert wide class counts to C/B/U proportions with regional mixed split.
    """
    if counts_wide.empty:
        out = era5_grid.copy()
        out["lon"] = np.where(out["lon360"] > 180.0, out["lon360"] - 360.0, out["lon360"])
        out[["C", "B", "U"]] = np.nan
        return out[["lon", "lat", "C", "B", "U"]]

    # Ensure expected class columns exist (1..19, plus 6, etc.)
    for cls in range(0, 20):
        if cls not in counts_wide.columns:
            counts_wide[cls] = 0

    # Pure counts
    counts_wide["n_con"] = counts_wide[1] + counts_wide[2]
    counts_wide["n_brd"] = counts_wide[3] + counts_wide[4] + counts_wide[5]
    counts_wide["n_mix"] = counts_wide[6]
    counts_wide["n_open"] = sum(counts_wide[c] for c in OPEN_CLASSES)

    # Within-cell ratio
    counts_wide["cell_ratio"] = [
        compute_ratio(c, b) for c, b in zip(counts_wide["n_con"].to_numpy(), counts_wide["n_brd"].to_numpy())
    ]

    # Latitude band ratio
    counts_wide["lat_band"] = lat_band(counts_wide["lat"].to_numpy(), LAT_BAND_DEG)
    band = counts_wide.groupby("lat_band", as_index=False)[["n_con", "n_brd"]].sum()
    band["band_ratio"] = [compute_ratio(c, b) for c, b in zip(band["n_con"].to_numpy(), band["n_brd"].to_numpy())]
    counts_wide = counts_wide.merge(band[["lat_band", "band_ratio"]], on="lat_band", how="left")

    # Global ratio
    global_ratio = compute_ratio(float(counts_wide["n_con"].sum()), float(counts_wide["n_brd"].sum()))

    # Choose mix ratio: cell -> band -> global -> 0.5
    mix_ratio = pd.Series(counts_wide["cell_ratio"]).copy()
    mix_ratio = mix_ratio.fillna(counts_wide["band_ratio"])
    mix_ratio = mix_ratio.fillna(global_ratio)
    mix_ratio = mix_ratio.fillna(0.5)

    # Allocate mixed
    C_eff = counts_wide["n_con"] + counts_wide["n_mix"] * mix_ratio
    B_eff = counts_wide["n_brd"] + counts_wide["n_mix"] * (1.0 - mix_ratio)
    U_eff = counts_wide["n_open"]

    denom = C_eff + B_eff + U_eff
    counts_wide["C"] = np.where(denom > 0, C_eff / denom, np.nan)
    counts_wide["B"] = np.where(denom > 0, B_eff / denom, np.nan)
    counts_wide["U"] = np.where(denom > 0, U_eff / denom, np.nan)

    out = counts_wide[["lon360", "lat", "C", "B", "U"]].copy()
    out = era5_grid.merge(out, on=["lon360", "lat"], how="left")

    out["lon"] = np.where(out["lon360"] > 180.0, out["lon360"] - 360.0, out["lon360"])
    out = out[["lon", "lat", "C", "B", "U"]].sort_values(["lat", "lon"]).reset_index(drop=True)
    return out


def main() -> int:
    tif = Path(TIF_2020)
    if not tif.exists():
        print(f"ERROR: TIFF not found: {tif}")
        return 2

    # Load ERA5 grid
    era5 = pd.read_csv(ERA5_GRID_CSV)
    if "X" in era5.columns and "Y" in era5.columns:
        era5 = era5.rename(columns={"X": "lon", "Y": "lat"})
    if not {"lon", "lat"}.issubset(era5.columns):
        print(f"ERROR: ERA5 grid must have lon/lat or X/Y columns. Got: {era5.columns.tolist()}")
        return 2

    era5["lon"] = era5["lon"].astype(float)
    era5["lat"] = era5["lat"].astype(float)
    era5["lon360"] = to_360(era5["lon"].to_numpy())
    era5 = era5[["lon360", "lat"]].drop_duplicates().reset_index(drop=True)

    # Infer ERA5 resolution
    res_lon = infer_resolution(era5["lon360"].to_numpy())
    res_lat = infer_resolution(era5["lat"].to_numpy())
    print(f"ERA5 resolution inferred: lon={res_lon}°, lat={res_lat}°")
    print(f"Raster: {tif}")

    # Quick sanity check: sample values (small center window)
    with rasterio.open(tif) as src:
        nodata = src.nodata if src.nodata is not None else RASTER_NODATA_FALLBACK
        r0 = src.height // 2
        c0 = src.width // 2
        h = min(2000, src.height)
        w = min(2000, src.width)
        arr = src.read(1, window=((r0, min(r0 + h, src.height)), (c0, min(c0 + w, src.width))))
        vals = np.unique(arr)
        print(f"Sanity subset unique values (first 30): {vals[:30]}  (n={len(vals)}) | nodata={nodata}")
        if len(vals) == 1 and vals[0] == nodata:
            print("WARNING: subset is all NoData. The raster might still be valid elsewhere, but verify clipping/extent.")

    t0 = time.time()
    counts = aggregate_tif_to_class_counts(str(tif), era5, res_lon=res_lon, res_lat=res_lat)
    print(f"Aggregation done. Counts rows: {len(counts):,}. Elapsed: {(time.time()-t0)/60:.1f} min")

    out = counts_to_CBU(counts, era5)
    out.to_csv(OUT_CSV, index=False)
    print(f"Wrote: {OUT_CSV}  (rows={len(out):,})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
