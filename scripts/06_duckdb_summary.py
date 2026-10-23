"""Step 6. Join all layers into one GeoParquet table and query it with DuckDB.

Each row is a 20 m pixel inside the basin, with NDVI, NDMI, elevation, canopy height
(only where there is LiDAR) and elevation band.

Why: having everything in one table lets you ask questions in SQL, for example
"how does NDVI change as the band goes up?", instead of crossing rasters by hand.

In short:
  1. build_table: bring everything onto the 20 m NDVI grid and build one row per pixel.
  2. Save the table as GeoParquet (pixels.parquet): a compact format that also stores the geometry.
  3. DuckDB (a SQL engine that reads the file directly, no server) computes three summaries:
     per band, canopy height per band (only where there is LiDAR) and correlation with elevation.
  4. Save the result to outputs/band_summary.csv and docs/data/band_summary.json (Step 7).

Things to keep in mind when interpreting:
  - LiDAR only covers the bbox_lidar area: outside it there is no canopy height.
  - Bands with very few pixels give unreliable means.
  - elev comes from the Copernicus DEM, which measures the surface (it includes crowns), not the ground.
  - The LiDAR is from 2021 and Sentinel-2 from 2025: anything cut between the two dates still shows tall canopy.

Usage:  python scripts/06_duckdb_summary.py
"""
import json

import common
import duckdb
import geopandas as gpd
import numpy as np
import rasterio
import rioxarray
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.features import rasterize


def _open(name):
    """Open a raster from data/processed/ as a 2D array (without the band axis)."""
    p = common.project_path("data", "processed", name)
    if not p.exists():
        raise SystemExit(f"{p} is missing. Check the previous steps.")
    return rioxarray.open_rasterio(p).squeeze("band", drop=True)


def build_table(cfg):
    ndvi = _open("ndvi_summer.tif")
    ndmi = _open("ndmi_summer.tif")
    # The reference grid is the NDVI's (20 m pixels). The DEM (30 m) and the CHM (finer) are
    # resampled onto it so every row has one value from each layer at the same place.
    # CHM: averaged (average) because it goes from small cells to a much larger one.
    dem = _open("dem_utm.tif").rio.reproject_match(ndvi, resampling=Resampling.bilinear)
    chm_path = common.project_path("data", "processed", "chm.tif")
    chm = (_open("chm.tif").rio.reproject_match(ndvi, resampling=Resampling.average)
           if chm_path.exists() else None)

    basins = gpd.read_file(common.project_path("data", "processed", "basins.gpkg"))
    bands = gpd.read_file(common.project_path("data", "processed", "elevation_bands.gpkg"))
    transform = ndvi.rio.transform()
    shape = ndvi.shape

    # "inside" marks with True the pixels that fall inside the basin(s).
    inside = rasterize(basins.geometry, out_shape=shape, transform=transform,
                       fill=0, default_value=1, dtype="uint8").astype(bool)
    # band_id holds each pixel's band number (-1 = outside every band).
    # 1 is added when rasterising because 0 means "empty"; it is subtracted afterwards.
    band_id = rasterize(zip(bands.geometry, bands["band"].astype(int) + 1), out_shape=shape,
                        transform=transform, fill=0, dtype="int32") - 1

    # One row per valid pixel: centre coordinates and the value of each layer.
    rows, cols = np.where(inside & (band_id >= 0))
    x, y = rasterio.transform.xy(transform, rows, cols)
    labels = dict(zip(bands["band"].astype(int), bands["label"]))
    table = gpd.GeoDataFrame({
        "x": np.asarray(x), "y": np.asarray(y),
        "elev": dem.values[rows, cols],
        "ndvi": ndvi.values[rows, cols],
        "ndmi": ndmi.values[rows, cols],
        "chm": chm.values[rows, cols] if chm is not None else np.nan,
        "band": band_id[rows, cols],
    }, geometry=gpd.points_from_xy(np.asarray(x), np.asarray(y)), crs=ndvi.rio.crs)
    table["label"] = table["band"].map(labels)
    return table


def main():
    cfg = common.load_config()
    table = build_table(cfg)
    dest = common.project_path("data", "processed", "pixels.parquet")
    table.to_parquet(dest)
    print(f"GeoParquet: {len(table):,} pixels in {dest}")

    # LiDAR-area box in UTM coordinates, to filter with a spatial query.
    lon0, lat0, lon1, lat1 = cfg["area"]["bbox_lidar"]
    to_utm = Transformer.from_crs(4326, cfg["area"]["epsg_utm"], always_xy=True)
    x0, y0 = to_utm.transform(lon0, lat0)
    x1, y1 = to_utm.transform(lon1, lat1)

    # --- DuckDB: in-memory connection; the table is read straight from the .parquet file ---
    con = duckdb.connect()
    spatial = False
    try:
        con.sql("INSTALL spatial; LOAD spatial;")
        spatial = True
    except Exception as e:  # no internet to fetch the extension
        print(f"  (DuckDB spatial extension unavailable: {e}. Falling back to a coordinate filter.)")
    con.sql(f"CREATE VIEW pix AS SELECT * FROM read_parquet('{dest.as_posix()}')")

    # Same LiDAR-box filter, with the spatial extension or with plain comparisons.
    where_lidar = (f"ST_Intersects(ST_Point(x, y), ST_MakeEnvelope({x0}, {y0}, {x1}, {y1}))"
                   if spatial else f"x BETWEEN {x0} AND {x1} AND y BETWEEN {y0} AND {y1}")

    # Query 1: summary per band (all pixels in the basin).
    summary = con.sql("""
        SELECT band, any_value(label) AS label, count(*) AS pixels,
               round(avg(elev), 0) AS elev_mean,
               round(avg(ndvi), 3) AS ndvi_mean,
               round(avg(ndmi), 3) AS ndmi_mean
        FROM pix GROUP BY band ORDER BY band
    """).df()
    # Query 2: canopy height per band, only with pixels inside the LiDAR box.
    # p95 = height below which 95 % of pixels lie (the tall crown, without extremes).
    canopy = con.sql(f"""
        SELECT band, count(chm) AS chm_pixels,
               round(avg(chm::DOUBLE), 1) AS chm_mean,
               round(quantile_cont(chm::DOUBLE, 0.95), 1) AS chm_p95
        FROM pix WHERE {where_lidar} AND chm IS NOT NULL
        GROUP BY band ORDER BY band
    """).df()
    # Query 3: Pearson correlation (-1 to 1) between elevation and each index.
    correlation = con.sql("""
        SELECT round(corr(elev, ndvi), 3) AS elev_vs_ndvi, round(corr(elev, ndmi), 3) AS elev_vs_ndmi
        FROM pix
    """).df()

    # Bands without LiDAR end up empty in the canopy columns (how="left").
    final = summary.merge(canopy, on="band", how="left")
    print("\nSummary per elevation band:")
    print(final.to_string(index=False))
    print("\nCorrelation with elevation:")
    print(correlation.to_string(index=False))

    out_csv = common.project_path("outputs", "band_summary.csv")
    final.to_csv(out_csv, index=False)
    out_json = common.project_path("docs", "data", "band_summary.json")
    out_json.write_text(json.dumps(json.loads(final.to_json(orient="records")), ensure_ascii=False, indent=2),
                        encoding="utf-8")
    print(f"\nSaved to {out_csv} and {out_json}")


if __name__ == "__main__":
    main()
