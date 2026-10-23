"""Step 1. Sentinel-2 L2A: monthly NDVI and NDMI (Zarr) and a summer mean (COG).

Source: the public Earth Search STAC catalogue (no account). The images are COGs on AWS;
only the pixels of your area are read, not the whole scenes.

In short:
  1. Ask the catalogue which satellite scenes cover your area in 2025 (with few clouds).
  2. Download only 4 bands of those scenes, clipped to your area.
  3. Compute two indices per date: NDVI (greenness) and NDMI (leaf water).
  4. Take the median of each month and the July-September mean.
  5. Save the results and check one pixel by hand to confirm the calculation is right.

Usage:  python scripts/01_sentinel2.py
"""
import common
from odc.stac import load
from pystac_client import Client

# STAC is a standard for searching satellite imagery. This is the catalogue's address.
STAC_URL = "https://earth-search.aws.element84.com/v1"


def main():
    # --- Read parameters from config.yml (area, dates, clouds, resolution) ---
    cfg = common.load_config()
    bbox = cfg["area"]["bbox"]                    # [lon_min, lat_min, lon_max, lat_max]
    epsg = f"EPSG:{cfg['area']['epsg_utm']}"      # a metric system (UTM 10N) to work in
    s2 = cfg["sentinel2"]
    period = f"{cfg['dates']['start']}/{cfg['dates']['end']}"

    # --- 1. Search: only asks for the list of scenes, nothing is downloaded yet ---
    print("Searching Earth Search for scenes...")
    catalogue = Client.open(STAC_URL)
    items = list(catalogue.search(
        collections=["sentinel-2-l2a"], bbox=bbox, datetime=period,
        query={"eo:cloud_cover": {"lt": s2["max_cloud"]}},   # discard very cloudy scenes
    ).items())
    print(f"  {len(items)} scenes with less than {s2['max_cloud']} % cloud")
    if not items:
        raise SystemExit("No scenes found. Check bbox and dates, or raise max_cloud in config.yml.")

    # --- 2. Load: build a data "cube" with time, y, x axes ---
    # Bands: red, nir (near infrared), swir16 (shortwave infrared, sensitive to water) and
    # scl (a per-pixel classification: cloud, shadow, vegetation...).
    # Lazy: pixels are read from the internet only when computed (.compute()).
    print("Loading bands (this can take a few minutes)...")
    cube = load(
        items, bands=["red", "nir", "swir16", "scl"], bbox=bbox, crs=epsg,
        resolution=s2["resolution"], groupby="solar_day",    # one layer per day
        # scl holds categories (not averaged): "nearest". The rest: "bilinear".
        resampling={"scl": "nearest", "*": "bilinear"},
        chunks={"x": 512, "y": 512},                         # chunks, to keep memory in check
    )

    # --- 3. Indices per date (see common.py: converts to reflectance and drops clouds) ---
    indices = common.spectral_indices(cube)
    # Monthly median: summarises several dates into one and reduces the effect of bad pixels.
    # This is when the data is actually downloaded.
    monthly = common.aggregate_monthly(indices).compute()

    # --- 4. Quality control: recompute the NDVI of the central pixel by hand ---
    # If the program's value and the manual one agree, the scale conversion is correct.
    lon_c = (bbox[0] + bbox[2]) / 2
    lat_c = (bbox[1] + bbox[3]) / 2
    print("Check with a known pixel:")
    common.check_pixel(cube, indices, lon_c, lat_c, epsg)

    # --- 5. Save ---
    # Zarr: keeps the whole monthly series (one month per layer).
    monthly.to_zarr(common.project_path("data", "processed", "monthly_indices.zarr"), mode="w")
    # Mean of the summer months -> one map per index, as a GeoTIFF (COG) that QGIS opens.
    summer = common.season_mean(monthly, s2["summer_months"])
    common.write_cog(summer["ndvi"], common.project_path("data", "processed", "ndvi_summer.tif"), epsg)
    common.write_cog(summer["ndmi"], common.project_path("data", "processed", "ndmi_summer.tif"), epsg)

    print("\nDone. Files in data/processed/:")
    print("  monthly_indices.zarr  monthly NDVI and NDMI series")
    print("  ndvi_summer.tif       mean summer NDVI")
    print("  ndmi_summer.tif       mean summer NDMI")
    print(f"  Summer NDVI: mean {float(summer['ndvi'].mean()):.3f}")
    print(f"  Summer NDMI: mean {float(summer['ndmi'].mean()):.3f}")


if __name__ == "__main__":
    main()
