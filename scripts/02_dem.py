"""Step 2. Copernicus GLO-30 digital elevation model (30 m), reprojected to UTM.

Source: Earth Search, collection cop-dem-glo-30 (no account).

What a DEM is: a raster where each pixel stores the terrain height above sea level, in metres.
Here each pixel is 30 x 30 m. It is used to draw the relief and, in Step 4, to work out where
water flows, delineate the basin and create elevation bands.

Note: the Copernicus DEM is a surface model (it includes tree crowns), not bare ground.

In short:
  1. Find the elevation "tiles" that cover your area in the catalogue.
  2. Download only your clip and reproject it to UTM 10N (metres), like Sentinel-2.
  3. Save it as a GeoTIFF to open in QGIS.

Usage:  python scripts/02_dem.py
"""
import common
from odc.stac import configure_rio, load
from pystac_client import Client

STAC_URL = "https://earth-search.aws.element84.com/v1"


def main():
    # The DEM tiles sit in a public AWS S3 bucket (s3://...). Without this line GDAL tries to
    # sign the request with AWS credentials and fails with "InvalidCredentials".
    # "aws_unsigned" tells it to request the file unsigned, like any anonymous visitor.
    configure_rio(cloud_defaults=True, aws={"aws_unsigned": True})

    cfg = common.load_config()
    bbox = cfg["area"]["bbox"]
    epsg = f"EPSG:{cfg['area']['epsg_utm']}"

    # --- 1. Search. No dates or cloud filter: relief does not change over time ---
    catalogue = Client.open(STAC_URL)
    items = list(catalogue.search(collections=["cop-dem-glo-30"], bbox=bbox).items())
    if not items:
        raise SystemExit("No DEM found for the area.")
    # Each tile has a file ("asset") called "data" holding the elevations.
    asset = "data" if "data" in items[0].assets else list(items[0].assets)[0]

    # --- 2. Load, clip to your bbox and reproject to metres (UTM 10N) ---
    # 30 m resolution = the DEM's native one. "bilinear" smooths while reprojecting.
    cube = load(items, bands=[asset], bbox=bbox, crs=epsg, resolution=30,
                resampling="bilinear", chunks={})
    dem = cube[asset]
    # If the area falls on several tiles, keep the highest value where they overlap.
    if "time" in dem.dims:
        dem = dem.max("time", skipna=True)
    dem = dem.astype("float32").compute()   # the data is downloaded here

    # --- 3. Save to data/processed/dem_utm.tif ---
    dest = common.project_path("data", "processed", "dem_utm.tif")
    common.write_cog(dem, dest, epsg)
    print(f"DEM saved to {dest}")
    print(f"  Minimum elevation {float(dem.min()):.0f} m, maximum {float(dem.max()):.0f} m")


if __name__ == "__main__":
    main()
