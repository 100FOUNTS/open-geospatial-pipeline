"""Step 7. Prepare the data for the web report: Web Mercator PNGs, GeoJSON and JSON.

A browser cannot open a GeoTIFF or a GeoPackage, so this script converts them:
  - Each raster (NDVI, NDMI, CHM) becomes a PNG with a colour ramp; no-data pixels are
    transparent. Their corners are recorded in lat/lon so Leaflet can place them, along with
    the ramp's colours to draw the legend.
  - The basins, drainage network and elevation bands become GeoJSON.
  - The DEM contour lines become an SVG (the background of the page header).
  - The interpolation metrics (CSV) become JSON and the variogram is copied to docs/img.
  - layers.json is the list of layers the page reads.

The page (docs/index.html) uses Leaflet for the map and hand-made SVG for the charts. View it
with a local server (not by double-clicking: the browser blocks reading the JSON files when you
open the file directly):
  python -m http.server 8000 --directory docs

Usage:  python scripts/07_web_report.py
"""
import json
import shutil

import common
import geopandas as gpd
import matplotlib
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject, transform_bounds
from scipy.ndimage import gaussian_filter, zoom

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LAYERS = [
    # Own basemap: hillshade of the DEM (does not depend on external services).
    ("relief", "Relief", "hillshade.tif", "gray", 0, 255, "", None),
    # id, name, file, ramp (from matplotlib), minimum and maximum painted, unit,
    # cell size in Web Mercator metres for the image (None = the raster's own)
    ("ndvi", "Summer NDVI", "ndvi_summer.tif", "RdYlGn", 0.2, 0.95, "", None),
    ("ndmi", "Summer NDMI", "ndmi_summer.tif", "BrBG", -0.1, 0.6, "", None),
    # The CHM (3 m) is heavy as a PNG: it is averaged to ~8 m cells for the web.
    ("chm", "Canopy height", "chm.tif", "YlGn", 0, 60, "m", 8),
]


def raster_to_png(tif, png, ramp, vmin, vmax, resolution=None):
    """Reproject to EPSG:3857, apply a colour ramp and return the bounds in lat/lon.

    EPSG:3857 (Web Mercator) is the web-map system: the image must be in it to fit on the
    basemap. Values outside [vmin, vmax] are painted with the extreme colour.
    With 'resolution' (Web Mercator metres) the image is resampled by averaging to that cell: lighter.
    """
    with rasterio.open(tif) as src:
        transform, width, height = calculate_default_transform(
            src.crs, "EPSG:3857", src.width, src.height, *src.bounds, resolution=resolution)
        dest = np.full((height, width), np.nan, dtype="float32")
        reproject(source=rasterio.band(src, 1), destination=dest,
                  src_transform=src.transform, src_crs=src.crs,
                  src_nodata=src.nodata if src.nodata is not None else np.nan,
                  dst_transform=transform, dst_crs="EPSG:3857", dst_nodata=np.nan,
                  resampling=Resampling.average if resolution else Resampling.bilinear)
    left, top = transform.c, transform.f
    right, bottom = left + width * transform.a, top + height * transform.e
    w, s, e, n = transform_bounds("EPSG:3857", "EPSG:4326", left, bottom, right, top)
    # NaN (no data) stays transparent in the PNG.
    plt.imsave(png, dest, cmap=ramp, vmin=vmin, vmax=vmax)
    return [[s, w], [n, e]]


def hillshade(dem_path, dest, azimuth=315, altitude=45):
    """Hillshade of the DEM, as a 0-255 GeoTIFF.

    Same method as QGIS/GDAL 'Hillshade': light from the northwest (azimuth 315) at 45 degrees.
    The surface DEM is noisy (crowns), so it is smoothed a little before computing.
    """
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype("float64")
        profile = src.profile.copy()
        step = src.res[0]
    dem = gaussian_filter(np.where(np.isfinite(dem), dem, np.nanmean(dem)), 1)
    dz_row, dz_x = np.gradient(dem, step)
    dz_y = -dz_row                        # rows grow southwards; y grows northwards
    slope = np.arctan(np.hypot(dz_x, dz_y))
    aspect = np.arctan2(dz_y, -dz_x)
    zenith = np.radians(90 - altitude)
    az = np.radians(360 - azimuth + 90)
    hs = 255 * (np.cos(zenith) * np.cos(slope) + np.sin(zenith) * np.sin(slope) * np.cos(az - aspect))
    profile.update(dtype="float32", count=1, nodata=np.nan, compress="deflate")
    with rasterio.open(dest, "w", **profile) as dst:
        dst.write(np.clip(hs, 0, 255).astype("float32"), 1)


def ramp_colors(ramp, n=9):
    """n hexadecimal colours from the ramp, low to high, to draw the legend."""
    cm = plt.get_cmap(ramp)
    return [mcolors.to_hex(cm(i / (n - 1))) for i in range(n)]


def contours_svg(dem_path, dest, interval=40, every=5):
    """Write an SVG with the DEM's contour lines (the page header background).

    The DEM (30 m) is coarse; it is enlarged and smoothed so the lines look like a topographic
    plan. Every 'every' contours a thicker one is drawn (index contour).
    """
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype("float64")
    dem = np.where(np.isfinite(dem), dem, np.nanmean(dem))
    f = 4
    fine = gaussian_filter(zoom(dem, f, order=1), sigma=f * 0.8)
    levels = np.arange(np.floor(dem.min() / interval) * interval, dem.max() + interval, interval)
    cs = plt.figure().add_subplot().contour(fine, levels=levels)
    plt.close("all")
    minor, major = [], []
    for i, segments in enumerate(cs.allsegs):
        target = major if round(levels[i] / interval) % every == 0 else minor
        for seg in segments:
            if len(seg) < 4:
                continue
            pts = " L".join(f"{x:.1f} {y:.1f}" for x, y in seg[::2])
            target.append(f"M{pts}")
    h, w = fine.shape
    svg = (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" preserveAspectRatio="xMidYMid slice">'
           f'<path d="{" ".join(minor)}" fill="none" stroke="#7F9A8B" stroke-width="0.9"/>'
           f'<path d="{" ".join(major)}" fill="none" stroke="#5E7D6D" stroke-width="2"/></svg>')
    dest.write_text(svg, encoding="utf-8")


def main():
    hillshade(common.project_path("data", "processed", "dem_utm.tif"),
              common.project_path("data", "processed", "hillshade.tif"))
    layers = []
    for id_, name, file, ramp, vmin, vmax, unit, res in LAYERS:
        tif = common.project_path("data", "processed", file)
        if not tif.exists():
            print(f"  (skipping {name}: {file} is missing)")
            continue
        png = common.project_path("docs", "img", f"{id_}.png")
        bounds = raster_to_png(tif, png, ramp, vmin, vmax, res)
        layers.append({"id": id_, "name": name, "file": f"img/{id_}.png", "bounds": bounds,
                       "ramp": ramp, "min": vmin, "max": vmax, "unit": unit,
                       "legend": ramp_colors(ramp)})
        print(f"  {name}: {png}")

    # Vectors to GeoJSON in lat/lon (EPSG:4326), which is what Leaflet expects.
    for name, file, layer in [("basins", "basins.gpkg", "basins"), ("drainage", "drainage.gpkg", "drainage"),
                              ("bands", "elevation_bands.gpkg", "bands")]:
        gdf = gpd.read_file(common.project_path("data", "processed", file), layer=layer).to_crs(4326)
        gdf.to_file(common.project_path("docs", "data", f"{name}.geojson"), driver="GeoJSON")
        print(f"  {name}.geojson")

    # DEM contour lines for the header.
    contours_svg(common.project_path("data", "processed", "dem_utm.tif"), common.project_path("docs", "img", "contours.svg"))
    print("  contours.svg")

    # Kriging results: metrics to JSON and the variogram to docs/img (if Step 5 has already run).
    metrics = common.project_path("outputs", "interpolation_metrics.csv")
    if metrics.exists():
        df = pd.read_csv(metrics).rename(columns={"Unnamed: 0": "method"})
        common.project_path("docs", "data", "metrics.json").write_text(
            df.to_json(orient="records", force_ascii=False, indent=2), encoding="utf-8")
        print("  metrics.json")
    variogram = common.project_path("outputs", "variogram.png")
    if variogram.exists():
        shutil.copy(variogram, common.project_path("docs", "img", "variogram.png"))
        print("  variogram.png")

    common.project_path("docs", "data", "layers.json").write_text(
        json.dumps({"layers": layers}, ensure_ascii=False, indent=2), encoding="utf-8")
    print("\nDone. Run:  python -m http.server 8000 --directory docs")
    print("and open http://localhost:8000")


if __name__ == "__main__":
    main()
