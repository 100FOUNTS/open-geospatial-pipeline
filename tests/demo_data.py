"""Generate SYNTHETIC demo data in data/processed/ to try steps 4 to 7 without internet.

This is not real data. It is only there to check that your installation works.
If you have already downloaded the real data, do not run this file: it would overwrite your results.

Usage:  python tests/demo_data.py
"""
import sys
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import common  # noqa: E402


def _save(name, a, res, x0, y0, crs="EPSG:32610"):
    dest = common.project_path("data", "processed", name)
    with rasterio.open(dest, "w", driver="GTiff", height=a.shape[0], width=a.shape[1], count=1,
                       dtype="float32", crs=crs, transform=from_origin(x0, y0, res, res),
                       nodata=np.nan) as dst:
        dst.write(a.astype("float32"), 1)


def generate():
    x0, y0 = 560000.0, 4900000.0
    n30 = 150
    yy, xx = np.mgrid[0:n30, 0:n30].astype("float32")
    dem = 800 + 5.0 * (n30 - yy) + 3.0 * np.abs(xx - n30 / 2) + 8 * np.sin(xx / 9.0)
    _save("dem_utm.tif", dem, 30, x0, y0)

    n20 = 225
    y2, x2 = np.mgrid[0:n20, 0:n20].astype("float32")
    elev = 800 + 5.0 * (n20 - y2) * 1.5 + 3.0 * np.abs(x2 - n20 / 2) * 1.5
    rng = np.random.default_rng(1)
    noise = rng.normal(0, 0.03, elev.shape)
    ndvi = np.clip(0.9 - (elev - 800) / 4000 + noise, 0.2, 0.95)
    ndmi = np.clip(0.45 - (elev - 800) / 5000 + rng.normal(0, 0.03, elev.shape), -0.1, 0.6)
    _save("ndvi_summer.tif", ndvi, 20, x0, y0)
    _save("ndmi_summer.tif", ndmi, 20, x0, y0)

    n1 = 400
    chm = np.clip(35 - (np.arange(n1)[:, None] * 0.03) + rng.normal(0, 3, (n1, n1)), 0, 60)
    _save("chm.tif", chm, 1.5, x0 + 1500, y0 - 1500)

    basins, drainage = common.delineate_basins(common.project_path("data", "processed", "dem_utm.tif"), 100)
    bands = common.elevation_bands(common.project_path("data", "processed", "dem_utm.tif"), basins, 150)
    basins.to_file(common.project_path("data", "processed", "basins.gpkg"), layer="basins", driver="GPKG")
    drainage.to_file(common.project_path("data", "processed", "drainage.gpkg"), layer="drainage", driver="GPKG")
    bands.to_file(common.project_path("data", "processed", "elevation_bands.gpkg"), layer="bands", driver="GPKG")
    print("Demo data created in data/processed/")


if __name__ == "__main__":
    generate()
