"""Offline tests on synthetic data.

They confirm your installation works before you download anything:
  python -m pytest tests -q
"""
import json
import sys
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
import rasterio
import rioxarray  # noqa: F401
import xarray as xr
from rasterio.transform import from_origin

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import common  # noqa: E402


def _synthetic_cube(n_months=4, ny=6, nx=7):
    times = pd.date_range("2025-06-15", periods=n_months, freq="MS") + pd.Timedelta(days=14)
    x = 500000 + 20 * np.arange(nx)
    y = 4900000 - 20 * np.arange(ny)
    shape = (n_months, ny, nx)
    coords = {"time": times, "y": y, "x": x}
    dims = ("time", "y", "x")
    red = np.full(shape, 500, dtype="uint16")      # reflectance 0.05
    nir = np.full(shape, 3500, dtype="uint16")     # 0.35
    swir = np.full(shape, 1500, dtype="uint16")    # 0.15
    scl = np.full(shape, 4, dtype="uint8")
    scl[0, 0, 0] = 9                               # cloud in one pixel
    red[1, 1, 1] = 0                               # no data
    return xr.Dataset({
        "red": (dims, red), "nir": (dims, nir), "swir16": (dims, swir), "scl": (dims, scl),
    }, coords=coords)


def test_reflectance_scale_and_no_data():
    dn = xr.DataArray(np.array([0, 1000, 5000, 10000, 20000], dtype="uint16"))
    r = common.to_reflectance(dn).values
    assert np.isnan(r[0])
    assert r[1] == pytest.approx(0.1)
    assert r[2] == pytest.approx(0.5)
    assert r[3] == pytest.approx(1.0)
    assert r[4] == pytest.approx(1.0)   # clipped to 1


def test_indices_and_masking():
    cube = _synthetic_cube()
    ind = common.spectral_indices(cube)
    # reflectances: red 0.05 ; NIR 0.35 ; SWIR 0.15
    assert float(ind["ndvi"][2, 3, 3]) == pytest.approx((0.35 - 0.05) / (0.35 + 0.05), abs=1e-4)
    assert float(ind["ndmi"][2, 3, 3]) == pytest.approx((0.35 - 0.15) / (0.35 + 0.15), abs=1e-4)
    assert np.isnan(float(ind["ndvi"][0, 0, 0]))   # cloud
    assert np.isnan(float(ind["ndvi"][1, 1, 1]))   # no data


def test_monthly_and_summer_mean_to_cog(tmp_path):
    ind = common.spectral_indices(_synthetic_cube(n_months=5))
    monthly = common.aggregate_monthly(ind)
    assert monthly.sizes["time"] == 5
    summer = common.season_mean(monthly, [7, 8, 9])
    assert summer["ndvi"].dims == ("y", "x")
    dest = tmp_path / "ndvi.tif"
    common.write_cog(summer["ndvi"], dest, "EPSG:32610")
    with rasterio.open(dest) as src:
        assert src.crs.to_epsg() == 32610
        assert src.count == 1
        assert src.read(1).dtype == np.float32
    zarr_dir = tmp_path / "m.zarr"
    monthly.to_zarr(zarr_dir, mode="w")
    assert zarr_dir.exists()


def test_season_mean_without_months_fails():
    ind = common.spectral_indices(_synthetic_cube(n_months=2))
    monthly = common.aggregate_monthly(ind)
    with pytest.raises(ValueError):
        common.season_mean(monthly, [1, 2])


def test_prepare_lidar_raster_fills_and_cleans(tmp_path):
    a = np.full((40, 40), 20.0, dtype="float32")
    a[10:12, 10:12] = -9999          # small gap
    a[30, 30] = 120.0                # impossible spike
    a[5, 5] = -2.0                   # negative
    src_path = tmp_path / "raw.tif"
    with rasterio.open(src_path, "w", driver="GTiff", height=40, width=40, count=1,
                       dtype="float32", crs="EPSG:32610", transform=from_origin(0, 40, 1, 1),
                       nodata=-9999) as dst:
        dst.write(a, 1)
    out = tmp_path / "chm.tif"
    r = common.prepare_lidar_raster(src_path, out, maximum=70, negatives_to_zero=True, fill_px=10)
    assert r[10, 10] == pytest.approx(20.0)      # gap filled
    assert r[5, 5] == 0.0                         # negative set to zero
    assert r[30, 30] == pytest.approx(20.0)      # spike removed and filled
    with rasterio.open(out) as src:
        assert src.profile["driver"] == "GTiff"
        assert src.read(1).shape == (40, 40)


def test_pdal_pipelines_are_consistent():
    cfg = common.load_config()
    pipes = common.build_pdal_pipelines(cfg, "/tmp/lidar")
    assert list(pipes) == ["download", "classify", "sample", "chm", "dtm"]
    json.dumps(pipes)                               # serialisable
    for name, pipe in pipes.items():
        stages = pipe["pipeline"]
        tags = {s["tag"] for s in stages if "tag" in s}
        for s in stages:
            for inp in s.get("inputs", []):
                assert inp in tags, f"{name}: {inp} is not defined"
        # PDAL runs a single output per pipeline: exactly one writer, at the end
        assert sum(s["type"].startswith("writers.") for s in stages) == 1
        assert stages[-1]["type"].startswith("writers.")
    assert pipes["download"]["pipeline"][0]["type"] == "readers.ept"
    xs = pipes["download"]["pipeline"][0]["bounds"]
    assert xs.startswith("([") and "],[" in xs


def _synthetic_dem(tif_path, n=90):
    """V-shaped valley sloping south, with a small undulation."""
    yy, xx = np.mgrid[0:n, 0:n].astype("float32")
    dem = 1200 - 6.0 * yy + 4.0 * np.abs(xx - n / 2) + 3.0 * np.sin(xx / 5.0)
    dem = (dem - dem.min() + 800).astype("float32")
    with rasterio.open(tif_path, "w", driver="GTiff", height=n, width=n, count=1, dtype="float32",
                       crs="EPSG:32610", transform=from_origin(500000, 4900000, 30, 30)) as dst:
        dst.write(dem, 1)


def test_synthetic_hydrology(tmp_path):
    dem = tmp_path / "dem.tif"
    _synthetic_dem(dem)
    basins, drainage = common.delineate_basins(dem, threshold_cells=30)
    assert len(basins) == 1
    assert float(basins["area_km2"].iloc[0]) > 0.5
    assert len(drainage) > 0
    bands = common.elevation_bands(dem, basins, step=150)
    assert len(bands) >= 2
    assert {"band", "label", "z_min", "z_max"} <= set(bands.columns)
    assert isinstance(bands, gpd.GeoDataFrame)
    # With two basins: they do not overlap and each has its own area
    two, _ = common.delineate_basins(dem, threshold_cells=30, n_basins=2)
    assert len(two) == 2
    assert two.geometry.iloc[0].intersection(two.geometry.iloc[1]).area < 1.0
