"""Functions shared by the project scripts.

The logic that can be tested without internet lives here (see tests/).
The numbered scripts only download data, call these functions and save results.
"""
from pathlib import Path

import geopandas as gpd
import numpy as np
import rasterio
import rioxarray  # noqa: F401  (enables the .rio accessor on xarray)
import xarray as xr
import yaml
from pyproj import Transformer
from rasterio.features import rasterize, shapes
from rasterio.fill import fillnodata
from rasterio.shutil import copy as rio_copy
from shapely.geometry import shape
from shapely.ops import unary_union

ROOT = Path(__file__).resolve().parents[1]

# Sentinel-2 L2A via Earth Search: reflectance = (DN + OFFSET_DN) * SCALE. DN = 0 means "no data".
# The catalogue declares an offset of -0.1 (-1000 DN), but the values that actually arrive no longer
# carry it: red DNs below 1000 were measured (minimum 196), impossible if the offset applied.
# With -1000 the forest NDVI came out as exactly 1.000, because negative reflectances were clipped to 0.
# Hence OFFSET_DN = 0. If you change data source, check it again with the control pixel.
SCALE = 1e-4
OFFSET_DN = 0


# ---------------------------------------------------------------- configuration
def load_config():
    with open(ROOT / "config.yml", encoding="utf-8") as f:
        return yaml.safe_load(f)


def project_path(*parts):
    """Path inside the project; creates the parent folder if it does not exist."""
    p = ROOT.joinpath(*parts)
    p.parent.mkdir(parents=True, exist_ok=True)
    return p


# -------------------------------------------------------------------- Sentinel-2
def to_reflectance(dn):
    """Sentinel-2 L2A DN to reflectance between 0 and 1."""
    dn = dn.astype("float32")
    refl = (dn + OFFSET_DN) * SCALE
    return refl.where(dn > 0).clip(0.0, 1.0)


def spectral_indices(cube, valid_classes=(4, 5)):
    """NDVI and NDMI per date, only on pixels without clouds or shadows.

    SCL classes used: 4 = vegetation, 5 = bare soil.
    NDVI = (NIR - red) / (NIR + red)
    NDMI = (NIR - SWIR1) / (NIR + SWIR1), associated with the water content of vegetation.
    """
    valid = cube["scl"].isin(list(valid_classes))
    red = to_reflectance(cube["red"]).where(valid)
    nir = to_reflectance(cube["nir"]).where(valid)
    swir = to_reflectance(cube["swir16"]).where(valid)
    ndvi = (nir - red) / (nir + red)
    ndmi = (nir - swir) / (nir + swir)
    ndvi = ndvi.where(np.isfinite(ndvi))
    ndmi = ndmi.where(np.isfinite(ndmi))
    return xr.Dataset({"ndvi": ndvi, "ndmi": ndmi})


def aggregate_monthly(indices):
    return indices.resample(time="1MS").median(skipna=True)


def season_mean(monthly, months):
    sel = monthly.sel(time=monthly["time"].dt.month.isin(list(months)))
    if sel.sizes["time"] == 0:
        raise ValueError(f"Months {list(months)} are not in the downloaded period.")
    return sel.mean("time", skipna=True)


def check_pixel(cube, indices, lon, lat, epsg):
    """Recompute the NDVI of one pixel by hand and compare it with the cube's value."""
    x, y = Transformer.from_crs(4326, epsg, always_xy=True).transform(lon, lat)
    px = cube.sel(x=x, y=y, method="nearest")
    ix = indices.sel(x=x, y=y, method="nearest")
    ndvi_series = ix["ndvi"].values
    valid = np.where(np.isfinite(ndvi_series))[0]
    if len(valid) == 0:
        print("  Control pixel: no clean observations at that point.")
        return None
    t = int(valid[0])
    dn = {b: int(px[b].isel(time=t).values) for b in ("red", "nir", "swir16")}
    refl = {b: (v + OFFSET_DN) * SCALE for b, v in dn.items()}
    ndvi_by_hand = (refl["nir"] - refl["red"]) / (refl["nir"] + refl["red"])
    date = str(px["time"].isel(time=t).values)[:10]
    print(f"  Control pixel ({lon}, {lat}) on {date}")
    print(f"    DN red={dn['red']}  NIR={dn['nir']}  SWIR1={dn['swir16']}")
    print(f"    NDVI by hand = {ndvi_by_hand:.3f}   NDVI from the cube = {float(ndvi_series[t]):.3f}")
    print("    Dense summer forest should give values between 0.75 and 0.90.")
    return ndvi_by_hand, float(ndvi_series[t])


# ----------------------------------------------------------------------- rasters
def write_cog(da, dest, crs):
    """Save a (y, x) DataArray as a float32 Cloud Optimized GeoTIFF."""
    da = da.astype("float32").rio.write_crs(crs).rio.write_nodata(np.nan)
    da.rio.to_raster(str(dest), driver="COG", compress="deflate")


def _to_cog(src_path, dest):
    rio_copy(str(src_path), str(dest), driver="COG", compress="deflate")


def prepare_lidar_raster(src_path, dest, maximum=None, negatives_to_zero=False, fill_px=0):
    """Clean a PDAL raster: nodata to NaN, remove spikes, fill small gaps and save as COG."""
    with rasterio.open(src_path) as src:
        a = src.read(1).astype("float32")
        profile = src.profile.copy()
        nodata = src.nodata
    valid = np.isfinite(a)
    if nodata is not None and np.isfinite(nodata):
        valid &= a != nodata
    a = np.where(valid, a, np.nan).astype("float32")
    if negatives_to_zero:
        a = np.where(a < 0, 0.0, a)
    if maximum is not None:
        a = np.where(a > maximum, np.nan, a)
    valid = np.isfinite(a)
    if fill_px > 0 and valid.any() and not valid.all():
        mask = valid.astype("uint8")
        filled = fillnodata(np.where(valid, a, 0).astype("float32"), mask=mask,
                            max_search_distance=fill_px, smoothing_iterations=0)
        covered = fillnodata(valid.astype("float32"), mask=mask,
                             max_search_distance=fill_px, smoothing_iterations=0) > 0.5
        a = np.where(valid, a, np.where(covered, filled, np.nan)).astype("float32")
    profile.update(dtype="float32", nodata=np.nan, count=1, driver="GTiff", compress="deflate")
    tmp = Path(str(dest) + ".tmp.tif")
    with rasterio.open(tmp, "w", **profile) as dst:
        dst.write(a, 1)
    _to_cog(tmp, dest)
    tmp.unlink()
    return a


# ---------------------------------------------------------------------- LiDAR
def build_pdal_pipelines(cfg, out_dir):
    """Build the PDAL LiDAR pipelines, one per stage, in the order they must run.

    They are split because PDAL runs only one output (writer) per pipeline, and so that the
    download, the slow part, happens once: everything else is computed from local files.

      download : reads 3DEP (EPT) -> reprojects to UTM -> saves raw_points.laz
      classify : removes noise -> classifies ground -> height above ground -> saves classified.laz
      sample   : keeps 1 of every 4 points of classified.laz -> sample_lidar.laz (for CloudCompare)
      chm      : maximum height above ground per cell -> chm_raw.tif
      dtm      : ground points only, interpolated per cell -> dtm_raw.tif
    """
    folder = Path(out_dir)
    raw = str(folder / "raw_points.laz")
    classified = str(folder / "classified.laz")
    lon0, lat0, lon1, lat1 = cfg["area"]["bbox_lidar"]
    to_3857 = Transformer.from_crs(4326, 3857, always_xy=True)
    x0, y0 = to_3857.transform(lon0, lat0)
    x1, y1 = to_3857.transform(lon1, lat1)
    lid = cfg["lidar"]
    utm = f"EPSG:{cfg['area']['epsg_utm']}"
    cell = lid["raster_resolution"]
    # LAZ in metres: 1 cm precision and automatic offset.
    scale = {"scale_x": 0.01, "scale_y": 0.01, "scale_z": 0.01,
             "offset_x": "auto", "offset_y": "auto", "offset_z": "auto"}

    download = [
        # 1. Read the points from the public collection, only inside bounds (Web Mercator metres).
        #    "resolution" limits the density: the larger it is, the fewer points are downloaded.
        {"type": "readers.ept", "tag": "read", "filename": lid["ept"],
         "bounds": f"([{x0:.1f},{x1:.1f}],[{y0:.1f},{y1:.1f}])",
         "resolution": lid["ept_resolution"]},
        # 2. Reproject to UTM 10N so it matches the other layers (everything in metres).
        {"type": "filters.reprojection", "tag": "utm", "inputs": ["read"], "out_srs": utm},
        {"type": "writers.las", "inputs": ["utm"], "filename": raw,
         "compression": "laszip", **scale},
    ]
    classify = [
        {"type": "readers.las", "tag": "read", "filename": raw},
        # 3. Flag isolated points (few neighbours nearby) as noise: birds, sensor errors.
        {"type": "filters.outlier", "tag": "denoised", "inputs": ["read"],
         "method": "statistical", "mean_k": 8, "multiplier": 3.0},
        # 4. Drop noise (classes 7 and 18 onward). Ranges on the same dimension are combined with OR.
        {"type": "filters.range", "tag": "clean", "inputs": ["denoised"],
         "limits": "Classification[0:6],Classification[8:17]"},
        # 5. SMRF stops if some points have ReturnNumber/NumberOfReturns = 0 and others do not.
        #    Setting them to 1 (a single return) lets it classify.
        {"type": "filters.assign", "tag": "returns", "inputs": ["clean"],
         "assignment": ["ReturnNumber[0:0]=1", "NumberOfReturns[0:0]=1"]},
        # 6. SMRF classifies the ground (class 2). "slope" is the maximum expected terrain slope.
        {"type": "filters.smrf", "tag": "ground", "inputs": ["returns"],
         "slope": lid["smrf_slope"], "window": 18, "threshold": 0.5, "scalar": 1.25},
        # 7. Compute HeightAboveGround: the height of each point above the nearest ground.
        {"type": "filters.hag_nn", "tag": "with_hag", "inputs": ["ground"]},
        {"type": "writers.las", "inputs": ["with_hag"], "filename": classified,
         "compression": "laszip", "extra_dims": "all", **scale},
    ]
    sample = [
        {"type": "readers.las", "tag": "read", "filename": classified},
        # 8. One point in every 4: the full cloud is too heavy for CloudCompare.
        {"type": "filters.decimation", "tag": "thinned", "inputs": ["read"], "step": 4},
        {"type": "writers.las", "inputs": ["thinned"], "filename": str(folder / "sample_lidar.laz"),
         "compression": "laszip", "extra_dims": "all", **scale},
    ]
    chm = [
        {"type": "readers.las", "tag": "read", "filename": classified},
        # 9. CHM: for each cell, the maximum height above ground (the tallest crown).
        {"type": "writers.gdal", "inputs": ["read"],
         "filename": str(folder / "chm_raw.tif"), "resolution": cell,
         "output_type": "max", "dimension": "HeightAboveGround",
         "nodata": -9999, "gdaldriver": "GTiff", "data_type": "float32"},
    ]
    dtm = [
        {"type": "readers.las", "tag": "read", "filename": classified},
        # 10. DTM: ground only (class 2), interpolated (IDW) in each cell.
        {"type": "filters.range", "tag": "ground_only", "inputs": ["read"], "limits": "Classification[2:2]"},
        {"type": "writers.gdal", "inputs": ["ground_only"],
         "filename": str(folder / "dtm_raw.tif"), "resolution": cell,
         "output_type": "idw", "window_size": 6,
         "nodata": -9999, "gdaldriver": "GTiff", "data_type": "float32"},
    ]
    return {"download": {"pipeline": download}, "classify": {"pipeline": classify},
            "sample": {"pipeline": sample}, "chm": {"pipeline": chm}, "dtm": {"pipeline": dtm}}


# ------------------------------------------------------------------ hydrology
def delineate_basins(dem_path, threshold_cells, n_basins=1):
    """Fill depressions, compute flow direction and accumulation, and delineate the basins.

    The outlet of the first basin is the cell with the highest accumulation in the DEM, so no
    coordinates are needed. With n_basins > 1 the process repeats, each time choosing the cell
    with the highest accumulation among those that do not yet belong to a delineated basin.
    Basins are clipped by the DEM edge.
    """
    # pysheds 0.5 uses np.in1d, which NumPy 2.4 removed. This alias avoids the error.
    if not hasattr(np, "in1d"):
        np.in1d = lambda a, b, **kw: np.isin(np.ravel(a), b, **kw)
    from pysheds.grid import Grid

    grid = Grid.from_raster(str(dem_path))
    dem = grid.read_raster(str(dem_path))
    no_pits = grid.fill_pits(dem)
    no_depressions = grid.fill_depressions(no_pits)
    flat = grid.resolve_flats(no_depressions)
    fdir = grid.flowdir(flat)
    acc = grid.accumulation(fdir)

    acc_np = np.asarray(acc)
    with rasterio.open(dem_path) as src:
        crs = src.crs

    taken = np.zeros(acc_np.shape, dtype=bool)   # cells that already belong to some basin
    rows = []
    for k in range(1, n_basins + 1):
        candidates = np.where(taken, -1, acc_np)
        iy, ix = np.unravel_index(int(np.argmax(candidates)), candidates.shape)
        # With xytype="index" (column, row) the outlet cell is not lost at the DEM edge.
        catchment = grid.catchment(x=int(ix), y=int(iy), fdir=fdir, xytype="index")
        mask = np.asarray(catchment).astype(bool)
        taken |= mask
        polygons = [shape(g) for g, v in shapes(mask.astype("uint8"), mask=mask,
                                                transform=grid.affine) if v == 1]
        rows.append({"name": f"basin_{k}" if n_basins > 1 else "basin",
                     "geometry": unary_union(polygons)})
    basins = gpd.GeoDataFrame(rows, geometry="geometry", crs=crs)
    basins["area_km2"] = basins.geometry.area / 1e6

    # The network is computed once over the whole DEM and clipped to the delineated basins.
    branches = grid.extract_river_network(fdir, acc > threshold_cells)
    drainage = gpd.GeoDataFrame.from_features(branches["features"], crs=crs)
    drainage = gpd.clip(drainage, unary_union(basins.geometry))
    drainage = drainage[~drainage.geometry.is_empty].reset_index(drop=True)
    return basins, drainage


def elevation_bands(dem_path, basins, step):
    """Elevation-band polygons inside the basin (the analysis zones)."""
    with rasterio.open(dem_path) as src:
        dem = src.read(1).astype("float32")
        trans, crs = src.transform, src.crs
    inside = rasterize(basins.geometry, out_shape=dem.shape, transform=trans,
                       fill=0, default_value=1, dtype="uint8").astype(bool)
    dem_b = np.where(inside & np.isfinite(dem), dem, np.nan)
    base = np.floor(np.nanmin(dem_b) / step) * step
    ids = np.floor((dem_b - base) / step)
    ids_i = np.where(np.isfinite(ids), ids, -1).astype("int32")
    geoms, vals = [], []
    for g, v in shapes(ids_i, mask=ids_i >= 0, transform=trans):
        geoms.append(shape(g))
        vals.append(int(v))
    gdf = gpd.GeoDataFrame({"band": vals}, geometry=geoms, crs=crs).dissolve(by="band").reset_index()
    gdf["z_min"] = base + gdf["band"] * step
    gdf["z_max"] = gdf["z_min"] + step
    gdf["label"] = gdf.apply(lambda r: f"{int(r['z_min'])}-{int(r['z_max'])} m", axis=1)
    return gdf
