"""Step 8. Static figures for the README (figures/*.png).

  overview.png          hillshade with NDVI, NDMI and canopy height, and the basins outlined
  ndvi_vs_elevation.png NDVI and NDMI against elevation (hexbin density + median per 100 m)
  lidar_profile.png     a 1000 m x 20 m slice of the point cloud through the centre of the area
  variogram.png         copy of the R variogram

Needs the outputs of steps 1 to 6. The LiDAR profile also needs data/raw/lidar/classified.laz and
the `pdal` command; if either is missing that figure is skipped.

Usage:  python scripts/08_figures.py
"""
import shutil
import subprocess

import common
import geopandas as gpd
import matplotlib
import numpy as np
import pandas as pd
import rasterio

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

OUT = common.project_path("figures", "placeholder").parent


def _read(name):
    """Raster from data/processed/ with its extent in the (left, right, bottom, top) order imshow expects."""
    with rasterio.open(common.project_path("data", "processed", name)) as src:
        a = src.read(1).astype("float32")
        if src.nodata is not None and np.isfinite(src.nodata):
            a[a == src.nodata] = np.nan
        b = src.bounds
    return a, (b.left, b.right, b.bottom, b.top)


def overview():
    """Three panels (NDVI, NDMI, canopy height) over the hillshade, with the basins outlined."""
    shade, shade_ext = _read("hillshade.tif")
    basins = gpd.read_file(common.project_path("data", "processed", "basins.gpkg"))
    panels = [("Summer NDVI", "ndvi_summer.tif", "RdYlGn", 0.2, 0.95, "NDVI"),
              ("Summer NDMI", "ndmi_summer.tif", "BrBG", -0.1, 0.6, "NDMI"),
              ("Canopy height (LiDAR 2021)", "chm.tif", "YlGn", 0, 60, "metres")]
    # Skip panels whose raster does not exist (for example, no LiDAR for the area).
    panels = [p for p in panels if common.project_path("data", "processed", p[1]).exists()]
    fig, axes = plt.subplots(1, len(panels), figsize=(5 * len(panels), 5.2), constrained_layout=True, squeeze=False)
    for ax, (title, file, ramp, vmin, vmax, label) in zip(axes[0], panels):
        a, ext = _read(file)
        ax.imshow(shade, extent=shade_ext, cmap="gray", vmin=0, vmax=255)
        im = ax.imshow(a, extent=ext, cmap=ramp, vmin=vmin, vmax=vmax, alpha=0.85, interpolation="nearest")
        basins.boundary.plot(ax=ax, color="#14231D", linewidth=1.2)
        ax.set_xlim(shade_ext[0], shade_ext[1])
        ax.set_ylim(shade_ext[2], shade_ext[3])
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(im, ax=ax, shrink=0.8, label=label)
    fig.savefig(OUT / "overview.png", dpi=110)
    plt.close(fig)
    print("  overview.png")


def ndvi_vs_elevation():
    """Hexbin density of NDVI and NDMI against elevation, with the median per 100 m of elevation."""
    t = pd.read_parquet(common.project_path("data", "processed", "pixels.parquet"),
                        columns=["elev", "ndvi", "ndmi"]).dropna()
    z = t["elev"].to_numpy()
    edges = np.arange(np.floor(z.min() / 100) * 100, z.max() + 100, 100)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
    for ax, col, label in zip(axes, ["ndvi", "ndmi"], ["NDVI", "NDMI"]):
        v = t[col].to_numpy()
        ax.hexbin(z, v, gridsize=50, bins="log", cmap="viridis")
        centres, medians = [], []
        for a in edges:
            s = (z >= a) & (z < a + 100)
            if s.sum() >= 20:
                centres.append(a + 50)
                medians.append(np.median(v[s]))
        ax.plot(centres, medians, "r-o", label="median per 100 m")
        ax.set_xlabel("Elevation (m)")
        ax.set_ylabel(f"Summer {label}")
        ax.set_title(f"{label} vs elevation")
        ax.legend(loc="lower left")
    fig.savefig(OUT / "ndvi_vs_elevation.png", dpi=110)
    plt.close(fig)
    print("  ndvi_vs_elevation.png")


def lidar_profile(length=1000, width=20):
    """Slice of the classified cloud: ground vs vegetation (top) and height above ground (bottom)."""
    las = common.project_path("data", "raw", "lidar", "classified.laz")
    chm = common.project_path("data", "processed", "chm.tif")
    if not las.exists() or not shutil.which("pdal"):
        print("  (skipping lidar_profile.png: needs classified.laz and pdal)")
        return
    with rasterio.open(chm) as src:
        b = src.bounds
    cx, cy = (b.left + b.right) / 2, (b.bottom + b.top) / 2
    bounds = f"([{cx - length / 2:.1f},{cx + length / 2:.1f}],[{cy - width / 2:.1f},{cy + width / 2:.1f}])"
    csv = common.project_path("data", "raw", "lidar", "profile_strip.csv")
    subprocess.run(["pdal", "translate", str(las), str(csv), f"--filters.crop.bounds={bounds}",
                    "--writers.text.order=X,Y,Z,Classification,HeightAboveGround",
                    "--writers.text.keep_unspecified=false", "-f", "crop", "-w", "text"],
                   check=True, capture_output=True)
    d = np.genfromtxt(csv, delimiter=",", names=True)
    x = d["X"] - d["X"].min()
    z, c, h = d["Z"], d["Classification"], d["HeightAboveGround"]
    ground = c == 2
    fig, ax = plt.subplots(2, 1, figsize=(14, 8), sharex=True, constrained_layout=True)
    ax[0].scatter(x[~ground], z[~ground], s=1, c="#2e7d32", label="non-ground (crowns, trunks)")
    ax[0].scatter(x[ground], z[ground], s=2.5, c="#8d5524", label="ground (class 2)")
    ax[0].set_ylabel("Elevation Z (m)")
    ax[0].legend(markerscale=8, loc="upper right")
    ax[0].set_title(f"{length} m x {width} m slice through the centre of the area (LiDAR, EPT resolution 5)")
    s = ax[1].scatter(x, z, s=1, c=h, cmap="viridis")
    fig.colorbar(s, ax=ax[1], label="HeightAboveGround (m)")
    ax[1].set_xlabel("East-west distance (m)")
    ax[1].set_ylabel("Elevation Z (m)")
    fig.savefig(OUT / "lidar_profile.png", dpi=100)
    plt.close(fig)
    csv.unlink()
    print("  lidar_profile.png")


def main():
    overview()
    ndvi_vs_elevation()
    lidar_profile()
    variogram = common.project_path("outputs", "variogram.png")
    if variogram.exists():
        shutil.copy(variogram, OUT / "variogram.png")
        print("  variogram.png")


if __name__ == "__main__":
    main()
