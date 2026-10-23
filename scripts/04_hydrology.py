"""Step 4. Hydrology with pysheds: basins, drainage network and elevation bands (analysis zones).

Everything comes from the DEM (data/processed/dem_utm.tif); it uses neither Sentinel-2 nor LiDAR.

In short (details are in common.py: delineate_basins and elevation_bands):
  1. Fill the DEM's pits so water always has a way out.
  2. Work out where water flows from each pixel and how many pixels drain through it (accumulation).
  3. Take the pixel with the highest accumulation as the outlet and delineate everything that drains to it.
     If n_basins > 1, repeat with the highest accumulation among cells not yet in a basin.
  4. Mark as stream the pixels whose accumulation exceeds accumulation_threshold (config.yml).
  5. Split the basin into elevation strips of band_step_m metres: these are the analysis zones
     used in Step 6 to summarise NDVI, NDMI and canopy height.

Limits to keep in mind:
  - The basin is clipped by the edge of the downloaded area; the real one continues beyond it.
  - The Copernicus DEM measures the surface (it includes crowns), not the ground, so water may
    be somewhat misrouted in narrow forested valleys.

Usage:  python scripts/04_hydrology.py
"""
import common


def main():
    cfg = common.load_config()
    dem = common.project_path("data", "processed", "dem_utm.tif")
    if not dem.exists():
        raise SystemExit("data/processed/dem_utm.tif is missing. Run scripts/02_dem.py first")

    # --- Calculation: basins and drainage network (polygons and lines), and elevation bands ---
    h = cfg["hydrology"]
    basins, drainage = common.delineate_basins(dem, h["accumulation_threshold"], h.get("n_basins", 1))
    bands = common.elevation_bands(dem, basins, h["band_step_m"])

    # --- Save as GeoPackage (.gpkg), which QGIS opens directly ---
    basins.to_file(common.project_path("data", "processed", "basins.gpkg"), layer="basins", driver="GPKG")
    drainage.to_file(common.project_path("data", "processed", "drainage.gpkg"), layer="drainage", driver="GPKG")
    bands.to_file(common.project_path("data", "processed", "elevation_bands.gpkg"), layer="bands", driver="GPKG")

    for _, b in basins.iterrows():
        print(f"{b['name']}: {b['area_km2']:.2f} km2 inside the downloaded area")
    print(f"Drainage network: {len(drainage)} segments, {drainage.geometry.length.sum() / 1000:.1f} km")
    print(f"Elevation bands: {len(bands)}")
    print(bands[["band", "label"]].to_string(index=False))


if __name__ == "__main__":
    main()
