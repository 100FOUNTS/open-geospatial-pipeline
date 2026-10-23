"""Step 3. 3DEP LiDAR with PDAL: canopy height model (CHM), terrain (DTM) and a LAZ sample.

PDAL reads the points straight from USGS's EPT collection on AWS, with no account and
without downloading the whole collection.

What LiDAR is: an aircraft fires laser pulses at the ground and times how long they take to
bounce back. Each return is a 3D point (x, y, z). In a forest some pulses bounce off the crowns
and others slip through the leaves to the ground, so the cloud holds both crowns and terrain.

In short:
  1. Build several PDAL "pipelines", one per stage (see build_pdal_pipelines in common.py).
  2. Run them in order: download the points for your area (the slow part, done once), clean them,
     separate the ground, compute each point's height above ground and produce the CHM, DTM and sample.
  3. Clean the two resulting rasters (CHM and DTM) and save them as COGs.
If a stage already has its output file it is skipped (delete data/raw/lidar/ to redo everything).

Outputs:
  CHM (data/processed/chm.tif): canopy height above ground, in metres.
  DTM (data/processed/dtm.tif): terrain height without trees, in metres above sea level.
  sample_lidar.laz: the point cloud, to view in 3D in CloudCompare.

Usage:
  python scripts/03_lidar.py                # download and process
  python scripts/03_lidar.py --json-only    # only write the pipelines so you can read them
"""
import argparse
import json
import shutil
import subprocess

import common


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json-only", action="store_true", help="only write the pipeline JSON files")
    args = ap.parse_args()

    # --- 1. Write each pipeline to a JSON file (the "recipe" PDAL reads) ---
    cfg = common.load_config()
    folder = common.project_path("data", "raw", "lidar", "pipeline.json").parent
    pipelines = common.build_pdal_pipelines(cfg, folder)
    for name, pipeline in pipelines.items():
        (folder / f"pipeline_{name}.json").write_text(json.dumps(pipeline, indent=2), encoding="utf-8")
    print(f"Pipelines written to {folder}")
    if args.json_only:
        return

    # --- 2. Run PDAL in order. The download is the slow part ---
    if shutil.which("pdal") is None:
        raise SystemExit("'pdal' not found. Activate the environment: mamba activate andrews")

    # File that proves each stage has finished
    outputs = {"download": "raw_points.laz", "classify": "classified.laz",
               "sample": "sample_lidar.laz", "chm": "chm_raw.tif", "dtm": "dtm_raw.tif"}
    for name, output in outputs.items():
        if (folder / output).exists():
            print(f"[{name}] {output} already exists, skipping")
            continue
        print(f"[{name}] running PDAL...")
        subprocess.run(["pdal", "pipeline", str(folder / f"pipeline_{name}.json")], check=True)

    # --- 3. Clean the rasters PDAL left (see prepare_lidar_raster in common.py) ---
    # CHM: negative values (measurement error) and heights above the plausible maximum
    # (noise) are cut; small gaps are filled up to 10 pixels away.
    chm = common.prepare_lidar_raster(
        folder / "chm_raw.tif", common.project_path("data", "processed", "chm.tif"),
        maximum=cfg["lidar"]["max_canopy_height"], negatives_to_zero=True, fill_px=10)
    # DTM: only fill gaps (where there were no ground points, for example under dense crowns).
    common.prepare_lidar_raster(folder / "dtm_raw.tif", common.project_path("data", "processed", "dtm.tif"),
                                fill_px=10)

    print("\nDone.")
    print("  data/processed/chm.tif  canopy height (m)")
    print("  data/processed/dtm.tif  terrain (m)")
    print("  data/raw/lidar/sample_lidar.laz  point cloud for CloudCompare")
    print(f"  Canopy height: mean {float(chm[chm == chm].mean()):.1f} m, "
          f"maximum {float(chm[chm == chm].max()):.1f} m")


if __name__ == "__main__":
    main()
