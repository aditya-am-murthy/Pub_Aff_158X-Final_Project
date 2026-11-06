"""Map MODIS MCD64A1 burned area over Altadena from Planetary Computer.

The public STAC catalog signs blob URLs with a SAS token. DefaultAzureCredential
is for a private Planetary Computer Pro GeoCatalog, not this public collection.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import planetary_computer
import rasterio
import requests
from matplotlib.colors import BoundaryNorm, ListedColormap
from pyproj import Transformer
from rasterio.windows import from_bounds

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "modis-64A1-061"

# Monthly tiles that bracket the Eaton Fire (started 7 Jan 2025, contained 31 Jan).
ITEMS = {
    "before": "MCD64A1.A2024336.h08v05.061.2025035195913",  # December 2024
    "fire": "MCD64A1.A2025001.h08v05.061.2025070160807",  # January 2025
    "after": "MCD64A1.A2025032.h08v05.061.2025098201953",  # February 2025
}

# Altadena plus the Eaton Canyon foothills. South edge is a control band.
WEST, SOUTH, EAST, NORTH = -118.20, 34.155, -118.03, 34.265

PLACES = {
    "JPL": (-118.171, 34.201),
    "Christmas Tree Lane": (-118.131, 34.189),
    "Eaton Canyon": (-118.096, 34.190),
    "Kinneloa Mesa": (-118.078, 34.185),
    "Farnsworth Park": (-118.128, 34.203),
    "Fair Oaks / Loma Alta": (-118.151, 34.193),
}

OUT = Path(__file__).resolve().parent / "output"


def load_item(item_id: str):
    url = f"{STAC}/collections/{COLLECTION}/items/{item_id}"
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    item = resp.json()
    signed = planetary_computer.sign(item)
    return signed


def read_window(asset_href: str):
    with rasterio.open(asset_href) as src:
        transformer = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        x0, y0 = transformer.transform(WEST, SOUTH)
        x1, y1 = transformer.transform(EAST, NORTH)
        left, right = sorted((x0, x1))
        bottom, top = sorted((y0, y1))
        window = from_bounds(left, bottom, right, top, transform=src.transform)
        window = window.round_offsets().round_lengths()
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        crs = src.crs
    return data, transform, crs


def pixel_lonlat(transform, crs, rows, cols):
    xs, ys = rasterio.transform.xy(transform, rows, cols, offset="center")
    to_wgs = Transformer.from_crs(crs, "EPSG:4326", always_xy=True)
    lons, lats = to_wgs.transform(xs, ys)
    return np.asarray(lons), np.asarray(lats)


def month_label(name: str) -> str:
    return {
        "before": "December 2024 (before)",
        "fire": "January 2025 (Eaton Fire)",
        "after": "February 2025 (following month)",
    }[name]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    panels = {}
    summaries = {}

    for name, item_id in ITEMS.items():
        item = load_item(item_id)
        href = item["assets"]["Burn_Date"]["href"]
        data, transform, crs = read_window(href)
        burned = data > 0
        rows, cols = np.where(burned)
        lons, lats = pixel_lonlat(transform, crs, rows, cols) if len(rows) else (np.array([]), np.array([]))
        days = data[rows, cols] if len(rows) else np.array([])
        panels[name] = {
            "data": data,
            "transform": transform,
            "crs": crs,
            "lons": lons,
            "lats": lats,
            "days": days.astype(int),
        }
        summaries[name] = {
            "item": item_id,
            "start": item["properties"].get("start_datetime"),
            "end": item["properties"].get("end_datetime"),
            "burned_pixels": int(burned.sum()),
            "burn_days": sorted({int(d) for d in days.tolist()}),
            "pixels": [
                {"lon": round(float(lon), 5), "lat": round(float(lat), 5), "day_of_year": int(day)}
                for lon, lat, day in zip(lons, lats, days)
            ],
        }
        print(
            f"{name}: {summaries[name]['burned_pixels']} burned 500 m pixels, "
            f"days {summaries[name]['burn_days']}"
        )

    fig, axes = plt.subplots(1, 3, figsize=(14.8, 5.6), sharex=True, sharey=True)
    day_edges = [7, 10, 13, 16, 19, 22, 25, 32]
    colors = ["#fee08b", "#fdae61", "#f46d43", "#d73027", "#a50026", "#67001f", "#2d004b"]
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(day_edges, cmap.N)

    for ax, name in zip(axes, ("before", "fire", "after")):
        panel = panels[name]
        ax.set_xlim(WEST, EAST)
        ax.set_ylim(SOUTH, NORTH)
        ax.set_aspect("equal", adjustable="box")
        ax.set_facecolor("#f4f1ea")
        if len(panel["lons"]):
            sc = ax.scatter(
                panel["lons"],
                panel["lats"],
                c=panel["days"],
                cmap=cmap,
                norm=norm,
                s=42,
                marker="s",
                linewidths=0.2,
                edgecolors="#3d2914",
                zorder=3,
            )
        for label, (lon, lat) in PLACES.items():
            ax.scatter([lon], [lat], c="black", s=12, zorder=4)
            ax.annotate(
                label,
                (lon, lat),
                textcoords="offset points",
                xytext=(4, 4),
                fontsize=7,
                color="#1a1a1a",
            )
        ax.set_title(f"{month_label(name)}\n{summaries[name]['burned_pixels']} burned pixels")
        ax.set_xlabel("Longitude")
        ax.tick_params(labelsize=8)
        ax.grid(True, linestyle=":", linewidth=0.4, color="#b8b1a4")

    axes[0].set_ylabel("Latitude")
    fig.subplots_adjust(left=0.05, right=0.90, bottom=0.16, top=0.82, wspace=0.18)
    if any(len(panels[n]["lons"]) for n in panels):
        cax = fig.add_axes((0.915, 0.22, 0.015, 0.52))
        cbar = fig.colorbar(sc, cax=cax)
        cbar.set_label("Burn day of 2025")

    fig.suptitle(
        "Altadena / Eaton Fire — MODIS MCD64A1 burned area (500 m), Planetary Computer",
        fontsize=13,
    )
    fig.text(
        0.5,
        0.02,
        "Each square is one 500 m pixel with Burn Date > 0. This product maps burned land, not individual trees.\n"
        "December 2024 is the pre-fire control. January 2025 is the Eaton Fire month. Tile h08v05.",
        ha="center",
        fontsize=8,
        color="#333",
    )
    png = OUT / "altadena_mcd64a1_burn.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    (OUT / "altadena_mcd64a1_burn.json").write_text(json.dumps(summaries, indent=2))
    print(f"wrote {png}")


if __name__ == "__main__":
    from visualize_altadena_maps import main as render_maps

    render_maps()
