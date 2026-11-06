"""Overlay MODIS burned area and leaf area on a map of Altadena.

Public Planetary Computer assets are signed with a SAS token. The
DefaultAzureCredential snippet is for a private Planetary Computer Pro
catalog, not this public STAC API.
"""

from __future__ import annotations

import json
from pathlib import Path

import contextily as cx
import matplotlib.pyplot as plt
import numpy as np
import planetary_computer
import rasterio
import requests
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.transform import array_bounds, from_bounds as transform_from_bounds
from rasterio.warp import reproject
from rasterio.windows import from_bounds

STAC = "https://planetarycomputer.microsoft.com/api/stac/v1"
OUT = Path(__file__).resolve().parent / "output"

WEST, SOUTH, EAST, NORTH = -118.20, 34.155, -118.03, 34.265
DISPLAY_RES_M = 90

PLACES = {
    "JPL": (-118.171, 34.201),
    "Christmas Tree Lane": (-118.131, 34.189),
    "Eaton Canyon": (-118.096, 34.190),
    "Kinneloa Mesa": (-118.078, 34.185),
    "Farnsworth Park": (-118.128, 34.203),
    "Fair Oaks / Loma Alta": (-118.151, 34.193),
}

BURN_ITEMS = {
    "before": "MCD64A1.A2024336.h08v05.061.2025035195913",
    "fire": "MCD64A1.A2025001.h08v05.061.2025070160807",
    "after": "MCD64A1.A2025032.h08v05.061.2025098201953",
}
LAI_ITEMS = {
    "before": "MYD15A2H.A2024281.h08v05.061.2024290063654",
    "after": "MYD15A2H.A2026281.h08v05.061.2026290035710",
}

TO_MERCATOR = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)
MAP_LEFT, MAP_BOTTOM = TO_MERCATOR.transform(WEST, SOUTH)
MAP_RIGHT, MAP_TOP = TO_MERCATOR.transform(EAST, NORTH)


def load_signed_item(collection: str, item_id: str) -> dict:
    url = f"{STAC}/collections/{collection}/items/{item_id}"
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return planetary_computer.sign(response.json())


def read_window(collection: str, item_id: str, asset: str):
    item = load_signed_item(collection, item_id)
    href = item["assets"][asset]["href"]
    with rasterio.open(href) as src:
        to_src = Transformer.from_crs("EPSG:4326", src.crs, always_xy=True)
        x0, y0 = to_src.transform(WEST, SOUTH)
        x1, y1 = to_src.transform(EAST, NORTH)
        left, right = sorted((x0, x1))
        bottom, top = sorted((y0, y1))
        window = from_bounds(left, bottom, right, top, transform=src.transform)
        window = window.round_offsets().round_lengths()
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        meta = src.meta.copy()
    return data, transform, meta, item


def mercator_grid(src, transform, crs, resampling):
    width = max(2, int(round((MAP_RIGHT - MAP_LEFT) / DISPLAY_RES_M)))
    height = max(2, int(round((MAP_TOP - MAP_BOTTOM) / DISPLAY_RES_M)))
    dst_transform = transform_from_bounds(MAP_LEFT, MAP_BOTTOM, MAP_RIGHT, MAP_TOP, width, height)
    destination = np.full((height, width), np.nan, dtype="float32")
    source = src.astype("float32")
    reproject(
        source,
        destination,
        src_transform=transform,
        src_crs=crs,
        src_nodata=np.nan,
        dst_transform=dst_transform,
        dst_crs="EPSG:3857",
        dst_nodata=np.nan,
        resampling=resampling,
    )
    return destination, dst_transform


def style_map(ax, title: str):
    ax.set_xlim(MAP_LEFT, MAP_RIGHT)
    ax.set_ylim(MAP_BOTTOM, MAP_TOP)
    ax.set_aspect("equal")
    cx.add_basemap(
        ax,
        crs="EPSG:3857",
        source=cx.providers.Esri.WorldTopoMap,
        zoom=14,
        attribution_size=5,
    )
    ax.set_xlim(MAP_LEFT, MAP_RIGHT)
    ax.set_ylim(MAP_BOTTOM, MAP_TOP)
    ax.set_title(title, fontsize=11)
    ax.set_axis_off()
    for label, (lon, lat) in PLACES.items():
        x, y = TO_MERCATOR.transform(lon, lat)
        ax.annotate(
            label,
            (x, y),
            textcoords="offset points",
            xytext=(3, 3),
            fontsize=6.5,
            color="#1a1a1a",
            zorder=5,
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.78},
        )


def imshow_grid(ax, grid, cmap, norm):
    left, bottom, right, top = array_bounds(grid.shape[0], grid.shape[1], 
        transform_from_bounds(MAP_LEFT, MAP_BOTTOM, MAP_RIGHT, MAP_TOP, grid.shape[1], grid.shape[0])
    )
    return ax.imshow(
        grid,
        extent=(left, right, bottom, top),
        origin="upper",
        cmap=cmap,
        norm=norm,
        interpolation="nearest",
        alpha=0.78,
        zorder=2,
    )


def burn_cmap():
    colors = ["#fee08b", "#fdae61", "#f46d43", "#d73027", "#a50026", "#67001f"]
    cmap = ListedColormap(colors)
    cmap.set_bad((0, 0, 0, 0))
    norm = BoundaryNorm([6, 8, 9, 10, 11, 12, 14], cmap.N)
    return cmap, norm


def lai_array(raw):
    lai = raw.astype("float32") * 0.1
    lai[raw > 100] = np.nan
    return lai


def summarize_lai(before, after, burned):
    both = np.isfinite(before) & np.isfinite(after)
    inside = both & burned
    outside = both & ~burned

    def block(mask):
        delta = after[mask] - before[mask]
        return {
            "pixels": int(mask.sum()),
            "lai_2024": round(float(before[mask].mean()), 3),
            "lai_2026": round(float(after[mask].mean()), 3),
            "change": round(float(delta.mean()), 3),
            "pixels_down_by_0_3": int((delta <= -0.3).sum()),
        }

    return {
        "valid_leaf_area_pixels": int(both.sum()),
        "burned_pixels": int(burned.sum()),
        "burned_pixels_without_lai": int((burned & ~np.isfinite(before)).sum()),
        "inside_burn": block(inside),
        "outside_burn": block(outside),
    }


def draw_burn_outline(ax, burned, transform, crs):
    mask = np.where(burned, 1.0, np.nan).astype("float32")
    grid, dst_transform = mercator_grid(mask, transform, crs, Resampling.nearest)
    height, width = grid.shape
    cols = np.arange(width)
    rows = np.arange(height)
    xs = np.asarray(rasterio.transform.xy(dst_transform, np.zeros(width), cols, offset="center")[0])
    ys = np.asarray(rasterio.transform.xy(dst_transform, rows, np.zeros(height), offset="center")[0])
    filled = np.nan_to_num(grid, nan=0.0)
    ax.contour(xs, ys[::-1], filled[::-1], levels=[0.5], colors="white", linewidths=1.8, zorder=4)
    ax.contour(xs, ys[::-1], filled[::-1], levels=[0.5], colors="#1a1a1a", linewidths=0.7, zorder=5)


def plot_burn(panels):
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 6.2))
    cmap, norm = burn_cmap()
    image = None
    titles = {
        "before": "December 2024\nno burned pixels",
        "fire": "January 2025, Eaton Fire\nburn day on the Altadena map",
        "after": "February 2025\nno burned pixels",
    }
    for ax, name in zip(axes, ("before", "fire", "after")):
        panel = panels[name]
        style_map(ax, titles[name])
        raw = panel["data"].astype("float32")
        raw[raw <= 0] = np.nan
        if np.isfinite(raw).any():
            grid, _ = mercator_grid(raw, panel["transform"], panel["crs"], Resampling.nearest)
            image = imshow_grid(ax, grid, cmap, norm)
    fig.suptitle("MODIS burned area (MCD64A1, 500 m) on a map of Altadena", fontsize=14)
    if image is not None:
        cax = fig.add_axes((0.90, 0.22, 0.014, 0.52))
        cbar = fig.colorbar(image, cax=cax)
        cbar.set_label("Burn day of 2025")
    fig.text(
        0.44,
        0.03,
        "Basemap: Esri World Topographic Map. Colored cells are 500 m burned-land pixels, not individual trees.",
        ha="center",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.02, right=0.88, bottom=0.08, top=0.84, wspace=0.04)
    path = OUT / "altadena_mcd64a1_burn.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def plot_lai(before, after, burned, burn_transform, burn_crs, summary):
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 6.2))
    lai_cmap = plt.cm.YlGn.copy()
    lai_cmap.set_bad((0, 0, 0, 0))
    lai_norm = plt.Normalize(vmin=0, vmax=2.5)
    change = after - before
    change_cmap = plt.cm.RdYlGn.copy()
    change_cmap.set_bad((0, 0, 0, 0))
    change_norm = TwoSlopeNorm(vmin=-1.5, vcenter=0, vmax=1.5)

    panels = (
        (before, lai_cmap, lai_norm, "7–14 October 2024\nleaf area before the fire"),
        (after, lai_cmap, lai_norm, "8–15 October 2026\nleaf area, same 8-day period"),
        (change, change_cmap, change_norm, "Change, 2026 minus 2024\nred is canopy loss"),
    )
    images = []
    for ax, (grid_src, cmap, norm, title) in zip(axes, panels):
        style_map(ax, title)
        grid, _ = mercator_grid(grid_src, before_transform_holder[0], before_transform_holder[1], Resampling.nearest)
        images.append(imshow_grid(ax, grid, cmap, norm))
    draw_burn_outline(axes[2], burned, burn_transform, burn_crs)

    cax_lai = fig.add_axes((0.885, 0.22, 0.012, 0.52))
    fig.colorbar(images[1], cax=cax_lai).set_label("Leaf area index")
    cax_change = fig.add_axes((0.945, 0.22, 0.012, 0.52))
    fig.colorbar(images[2], cax=cax_change).set_label("Change")

    inside = summary["inside_burn"]
    outside = summary["outside_burn"]
    fig.suptitle("MODIS leaf area (MYD15A2H Lai_500m, 500 m) on a map of Altadena", fontsize=14)
    fig.text(
        0.46,
        0.02,
        "Open areas are built-up land, where this product stores 250 and does not report leaf area. "
        "The black line is the January 2025 burn.\n"
        f"Inside the burn, leaf area {inside['lai_2024']:.2f} → {inside['lai_2026']:.2f}. "
        f"Outside the burn, {outside['lai_2024']:.2f} → {outside['lai_2026']:.2f}. "
        "Basemap: Esri World Topographic Map.",
        ha="center",
        fontsize=8,
    )
    fig.subplots_adjust(left=0.02, right=0.87, bottom=0.10, top=0.84, wspace=0.05)
    path = OUT / "altadena_lai.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


# Filled in main so plot_lai can reproject with the pre-fire grid.
before_transform_holder = [None, None]


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    burn_panels = {}
    for name, item_id in BURN_ITEMS.items():
        data, transform, meta, item = read_window("modis-64A1-061", item_id, "Burn_Date")
        burn_panels[name] = {
            "data": data,
            "transform": transform,
            "crs": meta["crs"],
            "item": item_id,
            "burned_pixels": int((data > 0).sum()),
        }
        print(f"burn {name}: {burn_panels[name]['burned_pixels']} pixels  {item['properties']['start_datetime']}")

    before_raw, before_transform, before_meta, before_item = read_window(
        "modis-15A2H-061", LAI_ITEMS["before"], "Lai_500m"
    )
    after_raw, after_transform, after_meta, after_item = read_window(
        "modis-15A2H-061", LAI_ITEMS["after"], "Lai_500m"
    )
    print("LAI 2026 meta", {k: after_meta[k] for k in ("driver", "dtype", "width", "height", "crs", "nodata")})
    before = lai_array(before_raw)
    after = lai_array(after_raw)
    before_transform_holder[0] = before_transform
    before_transform_holder[1] = before_meta["crs"]

    fire = burn_panels["fire"]
    burned_on_lai = np.zeros(before.shape, dtype="float32")
    reproject(
        np.where(fire["data"] > 0, 1, 0).astype("float32"),
        burned_on_lai,
        src_transform=fire["transform"],
        src_crs=fire["crs"],
        dst_transform=before_transform,
        dst_crs=before_meta["crs"],
        resampling=Resampling.nearest,
    )
    burned = burned_on_lai > 0.5
    summary = summarize_lai(before, after, burned)
    summary["items"] = {
        "lai_before": LAI_ITEMS["before"],
        "lai_after": LAI_ITEMS["after"],
        "lai_before_dates": [before_item["properties"]["start_datetime"], before_item["properties"]["end_datetime"]],
        "lai_after_dates": [after_item["properties"]["start_datetime"], after_item["properties"]["end_datetime"]],
        "burn": fire["item"],
    }
    print(json.dumps(summary, indent=2))

    burn_path = plot_burn(burn_panels)
    lai_path = plot_lai(before, after, burned, before_transform, before_meta["crs"], summary)
    (OUT / "altadena_lai_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {burn_path}")
    print(f"wrote {lai_path}")


if __name__ == "__main__":
    main()
