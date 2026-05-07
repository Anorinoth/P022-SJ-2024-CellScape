#!/usr/bin/env python3
"""
test_nuclear_enhancement.py — Compare nuclear preprocessing strategies for
Cellpose segmentation on a 2000×2000 center crop.

Tests multiple conditions (preprocessing + Cellpose normalize params) and
generates side-by-side overlay PNGs + a summary CSV for comparison.

Usage:
  source bioio_env/bin/activate
  python3 -u MAIN/PROG/Scripts/test_nuclear_enhancement.py
"""

import os
import sys
import time

os.environ["PYTHONUNBUFFERED"] = "1"

import numpy as np
import pandas as pd
import tifffile
import torch
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from skimage import exposure, filters, morphology, measure
from skimage.segmentation import find_boundaries

# ─── CONFIGURATION ──────────────────────────────────────────────────────────
OME_TIFF = "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/data/input/EGFR_TREG_2.ome.tif"
NUCLEAR_CHANNEL_INDEX = None  # Will be discovered from OME-XML
NUCLEAR_CHANNEL_NAME = "DNA-Brilliant Violet 421"
CROP_SIZE = 2000
OUTPUT_DIR = "/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2/MAIN/PROG/SEGMENTATION/EGFR_TREG_2/enhancement_test"

# ─── TEST CONDITIONS ────────────────────────────────────────────────────────
# Each condition: (label, preprocess_kwargs, cellpose_eval_kwargs)
#
# Levers:
#   1. denoise (Gaussian blur σ) — smooths noise, may help diffuse nuclei
#   2. sharpen (unsharp mask)    — boosts edges, but hurts diffuse nuclei
#   3. CLAHE                     — local histogram EQ, brightens dim regions
#   4. cellpose normalize dict:
#      - tile_norm_blocksize: local normalization (brightens dark areas)
#      - sharpen_radius: built-in high-pass filter
#      - percentile: [lo, hi] for normalization
#   5. cellprob_threshold: going negative captures borderline detections

CONDITIONS = [
    # --- Condition A: Current baseline (denoise=0.5, sharpen=1.0) ---
    {
        "label": "A_baseline",
        "desc": "Current: denoise=0.5, sharpen=1.0",
        "denoise_sigma": 0.5,
        "sharpen_amount": 1.0,
        "clahe": False,
        "normalize": True,
        "cellprob_threshold": 0.0,
    },
    # --- Condition B: No preprocessing (let Cellpose handle everything) ---
    {
        "label": "B_raw",
        "desc": "Raw image, Cellpose default normalize",
        "denoise_sigma": 0,
        "sharpen_amount": 0,
        "clahe": False,
        "normalize": True,
        "cellprob_threshold": 0.0,
    },
    # --- Condition C: Tile normalization (brightens dim regions) ---
    {
        "label": "C_tilenorm",
        "desc": "Tile norm (100px blocks), no preproc",
        "denoise_sigma": 0,
        "sharpen_amount": 0,
        "clahe": False,
        "normalize": {
            "normalize": True,
            "tile_norm_blocksize": 100,
        },
        "cellprob_threshold": 0.0,
    },
    # --- Condition D: CLAHE + lower cellprob ---
    {
        "label": "D_clahe",
        "desc": "CLAHE preproc, cellprob=-1.0",
        "denoise_sigma": 0.5,
        "sharpen_amount": 0,
        "clahe": True,
        "normalize": True,
        "cellprob_threshold": -1.0,
    },
    # --- Condition E: Tile normalization + lower cellprob ---
    {
        "label": "E_tilenorm_lowprob",
        "desc": "Tile norm (100px) + cellprob=-1.0",
        "denoise_sigma": 0,
        "sharpen_amount": 0,
        "clahe": False,
        "normalize": {
            "normalize": True,
            "tile_norm_blocksize": 100,
        },
        "cellprob_threshold": -1.0,
    },
    # --- Condition F: Gentle denoise + tile norm + lowered cellprob ---
    {
        "label": "F_denoise_tilenorm",
        "desc": "Denoise σ=0.5 + tile norm (100px) + cellprob=-0.5",
        "denoise_sigma": 0.5,
        "sharpen_amount": 0,
        "clahe": False,
        "normalize": {
            "normalize": True,
            "tile_norm_blocksize": 100,
        },
        "cellprob_threshold": -0.5,
    },
]


# ─── HELPERS ────────────────────────────────────────────────────────────────

def discover_nuclear_channel(ome_path, channel_name):
    """Find nuclear channel page index from OME-XML."""
    import xml.etree.ElementTree as ET
    with tifffile.TiffFile(ome_path) as tif:
        root = ET.fromstring(tif.ome_metadata)
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        for i, ch in enumerate(root.findall('.//ome:Channel', ns)):
            name = ch.get('Name', '')
            if channel_name.lower() in name.lower():
                return i
    raise ValueError(f"Nuclear channel '{channel_name}' not found")


def load_center_crop(ome_path, page_index, crop_size):
    """Load a center crop from the OME-TIFF."""
    with tifffile.TiffFile(ome_path) as tif:
        page = tif.series[0].levels[0].pages[page_index]
        full = page.asarray()
    h, w = full.shape
    cy, cx = h // 2, w // 2
    hs = crop_size // 2
    y0, x0 = max(0, cy - hs), max(0, cx - hs)
    crop = full[y0:y0+crop_size, x0:x0+crop_size]
    print(f"  Loaded crop {crop.shape} from ({y0},{x0}), range [{crop.min():.1f}, {crop.max():.1f}]")
    return crop


def preprocess(image, denoise_sigma=0, sharpen_amount=0, clahe=False):
    """Apply optional preprocessing to nuclear image."""
    img = image.astype(np.float64)
    if denoise_sigma > 0:
        img = filters.gaussian(img, sigma=denoise_sigma, preserve_range=True)
    if sharpen_amount > 0:
        img = filters.unsharp_mask(img, radius=1.0, amount=sharpen_amount,
                                   preserve_range=True)
    if clahe:
        # Normalize to [0,1] for CLAHE
        nz = img[img > 0]
        if len(nz) > 0:
            p_low, p_high = np.percentile(nz, [0.5, 99.5])
            if p_high > p_low:
                img_norm = np.clip((img - p_low) / (p_high - p_low), 0, 1)
                img_clahe = exposure.equalize_adapthist(img_norm, clip_limit=0.01)
                # Blend: 50% original + 50% CLAHE for diffuse nuclei
                img = (0.5 * img_norm + 0.5 * img_clahe)
                # Scale back to original range
                img = img * (p_high - p_low) + p_low
    img = np.clip(img, image.min(), max(image.max(), 1.0))
    return img.astype(image.dtype)


def segment(image, normalize=True, cellprob_threshold=0.0):
    """Run Cellpose segmentation."""
    from cellpose import models
    model = models.CellposeModel(gpu=True)
    masks, flows, styles = model.eval(
        image,
        diameter=None,
        flow_threshold=0.4,
        cellprob_threshold=cellprob_threshold,
        min_size=15,
        normalize=normalize,
        augment=False,
        tile_overlap=0.1,
        batch_size=256,
    )
    return masks


def make_overlay(nuclear, masks, label_text="", crop_size=2000):
    """Generate overlay image: green=nuclear, red=boundaries, yellow=centroids."""
    nuc = nuclear.astype(np.float64)
    nz = nuc[nuc > 0]
    if len(nz) > 0:
        p995 = np.percentile(nz, 99.5)
        nuc = np.clip(nuc / p995, 0, 1)
    nuc_u8 = (nuc * 255).astype(np.uint8)

    h, w = nuc_u8.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[:, :, 1] = nuc_u8  # green

    # Red boundaries
    if masks.max() > 0:
        boundaries = find_boundaries(masks, mode='outer')
        rgb[boundaries, 0] = 255
        rgb[boundaries, 1] = 0
        rgb[boundaries, 2] = 0

    img = Image.fromarray(rgb)
    draw = ImageDraw.Draw(img)

    # Yellow centroids
    if masks.max() > 0:
        props = measure.regionprops(masks)
        for p in props:
            cy, cx = p.centroid
            r = 2
            draw.ellipse([cx-r, cy-r, cx+r, cy+r],
                         fill=(255, 255, 0), outline=(255, 255, 0))

    # Label
    try:
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 28)
    except:
        font = ImageFont.load_default()

    n_cells = masks.max()
    text = f"{label_text}\n{n_cells:,} cells"
    # Black background for text
    draw.rectangle([5, 5, 500, 70], fill=(0, 0, 0, 180))
    draw.text((10, 10), text, fill=(255, 255, 255), font=font)

    return img


# ─── MAIN ───────────────────────────────────────────────────────────────────

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 70)
    print(" Nuclear Enhancement Test — EGFR_TREG_2 (2000×2000 crop)")
    print("=" * 70)

    # 1. Discover nuclear channel
    print("\n--- Discovering nuclear channel ---")
    nuc_idx = discover_nuclear_channel(OME_TIFF, NUCLEAR_CHANNEL_NAME)
    print(f"  Nuclear channel: index {nuc_idx} ({NUCLEAR_CHANNEL_NAME})")

    # 2. Load center crop
    print("\n--- Loading center crop ---")
    raw_crop = load_center_crop(OME_TIFF, nuc_idx, CROP_SIZE)

    # 3. Run each condition
    results = []
    overlay_images = []

    for i, cond in enumerate(CONDITIONS):
        label = cond["label"]
        desc = cond["desc"]
        print(f"\n{'─' * 60}")
        print(f"  [{i+1}/{len(CONDITIONS)}] {label}: {desc}")
        print(f"{'─' * 60}")

        # Preprocess
        t0 = time.time()
        processed = preprocess(
            raw_crop,
            denoise_sigma=cond.get("denoise_sigma", 0),
            sharpen_amount=cond.get("sharpen_amount", 0),
            clahe=cond.get("clahe", False),
        )
        t_preproc = time.time() - t0
        print(f"  Preprocessing: {t_preproc:.1f}s")

        # Segment
        t0 = time.time()
        masks = segment(
            processed,
            normalize=cond.get("normalize", True),
            cellprob_threshold=cond.get("cellprob_threshold", 0.0),
        )
        t_seg = time.time() - t0
        n_cells = masks.max()
        print(f"  Segmentation: {n_cells:,} cells in {t_seg:.1f}s")

        # Compute stats
        if n_cells > 0:
            props = measure.regionprops(masks)
            areas = np.array([p.area for p in props])
            median_area = np.median(areas)
            median_diam = 2 * np.sqrt(median_area / np.pi)
        else:
            median_area = 0
            median_diam = 0

        results.append({
            "condition": label,
            "description": desc,
            "n_cells": n_cells,
            "median_area_px": round(median_area, 1),
            "median_diameter_px": round(median_diam, 1),
            "preproc_time_s": round(t_preproc, 1),
            "seg_time_s": round(t_seg, 1),
        })

        # Generate overlay
        overlay = make_overlay(raw_crop, masks, label_text=f"{label}: {desc}")
        overlay_path = os.path.join(OUTPUT_DIR, f"{label}_overlay.png")
        overlay.save(overlay_path)
        overlay_images.append((label, overlay))
        print(f"  Saved: {overlay_path}")

        # GPU cleanup between runs
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    # 4. Summary table
    print(f"\n{'=' * 70}")
    print(" RESULTS SUMMARY")
    print(f"{'=' * 70}")
    df = pd.DataFrame(results)
    print(df.to_string(index=False))

    csv_path = os.path.join(OUTPUT_DIR, "enhancement_test_results.csv")
    df.to_csv(csv_path, index=False)
    print(f"\n  Summary saved: {csv_path}")

    # 5. Composite comparison image (2×3 grid)
    print("\n--- Generating comparison grid ---")
    n_conds = len(overlay_images)
    cols = 3
    rows = (n_conds + cols - 1) // cols
    thumb_size = 800

    grid_w = cols * thumb_size
    grid_h = rows * thumb_size
    grid = Image.new("RGB", (grid_w, grid_h), "black")

    for idx, (label, img) in enumerate(overlay_images):
        r, c = divmod(idx, cols)
        thumb = img.resize((thumb_size, thumb_size), Image.LANCZOS)
        grid.paste(thumb, (c * thumb_size, r * thumb_size))

    grid_path = os.path.join(OUTPUT_DIR, "comparison_grid.png")
    grid.save(grid_path)
    print(f"  Grid saved: {grid_path}")
    print(f"\n✓ All tests complete. Review overlays in:\n  {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
