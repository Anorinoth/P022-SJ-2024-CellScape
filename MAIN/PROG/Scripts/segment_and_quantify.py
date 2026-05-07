#!/usr/bin/env python3
"""
segment_and_quantify.py — GPU-Accelerated Cell Segmentation & Quantification
Author: Dr. Michael Zaiken

Complete Python pipeline replacing R's split_ome_tiff.R + robust_segment_quantify.R:
  1. Reads channels directly from pyramidal OME-TIFF via tifffile (no R/Java split)
  2. Parses channel names from OME-XML metadata
  3. Deduplicates channels and performs channel operations (e.g., Corrected-FoxP3)
  4. Generates foreground mask to skip empty background regions
  5. Segments nuclei using Cellpose v4 CP-SAM GPU (with native tiling)
  6. Applies per-channel adaptive enhancement (matching R pipeline)
  7. Quantifies mean intensity per cell per channel via regionprops
  8. Exports CSV ready for R SpatialExperiment import

Performance optimizations:
  - Direct OME-TIFF access: eliminates 30-60 min R split step
  - Float32 native: preserves original data quality (no uint16 quantization)
  - augment=False: 4x faster than test-time augmentation
  - batch_size=256: saturates 48GB A6000 GPU VRAM
  - torch.set_num_threads: all CPU cores for dynamics post-processing
  - foreground masking: skips background regions during quantification
  - --test mode: 2000x2000 center crop for parameter validation

Usage:
  python3 -u segment_and_quantify.py \\
    --image-name CAR_TREG_2 \\
    --ome-tiff data/input/CAR_TREG_2.ome.tif \\
    --nuclear-channel "DNA-Brilliant Violet 421" \\
    --output-dir MAIN/PROG \\
    --gpu

  # Test mode:
  python3 -u segment_and_quantify.py \\
    --image-name CAR_TREG_2 \\
    --ome-tiff data/input/CAR_TREG_2.ome.tif \\
    --nuclear-channel "DNA-Brilliant Violet 421" \\
    --output-dir MAIN/PROG \\
    --gpu --test
"""

import argparse
import json
import multiprocessing
import os
import re
import sys
import time
import uuid
import xml.etree.ElementTree as ET
from pathlib import Path

os.environ["PYTHONUNBUFFERED"] = "1"

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage import exposure, filters, measure, morphology


class Logger(object):
    """Reflects stdout/stderr to a file while maintaining console output."""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()

    def flush(self):
        self.terminal.flush()
        self.log.flush()


# =============================================================================
# OME-TIFF Channel Discovery
# =============================================================================
def discover_channels_ome(ome_path: str) -> dict:
    """
    Parse OME-XML metadata to discover channels, pixel size, and image dimensions.
    Returns a dict with channel info and metadata.
    """
    with tifffile.TiffFile(ome_path) as tif:
        if not tif.is_ome:
            raise ValueError(f"Not an OME-TIFF: {ome_path}")

        series = tif.series[0]
        shape = series.shape  # (C, Y, X)
        dtype = series.dtype
        axes = series.axes
        n_levels = len(series.levels)

        # Parse OME-XML
        xml_str = tif.ome_metadata
        root = ET.fromstring(xml_str)
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

        pixels = root.find('.//ome:Pixels', ns)
        px_size_x = float(pixels.get('PhysicalSizeX', '0'))
        px_size_y = float(pixels.get('PhysicalSizeY', '0'))
        px_unit = pixels.get('PhysicalSizeXUnit', 'um')

        # Extract channel names
        xml_channels = root.findall('.//ome:Channel', ns)
        channels = []
        for i, ch in enumerate(xml_channels):
            name = ch.get('Name', f'Channel_{i:02d}')
            channels.append({
                "index": i,
                "name": name,
                "page_index": i,  # page index in the TIFF
            })

    return {
        "channels": channels,
        "shape": shape,
        "dtype": str(dtype),
        "axes": axes,
        "n_levels": n_levels,
        "pixel_size_x": px_size_x,
        "pixel_size_y": px_size_y,
        "pixel_unit": px_unit,
    }


def deduplicate_channels(channels: list[dict]) -> list[dict]:
    """Keep only the first occurrence of each unique channel name."""
    seen = set()
    unique = []
    duplicates = []
    for ch in channels:
        if ch["name"] not in seen:
            seen.add(ch["name"])
            unique.append(ch)
        else:
            duplicates.append(ch)

    if duplicates:
        print(f"  Deduplication: {len(channels)} → {len(unique)} channels")
        for d in duplicates:
            print(f"    Removed duplicate: [{d['index']:2d}] {d['name']}")
    else:
        print(f"  No duplicates found ({len(channels)} channels)")
    return unique


def read_channel(ome_path: str, page_index: int,
                 crop: tuple = None) -> np.ndarray:
    """
    Read a single channel from the OME-TIFF.

    Args:
        ome_path: Path to the OME-TIFF file
        page_index: Page index (0-based channel index)
        crop: Optional (y0, x0, height, width) tuple for cropping
    """
    with tifffile.TiffFile(ome_path) as tif:
        page = tif.series[0].levels[0].pages[page_index]
        data = page.asarray()  # full channel, float32

    if crop is not None:
        y0, x0, h, w = crop
        data = data[y0:y0+h, x0:x0+w]

    return data


# =============================================================================
# Channel Operations
# =============================================================================
def normalize_to_01(data: np.ndarray, percentile: float = 99.5) -> np.ndarray:
    """
    Normalize float32 data to [0, 1] range using percentile of nonzero pixels.
    This matches the R pipeline's /65535 normalization for uint16 data.
    """
    nonzero = data[data > 0]
    if len(nonzero) == 0:
        return np.zeros_like(data, dtype=np.float64)

    p_high = np.percentile(nonzero, percentile)
    if p_high <= 0:
        return np.zeros_like(data, dtype=np.float64)

    normalized = data.astype(np.float64) / p_high
    return np.clip(normalized, 0.0, 1.0)


def compute_channel_operations(ome_path: str, channels: list[dict],
                                channel_ops: list[dict],
                                crop: tuple = None) -> list[dict]:
    """
    Perform channel operations (e.g., subtraction) and return new derived channels.

    Each op dict contains:
      - name: Name of the derived channel
      - op: Operation type ("-" for subtraction)
      - ch1: Name of the first source channel
      - ch2: Name of the second source channel
      - ch2_blur_sigma: Gaussian blur sigma for ch2 (in pixels)
      - ch2_brightness: Brightness offset for ch2 (in [0,1] space after normalization)
    """
    derived = []
    ch_name_to_page = {ch["name"]: ch["page_index"] for ch in channels}

    for op in channel_ops:
        op_name = op["name"]
        op_type = op["op"]
        ch1_name = op["ch1"]
        ch2_name = op["ch2"]

        print(f"  Computing: {op_name} = {ch1_name} {op_type} {ch2_name}")

        # Find source channels
        if ch1_name not in ch_name_to_page:
            print(f"    ⚠ Source channel '{ch1_name}' not found, skipping")
            continue
        if ch2_name not in ch_name_to_page:
            print(f"    ⚠ Source channel '{ch2_name}' not found, skipping")
            continue

        t0 = time.time()

        # Read source channels
        ch1_data = read_channel(ome_path, ch_name_to_page[ch1_name], crop=crop)
        ch2_data = read_channel(ome_path, ch_name_to_page[ch2_name], crop=crop)

        # Normalize both channels to [0, 1]
        ch1_norm = normalize_to_01(ch1_data)
        ch2_norm = normalize_to_01(ch2_data)
        del ch1_data, ch2_data

        # Apply preprocessing to ch2
        blur_sigma = op.get("ch2_blur_sigma", 0)
        brightness = op.get("ch2_brightness", 0)

        if blur_sigma > 0:
            print(f"    Blurring {ch2_name} (σ={blur_sigma} px)")
            ch2_norm = ndi.gaussian_filter(ch2_norm, sigma=blur_sigma)

        if brightness != 0:
            print(f"    Brightness offset: +{brightness}")
            ch2_norm = ch2_norm + brightness

        # Perform operation
        if op_type == "-":
            result = ch1_norm - ch2_norm
            result = np.clip(result, 0.0, 1.0)
        else:
            print(f"    ⚠ Unknown operation '{op_type}', skipping")
            continue

        elapsed = time.time() - t0
        print(f"    ✓ {op_name} computed in {elapsed:.1f}s "
              f"(range: [{result.min():.4f}, {result.max():.4f}])", flush=True)

        derived.append({
            "name": op_name,
            "data": result,
        })

    return derived


# =============================================================================
# Foreground Masking
# =============================================================================
def create_foreground_mask(nuclear_image: np.ndarray,
                           downsample: int = 8,
                           dilation_radius: int = 10) -> np.ndarray:
    """
    Create a foreground mask from the nuclear channel using Otsu thresholding.
    Returns a boolean mask at full resolution where True = tissue/foreground.
    """
    print("  Creating foreground mask...", flush=True)
    t0 = time.time()

    small = nuclear_image[::downsample, ::downsample]
    nonzero = small[small > 0]

    if len(nonzero) == 0:
        print("    ⚠ No nonzero pixels, treating entire image as foreground")
        return np.ones(nuclear_image.shape, dtype=bool)

    thresh = filters.threshold_otsu(nonzero)

    fg_small = small > thresh
    selem = morphology.disk(dilation_radius)
    fg_small = morphology.binary_dilation(fg_small, selem)
    fg_small = morphology.remove_small_holes(fg_small, area_threshold=500)

    fg_full = np.repeat(np.repeat(fg_small, downsample, axis=0),
                        downsample, axis=1)
    fg_full = fg_full[:nuclear_image.shape[0], :nuclear_image.shape[1]]

    fg_fraction = fg_small.sum() / fg_small.size
    print(f"    Tissue coverage: {fg_fraction*100:.1f}% of image")
    print(f"    Otsu threshold: {thresh:.1f}")
    print(f"    ✓ Foreground mask in {time.time()-t0:.1f}s", flush=True)

    return fg_full


# =============================================================================
# Adaptive Channel Enhancement (matches R pipeline)
# =============================================================================
def enhance_channel(channel_data: np.ndarray, gamma: float = 0.5,
                    percentile_range: tuple = (0.005, 0.995),
                    exclude_outliers: bool = True,
                    outlier_percentile: float = 0.98,
                    adaptive: bool = True) -> np.ndarray:
    """
    Adaptive channel enhancement matching the R enhance_channel function:
      1. Percentile normalization to [0, 1]
      2. Outlier masking (replace top percentile with median)
      3. CLAHE blend (70% original + 30% equalized)
      4. Gamma correction
    """
    data = channel_data.astype(np.float64)

    nonzero_mask = data > 0
    nonzero_pixels = data[nonzero_mask]

    if len(nonzero_pixels) == 0:
        return data

    # 1. Percentile normalization
    q_lower = np.percentile(nonzero_pixels, percentile_range[0] * 100)
    q_upper = np.percentile(nonzero_pixels, percentile_range[1] * 100)

    # 2. Handle bright outliers
    outlier_mask = None
    if exclude_outliers:
        outlier_threshold = np.percentile(nonzero_pixels, outlier_percentile * 100)
        outlier_mask = data > outlier_threshold
        if np.any(outlier_mask):
            q_upper = outlier_threshold

    denom = q_upper - q_lower
    if denom <= 0:
        denom = 1.0
    enhanced = (data - q_lower) / denom
    enhanced = np.clip(enhanced, 0.0, 1.0)

    if outlier_mask is not None and np.any(outlier_mask):
        median_val = np.median(enhanced[enhanced > 0]) if np.any(enhanced > 0) else 0.0
        enhanced[outlier_mask] = median_val

    # 3. CLAHE blend
    if adaptive:
        enhanced_clahe = exposure.equalize_adapthist(enhanced, clip_limit=0.005)
        enhanced = 0.7 * enhanced + 0.3 * enhanced_clahe

    # 4. Gamma correction
    enhanced = np.power(enhanced, gamma)

    return enhanced


# =============================================================================
# Cellpose Segmentation
# =============================================================================
def segment_nuclei(nuclear_image: np.ndarray, use_gpu: bool = True,
                   diameter: float = None, flow_threshold: float = 0.4,
                   cellprob_threshold: float = 0.0,
                   min_size: int = 15,
                   batch_size: int = 256,
                   n_threads: int = None,
                   verbose: bool = False,
                   tile_norm_blocksize: int = 0) -> np.ndarray:
    """
    Segment nuclei using Cellpose v4 (CP-SAM) with GPU acceleration.
    augment=False for fast mode (no test-time augmentation).

    Args:
        tile_norm_blocksize: If > 0, enable Cellpose tile normalization with this
            block size (in pixels). This locally normalizes intensity within each
            tile, brightening dim/diffuse nuclei that global normalization misses.
            Recommended value: 100. Set to 0 to use default global normalization.
    """
    import torch
    from cellpose import models, io

    if verbose:
        io.logger_setup()

    # --- Auto-disable flow QC for very large images (>100M pixels) -----------
    # Cellpose's remove_bad_flow_masks allocates a neighbors tensor of shape
    # (2, 9, N_mask_pixels) on GPU that can exceed VRAM for large images.
    # Cellpose's own internal threshold is 100M pixels (dynamics.py:419) and
    # its source recommends: "turn off QC step with flow_threshold=0 if too slow"
    LARGE_IMAGE_THRESHOLD = 100_000_000  # 100M pixels
    n_pixels = nuclear_image.shape[0] * nuclear_image.shape[1]
    if n_pixels > LARGE_IMAGE_THRESHOLD and flow_threshold > 0:
        print(f"  ⚠️  Large image ({n_pixels/1e6:.0f}M px > 100M threshold): "
              f"disabling flow_threshold QC", flush=True)
        print(f"       flow_threshold: {flow_threshold} → 0 "
              f"(GPU OOM prevention for flow validation step)", flush=True)
        flow_threshold = 0

    if n_threads is None:
        n_threads = min(os.cpu_count(), 64)
    torch.set_num_threads(n_threads)
    print(f"  CPU threads for dynamics: {n_threads}", flush=True)

    print(f"  Initializing Cellpose v4 CellposeModel (GPU={use_gpu})...", flush=True)
    model = models.CellposeModel(gpu=use_gpu)

    # Build normalize parameter
    if tile_norm_blocksize > 0:
        normalize_param = {
            "normalize": True,
            "tile_norm_blocksize": tile_norm_blocksize,
        }
        print(f"  Tile normalization: enabled (blocksize={tile_norm_blocksize}px)", flush=True)
    else:
        normalize_param = True

    print(f"  Running segmentation on {nuclear_image.shape} image...", flush=True)
    print(f"    diameter={diameter}, flow_threshold={flow_threshold}, "
          f"cellprob_threshold={cellprob_threshold}, min_size={min_size}", flush=True)
    print(f"    batch_size={batch_size}, augment=False (fast mode)", flush=True)

    t0 = time.time()

    masks, flows, styles = model.eval(
        nuclear_image,
        diameter=diameter,
        flow_threshold=flow_threshold,
        cellprob_threshold=cellprob_threshold,
        min_size=min_size,
        normalize=normalize_param,
        augment=False,
        tile_overlap=0.1,
        batch_size=batch_size,
    )

    elapsed = time.time() - t0
    n_cells = masks.max()
    print(f"  ✓ Segmentation: {n_cells:,} cells in {elapsed:.1f}s "
          f"({elapsed/60:.1f} min)", flush=True)

    # --- GPU memory diagnostics & cleanup ------------------------------------
    if use_gpu and torch.cuda.is_available():
        peak_mb = torch.cuda.max_memory_allocated() / 1e6
        curr_mb = torch.cuda.memory_allocated() / 1e6
        free_mb, total_mb = torch.cuda.mem_get_info()
        free_mb /= 1e6; total_mb /= 1e6
        print(f"  GPU VRAM: {curr_mb:.0f} MB allocated, {peak_mb:.0f} MB peak, "
              f"{free_mb:.0f}/{total_mb:.0f} MB free/total", flush=True)
        torch.cuda.empty_cache()
        torch.cuda.reset_peak_memory_stats()
        print(f"  GPU cache cleared for downstream steps", flush=True)

    return masks


# =============================================================================
# Test Mode: Segmentation Overlay
# =============================================================================
def generate_overlay(nuclear_image: np.ndarray, masks: np.ndarray,
                     output_path: str) -> None:
    """
    Generate a nuclear channel + mask overlay PNG for visual QC.
    Nuclear signal = green, cell boundaries = red, centroids = yellow.
    """
    from PIL import Image, ImageDraw
    from skimage.segmentation import find_boundaries

    print(f"  Generating segmentation overlay...", flush=True)
    t0 = time.time()

    # Normalize nuclear image to [0, 255] for display
    nuc = nuclear_image.astype(np.float64)
    nz = nuc[nuc > 0]
    if len(nz) > 0:
        p995 = np.percentile(nz, 99.5)
        nuc = np.clip(nuc / p995, 0, 1)
    nuc_u8 = (nuc * 255).astype(np.uint8)

    # RGB: nuclear = green channel
    h, w = nuc_u8.shape
    rgb = np.zeros((h, w, 3), dtype=np.uint8)
    rgb[:, :, 1] = nuc_u8  # green

    # Red boundaries
    boundaries = find_boundaries(masks, mode='outer')
    rgb[boundaries, 0] = 255
    rgb[boundaries, 1] = 0
    rgb[boundaries, 2] = 0

    # Yellow centroids
    img = Image.fromarray(rgb)
    draw = ImageDraw.Draw(img)
    props = measure.regionprops(masks)
    for p in props:
        cy, cx = p.centroid
        r = 2
        draw.ellipse([cx - r, cy - r, cx + r, cy + r],
                     fill=(255, 255, 0), outline=(255, 255, 0))

    img.save(output_path)
    print(f"    ✓ Overlay saved: {output_path} ({time.time()-t0:.1f}s)", flush=True)


# =============================================================================
# Quantification
# =============================================================================
def _quantify_single_channel_worker(args):
    """Worker function for parallel quantification."""
    (ch_name, ch_data, masks, labels, enhance,
     enhancement_params, return_enhanced) = args

    if enhance:
        ch_enhanced = enhance_channel(ch_data, **enhancement_params)
    else:
        # Use a copy to avoid modifying original or issues with shared mem
        ch_enhanced = ch_data.astype(np.float64, copy=True)

    mean_intensities = ndi.mean(ch_enhanced, labels=masks, index=labels)

    enhanced_result = None
    if return_enhanced:
        enhanced_result = ch_enhanced.astype(np.float32)

    return ch_name, mean_intensities, enhanced_result


def quantify_cells(masks: np.ndarray,
                   channel_data_list: list[tuple[str, np.ndarray]],
                   enhance: bool = True,
                   enhancement_params: dict = None,
                   foreground_mask: np.ndarray = None,
                   return_enhanced: bool = False,
                   n_workers: int = 1) -> tuple[pd.DataFrame, dict]:
    """
    Measure per-cell mean intensity for each channel.
    channel_data_list: list of (channel_name, image_array) tuples.
    """
    if enhancement_params is None:
        enhancement_params = {
            "gamma": 0.5,
            "percentile_range": (0.005, 0.995),
            "exclude_outliers": True,
            "outlier_percentile": 0.98,
            "adaptive": True,
        }

    n_cells = masks.max()
    if n_cells == 0:
        return pd.DataFrame()

    print(f"\n=== Quantification: {n_cells:,} cells × "
          f"{len(channel_data_list)} channels ===", flush=True)

    enhanced_images = {}

    # Cell centroids
    print("  Computing cell centroids...", flush=True)
    t0 = time.time()
    props = measure.regionprops(masks)
    centroids = {p.label: p.centroid for p in props}
    areas = {p.label: p.area for p in props}
    cell_ids = sorted(centroids.keys())

    # Filter by foreground mask
    if foreground_mask is not None:
        fg_cells = []
        for cid in cell_ids:
            cy, cx = centroids[cid]
            iy, ix = int(round(cy)), int(round(cx))
            if 0 <= iy < foreground_mask.shape[0] and 0 <= ix < foreground_mask.shape[1]:
                if foreground_mask[iy, ix]:
                    fg_cells.append(cid)
            else:
                fg_cells.append(cid)
        n_filtered = len(cell_ids) - len(fg_cells)
        if n_filtered > 0:
            print(f"    Filtered {n_filtered:,} cells outside foreground", flush=True)
        cell_ids = fg_cells

    rows = [{"cell_id": cid, "x": centroids[cid][1], "y": centroids[cid][0],
             "area": areas[cid]} for cid in cell_ids]
    df = pd.DataFrame(rows)
    print(f"    ✓ {len(df):,} cells in {time.time()-t0:.1f}s", flush=True)

    # Quantify channels
    labels = np.array(cell_ids)
    enhanced_images = {}

    if n_workers > 1:
        print(f"  Parallel quantification with {n_workers} workers...", flush=True)
        t_start = time.time()

        # Prepare arguments for workers
        worker_args = []
        for ch_name, ch_data in channel_data_list:
            worker_args.append((
                ch_name, ch_data, masks, labels,
                enhance, enhancement_params, return_enhanced
            ))

        # Use multiprocessing Pool
        # On Linux, this defaults to 'fork' which efficiently shares 'masks' and 'labels'
        with multiprocessing.Pool(processes=n_workers) as pool:
            # imap_unordered for better progress visibility
            results = pool.imap_unordered(_quantify_single_channel_worker, worker_args)

            for i, (ch_name, mean_vals, enhanced_img) in enumerate(results):
                df[f"ch_{ch_name}"] = mean_vals
                if return_enhanced and enhanced_img is not None:
                    enhanced_images[ch_name] = enhanced_img
                print(f"    [{i+1:2d}/{len(channel_data_list)}] {ch_name} (parallel)", flush=True)

        print(f"  ✓ Parallel quantification complete in {time.time()-t_start:.1f}s", flush=True)

    else:
        # Serial quantification
        for i, (ch_name, ch_data) in enumerate(channel_data_list):
            t_ch = time.time()
            if enhance:
                ch_enhanced = enhance_channel(ch_data, **enhancement_params)
            else:
                ch_enhanced = ch_data.astype(np.float64)

            mean_intensities = ndi.mean(ch_enhanced, labels=masks, index=labels)
            df[f"ch_{ch_name}"] = mean_intensities

            if return_enhanced:
                enhanced_images[ch_name] = ch_enhanced.astype(np.float32)

            print(f"    [{i+1:2d}/{len(channel_data_list)}] {ch_name} {time.time()-t_ch:.1f}s", flush=True)

    return df, enhanced_images


# =============================================================================
# Run Diagnostics report
# =============================================================================
def generate_diagnostics(args, meta: dict, df: pd.DataFrame,
                         output_path: str) -> None:
    """
    Generate a human-friendly Markdown report of the pipeline run.
    """
    print(f"  Generating diagnostics report...", flush=True)
    t0 = time.time()

    n_cells = len(df)
    areas = df["area"].values
    mean_area = np.mean(areas)
    median_area = np.median(areas)
    
    # Calculate diameters from area = pi * (d/2)^2  => d = 2 * sqrt(area/pi)
    median_diam_px = 2 * np.sqrt(median_area / np.pi)
    px_size = meta.get("pixel_size_x", 0.182)
    median_diam_um = median_diam_px * px_size

    # Area distribution
    bins = [0, 100, 200, 500, 1000, 2000, 5000, 10000, 50000]
    dist_rows = []
    for i in range(len(bins)-1):
        count = int(np.sum((areas >= bins[i]) & (areas < bins[i+1])))
        pct = (count / n_cells * 100) if n_cells > 0 else 0
        dist_rows.append(f"| {bins[i]:6d}-{bins[i+1]:6d} px | {count:6d} | {pct:5.1f}% |")

    # Image stats
    h, w = meta["shape"][-2:]
    pixels_m = (h * w) / 1e6
    density = n_cells / pixels_m if not args.test else (n_cells / (args.test_size**2 / 1e6))

    full_est_str = ""
    if args.test:
        full_est = density * pixels_m
        full_est_str = f"| *Full Image Estimate* | *~{full_est:,.0f} cells* |"

    md = f"""# Pipeline Diagnostics: {args.image_name}
Run at: {time.strftime('%Y-%m-%d %H:%M:%S')}
Mode: {'TEST (center crop)' if args.test else 'FULL'}

## Cell Statistics
| Metric | Value |
|---|---|
| **Total Cells** | **{n_cells:,}** |
| Median Area | {median_area:.1f} px |
| Mean Area | {mean_area:.1f} px |
| **Median Diameter** | **{median_diam_px:.1f} px** |
| **Physical Diameter** | **{median_diam_um:.2f} µm** |
| Min / Max Area | {areas.min():.0f} / {areas.max():.0f} px |

## Area Distribution
| Area Range | Count | % |
|---|---|---|
{chr(10).join(dist_rows)}

## Image & Density
| Metric | Value |
|---|---|
| Image Dimensions | {h} x {w} px |
| Pixel Size | {px_size} {meta.get('pixel_unit', 'µm')}/px |
| Cell Density | {density:.1f} cells/M px |
{full_est_str}

## Channel Operations
- Corrected-FoxP3: GFP-FITC - (blurred CD19-FITC + 0.05) [Normalized float32]
"""
    with open(output_path, "w") as f:
        f.write(md)
    print(f"    ✓ Diagnostics saved: {output_path} ({time.time()-t0:.1f}s)", flush=True)


def export_enhanced_ome_tiff(ome_path: str,
                             channels: list[dict],
                             derived_channels: list[dict],
                             enhancement_params: dict,
                             output_path: str,
                             tile_size: int = 512,
                             pixel_size_um: float = 0.182,
                             enhanced_data_dict: dict = None,
                             crop: tuple = None) -> None:
    """
    Export enhanced channels to a tiled multi-channel OME-TIFF.
    Version 2.1: Pure Automated imwrite (Variant 4 logic).
    """
    print(f"\n--- Writing Multi-Channel OME-TIFF (v2.1/Pure Automated) ---", flush=True)
    t0 = time.time()
    
    # Ensure extension is .ome.tif and clean up old .tiff to avoid confusion
    if output_path.endswith(".ome.tiff"):
        old_path = output_path
        output_path = output_path[:-5] + ".tif"
        if os.path.exists(old_path) and old_path != output_path:
            print(f"  Cleaning up legacy file: {os.path.basename(old_path)}")
            try: os.remove(old_path)
            except: pass

    # 1. Identify and order all channels (Maintain Original Discovery Order)
    original_names = [ch["name"] for ch in channels]
    final_ordered_names = []
    
    # Add originals first (matching discovery order)
    for name in original_names:
        if (enhanced_data_dict and name in enhanced_data_dict) or any(ch["name"] == name for ch in channels):
             if name not in final_ordered_names:
                 final_ordered_names.append(name)

    # Add any derived channels (like Corrected-FoxP3) at the end
    if enhanced_data_dict:
        # Use sorted keys for deterministic derived ordering
        for name in sorted(enhanced_data_dict.keys()):
            if name not in final_ordered_names:
                final_ordered_names.append(name)
    else:
        for dch in derived_channels:
            if dch["name"] not in final_ordered_names:
                final_ordered_names.append(dch["name"])

    total_channels = len(final_ordered_names)
    print(f"  Total Channels to Export: {total_channels}")
    
    # 2. Assemble stack (Variant 4 requires a 3D float32 stack for imwrite)
    if enhanced_data_dict:
        sample = next(iter(enhanced_data_dict.values()))
    else:
        sample = read_channel(ome_path, channels[0]["page_index"], crop=crop)
    height, width = sample.shape
    dtype = np.float32

    print(f"  Assembling stack ({total_channels}, {height}, {width})...", flush=True)
    stack = np.empty((total_channels, height, width), dtype=dtype)
    
    for i, name in enumerate(final_ordered_names):
        t_ch = time.time()
        # Source from dict (Resumability path)
        if enhanced_data_dict and name in enhanced_data_dict:
            img_data = enhanced_data_dict[name]
        # Source from derived channels list
        elif any(dch["name"] == name for dch in derived_channels):
            dch_match = next(dch for dch in derived_channels if dch["name"] == name)
            img_data = dch_match["data"]
        # Source from original file (Quantification path)
        else:
            ch_match = next(ch for ch in channels if ch["name"] == name)
            img = read_channel(ome_path, ch_match["page_index"], crop=crop)
            img_data = enhance_channel(img, **enhancement_params).astype(dtype)
            del img
            
        stack[i] = img_data
        if (i + 1) % 10 == 0 or (i + 1) == total_channels:
            print(f"    Added {i+1}/{total_channels} [{name}] ({time.time()-t_ch:.1f}s)", flush=True)

    # 3. Build Metadata Dictionary (Matching Variant 4)
    # Note: Use simple numeric values for PhysicalSize to match working harness
    metadata = {
        'axes': 'CYX',
        'PhysicalSizeX': pixel_size_um,
        'PhysicalSizeY': pixel_size_um,
        'Channel': {'Name': final_ordered_names}
    }

    # 4. Write via imwrite (The most compatible way for multi-channel)
    print(f"  Writing OME-TIFF via imwrite...", flush=True)
    tifffile.imwrite(
        output_path,
        stack,
        metadata=metadata,
        ome=True,
        bigtiff=True,  # Enable BigTIFF for files > 4GB
        tile=(tile_size, tile_size),
        compression="lzw"
    )

    file_size_gb = os.path.getsize(output_path) / 1e9
    print(f"  ✓ Enhanced OME-TIFF (v2.1): {file_size_gb:.1f} GB in {time.time()-t0:.1f}s")
    print(f"    {output_path}", flush=True)
    
    del stack # Explicit cleanup




# =============================================================================
# Test Mode: Center Crop
# =============================================================================
def get_crop_params(height: int, width: int, crop_size: int) -> tuple:
    """Returns (y0, x0, h, w) for center crop."""
    cy, cx = height // 2, width // 2
    hs = crop_size // 2
    y0 = max(0, cy - hs)
    x0 = max(0, cx - hs)
    h = min(crop_size, height - y0)
    w = min(crop_size, width - x0)
    return (y0, x0, h, w)


# =============================================================================
# Main Pipeline
# =============================================================================
def run_pipeline(args):
    """Main pipeline: OME-TIFF → segmentation → quantification → CSV."""

    overall_start = time.time()
    mode_str = "TEST MODE" if args.test else "FULL"
    print(f"\n{'=' * 70}")
    print(f" Segment & Quantify: {args.image_name} [{mode_str}]")
    print(f" Source: {args.ome_tiff}")
    print(f"{'=' * 70}", flush=True)

    # -------------------------------------------------------------------------
    # Step 1: Parse OME-TIFF metadata and discover channels
    # -------------------------------------------------------------------------
    print(f"\n--- Step 1: OME-TIFF Channel Discovery ---")

    meta = discover_channels_ome(args.ome_tiff)
    all_channels = meta["channels"]
    shape = meta["shape"]
    n_total = len(all_channels)

    print(f"  Image: {shape} ({meta['dtype']}, {meta['axes']})")
    print(f"  Pixel size: {meta['pixel_size_x']} {meta['pixel_unit']}/px")
    print(f"  Pyramid levels: {meta['n_levels']}")
    print(f"  Channels: {n_total}")

    # Dimensions from shape (CYX)
    if meta["axes"] == "CYX":
        n_ch, height, width = shape
    elif meta["axes"] == "YXC":
        height, width, n_ch = shape
    else:
        raise ValueError(f"Unexpected axes order: {meta['axes']}")

    # Deduplicate
    channels = deduplicate_channels(all_channels)
    n_unique = len(channels)

    # Find nuclear channel
    nuclear_ch = None
    for ch in channels:
        if ch["name"] == args.nuclear_channel:
            nuclear_ch = ch
            break
    if nuclear_ch is None:
        # Try fuzzy match
        for ch in channels:
            if args.nuclear_channel.lower() in ch["name"].lower():
                nuclear_ch = ch
                break
    if nuclear_ch is None:
        print(f"\n  ERROR: Nuclear channel '{args.nuclear_channel}' not found.")
        print(f"  Available channels:")
        for ch in channels:
            print(f"    [{ch['index']:2d}] {ch['name']}")
        sys.exit(1)

    print(f"  Nuclear channel: [{nuclear_ch['index']}] {nuclear_ch['name']}")
    print(f"  Unique channels: {n_unique}", flush=True)

    # Compute crop parameters if test mode
    crop = None
    crop_offset = (0, 0)
    if args.test:
        crop = get_crop_params(height, width, args.test_size)
        crop_offset = (crop[0], crop[1])
        print(f"\n  🧪 TEST MODE: {args.test_size}×{args.test_size} center crop")
        print(f"     Offset: y={crop[0]}, x={crop[1]}", flush=True)

    # -------------------------------------------------------------------------
    # Step 1.5: Identify Cache Paths and Check for Early Loading
    # -------------------------------------------------------------------------
    suffix = "_test" if args.test else ""
    mask_file = os.path.join(args.output_dir, f"{args.image_name}{suffix}_masks.npz")
    enhanced_cache = os.path.join(args.output_dir, f"{args.image_name}{suffix}_enhanced.npz")
    
    enhanced_dict = None
    skip_preprocessing = False
    
    if args.cache_enhanced and os.path.exists(enhanced_cache):
        print(f"\n--- Early Cache Check ---")
        print(f"  ♻️  Loading enhanced channels from cache: {enhanced_cache}")
        try:
            with np.load(enhanced_cache) as loader:
                enhanced_dict = {k: loader[k] for k in loader.files}
            print(f"    ✓ Loaded {len(enhanced_dict)} channels")
            
            # Safety: If we have enhanced cache but NO masks, we still need to segment.
            # Segmentation REQUIRES the nuclear image.
            if not os.path.exists(mask_file):
                print(f"    ⚠️ Masks missing ({mask_file}). Re-loading nuclear for segmentation.")
                enhanced_dict = None # Reset to ensure full load
            else:
                skip_preprocessing = True
                print(f"    ✓ Resuming from enhanced cache. Skipping early steps.")
        except Exception as e:
            print(f"    ⚠️ Could not load cache: {e}. Proceeding with full processing.")

    # -------------------------------------------------------------------------
    # Step 2: Load Nuclear Channel (Skip if pre-enhanced cache found)
    # -------------------------------------------------------------------------
    nuclear_image = None
    if not skip_preprocessing:
        print(f"\n--- Step 2: Load Nuclear Channel ---")
        t0 = time.time()
        nuclear_image = read_channel(args.ome_tiff, nuclear_ch["index"], crop=crop)
        
        print(f"  Shape: {nuclear_image.shape}, dtype={nuclear_image.dtype}")
        print(f"  Range: [{nuclear_image.min():.2f}, {nuclear_image.max():.2f}]")
        print(f"  RAM: ~{nuclear_image.nbytes/1e9:.2f} GB")
        print(f"  Loaded in {time.time()-t0:.1f}s")

        # Optional: Nuclear Enhancement (before foreground/segmentation)
        if args.denoise_nuclear > 0 or args.sharpen_nuclear > 0:
            print(f"\n--- Step 2a: Nuclear Enhancement ---")
            t_enh = time.time()
            orig_max = nuclear_image.max()
            orig_min = nuclear_image.min()
            
            # Work in float for precision
            nuclear_f64 = nuclear_image.astype(np.float64)
            
            if args.denoise_nuclear > 0:
                print(f"  Denoising (Gaussian, σ={args.denoise_nuclear})...")
                nuclear_f64 = filters.gaussian(nuclear_f64, sigma=args.denoise_nuclear, 
                                            preserve_range=True)
                
            if args.sharpen_nuclear > 0:
                print(f"  Sharpening (Unsharp Mask, strength={args.sharpen_nuclear})...")
                nuclear_f64 = filters.unsharp_mask(nuclear_f64, radius=1.0, 
                                                amount=args.sharpen_nuclear,
                                                preserve_range=True)
            
            nuclear_f64 = np.clip(nuclear_f64, orig_min, max(orig_max, 1.0))
            nuclear_image = nuclear_f64.astype(nuclear_image.dtype)
            print(f"  ✓ Enhanced in {time.time()-t_enh:.1f}s")
    else:
        print(f"\n--- Step 2: Load Nuclear Channel (SKIPPED) ---")

    # -------------------------------------------------------------------------
    # Step 2b: Foreground Detection (Skip if pre-enhanced cache found)
    # -------------------------------------------------------------------------
    foreground_mask = None
    if not skip_preprocessing:
        print(f"\n--- Step 2b: Foreground Detection ---")
        foreground_mask = create_foreground_mask(nuclear_image)
    else:
        print(f"\n--- Step 2b: Foreground Detection (SKIPPED) ---")

    # -------------------------------------------------------------------------
    # Step 3: Cellpose Segmentation
    # -------------------------------------------------------------------------
    suffix = "_test" if args.test else ""
    mask_file = os.path.join(args.output_dir, f"{args.image_name}{suffix}_masks.npz")

    # Step 3: Check for existing mask (Resumability)
    if os.path.exists(mask_file):
        print(f"  ♻️  Loading existing masks from: {mask_file}")
        with np.load(mask_file) as loader:
            masks = loader['masks']
        n_cells = masks.max()
        print(f"    ✓ Loaded {n_cells:,} cells")
    else:
        masks = segment_nuclei(
            nuclear_image,
            use_gpu=args.gpu,
            diameter=args.diameter,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
            min_size=args.min_size,
            batch_size=args.batch_size,
            n_threads=args.threads,
            verbose=args.verbose,
            tile_norm_blocksize=args.tile_norm_blocksize,
        )
        np.savez_compressed(mask_file, masks=masks)
        print(f"  💾 Masks saved: {mask_file} ({os.path.getsize(mask_file)/1e6:.0f} MB)",
              flush=True)

    # Generate overlay in test mode (only if nuclear image was loaded)
    if args.test and nuclear_image is not None:
        overlay_path = os.path.join(
            args.output_dir, f"{args.image_name}{suffix}_overlay.png")
        generate_overlay(nuclear_image, masks, overlay_path)

    del nuclear_image

    # Channel operations definition (required for metadata/logging)
    channel_ops = [
        {
            "name": "Corrected-FoxP3",
            "op": "-",
            "ch1": "GFP-FITC",
            "ch2": "CD19-FITC",
            "ch2_blur_sigma": 1.5,
            "ch2_brightness": 0.05,
        },
    ]

    # -------------------------------------------------------------------------
    # Step 4: Channel Operations (Skip if pre-enhanced cache found)
    # -------------------------------------------------------------------------
    derived_channels = []
    if not skip_preprocessing:
        print(f"\n--- Step 4: Channel Operations ---")
        derived_channels = compute_channel_operations(
            args.ome_tiff, channels, channel_ops, crop=crop
        )
        if derived_channels:
            print(f"  {len(derived_channels)} derived channel(s) computed", flush=True)
    else:
        print(f"\n--- Step 4: Channel Operations (SKIPPED) ---")

    # -------------------------------------------------------------------------
    # Step 5: Load channels and quantify
    # -------------------------------------------------------------------------
    print(f"\n--- Step 5: Channel Enhancement & Quantification ---")
    
    channel_data_list = []
    
    if not skip_preprocessing:
        print(f"  Loading {n_unique} channels + {len(derived_channels)} derived...", flush=True)
        for i, ch in enumerate(channels):
            t_load = time.time()
            ch_data = read_channel(args.ome_tiff, ch["page_index"], crop=crop)
            channel_data_list.append((ch["name"], ch_data))
            if (i + 1) % 5 == 0 or i == n_unique - 1:
                print(f"    Loaded {i+1}/{n_unique} ({time.time()-t_load:.1f}s)", flush=True)
        for dc in derived_channels:
            channel_data_list.append((dc["name"], dc["data"]))
    else:
        print(f"  ♻️  Using {len(enhanced_dict)} channels from enhanced cache.", flush=True)
        # Sort keys to maintain deterministic ordering
        for name in sorted(enhanced_dict.keys()):
            channel_data_list.append((name, enhanced_dict[name]))
    
    total_channels = len(channel_data_list)

    # Enhancement parameters
    enhancement_params = {
        "gamma": args.enhancement_gamma,
        "percentile_range": (args.enhancement_perc_lower, args.enhancement_perc_upper),
        "exclude_outliers": args.enhancement_exclude_outliers,
        "outlier_percentile": args.enhancement_outlier_percentile,
        "adaptive": args.enhancement_adaptive,
    }

    # Determine if OME-TIFF export is needed (must be computed BEFORE quantify
    # to decide whether to retain enhanced images in memory)
    should_export = args.export_ome if args.export_ome is not None else (not args.test)

    # Quantify
    n_quant_workers = args.threads or min(os.cpu_count() // 2, 32)
    # If skip_preprocessing is True, enhance=False because images are already enhanced.
    cells_df, enhanced_dict_new = quantify_cells(
        masks, channel_data_list,
        enhance=(not skip_preprocessing),
        enhancement_params=enhancement_params,
        foreground_mask=foreground_mask,
        return_enhanced=(not args.test or args.cache_enhanced or should_export),
        n_workers=n_quant_workers
    )
    
    # If we didn't skip, use the new dict. If we did skip, enhanced_dict set at top is fine.
    if not skip_preprocessing:
        enhanced_dict = enhanced_dict_new

    # Free initial memory
    del channel_data_list, masks, foreground_mask

    if cells_df.empty:
        print("\n  ⚠️  No cells detected! Check parameters.", flush=True)
        sys.exit(1)

    # Adjust coordinates for test crop
    if args.test and crop_offset != (0, 0):
        cells_df["x"] += crop_offset[1]
        cells_df["y"] += crop_offset[0]

    # -------------------------------------------------------------------------
    # Step 6: Export CSV & Diagnostics (Before OME-TIFF to ensure results save)
    # -------------------------------------------------------------------------
    print(f"\n--- Step 6: Export (Core Results) ---", flush=True)

    cells_df["image_name"] = args.image_name

    csv_file = os.path.join(args.output_dir, f"{args.image_name}{suffix}_cells.csv")
    cells_df.to_csv(csv_file, index=False)
    print(f"  💾 CSV: {csv_file}")
    print(f"     {len(cells_df):,} cells × {len(cells_df.columns)} columns")
    print(f"     File size: {os.path.getsize(csv_file)/1e6:.1f} MB", flush=True)

    # Step 7: Diagnostics
    print(f"\n--- Step 7: Diagnostics ---")
    diag_file = os.path.join(args.output_dir, f"{args.image_name}{suffix}_diagnostics.md")
    generate_diagnostics(args, meta, cells_df, diag_file)

    # Step 8: Export Enhanced OME-TIFF (Graceful Failure)
    # should_export was computed earlier (before quantification) to determine
    # whether to retain enhanced images in memory
    
    if should_export:
        print(f"\n--- Step 8: OME-TIFF Export ---")
        
        # Enhanced Data Caching (Resumability for export - OPTIONAL)
        enhanced_cache = os.path.join(args.output_dir, f"{args.image_name}{suffix}_enhanced.npz")
        
        if args.cache_enhanced:
            if not skip_preprocessing and enhanced_dict:
                print(f"  💾 Saving enhanced channels to cache: {enhanced_cache}")
                try:
                    np.savez_compressed(enhanced_cache, **enhanced_dict)
                    print(f"    ✓ Cache saved ({os.path.getsize(enhanced_cache)/1e6:.1f} MB)")
                except Exception as e:
                    print(f"    ⚠️ Could not save cache: {e}")
            else:
                print(f"  ♻️  Using channels already loaded from cache.")

        try:
            ome_out_path = os.path.join(
                args.output_dir,
                f"{args.image_name}{suffix}_enhanced_multichannel.ome.tiff"
            )
            export_enhanced_ome_tiff(
                args.ome_tiff,
                channels,
                derived_channels,
                enhancement_params,
                ome_out_path,
                pixel_size_um=meta["pixel_size_x"],
                enhanced_data_dict=enhanced_dict,
                crop=crop
            )
        except Exception as e:
            print(f"\n  ⚠️ ERROR writing OME-TIFF: {str(e)}")
            print(f"     Quantification results are already saved to CSV.")

    # Free memory
    del enhanced_dict, derived_channels

    # Save parameters
    ch_names_final = [ch["name"] for ch in channels]
    ch_names_final += [dc for dc in
                       [op["name"] for op in channel_ops] if dc not in ch_names_final]

    params = {
        "image_name": args.image_name,
        "ome_tiff": args.ome_tiff,
        "pixel_size_um": meta["pixel_size_x"],
        "image_dims": [int(height), int(width)],
        "n_channels_original": n_total,
        "n_channels_unique": n_unique,
        "n_channels_total": total_channels,
        "channel_names": ch_names_final,
        "test_mode": args.test,
        "test_size": args.test_size if args.test else None,
        "segmentation": {
            "model": "cellpose_v4_cpsam",
            "gpu": args.gpu,
            "diameter": args.diameter,
            "flow_threshold": args.flow_threshold,
            "cellprob_threshold": args.cellprob_threshold,
            "min_size": args.min_size,
            "batch_size": args.batch_size,
            "augment": False,
            "cpu_threads": args.threads or min(os.cpu_count(), 64),
        },
        "channel_operations": channel_ops,
        "enhancement": enhancement_params,
        "n_cells": len(cells_df),
    }

    params_file = os.path.join(
        args.output_dir, f"{args.image_name}{suffix}_pipeline_params.json")
    with open(params_file, "w") as f:
        json.dump(params, f, indent=2, default=str)
    print(f"  💾 Params: {params_file}", flush=True)

    total_elapsed = time.time() - overall_start
    print(f"\n{'=' * 70}")
    print(f" ✓ PIPELINE COMPLETE: {args.image_name} [{mode_str}]")
    print(f"   Total time: {total_elapsed/60:.1f} min ({total_elapsed:.0f}s)")
    print(f"   Cells: {len(cells_df):,}")
    print(f"   Channels: {total_channels}")
    print(f"   Output: {csv_file}")
    print(f"{'=' * 70}\n", flush=True)


# =============================================================================
# CLI
# =============================================================================
def main():
    parser = argparse.ArgumentParser(
        description="GPU-accelerated cell segmentation & quantification for IMC",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    # Required
    parser.add_argument("--image-name", required=True,
                        help="Image identifier (e.g., CAR_TREG_2)")
    parser.add_argument("--ome-tiff", required=True,
                        help="Path to OME-TIFF file")
    parser.add_argument("--nuclear-channel", required=True,
                        help="Name of the nuclear channel "
                             '(e.g., "DNA-Brilliant Violet 421")')
    parser.add_argument("--output-dir", required=True,
                        help="Directory for output files")

    # Test mode
    parser.add_argument("--test", action="store_true", default=False,
                        help="Test mode: 2000x2000 center crop")
    parser.add_argument("--test-size", type=int, default=2000,
                        help="Size of test crop in pixels (default: 2000)")

    # GPU / Performance
    parser.add_argument("--gpu", action="store_true", default=True)
    parser.add_argument("--no-gpu", dest="gpu", action="store_false")
    parser.add_argument("--batch-size", type=int, default=256,
                        help="GPU batch size (default: 256 for 48GB VRAM)")
    parser.add_argument("--threads", type=int, default=None,
                        help="CPU threads for dynamics (default: auto)")
    parser.add_argument("--log-file", type=str, default=None,
                        help="Path to log file")
    parser.add_argument("--cache-enhanced", action="store_true",
                        help="Cache enhanced channels to .npz before OME-TIFF export")
    parser.add_argument("--export-ome", action="store_true", default=None,
                        help="Force OME-TIFF export (default: True in full, False in test)")
    parser.add_argument("--verbose", action="store_true", default=False,
                        help="Enable verbose output (Cellpose detailed logs)")

    # Segmentation
    parser.add_argument("--diameter", type=float, default=None,
                        help="Cellpose diameter (default: auto)")
    parser.add_argument("--sharpen-nuclear", type=float, default=0.0,
                        help="Nuclear sharpening strength (unsharp mask radius=1.0)")
    parser.add_argument("--denoise-nuclear", type=float, default=0.0,
                        help="Nuclear denoising (Gaussian blur sigma)")
    parser.add_argument("--flow-threshold", type=float, default=0.4)
    parser.add_argument("--cellprob-threshold", type=float, default=0.0)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--tile-norm-blocksize", type=int, default=0,
                        help="Cellpose tile normalization block size in pixels. "
                             "Locally normalizes intensity within tiles to brighten "
                             "dim/diffuse nuclei. Recommended: 100. 0=disabled.")

    # Enhancement
    parser.add_argument("--enhancement-gamma", type=float, default=0.5)
    parser.add_argument("--enhancement-perc-lower", type=float, default=0.005)
    parser.add_argument("--enhancement-perc-upper", type=float, default=0.995)
    parser.add_argument("--enhancement-exclude-outliers",
                        action="store_true", default=True)
    parser.add_argument("--enhancement-outlier-percentile",
                        type=float, default=0.98)
    parser.add_argument("--enhancement-adaptive",
                        action="store_true", default=True)
    parser.add_argument("--no-enhancement-adaptive",
                        dest="enhancement_adaptive", action="store_false")

    args = parser.parse_args()

    if not os.path.isfile(args.ome_tiff):
        print(f"ERROR: OME-TIFF not found: {args.ome_tiff}")
        sys.exit(1)
    os.makedirs(args.output_dir, exist_ok=True)

    # Setup logging
    if args.log_file:
        logger = Logger(args.log_file)
        sys.stdout = logger
        sys.stderr = logger

    run_pipeline(args)


if __name__ == "__main__":
    main()
