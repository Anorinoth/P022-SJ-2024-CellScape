#!/usr/bin/env python3
"""
Merge single-channel TIFFs into a multi-channel OME-TIFF for QuPath v6
Author: Dr. Michael Zaiken

This script reads single-channel TIFF files exported from the enhance_and_quantify
workflow and combines them into a single multi-channel OME-TIFF that can be 
loaded in QuPath v6.

Usage:
    python merge_channels_to_ome.py <input_dir> <output_file.ome.tiff>
    
Example:
    python merge_channels_to_ome.py \\
        /path/to/RES/image_name/enhanced_tiffs \\
        /path/to/RES/image_name/image_name_multichannel.ome.tiff

Requirements:
    pip install tifffile numpy

Date: 2025-10-29
"""

import tifffile
import numpy as np
from pathlib import Path
import sys
import re
from datetime import datetime


def extract_channel_info(filename):
    """
    Extract channel number and name from filename.
    
    Expected format: imagename_chXX_channelname.tif
    
    Returns:
        tuple: (channel_number, channel_name)
    """
    # Match pattern: _chXX_
    match = re.search(r'_ch(\d{2})_(.+)\.tif$', filename)
    if match:
        ch_num = int(match.group(1))
        ch_name = match.group(2).replace('_', ' ')
        return (ch_num, ch_name)
    else:
        # Fallback: try to extract just channel number
        match = re.search(r'ch(\d{2})', filename)
        if match:
            ch_num = int(match.group(1))
            return (ch_num, f"Channel {ch_num}")
        else:
            return (None, None)


def merge_channels_to_ome_tiff(input_dir, output_file, compression='lzw', verbose=True):
    """
    Merge single-channel TIFFs into multi-channel OME-TIFF.
    
    Args:
        input_dir (str): Directory containing single-channel TIFF files
        output_file (str): Path for output OME-TIFF file
        compression (str): Compression type ('lzw', 'zlib', or None)
        verbose (bool): Print progress messages
    
    Returns:
        None
    """
    
    input_path = Path(input_dir)
    
    if not input_path.exists():
        print(f"❌ Error: Input directory does not exist: {input_dir}")
        sys.exit(1)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Merging Single-Channel TIFFs to Multi-Channel OME-TIFF")
        print(f"{'='*70}")
        print(f"Input directory: {input_dir}")
        print(f"Output file: {output_file}")
        print(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"{'='*70}\n")
    
    # Get all TIFF files
    tiff_files = list(input_path.glob("*_ch*.tif"))
    
    if not tiff_files:
        print(f"❌ Error: No TIFF files matching pattern '*_ch*.tif' found in {input_dir}")
        sys.exit(1)
    
    # Extract channel info and sort by channel number
    channel_info = []
    for tiff_file in tiff_files:
        ch_num, ch_name = extract_channel_info(tiff_file.name)
        if ch_num is not None:
            channel_info.append((ch_num, ch_name, tiff_file))
        else:
            print(f"⚠️  Warning: Could not extract channel info from: {tiff_file.name}")
    
    # Sort by channel number
    channel_info.sort(key=lambda x: x[0])
    
    if verbose:
        print(f"Found {len(channel_info)} channel files:\n")
    
    # Read all channels
    channels = []
    channel_names = []
    
    for i, (ch_num, ch_name, tiff_path) in enumerate(channel_info):
        if verbose:
            print(f"  [{i+1:2d}/{len(channel_info):2d}] Ch{ch_num:02d}: {ch_name}")
            print(f"         Loading: {tiff_path.name}")
        
        # Read image
        img = tifffile.imread(tiff_path)
        
        # Check if 2D
        if img.ndim != 2:
            print(f"❌ Error: Expected 2D image, got {img.ndim}D: {tiff_path.name}")
            sys.exit(1)
        
        channels.append(img)
        channel_names.append(ch_name)
        
        if verbose and i == 0:
            print(f"         Dimensions: {img.shape[0]} x {img.shape[1]}")
            print(f"         Data type: {img.dtype}")
            print()
    
    # Verify all channels have same dimensions
    shapes = [ch.shape for ch in channels]
    if len(set(shapes)) > 1:
        print(f"❌ Error: Not all channels have the same dimensions:")
        for i, shape in enumerate(shapes):
            print(f"    Channel {i+1}: {shape}")
        sys.exit(1)
    
    if verbose:
        print(f"\n{'='*70}")
        print(f"Stacking {len(channels)} channels into multi-channel array...")
        print(f"{'='*70}\n")
    
    # Stack into multi-channel array (C, Y, X)
    # QuPath expects channel-first ordering
    multi_channel = np.stack(channels, axis=0)
    
    if verbose:
        print(f"Multi-channel array shape: {multi_channel.shape} (C, Y, X)")
        print(f"  Channels: {multi_channel.shape[0]}")
        print(f"  Height: {multi_channel.shape[1]}")
        print(f"  Width: {multi_channel.shape[2]}")
        print(f"  Data type: {multi_channel.dtype}")
        print(f"  Memory size: {multi_channel.nbytes / 1e9:.2f} GB\n")
    
    # Create OME-XML metadata for QuPath v6 compatibility
    # QuPath v6 uses Bio-Formats which expects proper OME-TIFF metadata
    metadata = {
        'axes': 'CYX',  # Channel, Y, X ordering
        'Channel': {'Name': channel_names},
        'SignificantBits': 16 if multi_channel.dtype == np.uint16 else 8
    }
    
    if verbose:
        print(f"{'='*70}")
        print(f"Writing OME-TIFF file...")
        print(f"{'='*70}\n")
        print(f"  Format: OME-TIFF (BigTIFF)")
        print(f"  Compression: {compression.upper() if compression else 'None'}")
        print(f"  Bits per sample: {metadata['SignificantBits']}")
        print(f"  Output: {output_file}\n")
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Write BigTIFF with OME-XML metadata
    # BigTIFF is required for files > 4 GB
    # LZW compression is lossless and widely supported
    tifffile.imwrite(
        output_file,
        multi_channel,
        bigtiff=True,
        compression=compression if compression else None,
        metadata=metadata,
        photometric='minisblack'  # Grayscale channels
    )
    
    if verbose:
        file_size = Path(output_file).stat().st_size / 1e9
        print(f"{'='*70}")
        print(f"✓ SUCCESS: Multi-channel OME-TIFF created")
        print(f"{'='*70}\n")
        print(f"  File: {output_file}")
        print(f"  Size: {file_size:.2f} GB")
        print(f"  Channels: {len(channel_names)}")
        print(f"  Dimensions: {multi_channel.shape[0]} x {multi_channel.shape[1]} x {multi_channel.shape[2]} (C x Y x X)")
        print(f"\n  Ready to load in QuPath v6!")
        print(f"\n{'='*70}\n")


def main():
    """Main entry point for command-line usage."""
    
    if len(sys.argv) < 3:
        print("Usage: python merge_channels_to_ome.py <input_dir> <output_file.ome.tiff>")
        print("\nExample:")
        print("  python merge_channels_to_ome.py \\")
        print("    /path/to/RES/image_name/enhanced_tiffs \\")
        print("    /path/to/RES/image_name/image_name_multichannel.ome.tiff")
        sys.exit(1)
    
    input_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    # Optional: compression argument
    compression = 'lzw'  # Default
    if len(sys.argv) > 3:
        compression = sys.argv[3].lower()
        if compression not in ['lzw', 'zlib', 'none']:
            print(f"Warning: Unknown compression '{compression}', using 'lzw'")
            compression = 'lzw'
        if compression == 'none':
            compression = None
    
    try:
        merge_channels_to_ome_tiff(input_dir, output_file, compression=compression, verbose=True)
    except Exception as e:
        print(f"\n❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()

