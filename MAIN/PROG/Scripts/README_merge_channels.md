# Merging Single-Channel TIFFs to Multi-Channel OME-TIFF

This directory contains scripts to merge single-channel TIFF files (exported from `enhance_and_quantify()`) into a multi-channel OME-TIFF that can be loaded in QuPath v6.

## Why Not Use bfconvert?

The Bio-Formats `bfconvert` command-line tool **cannot merge** multiple single-channel files into a multi-channel image. It only accepts one input file and is designed to convert or split existing multi-channel images, not combine separate files.

Reference: [Bio-Formats bfconvert documentation](https://docs.openmicroscopy.org/bio-formats/5.8.2/users/comlinetools/conversion.html)

## Installation

### 1. Install Python Dependencies

```bash
# From the Scripts directory
pip install -r requirements.txt

# Or install directly
pip install tifffile numpy
```

### 2. Verify Installation

```bash
python3 -c "import tifffile, numpy; print('✓ Dependencies installed')"
```

## Usage

### Option 1: Command Line (Direct Python)

```bash
# Basic usage
python3 merge_channels_to_ome.py \
  /path/to/enhanced_tiffs \
  /path/to/output_multichannel.ome.tiff

# With specific compression
python3 merge_channels_to_ome.py \
  /path/to/enhanced_tiffs \
  /path/to/output_multichannel.ome.tiff \
  lzw
```

**Arguments:**
- `input_dir`: Directory containing single-channel TIFF files (from `enhance_and_quantify()`)
- `output_file`: Path for output OME-TIFF file (should end in `.ome.tiff`)
- `compression`: Optional, compression type: `lzw` (default), `zlib`, or `none`

### Option 2: From R (Using Wrapper Function)

```r
# Load the wrapper function
source("MAIN/PROG/Scripts/merge_tiffs_wrapper.R")

# Basic usage
merge_tiffs_to_ome(
  tiff_dir = "MAIN/RES/sample_01/enhanced_tiffs",
  output_file = "MAIN/RES/sample_01/sample_01_multichannel.ome.tiff"
)

# Specify Python path (if needed)
merge_tiffs_to_ome(
  tiff_dir = "MAIN/RES/sample_01/enhanced_tiffs",
  output_file = "MAIN/RES/sample_01/sample_01_multichannel.ome.tiff",
  python_path = "/usr/bin/python3"  # or "python" on Windows
)

# Disable dependency checks (faster, if already verified)
merge_tiffs_to_ome(
  tiff_dir = "MAIN/RES/sample_01/enhanced_tiffs",
  output_file = "MAIN/RES/sample_01/sample_01_multichannel.ome.tiff",
  check_dependencies = FALSE
)
```

### Option 3: Integrated in SpicyFlow.Rmd

Add this after the `enhance_and_quantify()` call:

```r
# Merge single-channel TIFFs to multi-channel OME-TIFF for QuPath
if (!is.null(tiff_dir) && dir.exists(tiff_dir)) {
  message("\n=== Creating Multi-Channel OME-TIFF for QuPath ===")
  
  # Load wrapper function
  source(file.path(PROG_DIR, "Scripts", "merge_tiffs_wrapper.R"))
  
  # Define output file
  ome_tiff_file <- file.path(RES_DIR, image_name, 
                              paste0(image_name, "_multichannel.ome.tiff"))
  
  # Merge channels
  merge_result <- merge_tiffs_to_ome(
    tiff_dir = tiff_dir,
    output_file = ome_tiff_file,
    python_path = "python3",  # or "python" on Windows
    check_dependencies = TRUE
  )
  
  if (merge_result$success) {
    message("✓ Multi-channel OME-TIFF ready for QuPath v6")
    message("  File: ", ome_tiff_file)
  } else {
    warning("Failed to create multi-channel OME-TIFF")
  }
}
```

## Expected Input Files

The script expects single-channel TIFF files with this naming pattern:
```
imagename_ch01_channelname.tif
imagename_ch02_channelname.tif
...
imagename_ch31_channelname.tif
```

These are automatically created by `enhance_and_quantify()` when `export_tiffs = TRUE`.

## Output Format

**Multi-channel OME-TIFF** with:
- ✅ 31 color channels (or however many input TIFFs)
- ✅ Single timepoint (no time series)
- ✅ No Z-stacking (single plane)
- ✅ BigTIFF format (for files > 4 GB)
- ✅ LZW compression (lossless, widely supported)
- ✅ Proper OME-XML metadata for QuPath v6
- ✅ Channel names preserved from input filenames

## Loading in QuPath v6

1. Open QuPath v6
2. **File → Open → Choose image...**
3. Select the `*_multichannel.ome.tiff` file
4. QuPath will automatically recognize all 31 channels
5. Channels can be toggled on/off in the brightness/contrast panel

## Troubleshooting

### "Python not found"
```bash
# Check Python installation
python3 --version

# If not installed, install Python 3.8+
# Ubuntu/Debian: sudo apt install python3 python3-pip
# macOS: brew install python3
# Windows: Download from python.org
```

### "tifffile not found"
```bash
# Install Python dependencies
pip install tifffile numpy

# Or from requirements file
pip install -r requirements.txt
```

### "No TIFF files found"
- Ensure `export_tiffs = TRUE` in `enhance_and_quantify()` config
- Check that TIFF files follow naming pattern: `*_ch*.tif`
- Verify the input directory path is correct

### "Channel dimensions don't match"
- All channels must have identical dimensions (height x width)
- This should never happen if all TIFFs are from `enhance_and_quantify()`

## Technical Details

- **Input format**: Single-channel 16-bit TIFF files
- **Output format**: Multi-channel OME-TIFF (BigTIFF)
- **Channel order**: Determined by channel number in filename (ch01, ch02, ...)
- **Compression**: LZW (lossless, ~50% reduction)
- **Array shape**: (C, Y, X) - channels first
- **Photometric**: Minisblack (grayscale channels)

## Files

- `merge_channels_to_ome.py` - Python script (main implementation)
- `merge_tiffs_wrapper.R` - R wrapper function
- `requirements.txt` - Python dependencies
- `README_merge_channels.md` - This file

## References

- [Bio-Formats bfconvert documentation](https://docs.openmicroscopy.org/bio-formats/5.8.2/users/comlinetools/conversion.html)
- [tifffile Python package](https://pypi.org/project/tifffile/)
- [QuPath v6 documentation](https://qupath.readthedocs.io/)

