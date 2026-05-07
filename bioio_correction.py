#!/usr/bin/env python3
"""
MODULAR 4-PHASE OME-TIFF BLEEDTHROUGH CORRECTION
Memory-efficient, resumable pipeline using Bio-Formats CLI and ome-tiff-pyramid-tools

NEW ARCHITECTURE:
Phase 1: Extract individual channel TIFF files (resumable)
Phase 2: Create corrected channel TIFF (resumable) 
Phase 3: Assemble final OME-TIFF using pyramid_assemble.py (resumable)
Phase 4: Assign proper channel names for QuPath compatibility (resumable)

BENEFITS:
- Resumable: Skip completed phases
- Debuggable: Inspect intermediate results
- Efficient: Use dedicated tools for each task
- Memory-safe: Process channels individually
- QuPath-ready: Proper channel names for quantitative analysis

Documentation references:
- Bio-Formats patterns: https://bio-formats.readthedocs.io/en/stable/users/comlinetools/conversion.html
- ome-tiff-pyramid-tools: https://github.com/labsyspharm/ome-tiff-pyramid-tools
- NEUBIAS intensity guidelines: https://neubias.github.io/training-resources/measure_intensities/index.html
"""

import numpy as np
import os
import subprocess
import multiprocessing as mp
import psutil
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import gc
import tempfile

# Optional scipy import for efficient blurring
try:
    from scipy.ndimage import uniform_filter
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Channel mapping for QuPath compatibility - will be populated dynamically
CHANNEL_NAMES = {}

# Bleedthrough correction configuration
BLEEDTHROUGH_CONFIG = {
    'source_channel_idx': 8,      # CD19-FITC channel
    'target_channel_idx': 13,     # GFP-FITC channel (FoxP3-GFP)
    'scaling_factor': 5.0,       # Very aggressive scaling for complete bleedthrough removal
    'corrected_channel_name': 'CorrectedFOXP3-GFP'
}

def extract_channel_names_from_input(input_path):
    """Extract actual channel names from input OME-TIFF for accurate mapping."""
    
    print(f"🔍 EXTRACTING CHANNEL NAMES from input OME-TIFF...")
    
    try:
        # Get OME-XML metadata using Bio-Formats
        cmd = ['./bftools/showinf', '-nopix', '-omexml-only', input_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                              cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2')
        
        if result.returncode != 0:
            print(f"❌ Failed to read OME-XML from input: {result.stderr}")
            return None
        
        ome_xml = result.stdout
        print(f"✅ OME-XML extracted ({len(ome_xml)} characters)")
        
        # Parse XML to extract channel names
        import xml.etree.ElementTree as ET
        root = ET.fromstring(ome_xml)
        
        # Find namespace
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        if root.tag.startswith('{'):
            ns_uri = root.tag.split('}')[0][1:]
            ns = {'ome': ns_uri}
        
        # Find all Channel elements and extract names
        channels = root.findall('.//ome:Channel', ns)
        channel_names = {}
        
        print(f"📋 FOUND {len(channels)} CHANNELS:")
        for i, channel in enumerate(channels):
            channel_name = channel.get('Name', f'Channel_{i:02d}')
            channel_names[i] = channel_name
            
            # Highlight the bleedthrough channels
            if i == BLEEDTHROUGH_CONFIG['source_channel_idx']:
                print(f"  📍 Channel {i:2d}: {channel_name} (BLEEDTHROUGH SOURCE)")
            elif i == BLEEDTHROUGH_CONFIG['target_channel_idx']:
                print(f"  📍 Channel {i:2d}: {channel_name} (BLEEDTHROUGH TARGET)")
            else:
                print(f"     Channel {i:2d}: {channel_name}")
        
        # Add the corrected channel name
        corrected_idx = len(channels)  # Next available index
        channel_names[corrected_idx] = BLEEDTHROUGH_CONFIG['corrected_channel_name']
        print(f"  🧮 Channel {corrected_idx:2d}: {channel_names[corrected_idx]} (CORRECTED)")
        
        print(f"✅ CHANNEL NAME EXTRACTION COMPLETE: {len(channel_names)} total channels")
        return channel_names
        
    except Exception as e:
        print(f"❌ Error extracting channel names: {e}")
        import traceback
        traceback.print_exc()
        return None

def get_memory_usage():
    """Get current memory usage as percentage."""
    return psutil.virtual_memory().percent

def check_dependencies():
    """Check if required tools are available."""
    dependencies = {}
    
    # Check Bio-Formats CLI
    try:
        result = subprocess.run(['./bftools/bfconvert'], 
                              capture_output=True, text=True, timeout=10)
        # bfconvert shows help but returns exit code 1 when run without arguments
        # If it runs and shows output, it's available
        dependencies['bfconvert'] = (result.returncode in [0, 1] and 
                                   ('bfconvert' in result.stdout or 'To convert a file' in result.stdout))
    except:
        dependencies['bfconvert'] = False
    
    # Check pyramid_assemble.py (ome-tiff-pyramid-tools)
    try:
        # Only check for required dependencies: tifffile (required), scipy (optional)
        result = subprocess.run(['python3', '-c', 'import tifffile; import zarr'], 
                              capture_output=True, timeout=10)
        dependencies['pyramid_assemble_deps'] = result.returncode == 0
    except:
        dependencies['pyramid_assemble_deps'] = False
    
    # Check if pyramid_assemble.py exists
    pyramid_script = Path('pyramid_assemble.py')
    dependencies['pyramid_assemble_script'] = pyramid_script.exists()
    
    return dependencies

def setup_temp_directory(output_dir, test_mode=False):
    """Set up persistent temporary directory in output folder."""
    mode_suffix = "_TEST" if test_mode else ""
    temp_dir_name = f"bioformats_channels{mode_suffix}"
    temp_dir = os.path.join(output_dir, temp_dir_name)
    
    os.makedirs(temp_dir, exist_ok=True)
    return temp_dir

# ============================================================================
# PHASE 1: EXTRACT INDIVIDUAL CHANNEL TIFF FILES 
# ============================================================================

def extract_channel_bioformats(args):
    """Extract single channel using Bio-Formats CLI to individual TIFF file."""
    channel_idx, input_path, temp_dir, test_mode = args
    
    print(f"📥 Extracting channel {channel_idx} with Bio-Formats CLI...")
    
    try:
        # Generate meaningful filename using actual channel names
        channel_name = CHANNEL_NAMES.get(channel_idx, f"Channel_{channel_idx:02d}") if CHANNEL_NAMES else f"Channel_{channel_idx:02d}"
        output_filename = f"{channel_name}.tif"
        output_path = os.path.join(temp_dir, output_filename)
        
        # Skip if file already exists (resumable)
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024**3)  # GB
            print(f"✅ Channel {channel_idx} ({channel_name}) already exists: {file_size:.2f}GB")
            return channel_idx, output_path, True
        
        # Bio-Formats bfconvert command for single channel extraction
        cmd = [
            './bftools/bfconvert',
            '-channel', str(channel_idx),    # Extract only this channel
            '-tilex', '1024',                # Large tiles for efficiency
            '-tiley', '1024',                
            '-compression', 'LZW',           # LZW compression
            '-bigtiff',                      # BigTIFF support
            '-overwrite',                    # Overwrite if exists
            str(input_path),
            output_path
        ]
        
        # Run bfconvert with timeout
        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minutes per channel
            cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2'
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode == 0 and os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024**3)  # GB
            print(f"✅ Channel {channel_idx} ({channel_name}) extracted in {elapsed:.1f}s: {file_size:.2f}GB")
            return channel_idx, output_path, True
        else:
            print(f"❌ Channel {channel_idx} extraction failed: {result.stderr}")
            return channel_idx, None, False
            
    except subprocess.TimeoutExpired:
        print(f"❌ Channel {channel_idx} timed out after 600s")
        return channel_idx, None, False
    except Exception as e:
        print(f"❌ Channel {channel_idx} error: {e}")
        return channel_idx, None, False

def phase1_extract_channels(input_path, temp_dir, num_workers, test_mode=False):
    """Phase 1: Extract individual channel TIFF files (resumable)."""
    
    mode_info = "TEST MODE" if test_mode else "FULL MODE"
    print(f"\n🔥 PHASE 1: {mode_info} - Extract Individual Channel TIFFs")
    print("=" * 60)
    print(f"📁 Temp directory: {temp_dir}")
    print(f"💡 RESUMABLE: Will skip channels that already exist")
    
    # Determine which channels to extract
    if test_mode:
        # TEST MODE: First 5 channels + bleedthrough channels (source and target)
        source_idx = BLEEDTHROUGH_CONFIG['source_channel_idx']
        target_idx = BLEEDTHROUGH_CONFIG['target_channel_idx']
        channel_list = [0, 1, 2, 3, 4, source_idx, target_idx]
        # Remove duplicates and sort
        channel_list = sorted(list(set(channel_list)))
        print(f"🧪 TEST SCOPE: Extracting {len(channel_list)} channels: {channel_list}")
        print(f"  📍 Includes bleedthrough source (ch{source_idx}) and target (ch{target_idx})")
    else:
        # FULL MODE: All 30 original channels
        channel_list = list(range(30))
        print(f"🔄 FULL SCOPE: Extracting {len(channel_list)} original channels")
    
    # Prepare extraction tasks
    tasks = [(i, input_path, temp_dir, test_mode) for i in channel_list]
    
    print(f"👥 Using {num_workers} parallel workers")
    
    start_time = time.time()
    extracted_files = {}
    
    # Process channels in parallel
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        future_to_channel = {executor.submit(extract_channel_bioformats, task): task[0] 
                           for task in tasks}
        
        completed_count = 0
        for future in as_completed(future_to_channel):
            channel_idx = future_to_channel[future]
            try:
                idx, temp_file, success = future.result()
                if success and temp_file:
                    extracted_files[idx] = temp_file
                    completed_count += 1
                    print(f"✅ Channel {idx} completed ({completed_count}/{len(tasks)}) - Memory: {get_memory_usage():.1f}%")
                else:
                    print(f"❌ Failed to extract channel {idx}")
                    return None
            except Exception as exc:
                print(f"❌ Channel {channel_idx} generated exception: {exc}")
                return None
    
    elapsed_time = time.time() - start_time
    
    print(f"\n✅ PHASE 1 COMPLETE: Extracted {len(extracted_files)} channels in {elapsed_time:.1f} seconds")
    print(f"⚡ Average: {elapsed_time/len(extracted_files):.2f} seconds per channel")
    print(f"📊 Memory usage: {get_memory_usage():.1f}%")
    
    return extracted_files

def re_extract_corrupted_channel(channel_idx, input_path, temp_dir):
    """Re-extract a single corrupted channel using Bio-Formats CLI."""
    
    print(f"🔄 ATTEMPTING RE-EXTRACTION: Channel {channel_idx}")
    
    try:
        # Generate meaningful filename using actual channel names
        channel_name = CHANNEL_NAMES.get(channel_idx, f"Channel_{channel_idx:02d}") if CHANNEL_NAMES else f"Channel_{channel_idx:02d}"
        output_filename = f"{channel_name}.tif"
        output_path = os.path.join(temp_dir, output_filename)
        
        # Remove corrupted file first
        if os.path.exists(output_path):
            print(f"  🗑️  Removing corrupted file: {output_filename}")
            os.remove(output_path)
        
        # Bio-Formats bfconvert command for single channel extraction
        cmd = [
            './bftools/bfconvert',
            '-channel', str(channel_idx),    # Extract only this channel
            '-tilex', '1024',                # Large tiles for efficiency
            '-tiley', '1024',                
            '-compression', 'LZW',           # LZW compression
            '-bigtiff',                      # BigTIFF support
            '-overwrite',                    # Overwrite if exists
            str(input_path),
            output_path
        ]
        
        # Run bfconvert with timeout
        print(f"  ⚙️  Running: bfconvert -channel {channel_idx} ...")
        start_time = time.time()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600,  # 10 minutes per channel
            cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2'
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode == 0 and os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024**3)  # GB
            print(f"  ✅ RE-EXTRACTION SUCCESS: Channel {channel_idx} ({channel_name}) in {elapsed:.1f}s: {file_size:.2f}GB")
            return True
        else:
            print(f"  ❌ RE-EXTRACTION FAILED: {result.stderr}")
            return False
            
    except subprocess.TimeoutExpired:
        print(f"  ❌ RE-EXTRACTION TIMEOUT: Channel {channel_idx} after 600s")
        return False
    except Exception as e:
        print(f"  ❌ RE-EXTRACTION ERROR: {e}")
        return False

# ============================================================================
# PHASE 2: CREATE CORRECTED CHANNEL TIFF
# ============================================================================

def phase2_create_corrected_channel(temp_dir, test_mode=False):
    """Phase 2: Create corrected channel TIFF (resumable)."""
    
    print(f"\n🧮 PHASE 2: Create Corrected Channel TIFF")
    print("=" * 60)
    
    # Define corrected channel filename
    corrected_filename = BLEEDTHROUGH_CONFIG['corrected_channel_name'] + ".tif"
    corrected_path = os.path.join(temp_dir, corrected_filename)
    
    # Skip if corrected file already exists (resumable)
    if os.path.exists(corrected_path):
        file_size = os.path.getsize(corrected_path) / (1024**3)  # GB
        print(f"✅ Corrected channel already exists: {file_size:.2f}GB")
        print(f"📁 Path: {corrected_path}")
        return corrected_path
    
    # Check if required channels exist
    source_idx = BLEEDTHROUGH_CONFIG['source_channel_idx']
    target_idx = BLEEDTHROUGH_CONFIG['target_channel_idx']
    
    source_name = CHANNEL_NAMES.get(source_idx, f"Channel_{source_idx:02d}") if CHANNEL_NAMES else f"Channel_{source_idx:02d}"
    target_name = CHANNEL_NAMES.get(target_idx, f"Channel_{target_idx:02d}") if CHANNEL_NAMES else f"Channel_{target_idx:02d}"
    
    cd19_path = os.path.join(temp_dir, f"{source_name}.tif")
    gfp_path = os.path.join(temp_dir, f"{target_name}.tif")
    
    if not os.path.exists(cd19_path):
        print(f"❌ CD19-FITC channel not found: {cd19_path}")
        return None
    
    if not os.path.exists(gfp_path):
        print(f"❌ FoxP3-GFP channel not found: {gfp_path}")
        return None
    
    print(f"📥 Using channels:")
    print(f"   CD19-FITC: {cd19_path}")
    print(f"   FoxP3-GFP: {gfp_path}")
    
    # Use existing corrected channel creation logic (memory-safe approach)
    print(f"🧮 Creating corrected channel using memory-safe Bio-Formats approach...")
    
    try:
        # Reuse the existing create_corrected_channel_pure_bioformats function
        # but modify to save directly to temp_dir with proper name
        result = create_corrected_channel_simple_python(cd19_path, gfp_path, corrected_path)
        
        if result and os.path.exists(corrected_path):
            file_size = os.path.getsize(corrected_path) / (1024**3)  # GB
            print(f"✅ PHASE 2 COMPLETE: Corrected channel created: {file_size:.2f}GB")
            print(f"📁 Saved: {corrected_path}")
            return corrected_path
        else:
            print(f"❌ PHASE 2 FAILED: Could not create corrected channel")
            return None
    
    except Exception as e:
        print(f"❌ PHASE 2 ERROR: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_corrected_channel_simple_python(cd19_path, gfp_path, output_path):
    """Create corrected channel using enhanced Python raster subtraction with scaling and blur."""
    
    try:
        import numpy as np
        import tifffile
        
        scaling_factor = BLEEDTHROUGH_CONFIG['scaling_factor']
        blur_kernel_size = 15  # Large blur kernel for extensive peripheral bleedthrough capture
        
        print(f"  🧮 Enhanced bleedthrough correction:")
        print(f"    📈 Scaling factor: {scaling_factor}x")
        print(f"    🌀 Large blur kernel: {blur_kernel_size}x{blur_kernel_size} (for extensive peripheral capture)")
        print(f"  💡 Formula: Corrected = max(GFP - {scaling_factor} × blur(CD19), 0)")
        
        def simple_blur(image, kernel_size=5):
            """Apply simple averaging blur to capture peripheral bleedthrough using efficient numpy operations."""
            # For large images, use a more efficient approach
            try:
                # Try using scipy if available
                if HAS_SCIPY:
                    return uniform_filter(image, size=kernel_size)
                else:
                    # Fallback to simpler approach for smaller kernel
                    # Use separable filter for efficiency: blur_x then blur_y
                    kernel_1d = np.ones(kernel_size) / kernel_size
                    
                    # Apply 1D blur in x direction
                    blurred_x = np.zeros_like(image)
                    pad_size = kernel_size // 2
                    padded = np.pad(image, ((0, 0), (pad_size, pad_size)), mode='reflect')
                    
                    for i in range(image.shape[0]):
                        for j in range(image.shape[1]):
                            blurred_x[i, j] = np.sum(padded[i, j:j+kernel_size] * kernel_1d)
                    
                    # Apply 1D blur in y direction  
                    blurred = np.zeros_like(image)
                    padded_y = np.pad(blurred_x, ((pad_size, pad_size), (0, 0)), mode='reflect')
                    
                    for i in range(image.shape[0]):
                        for j in range(image.shape[1]):
                            blurred[i, j] = np.sum(padded_y[i:i+kernel_size, j] * kernel_1d)
                    
                    return blurred
            except ImportError:
                # Fallback to simpler approach if scipy is not available
                print("  ⚠️  scipy.ndimage not available, falling back to slower blurring.")
                # Use separable filter for efficiency: blur_x then blur_y
                kernel_1d = np.ones(kernel_size) / kernel_size
                
                # Apply 1D blur in x direction
                blurred_x = np.zeros_like(image)
                pad_size = kernel_size // 2
                padded = np.pad(image, ((0, 0), (pad_size, pad_size)), mode='reflect')
                
                for i in range(image.shape[0]):
                    for j in range(image.shape[1]):
                        blurred_x[i, j] = np.sum(padded[i, j:j+kernel_size] * kernel_1d)
                
                # Apply 1D blur in y direction  
                blurred = np.zeros_like(image)
                padded_y = np.pad(blurred_x, ((pad_size, pad_size), (0, 0)), mode='reflect')
                
                for i in range(image.shape[0]):
                    for j in range(image.shape[1]):
                        blurred[i, j] = np.sum(padded_y[i:i+kernel_size, j] * kernel_1d)
                
                return blurred
        
        def read_bioformats_tiff(filepath):
            """Read Bio-Formats generated TIFF with multiple strategies."""
            try:
                # Strategy 1: Default tifffile
                return tifffile.imread(filepath)
            except Exception as e1:
                print(f"    Strategy 1 failed ({e1}), trying alternative...")
                try:
                    # Strategy 2: Read without series detection
                    with tifffile.TiffFile(filepath) as tif:
                        # Read first page/series directly
                        return tif.pages[0].asarray()
                except Exception as e2:
                    print(f"    Strategy 2 failed ({e2}), trying raw page access...")
                    try:
                        # Strategy 3: Read all pages and combine if needed
                        with tifffile.TiffFile(filepath) as tif:
                            if len(tif.pages) == 1:
                                return tif.pages[0].asarray()
                            else:
                                # Multiple pages - Bio-Formats might have created multi-page
                                print(f"    Found {len(tif.pages)} pages, reading first page...")
                                return tif.pages[0].asarray()
                    except Exception as e3:
                        print(f"    All strategies failed: {e1}, {e2}, {e3}")
                        raise e3
        
        print(f"  📖 Loading CD19-FITC channel...")
        cd19_image = read_bioformats_tiff(cd19_path)
        print(f"    CD19 shape: {cd19_image.shape}, dtype: {cd19_image.dtype}")
        print(f"    CD19 stats: min={cd19_image.min():.2f}, max={cd19_image.max():.2f}, mean={cd19_image.mean():.2f}")
        
        print(f"  📖 Loading GFP-FITC channel...")
        gfp_image = read_bioformats_tiff(gfp_path)
        print(f"    GFP shape: {gfp_image.shape}, dtype: {gfp_image.dtype}")
        print(f"    GFP stats: min={gfp_image.min():.2f}, max={gfp_image.max():.2f}, mean={gfp_image.mean():.2f}")
        
        # Verify dimensions match
        if cd19_image.shape != gfp_image.shape:
            print(f"  ❌ Dimension mismatch: CD19 {cd19_image.shape} vs GFP {gfp_image.shape}")
            return False
        
        print(f"  🧮 Performing enhanced bleedthrough correction with blur and {scaling_factor}x scaling...")
        
        # Ensure we're working with the same data type
        if cd19_image.dtype != gfp_image.dtype:
            print(f"  🔄 Converting data types to match...")
            # Convert to the higher precision type
            if cd19_image.dtype == np.uint8 or gfp_image.dtype == np.uint8:
                target_dtype = np.uint16
            else:
                target_dtype = np.float32
            cd19_image = cd19_image.astype(target_dtype)
            gfp_image = gfp_image.astype(target_dtype)
        
        # Apply simple averaging blur to CD19 channel for peripheral bleedthrough capture
        print(f"  🌀 Applying large averaging blur ({blur_kernel_size}x{blur_kernel_size}) to CD19 channel...")
        cd19_blurred = simple_blur(cd19_image.astype(np.float32), blur_kernel_size)
        print(f"    Blurred CD19 stats: min={cd19_blurred.min():.2f}, max={cd19_blurred.max():.2f}, mean={cd19_blurred.mean():.2f}")
        
        # Apply scaling factor to blurred CD19 channel
        print(f"  📈 Scaling blurred CD19 channel by factor {scaling_factor}...")
        scaled_blurred_cd19 = cd19_blurred * scaling_factor
        print(f"    Scaled+blurred CD19 stats: min={scaled_blurred_cd19.min():.2f}, max={scaled_blurred_cd19.max():.2f}, mean={scaled_blurred_cd19.mean():.2f}")
        
        # Convert back to original data type for subtraction
        if np.issubdtype(gfp_image.dtype, np.unsignedinteger):
            # Convert to signed integer for subtraction, then clip
            scaled_blurred_cd19_signed = scaled_blurred_cd19.astype(np.int32)
            gfp_signed = gfp_image.astype(np.int32)
            corrected_signed = gfp_signed - scaled_blurred_cd19_signed
            # Clip negative values to 0
            corrected_signed = np.maximum(corrected_signed, 0)
            # Convert back to original unsigned type
            corrected_image = corrected_signed.astype(gfp_image.dtype)
        else:
            # For float types, simple subtraction with clipping
            corrected_image = np.maximum(gfp_image.astype(np.float32) - scaled_blurred_cd19, 0)
            # Convert back to original dtype if needed
            corrected_image = corrected_image.astype(gfp_image.dtype)
        
        print(f"  📊 Corrected image stats:")
        print(f"    Shape: {corrected_image.shape}")
        print(f"    Dtype: {corrected_image.dtype}")
        print(f"    Min: {np.min(corrected_image):.2f}, Max: {np.max(corrected_image):.2f}")
        print(f"    Mean: {corrected_image.mean():.2f}")
        
        # Calculate improvement metrics
        original_mean = gfp_image.mean()
        corrected_mean = corrected_image.mean()
        reduction_percent = ((original_mean - corrected_mean) / original_mean * 100) if original_mean > 0 else 0
        print(f"  📉 Signal reduction: {reduction_percent:.1f}% (indicates bleedthrough removal)")
        print(f"  🎯 Peripheral capture: Blur helps remove misaligned bleedthrough")
        
        print(f"  💾 Saving corrected channel...")
        
        # Save with same compression and tiling as input
        tifffile.imwrite(
            output_path,
            corrected_image,
            compression='lzw',
            photometric='minisblack',
            metadata={'axes': 'YX' if len(corrected_image.shape) == 2 else 'ZYX'},
            tile=(512, 512)
        )
        
        # Verify output file
        if os.path.exists(output_path):
            file_size = os.path.getsize(output_path) / (1024**3)  # GB
            print(f"  ✅ Enhanced corrected channel with blur saved: {file_size:.2f}GB")
            return True
        else:
            print(f"  ❌ Output file not created")
            return False
        
    except Exception as e:
        print(f"  ❌ Error in enhanced bleedthrough correction with blur: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# PHASE 3: ASSEMBLE FINAL OME-TIFF USING PYRAMID_ASSEMBLE.PY
# ============================================================================

def convert_bioformats_tiff_to_standard(input_path, output_path):
    """Convert Bio-Formats TIFF to standard TIFF for pyramid_assemble.py compatibility."""
    
    try:
        import numpy as np
        import tifffile
        
        def read_bioformats_tiff(filepath):
            """Read Bio-Formats generated TIFF with multiple strategies."""
            try:
                # Strategy 1: Default tifffile
                return tifffile.imread(filepath)
            except Exception as e1:
                try:
                    # Strategy 2: Read without series detection
                    with tifffile.TiffFile(filepath) as tif:
                        return tif.pages[0].asarray()
                except Exception as e2:
                    try:
                        # Strategy 3: Read all pages and combine if needed
                        with tifffile.TiffFile(filepath) as tif:
                            if len(tif.pages) == 1:
                                return tif.pages[0].asarray()
                            else:
                                return tif.pages[0].asarray()
                    except Exception as e3:
                        print(f"    ❌ All reading strategies failed: {e1}, {e2}, {e3}")
                        raise e3
        
        # Read Bio-Formats TIFF with enhanced error handling
        image_data = read_bioformats_tiff(input_path)
        
        # Convert data type if needed for pyramid_assemble.py compatibility
        if image_data.dtype == np.float32:
            # Convert float32 to uint16 (standard for microscopy)
            # Clamp values to uint16 range and convert
            image_data = np.clip(image_data, 0, 65535).astype(np.uint16)
        elif image_data.dtype == np.float64:
            # Convert float64 to uint16  
            image_data = np.clip(image_data, 0, 65535).astype(np.uint16)
        # uint16, uint8, etc. are already compatible
        
        # Save as standard TIFF
        tifffile.imwrite(
            output_path,
            image_data,
            compression='lzw',
            photometric='minisblack',
            tile=(512, 512)
        )
        
        return True
        
    except Exception as e:
        print(f"    ❌ Failed to convert {os.path.basename(input_path)}: {e}")
        # Check if this is a corrupted tile error that might be fixable by re-extraction
        error_msg = str(e).lower()
        if any(keyword in error_msg for keyword in ['corrupted', 'tile', 'reshape', 'missing required tags']):
            return "corrupted"
        else:
            return False

def phase3_assemble_final_ometiff(temp_dir, output_path, input_path, test_mode=False):
    """Phase 3: Assemble final OME-TIFF using ome-tiff-pyramid-tools."""
    
    mode_info = "TEST MODE" if test_mode else "FULL MODE"
    print(f"\n🏗️ PHASE 3: {mode_info} - Assemble Final OME-TIFF")
    print("=" * 60)
    print(f"🛠️ Using: ome-tiff-pyramid-tools (pyramid_assemble.py)")
    print(f"📚 Research: https://github.com/labsyspharm/ome-tiff-pyramid-tools")
    
    # Check if pyramid_assemble.py is available
    pyramid_script = 'pyramid_assemble.py'
    if not os.path.exists(pyramid_script):
        print(f"❌ pyramid_assemble.py not found. Please download from:")
        print(f"   https://github.com/labsyspharm/ome-tiff-pyramid-tools")
        return False
    
    # Collect all TIFF files from temp directory
    tiff_files = []
    channel_order = []
    
    # Collect files in proper channel order
    if test_mode:
        # TEST MODE: First 5 + bleedthrough channels + corrected
        source_idx = BLEEDTHROUGH_CONFIG['source_channel_idx']
        target_idx = BLEEDTHROUGH_CONFIG['target_channel_idx']
        original_channels = [0, 1, 2, 3, 4, source_idx, target_idx]
        original_channels = sorted(list(set(original_channels)))  # Remove duplicates
    else:
        # FULL MODE: All 30 original channels + corrected
        original_channels = list(range(30))
    
    # First, collect all original channels
    for channel_idx in original_channels:
        channel_name = CHANNEL_NAMES.get(channel_idx, f"Channel_{channel_idx:02d}") if CHANNEL_NAMES else f"Channel_{channel_idx:02d}"
        tiff_file = os.path.join(temp_dir, f"{channel_name}.tif")
        
        if os.path.exists(tiff_file):
            tiff_files.append(tiff_file)
            channel_order.append(channel_name)
            file_size = os.path.getsize(tiff_file) / (1024**3)  # GB
            print(f"📎 Channel {channel_idx:2d}: {channel_name} ({file_size:.2f}GB)")
        else:
            print(f"❌ Missing channel {channel_idx}: {tiff_file}")
            return False
    
    # Then, add the corrected channel (special case)
    corrected_channel_name = BLEEDTHROUGH_CONFIG['corrected_channel_name']
    corrected_tiff_file = os.path.join(temp_dir, f"{corrected_channel_name}.tif")
    
    if os.path.exists(corrected_tiff_file):
        tiff_files.append(corrected_tiff_file)
        channel_order.append(corrected_channel_name)
        file_size = os.path.getsize(corrected_tiff_file) / (1024**3)  # GB
        corrected_idx = len(original_channels)  # Next index after original channels
        print(f"📎 Channel {corrected_idx:2d}: {corrected_channel_name} ({file_size:.2f}GB) [CORRECTED]")
    else:
        print(f"❌ Missing corrected channel: {corrected_tiff_file}")
        return False
    
    print(f"\n📊 ASSEMBLY INPUT: {len(tiff_files)} TIFF files")
    total_size = sum(os.path.getsize(f) for f in tiff_files) / (1024**3)
    print(f"📊 Total input size: {total_size:.2f}GB")
    
    # STEP 1: Convert Bio-Formats TIFFs to standard TIFFs
    print(f"\n🔄 STEP 1: Converting Bio-Formats TIFFs to standard TIFFs...")
    print(f"💡 REASON: pyramid_assemble.py can't read Bio-Formats TIFFs (incompatible keyframe)")
    print(f"🔬 STRICT MODE: Will re-extract corrupted channels and fail if unsuccessful")
    
    # Create subdirectory for standard TIFFs
    standard_dir = os.path.join(temp_dir, "standard_tiffs")
    os.makedirs(standard_dir, exist_ok=True)
    
    standard_tiff_files = []
    
    for i, (tiff_file, channel_name) in enumerate(zip(tiff_files, channel_order)):
        print(f"  🔄 Converting {i+1}/{len(tiff_files)}: {channel_name}...")
        
        # Standard TIFF filename
        standard_filename = f"standard_{channel_name}.tif"
        standard_path = os.path.join(standard_dir, standard_filename)
        
        # Skip if already converted
        if os.path.exists(standard_path):
            print(f"    ✅ Already converted: {standard_filename}")
            standard_tiff_files.append(standard_path)
            continue
        
        # Convert Bio-Formats TIFF to standard TIFF
        conversion_result = convert_bioformats_tiff_to_standard(tiff_file, standard_path)
        
        if conversion_result == "corrupted":
            print(f"    ⚠️  CORRUPTED CHANNEL DETECTED: {channel_name}")
            
            # Find channel index for re-extraction
            channel_idx = None
            if CHANNEL_NAMES:
                for idx, name in CHANNEL_NAMES.items():
                    if name == channel_name:
                        channel_idx = idx
                        break
            
            if channel_idx is not None:
                print(f"    🔄 ATTEMPTING AUTOMATIC RE-EXTRACTION...")
                
                # Attempt to re-extract the corrupted channel
                re_extract_success = re_extract_corrupted_channel(channel_idx, input_path, temp_dir)
                
                if re_extract_success:
                    print(f"    🔄 Retrying conversion after re-extraction...")
                    # Try conversion again after re-extraction
                    retry_result = convert_bioformats_tiff_to_standard(tiff_file, standard_path)
                    
                    if retry_result is True and os.path.exists(standard_path):
                        standard_size = os.path.getsize(standard_path) / (1024**3)  # GB
                        print(f"    ✅ RECOVERY SUCCESS: {standard_filename} ({standard_size:.2f}GB)")
                        standard_tiff_files.append(standard_path)
                        continue
                    else:
                        print(f"    ❌ CONVERSION STILL FAILED AFTER RE-EXTRACTION")
                        print(f"❌ CRITICAL ERROR: Cannot recover channel {channel_name}")
                        print(f"❌ ABORTING: Data integrity cannot be guaranteed with missing channels")
                        return False
                else:
                    print(f"    ❌ RE-EXTRACTION FAILED")
                    print(f"❌ CRITICAL ERROR: Cannot recover channel {channel_name}")
                    print(f"❌ ABORTING: Data integrity cannot be guaranteed with missing channels")
                    return False
            else:
                print(f"    ❌ Cannot find channel index for re-extraction")
                print(f"❌ CRITICAL ERROR: Cannot recover channel {channel_name}")
                print(f"❌ ABORTING: Data integrity cannot be guaranteed with missing channels")
                return False
        
        elif conversion_result is False:
            print(f"    ❌ Conversion failed with non-recoverable error")
            print(f"❌ CRITICAL ERROR: Cannot convert channel {channel_name}")
            print(f"❌ ABORTING: Data integrity cannot be guaranteed")
            return False
        
        else:  # conversion_result is True
            if os.path.exists(standard_path):
                standard_size = os.path.getsize(standard_path) / (1024**3)  # GB
                print(f"    ✅ Converted: {standard_filename} ({standard_size:.2f}GB)")
                standard_tiff_files.append(standard_path)
            else:
                print(f"    ❌ Conversion reported success but output not created")
                print(f"❌ CRITICAL ERROR: Cannot convert channel {channel_name}")
                print(f"❌ ABORTING: Data integrity cannot be guaranteed")
                return False
    
    print(f"✅ STEP 1 COMPLETE: All {len(standard_tiff_files)} TIFFs converted successfully")
    print(f"🔬 DATA INTEGRITY: Complete dataset ready for assembly")
    
    # STEP 2: Run pyramid_assemble.py with standard TIFFs
    print(f"\n🚀 STEP 2: Running pyramid_assemble.py with standard TIFFs...")
    
    # Build pyramid_assemble.py command
    cmd = ['python3', pyramid_script]
    cmd.extend(standard_tiff_files)  # Add all standard TIFF files
    cmd.append(output_path)  # Output OME-TIFF
    cmd.extend(['--pixel-size', '0.325'])  # Pixel size in micrometers
    
    print(f"🚀 ASSEMBLY COMMAND:")
    print(f"   python3 pyramid_assemble.py [inputs...] {os.path.basename(output_path)} --pixel-size 0.325")
    print(f"💡 PERFORMANCE: Will use all available CPU cores")
    
    try:
        start_time = time.time()
        print(f"⏰ Starting assembly... (Memory: {get_memory_usage():.1f}%)")
        
        # Run pyramid_assemble.py
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2'
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode == 0 and os.path.exists(output_path):
            output_size = os.path.getsize(output_path) / (1024**3)  # GB
            print(f"🎉 PHASE 3 COMPLETE: OME-TIFF assembled in {elapsed/60:.1f} minutes!")
            print(f"📁 Output: {output_path} ({output_size:.2f}GB)")
            print(f"⚡ Processing rate: {total_size/elapsed*60:.1f} GB/minute")
            
            if test_mode:
                print(f"🧪 TEST SUCCESS: {len(channel_order)} channels assembled")
                print(f"✅ Ready for QuPath v6.0 validation")
            else:
                print(f"🔬 FULL SUCCESS: {len(channel_order)} channels assembled")
                print(f"✅ Ready for QuPath v6.0 analysis")
            
            return True
        else:
            print(f"❌ PHASE 3 FAILED: pyramid_assemble.py failed")
            print(f"📄 Error: {result.stderr}")
            print(f"📄 Output: {result.stdout}")
            return False
            
    except Exception as e:
        print(f"❌ PHASE 3 ERROR: {e}")
        return False

# ============================================================================
# PHASE 4: ASSIGN PROPER CHANNEL NAMES TO OME-TIFF
# ============================================================================

def phase4_assign_channel_names(output_path, test_mode=False):
    """Phase 4: Assign proper channel names to the assembled OME-TIFF."""
    
    mode_info = "TEST MODE" if test_mode else "FULL MODE"
    print(f"\n🏷️ PHASE 4: {mode_info} - Assign Proper Channel Names")
    print("=" * 60)
    print(f"📚 RESEARCH: Proper channel labeling crucial for quantitative microscopy")
    print(f"🎯 GOAL: QuPath should show meaningful channel names instead of 'Channel 1', etc.")
    
    if not os.path.exists(output_path):
        print(f"❌ Output file not found: {output_path}")
        return False
    
    # Determine expected channel names based on mode
    if test_mode:
        # TEST MODE: First 5 + bleedthrough channels + corrected
        source_idx = BLEEDTHROUGH_CONFIG['source_channel_idx']
        target_idx = BLEEDTHROUGH_CONFIG['target_channel_idx']
        original_channel_indices = [0, 1, 2, 3, 4, source_idx, target_idx]
        original_channel_indices = sorted(list(set(original_channel_indices)))  # Remove duplicates
    else:
        # FULL MODE: All 30 original channels + corrected
        original_channel_indices = list(range(30))
    
    # Build expected names list: original channels + corrected channel
    expected_names = []
    
    # Add original channel names
    for idx in original_channel_indices:
        channel_name = CHANNEL_NAMES.get(idx, f"Channel_{idx:02d}") if CHANNEL_NAMES else f"Channel_{idx:02d}"
        expected_names.append(channel_name)
    
    # Add corrected channel name (special case)
    corrected_channel_name = BLEEDTHROUGH_CONFIG['corrected_channel_name']
    expected_names.append(corrected_channel_name)
    
    print(f"📝 Expected channel names:")
    for i, name in enumerate(expected_names):
        print(f"  Channel {i}: {name}")
    
    try:
        import xml.etree.ElementTree as ET
        
        # Skip backup creation to save disk space due to large file sizes
        print(f"💡 Skipping backup creation to conserve disk space")
        
        print(f"🔍 Reading OME-XML metadata from assembled file...")
        
        # Read the OME-XML using Bio-Formats
        cmd = ['./bftools/showinf', '-nopix', '-omexml-only', output_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60,
                              cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2')
        
        if result.returncode != 0:
            print(f"❌ Failed to read OME-XML: {result.stderr}")
            return False
        
        ome_xml = result.stdout
        print(f"✅ OME-XML extracted ({len(ome_xml)} characters)")
        
        # Parse XML and modify channel names
        print(f"🔧 Modifying channel names in OME-XML...")
        
        # Parse the XML
        root = ET.fromstring(ome_xml)
        
        # Find namespace
        ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
        if root.tag.startswith('{'):
            # Extract namespace from root tag
            ns_uri = root.tag.split('}')[0][1:]
            ns = {'ome': ns_uri}
        
        # Find all Channel elements
        channels = root.findall('.//ome:Channel', ns)
        print(f"  Found {len(channels)} channels in OME-XML")
        
        if len(channels) != len(expected_names):
            print(f"  ⚠️  Channel count mismatch: found {len(channels)}, expected {len(expected_names)}")
        
        # Update channel names
        channels_updated = 0
        for i, channel in enumerate(channels):
            if i < len(expected_names):
                old_name = channel.get('Name', f'Channel {i}')
                new_name = expected_names[i]
                channel.set('Name', new_name)
                print(f"  📝 Channel {i}: '{old_name}' → '{new_name}'")
                channels_updated += 1
            else:
                print(f"  ⚠️  No name mapping for channel {i}")
        
        print(f"✅ Updated {channels_updated} channel names")
        
        # Convert back to string
        updated_xml = ET.tostring(root, encoding='unicode')
        
        # Write updated OME-XML to temporary file
        temp_xml_file = output_path.replace('.ome.tif', '_updated.ome.xml')
        with open(temp_xml_file, 'w') as f:
            f.write(updated_xml)
        
        print(f"💾 Applying updated OME-XML to image file...")
        
        # Use Bio-Formats to update the OME-XML in the TIFF file
        cmd = ['./bftools/tiffcomment', '-set', temp_xml_file, output_path]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                              cwd='/mnt/b/Projects/Image_Analysis_Projects/SJ_04.2')
        
        if result.returncode == 0:
            print(f"✅ PHASE 4 COMPLETE: Channel names successfully assigned!")
            print(f"🏷️  QuPath should now show meaningful channel names")
            
            # Cleanup temporary files
            if os.path.exists(temp_xml_file):
                os.remove(temp_xml_file)
            
            return True
        else:
            print(f"❌ Failed to update OME-XML in TIFF: {result.stderr}")
            return False
    
    except Exception as e:
        print(f"❌ PHASE 4 ERROR: {e}")
        import traceback
        traceback.print_exc()
        return False

# ============================================================================
# MAIN ORCHESTRATOR
# ============================================================================

def modular_bioformats_correction(input_path, output_path, test_mode=False):
    """Main orchestrator for 4-phase modular OME-TIFF correction."""
    
    mode_title = "TEST MODE" if test_mode else "FULL MODE" 
    print(f"🏗️ MODULAR 4-PHASE OME-TIFF CORRECTION - {mode_title}")
    print("=" * 80)
    
    # System info
    memory_info = psutil.virtual_memory()
    print(f"💾 System Memory: {memory_info.total/(1024**3):.1f} GB total")
    print(f"💾 Available Memory: {memory_info.available/(1024**3):.1f} GB")
    print(f"💾 Current Usage: {memory_info.percent:.1f}%")
    
    # CRITICAL: Extract actual channel names from input OME-TIFF
    print(f"🔍 STEP 0: Extract channel names from input OME-TIFF...")
    global CHANNEL_NAMES
    CHANNEL_NAMES = extract_channel_names_from_input(input_path)
    
    if CHANNEL_NAMES is None:
        print(f"❌ FAILED: Could not extract channel names from input")
        return False
    
    # Display bleedthrough correction configuration
    source_idx = BLEEDTHROUGH_CONFIG['source_channel_idx']
    target_idx = BLEEDTHROUGH_CONFIG['target_channel_idx']
    scaling = BLEEDTHROUGH_CONFIG['scaling_factor']
    
    print(f"\n🧮 BLEEDTHROUGH CORRECTION CONFIGURATION:")
    print(f"  📍 Source Channel {source_idx}: {CHANNEL_NAMES.get(source_idx, 'Unknown')}")
    print(f"  📍 Target Channel {target_idx}: {CHANNEL_NAMES.get(target_idx, 'Unknown')}")  
    print(f"  📈 Scaling Factor: {scaling}x (to enhance bleedthrough removal)")
    print(f"  🧮 Formula: Corrected = max(Target - {scaling} × Source, 0)")
    
    # Check dependencies
    print(f"\n🔍 Checking dependencies...")
    deps = check_dependencies()
    for tool, available in deps.items():
        status = "✅" if available else "❌"
        print(f"  {status} {tool}")
    
    if not all(deps.values()):
        print(f"❌ Missing dependencies. Please install required tools.")
        return False
    
    # Calculate workers
    available_gb = memory_info.available / (1024**3)
    max_workers_memory = max(1, int(available_gb * 0.8 / 2.0))  # 2GB per worker
    max_workers_cpu = min(mp.cpu_count() - 2, 16)
    optimal_workers = min(max_workers_memory, max_workers_cpu)
    
    print(f"\n📖 Input: {input_path}")
    print(f"📝 Output: {output_path}")
    print(f"👥 Workers: {optimal_workers}")
    
    if test_mode:
        print(f"🧪 TEST STRATEGY: Limited channels for fast validation")
        expected_channels = len([i for i in [0, 1, 2, 3, 4, source_idx, target_idx]]) + 1  # +1 for corrected
    else:
        print(f"🔧 FULL STRATEGY: All original channels + corrected channel")
        expected_channels = len([i for i in range(30)]) + 1  # +1 for corrected
    
    print("=" * 80)
    
    overall_start = time.time()
    
    # Set up persistent temp directory
    output_dir = os.path.dirname(output_path)
    temp_dir = setup_temp_directory(output_dir, test_mode)
    print(f"📁 Persistent temp directory: {temp_dir}")
    
    try:
        # PHASE 1: Extract individual channel TIFF files
        print(f"🎯 RESUMABLE: Each phase will skip if outputs already exist")
        extracted_files = phase1_extract_channels(input_path, temp_dir, optimal_workers, test_mode)
        
        if extracted_files is None:
            print("❌ PHASE 1 FAILED")
            return False
        
        # PHASE 2: Create corrected channel TIFF
        corrected_path = phase2_create_corrected_channel(temp_dir, test_mode)
        
        if corrected_path is None:
            print("❌ PHASE 2 FAILED")
            return False
        
        # PHASE 3: Assemble final OME-TIFF
        assembly_success = phase3_assemble_final_ometiff(temp_dir, output_path, input_path, test_mode)
        
        if not assembly_success:
            print("❌ PHASE 3 FAILED")
            return False
        
        # PHASE 4: Assign proper channel names
        name_assignment_success = phase4_assign_channel_names(output_path, test_mode)
        
        if not name_assignment_success:
            print("❌ PHASE 4 FAILED")
            return False
        
        # SUCCESS!
        total_time = time.time() - overall_start
        final_memory = get_memory_usage()
        
        print(f"\n🎉 ALL 4 PHASES COMPLETE!")
        print(f"⏱️  Total processing time: {total_time/60:.1f} minutes")
        print(f"💾 Peak memory usage: {final_memory:.1f}%")
        print(f"📁 Final output: {output_path}")
        
        if test_mode:
            print(f"🧪 TEST MODE SUCCESS: Ready for QuPath v6.0 with proper channel names!")
            print(f"💡 NEXT: Set test_mode=False for full processing")
        else:
            print(f"🔬 FULL MODE SUCCESS: Ready for QuPath v6.0 analysis with proper channel names!")
            print(f"✅ All 31 channels processed with bleedthrough correction and proper naming")
        
        return True
        
    except Exception as e:
        mode_error = "TEST" if test_mode else "FULL"
        print(f"❌ {mode_error} MODE CORRECTION FAILED: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    # Configuration
    input_file = 'data/input/EGFR_TREG_2.ome.tif'
    
    # FULL MODE: Process all original channels + corrected channel
    test_mode = False  # Set to True for limited testing
    
    if test_mode:
        output_file = 'data/output/EGFR_TREG_2_modular_TEST.ome.tif'
        print("🧪 MODULAR 4-PHASE CORRECTION - TEST MODE")
        print("=" * 60)
        print("🎯 PHASE 1: Extract limited channels (first 5 + bleedthrough channels)")
        print("🎯 PHASE 2: Enhanced bleedthrough correction with scaling")
        print("🎯 PHASE 3: Assemble using pyramid_assemble.py")
        print("🎯 PHASE 4: Assign proper channel names for QuPath")
        print("⚡ RESUMABLE: Skip completed phases")
        print("⏰ EXPECTED TIME: 15-45 minutes")
    else:
        output_file = 'data/output/EGFR_TREG_2_modular.ome.tif'
        print("🏗️ MODULAR 4-PHASE CORRECTION - FULL MODE")
        print("=" * 60)
        print("🎯 PHASE 1: Extract all 30 original channels")
        print("🎯 PHASE 2: Enhanced bleedthrough correction with scaling")
        print("🎯 PHASE 3: Assemble using pyramid_assemble.py")
        print("🎯 PHASE 4: Assign proper channel names for QuPath")
        print("⚡ RESUMABLE: Skip completed phases")
        print("⏰ EXPECTED TIME: 1-3 hours")
    
    print("🔑 MODULAR BENEFITS:")
    print("  ✅ Resumable: Skip completed phases")
    print("  ✅ Debuggable: Inspect intermediate TIFFs")
    print("  ✅ Memory-efficient: Process channels individually")
    print("  ✅ Tool-optimized: Use best tool for each phase")
    print("=" * 60)
    
    # Run the modular correction
    success = modular_bioformats_correction(input_file, output_file, test_mode)
    
    if success:
        mode_result = "TEST" if test_mode else "FULL"
        print(f"\n🎯 SUCCESS: Modular {mode_result} MODE OME-TIFF created with proper channel names!")
        if test_mode:
            print("🧪 NEXT STEP: Validate channel names in QuPath v6.0, then set test_mode=False")
    else:
        mode_result = "TEST" if test_mode else "FULL"
        print(f"\n❌ FAILED: Modular {mode_result} MODE correction unsuccessful") 