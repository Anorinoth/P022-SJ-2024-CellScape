# Germinal Center and T Regulatory Cell Analysis Macros for QuPath

This directory contains Groovy macros for comprehensive analysis of germinal centers and regulatory T cells in corrected fluorescence images using QuPath v0.5.1.

## Contents

- `Practical_GC_Treg_Analysis.groovy` - Main working macro for production use
- `Germinal_Center_Treg_Analysis.groovy` - Advanced macro with InstaSeg integration (experimental)
- `README.md` - This documentation file

## Overview

These macros perform the following analysis pipeline:

1. **Germinal Center Core Identification** - Uses PNA+ staining to identify GC cores
2. **B Cell Zone Detection** - Uses CD19-FITC to identify B cell zones overlapping/extending GC cores
3. **T Cell Zone Detection** - Uses CD4-BV421 to identify T cell zones bordering B cell zones
4. **Cell Segmentation** - Segments individual cells using nuclear staining (C30 channel)
5. **FoxP3+ Cell Classification** - Identifies Treg cells using FoxP3-GFP (C31 channel)
6. **Proximity Analysis** - Finds cells within specified radius of FoxP3+ cells
7. **Quantification & Export** - Exports comprehensive measurements as TSV files

## Prerequisites

### QuPath Setup
- QuPath v0.5.1 or later
- Your corrected images loaded in QuPath
- Project created (recommended for saving classifiers and results)

### Image Requirements
- Multi-channel fluorescence images with the following channels:
  - **Peanut Agglutinin-FITC** (Channel 1) - For germinal center core identification
  - **CD19-FITC** (Channel 8) - For B cell zone identification  
  - **CD4-Brilliant Violet 421** (Channel 3) - For T cell zone identification
  - **CorrectedFOXP3-GFP** (Channel 30) - For Treg identification
  - **DNA-BUV395** (Channel 20) - Nuclear stain for cell segmentation

### Optional: InstaSeg Extension
For advanced cell segmentation, install the InstaSeg extension:
1. Download from the QuPath extension manager or GitHub
2. Follow installation instructions in QuPath documentation
3. Use the advanced macro version for enhanced segmentation

## Usage Instructions

### Quick Start (Recommended)

1. **Open Your Image**
   ```
   File → Open → Select your corrected image from data/output/
   ```

2. **Load the Practical Macro**
   ```
   Automate → Show script editor
   File → Open script → Select Practical_GC_Treg_Analysis.groovy
   ```

3. **Configure Parameters**
   Edit the configuration section at the top of the script:
   ```groovy
   // Threshold values (adjust based on your image intensities)
   def PNA_THRESHOLD = 0.3         // Threshold for PNA+ (Germinal Center cores)
   def CD19_THRESHOLD = 0.25       // Threshold for CD19-FITC (B cell zones)
   def CD4_THRESHOLD = 0.2         // Threshold for CD4-BV421 (T cell zones)
   def FOXP3_THRESHOLD = 0.35      // Threshold for FoxP3-GFP (Tregs)
   
   // Channel names - CONFIGURED FOR YOUR CORRECTED IMAGES
   def PNA_CHANNEL = "Peanut Agglutinin-FITC"    // Channel 1
   def CD19_CHANNEL = "CD19-FITC"                 // Channel 8
   def CD4_CHANNEL = "CD4-Brilliant Violet 421"   // Channel 3
   def FOXP3_CHANNEL = "CorrectedFOXP3-GFP"       // Channel 30
   def NUCLEAR_CHANNEL = "DNA-BUV395"             // Channel 20
   ```

4. **Run the Analysis**
   ```
   Run → Run (or Ctrl+R)
   ```

5. **Monitor Progress**
   - Watch the console output for progress updates
   - Analysis typically takes 5-30 minutes depending on image size
   - Check for any error messages or warnings

6. **Review Results**
   - Generated TSV files will be in your project directory
   - Annotations will be visible in QuPath for visual verification

### Advanced Configuration

#### Adjusting Detection Sensitivity

**If you're missing GC cores or zones:**
- Lower the thresholds (e.g., PNA_THRESHOLD = 0.2)
- Reduce minimum area constraints (MIN_GC_AREA = 500.0)

**If you're getting false positives:**
- Raise the thresholds (e.g., PNA_THRESHOLD = 0.4)
- Increase minimum area constraints (MIN_GC_AREA = 1500.0)

#### Cell Detection Parameters

```groovy
// Geometric constraints
def MIN_GC_AREA = 1000.0        // Minimum area for GC cores (µm²)
def MIN_ZONE_AREA = 500.0       // Minimum area for B/T cell zones (µm²)
def PROXIMITY_RADIUS = 50.0     // Radius for proximity analysis (µm)

// Cell detection parameters
def CELL_EXPANSION = 2.0        // Cell expansion from nucleus (µm)
def MIN_CELL_AREA = 20.0        // Minimum cell area (µm²)
def MAX_CELL_AREA = 500.0       // Maximum cell area (µm²)
```

#### Channel Name Mapping

Make sure channel names match your image exactly. You can check channel names in QuPath:
```
View → Brightness/Contrast → Channel names are listed there
```

## Output Files

The macro generates three TSV files:

### 1. FoxP3_Treg_Quantification.tsv
Contains data for each FoxP3+ cell:
- `Unique_ID` - Unique identifier for each Treg cell
- `GC_Location` - Location (GC_Core, B_Cell_Zone, T_Cell_Zone, Outside_GC)
- `Zone_ID` - ID number of the specific zone
- `FoxP3_Mean` - Average FoxP3 intensity in nucleus
- `Nuclear_Mean` - Average nuclear stain intensity
- `Cell_Area` - Cell area in µm²
- `Nucleus_Area` - Nucleus area in µm²
- `X_Centroid`, `Y_Centroid` - Cell coordinates

### 2. Proximity_Cell_Quantification.tsv
Contains data for cells within proximity radius of FoxP3+ cells:
- `Cell_ID` - Sequential cell identifier
- `Proximal_FoxP3_ID` - ID of the nearest FoxP3+ cell
- `Distance_to_FoxP3` - Distance to nearest FoxP3+ cell in µm
- `FoxP3_Mean` - FoxP3 intensity (may be low/negative)
- `Nuclear_Mean` - Nuclear stain intensity
- `Cell_Area`, `Nucleus_Area` - Morphometric measurements
- `X_Centroid`, `Y_Centroid` - Cell coordinates

### 3. Analysis_Summary.tsv
Contains overall statistics:
- Total cell counts
- FoxP3+ cell counts and percentages
- Proximity cell counts
- Number of detected GC cores, B cell zones, T cell zones

## Troubleshooting

### Common Issues

**"No image is currently open!"**
- Solution: Open an image in QuPath before running the script

**"Channel not found" errors**
- Solution: Check channel names in Brightness/Contrast dialog and update the configuration section

**No cells detected**
- Check nuclear channel name (C30)
- Adjust cell detection thresholds
- Verify image has proper pixel calibration

**No GC cores/zones detected**
- Lower the threshold values
- Check that the specified channels exist in your image
- Verify channel names match exactly (case-sensitive)

**Analysis very slow**
- Large images take time - be patient
- Consider reducing image resolution for initial testing
- Close other applications to free up memory

### Performance Optimization

**For large images (>2GB):**
- Increase QuPath memory allocation in preferences
- Consider analyzing smaller regions of interest first
- Use lower resolution for initial parameter testing

**For parameter optimization:**
- Start with a small test region
- Use manual thresholding in QuPath GUI first to determine good threshold values
- Test on a representative image before batch processing

### Getting Help

1. **Check Console Output** - Most errors will be reported in the script console
2. **QuPath Documentation** - https://qupath.readthedocs.io/
3. **QuPath Forum** - https://forum.image.sc/tag/qupath
4. **InstaSeg Documentation** - https://github.com/instanseg/instanseg

## Integration with InstaSeg

For enhanced cell segmentation, use the experimental macro:

1. Install InstaSeg extension in QuPath
2. Use `Germinal_Center_Treg_Analysis.groovy` instead
3. The macro will automatically attempt to use InstaSeg if available
4. Falls back to standard detection if InstaSeg is not found

## Batch Processing

To process multiple images:

1. Create a QuPath project with all your images
2. Open the first image and run the script
3. Use QuPath's "Run for project" feature:
   ```
   Run → Run for project
   ```
4. Select all images you want to process
5. Results will be generated for each image

## Data Analysis

The exported TSV files can be analyzed in:
- **R/RStudio** - For statistical analysis and visualization
- **Python/Pandas** - For data processing and machine learning
- **Excel/LibreOffice** - For basic analysis and visualization
- **GraphPad Prism** - For biostatistics and publication-quality plots
- **ImageJ/FIJI** - For spatial analysis and visualization

## Citation

If you use these macros in your research, please cite:
- QuPath: Bankhead, P. et al. QuPath: Open source software for digital pathology image analysis. Scientific Reports (2017)
- InstaSeg (if used): Goldsborough, T. et al. InstanSeg: an embedding-based instance segmentation algorithm optimized for accurate, efficient and portable cell segmentation. arXiv (2024)

## Version History

- **v1.0** - Initial release with basic functionality
- **v1.1** - Added InstaSeg integration and improved error handling
- **v1.2** - Enhanced proximity analysis and export functionality

## Contact

For questions about these specific macros, please refer to the project documentation or create an issue in the project repository. 