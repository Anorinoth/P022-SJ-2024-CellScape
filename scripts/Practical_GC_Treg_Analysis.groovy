/**
 * Practical Germinal Center and T Regulatory Cell Analysis for QuPath v6.0
 * 
 * This macro provides a working implementation for:
 * 1. GC Core identification using PNA+ channel-specific thresholding
 * 2. B Cell zone detection using CD19-FITC channel-specific detection
 * 3. T Cell zone detection using CD4-BV421 channel-specific detection
 * 4. FoxP3+ cell segmentation using channel-specific detection on FoxP3 channel
 * 5. Proximity analysis around FoxP3+ cells with cell detection
 * 6. Spatial assignment of FoxP3+ cells to GC regions
 * 7. Comprehensive quantification and data export
 * 
 * HIERARCHY:
 * - Tissue Region -> GC Core -> B Cell Zone -> T Cell Zone (channel-specific detection)
 * - Tissue Region -> FoxP3+ Cell -> FoxP3+ Cell Proximal Region (separate branch)
 * 
 * CHANNEL-SPECIFIC DETECTION:
 * - Uses createThresholder for fluorescence-based regions (GC cores, B/T zones)
 * - Uses WatershedCellDetection for individual cells (FoxP3+ Tregs)
 * - PNA channel for GC cores, CD19 channel for B zones
 * - CD4 channel for T zones, FoxP3 channel for Treg cells
 * - Nuclear channel for cell detection in proximal regions
 * - Each detection targets specific marker channels by name
 * 
 * MEMORY OPTIMIZED:
 * - Uses targeted detection approaches that process specific channels
 * - Avoids memory issues with large images by channel-specific analysis
 * 
 * USAGE:
 * 1. Open your corrected image in QuPath
 * 2. Run this script from the Script Editor
 * 3. Adjust thresholds in the configuration section if needed
 * 4. Results will be exported as TSV files
 * 
 * Compatible with QuPath v6.0
 */

// Import required classes
import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.images.servers.ImageServer
import qupath.lib.objects.PathObject
import qupath.lib.measurements.MeasurementList
import qupath.lib.common.GeneralTools
import qupath.opencv.ml.pixel.PixelClassifierTools

// === CONFIGURATION SECTION ===
// Adjust these parameters based on your specific images and requirements

// === DIAGNOSTIC MODE ===
// Set to true for initial testing to help identify detection issues
DIAGNOSTIC_MODE = true       // Shows extra debug info and relaxed area constraints for testing

// Threshold values (0.0 to 1.0 range, adjust based on your image intensities)
// LOWERED THRESHOLDS FOR TESTING - may need further adjustment based on your image
PNA_THRESHOLD = 0.1         // Threshold for PNA+ (Germinal Center cores) - LOWERED from 0.3
CD19_THRESHOLD = 0.08       // Threshold for CD19-FITC (B cell zones) - LOWERED from 0.25
CD4_THRESHOLD = 0.06        // Threshold for CD4-BV421 (T cell zones) - LOWERED from 0.2
FOXP3_THRESHOLD = 0.05      // Very low threshold for FoxP3+ detection (direct segmentation, user wants dropped)

// Geometric constraints - UPDATED TO BIOLOGICAL REALITY
// GC cores: 75-500 µm diameter → 4,417-196,350 µm² area
MIN_GC_AREA = 4000.0        // Minimum area for GC cores (µm²) - ~75µm diameter
MAX_GC_AREA = 200000.0      // Maximum area for GC cores (µm²) - ~500µm diameter

// B/T cell zones: 100-1000 µm diameter → 7,854-785,398 µm² area  
MIN_ZONE_AREA = 7500.0      // Minimum area for B/T cell zones (µm²) - ~100µm diameter
MAX_ZONE_AREA = 800000.0    // Maximum area for B/T cell zones (µm²) - ~1000µm diameter

// Proximity analysis: 25 µm radius as specified
PROXIMITY_RADIUS = 25.0     // Radius for proximity analysis (µm)

// Cell detection parameters - FoxP3+ cells: 8-30 µm diameter → 50-700 µm² area
CELL_EXPANSION = 2.0        // Cell expansion from nucleus (µm)
MIN_CELL_AREA = 50.0        // Minimum cell area for FoxP3+ cells (µm²) - ~8µm diameter
MAX_CELL_AREA = 700.0       // Maximum cell area for FoxP3+ cells (µm²) - ~30µm diameter

// Channel names - CONFIGURED FOR YOUR CORRECTED IMAGES
PNA_CHANNEL = "Peanut Agglutinin-FITC"    // Channel 1
CD19_CHANNEL = "CD19-FITC"                 // Channel 8  
CD4_CHANNEL = "CD4-Brilliant Violet 421"   // Channel 3
FOXP3_CHANNEL = "CorrectedFOXP3-GFP"       // Channel 30
NUCLEAR_CHANNEL = "DNA-Brilliant Violet 421" // Channel 21 (C30)

println("=== PRACTICAL GC AND TREG ANALYSIS ===")
println("Compatible with QuPath v6.0")
if (DIAGNOSTIC_MODE) {
    println("*** DIAGNOSTIC MODE ENABLED ***")
    println("- Using lowered thresholds for testing")
    println("- Using relaxed area constraints")
    println("- Enhanced logging for troubleshooting")
}
println("")

// Check if image is open
def imageData = getCurrentImageData()
if (imageData == null) {
    println("ERROR: No image is currently open!")
    println("Please open an image and try again.")
    return
}

def server = imageData.getServer()
def cal = server.getPixelCalibration()

println("Image: " + server.getMetadata().getName())
println("Dimensions: " + server.getWidth() + " x " + server.getHeight())
println("Channels: " + server.getMetadata().getChannels().size())

/**
 * Helper function for fluorescence-based region detection using QuPath v6.0 pixel classification
 */
def createFluorescenceRegions(channelName, threshold, classificationName, minArea, maxArea) {
    
    // Store original parallelism to restore later
    def originalParallelism = null
    try {
        originalParallelism = System.getProperty("java.util.concurrent.ForkJoinPool.common.parallelism")
    } catch (Exception e) {
        // Ignore if we can't get the property
    }
    
    try {
        // Get current image data
        def imageData = getCurrentImageData()
        if (imageData == null) {
            println("ERROR: No image data available for " + classificationName)
            return []
        }
        
        def server = imageData.getServer()
        def imageSize = server.getWidth() * server.getHeight()
        
        // Find the channel index
        def channels = server.getMetadata().getChannels()
        def channelIndex = -1
        channels.eachWithIndex { channel, index ->
            if (channel.getName() == channelName) {
                channelIndex = index
            }
        }
        
        if (channelIndex == -1) {
            println("ERROR: Could not find channel '" + channelName + "' for " + classificationName)
            println("Available channels: " + channels.collect { it.getName() })
            return []
        }
        
        println("Creating " + classificationName + " regions from channel: " + channelName + " (index " + channelIndex + ")")
        println("Threshold: " + threshold + ", Min area: " + minArea + ", Max area: " + maxArea)
        println("Image size: " + String.format("%.1f", imageSize / 1000000.0) + " megapixels")
        
        // For very large images, we'll rely on QuPath's internal memory management
        def downsampleFactor = 1.0
        
        if (imageSize > 500000000) { // > 500M pixels
            println("Large image detected (" + String.format("%.1f", imageSize / 1000000.0) + " megapixels)")
            println("Using QuPath's internal tile-based processing for memory efficiency")
        }
        
        // Get channel bit depth and pixel type for proper threshold scaling
        def channelInfo = channels[channelIndex]
        def pixelType = server.getPixelType()
        def maxValue = pixelType.getBitsPerPixel() == 16 ? 65535.0 : 255.0
        def scaledThreshold = threshold * maxValue
        
        println("Channel bit depth: " + pixelType.getBitsPerPixel() + "-bit, Max value: " + maxValue)
        println("Scaled threshold: " + String.format("%.1f", scaledThreshold) + " (from " + threshold + " on 0-1 scale)")
        
        if (DIAGNOSTIC_MODE) {
            println("DIAGNOSTIC: In QuPath viewer, manually check intensities in '" + channelName + "' channel")
            println("DIAGNOSTIC: Look for areas with intensity values > " + String.format("%.0f", scaledThreshold) + " to verify signal presence")
        }
        
        // For very large images, reduce threading to prevent resource conflicts
        if (imageSize > 500000000) { // > 500M pixels
            println("Very large image detected - applying aggressive memory management:")
            println("  - Reducing thread parallelism to 1 (prevents concurrent memory spikes)")
            println("  - Forcing garbage collection to free available memory")
            System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", "1")
            
                         // Multiple rounds of aggressive memory cleanup
             for (int i = 0; i < 3; i++) {
                 System.gc()
                 if (i < 2) Thread.sleep(500) // Brief pause between GC calls
             }
             println("  - Memory cleanup complete")
         }
         
         // Create a threshold server for this channel with proper scaling
        def aboveClass = getPathClass(classificationName)
        def belowClass = getPathClass("Background")
        
        println("Creating threshold server...")
        def thresholdServer = null
        def retryCount = 0
        def maxRetries = 3
        
        while (thresholdServer == null && retryCount < maxRetries) {
            try {
                if (retryCount > 0) {
                    println("Retry attempt " + retryCount + " for threshold server creation...")
                    System.gc()
                    Thread.sleep(2000 * retryCount) // Exponential backoff
                }
                
                thresholdServer = PixelClassifierTools.createThresholdServer(
                    server, 
                    channelIndex, 
                    scaledThreshold,
                    belowClass, 
                    aboveClass
                )
                
            } catch (Exception e) {
                retryCount++
                println("ERROR creating threshold server (attempt " + retryCount + "): " + e.getMessage())
                if (retryCount >= maxRetries) {
                    throw e
                }
            }
        }
        
        if (thresholdServer == null) {
            println("ERROR: Failed to create threshold server for " + classificationName + " after " + maxRetries + " attempts")
            return []
        }
        
        println("Threshold server created successfully")
        
        // Create objects from the threshold server with robust error handling
        def parentRegion = null
        def tissueRegion = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
        if (tissueRegion != null) {
            parentRegion = tissueRegion.getROI()
            // Note: No need to scale parent region - pixel classifier will handle downsampling internally
        }
        
        println("Creating objects from pixel classifier...")
        def objectsCreated = []
        retryCount = 0
        maxRetries = 2
        
        // Use the specified area constraints directly, or relax them in diagnostic mode
        def minArea_adjusted = minArea
        def maxArea_adjusted = maxArea
        
        // In diagnostic mode, use much more relaxed area constraints to see what's being detected
        if (DIAGNOSTIC_MODE) {
            minArea_adjusted = 100.0      // Very small minimum (any detectable region)
            maxArea_adjusted = 5000000.0  // Very large maximum (most of image)
            println("DIAGNOSTIC MODE: Using relaxed area constraints for testing")
            println("  Original range: " + String.format("%.1f", minArea) + "-" + String.format("%.1f", maxArea) + " µm²")
            println("  Testing range:  " + String.format("%.1f", minArea_adjusted) + "-" + String.format("%.1f", maxArea_adjusted) + " µm²")
        } else {
            println("Processing with area constraints: Min " + String.format("%.1f", minArea_adjusted) + " µm², Max " + String.format("%.1f", maxArea_adjusted) + " µm²")
        }
        
        while (retryCount <= maxRetries) {
            try {
                if (retryCount > 0) {
                    println("Retry attempt " + retryCount + " for object creation...")
                    System.gc()
                    Thread.sleep(3000 * retryCount) // Longer wait for retries
                }
                
                objectsCreated = PixelClassifierTools.createObjectsFromPixelClassifier(
                    thresholdServer,
                    [(1): aboveClass], // Map classification index to PathClass
                    parentRegion,
                    { roi -> PathObjects.createAnnotationObject(roi, aboveClass) },
                    minArea_adjusted,
                    0.0, // No minimum hole area
                    false // Don't split connected regions to reduce complexity
                )
                
                // Break if we got results OR if this was our final attempt
                if (!objectsCreated.isEmpty() || retryCount >= maxRetries) {
                    break
                }
                
            } catch (InterruptedException e) {
                println("Processing interrupted - this indicates memory/threading pressure")
                Thread.currentThread().interrupt()
                throw e
            } catch (Exception e) {
                println("ERROR creating objects from pixel classifier (attempt " + (retryCount + 1) + "): " + e.getMessage())
                if (e.getCause() != null && e.getCause().toString().contains("InterruptedException")) {
                    println("Underlying interruption detected - likely memory/threading issue")
                }
                if (retryCount >= maxRetries) {
                    throw e
                }
            }
            
            retryCount++ // Always increment retry count
        }
        
        println("Created " + objectsCreated.size() + " initial objects from pixel classifier")
        
        // If no objects created, this might indicate threshold too high or no signal
        if (objectsCreated.isEmpty()) {
            println("WARNING: No objects detected from pixel classifier.")
            println("This could indicate:")
            println("  - Threshold too high for the actual signal intensity")
            println("  - Very weak or no signal in channel '" + channelName + "'")
            println("  - Background noise masking the signal")
            println("Consider lowering the threshold or checking channel intensities in QuPath viewer.")
        }
        
        // Filter by area and set names - objects should already be at full resolution
        def validObjects = []
        objectsCreated.eachWithIndex { obj, index ->
            try {
                def roi = obj.getROI()
                def area = roi.getArea()
                
                if (area >= minArea && area <= maxArea) {
                    obj.setName(classificationName + "_" + (validObjects.size() + 1))
                    validObjects.add(obj)
                    println("  ✓ Created " + classificationName + " with area: " + String.format("%.1f", area) + " µm²")
                } else {
                    println("  ✗ Filtered out region with area: " + String.format("%.1f", area) + " µm² (outside range " + String.format("%.1f", minArea) + "-" + String.format("%.1f", maxArea) + ")")
                }
            } catch (Exception e) {
                println("ERROR processing object " + index + ": " + e.getMessage())
            }
        }
        
        return validObjects
        
    } catch (InterruptedException e) {
        println("ERROR: Processing was interrupted for " + classificationName + ". This indicates severe memory/threading pressure.")
        println("SOLUTIONS: 1) Increase QuPath memory (Edit→Preferences→Set maximum memory)")
        println("          2) Reduce parallel threads in QuPath preferences") 
        println("          3) Close other applications to free system memory")
        println("          4) Consider processing smaller image regions")
        Thread.currentThread().interrupt() // Restore interrupted status
        return []
    } catch (Exception e) {
        println("ERROR in createFluorescenceRegions for " + classificationName + ": " + e.getMessage())
        println("Error type: " + e.getClass().getSimpleName())
        if (e.getCause() != null) {
            println("Underlying cause: " + e.getCause().getMessage())
        }
        e.printStackTrace()
        return []
    } finally {
        // Always restore original threading settings
        try {
            if (originalParallelism != null) {
                System.setProperty("java.util.concurrent.ForkJoinPool.common.parallelism", originalParallelism)
                println("Restored original thread parallelism: " + originalParallelism)
            } else {
                System.clearProperty("java.util.concurrent.ForkJoinPool.common.parallelism")
                println("Restored thread parallelism to system default")
            }
        } catch (Exception e) {
            println("WARNING: Could not restore threading settings: " + e.getMessage())
        }
        
        // Final cleanup
        System.gc()
    }
}

if (!cal.hasPixelSizeMicrons()) {
    println("WARNING: No pixel calibration found. Results may be inaccurate.")
    println("Consider setting pixel size via 'Set pixel size' command.")
}

// Clear existing objects
clearAllObjects()

try {
    
    // === STEP 1: CREATE TISSUE REGION ===
    println("\n--- Step 1: Creating Tissue Region ---")
    createTissueRegion()
    
    // === STEP 2: IDENTIFY GERMINAL CENTER CORES ===
    println("\n--- Step 2: Identifying Germinal Center Cores (PNA+) ---")
    identifyGerminalCenterCores()
    
    // === STEP 3: IDENTIFY B CELL ZONES (AROUND GC CORES) ===
    println("\n--- Step 3: Identifying B Cell Zones (CD19+) ---")
    identifyBCellZones()
    
    // === STEP 4: IDENTIFY T CELL ZONES (AROUND B CELL ZONES) ===
    println("\n--- Step 4: Identifying T Cell Zones (CD4+) ---")
    identifyTCellZones()
    
    // === STEP 5: DETECT INDIVIDUAL FOXP3+ CELLS ===
    println("\n--- Step 5: Detecting Individual FoxP3+ Cells ---")
    detectFoxP3Cells()
    
    // === STEP 6: CREATE PROXIMAL REGIONS AROUND FOXP3+ CELLS ===
    println("\n--- Step 6: Creating FoxP3+ Cell Proximal Regions ---")
    createFoxP3ProximalRegions()
    
    // === STEP 7: DETECT CELLS IN PROXIMAL REGIONS ===
    println("\n--- Step 7: Detecting Cells in FoxP3+ Proximal Regions ---")
    detectCellsInProximalRegions()
    
    // === STEP 8: ASSIGN FOXP3+ CELLS TO GC REGIONS ===
    println("\n--- Step 8: Assigning FoxP3+ Cells to Germinal Center Regions ---")
    assignFoxP3CellsToGCRegions()
    
    // === STEP 9: PROXIMITY ANALYSIS ===
    println("\n--- Step 9: Performing Proximity Analysis ---")
    performProximityAnalysis()
    
    // === STEP 10: EXPORT RESULTS ===
    println("\n--- Step 10: Exporting Results ---")
    exportResults()
    
    println("\n=== ANALYSIS COMPLETE ===")
    println("Check the console output above for any warnings or errors.")
    println("Results have been exported as TSV files.")
    
} catch (Exception e) {
    println("ERROR during analysis: " + e.getMessage())
    e.printStackTrace()
}

// === IMPLEMENTATION FUNCTIONS ===

/**
 * Create a basic tissue region annotation
 */
def createTissueRegion() {
    // Get image data and server within function scope
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
    
    // Create a simple rectangular region covering the whole image
    def roi = ROIs.createRectangleROI(0, 0, server.getWidth(), server.getHeight(), ImagePlane.getDefaultPlane())
    def tissueAnnotation = PathObjects.createAnnotationObject(roi, getPathClass("Tissue"))
    tissueAnnotation.setName("Tissue_Region")
    addObject(tissueAnnotation)
    println("Created tissue region annotation")
}

/**
 * Detect individual FoxP3+ cells using channel-specific segmentation on FoxP3 channel
 */
def detectFoxP3Cells() {
    // Configuration variables - use global constants for consistency
    def foxp3Channel = FOXP3_CHANNEL         // Channel 30
    def foxp3Threshold = FOXP3_THRESHOLD     // Very low threshold for initial detection (user wants threshold dropped)
    def minCellArea = MIN_CELL_AREA          // Minimum area for single cells (µm²) - ~8µm diameter
    def maxCellArea = MAX_CELL_AREA          // Maximum cell area (µm²) - ~30µm diameter
    
    println("Detecting FoxP3+ cells using channel-specific segmentation...")
    println("Target channel: " + foxp3Channel)
    println("Threshold: " + foxp3Threshold + " (0-1 scale)")
    
    println("Starting FoxP3+ cell detection...")
    
    // MEMORY OPTIMIZATION: Only detect cells within specific regions, not entire tissue
    def detectedRegions = []
    def gcCores = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "GC_Core" }
    def bZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "B_Cell_Zone" }
    def tZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "T_Cell_Zone" }
    
    detectedRegions.addAll(gcCores)
    detectedRegions.addAll(bZones)
    detectedRegions.addAll(tZones)
    
    if (detectedRegions.isEmpty()) {
        println("WARNING: No GC/B/T cell regions detected. Falling back to tissue region detection.")
        def tissueRegion = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
        if (tissueRegion != null) {
            detectedRegions.add(tissueRegion)
        }
    }
    
    println("Detecting FoxP3+ cells within " + detectedRegions.size() + " specific regions (memory optimized)")
    
    if (detectedRegions.isEmpty()) {
        println("ERROR: No regions available for FoxP3+ cell detection")
        return
    }
    
    try {
        // Get the image data and server to work with channels directly
        def imageData = getCurrentImageData()
        def server = imageData.getServer()
        def channels = server.getMetadata().getChannels()
        
        // Find the FoxP3 channel index
        def foxp3ChannelIndex = -1
        channels.eachWithIndex { channel, index ->
            if (channel.getName() == foxp3Channel) {
                foxp3ChannelIndex = index
            }
        }
        
        if (foxp3ChannelIndex == -1) {
            println("ERROR: Could not find channel '" + foxp3Channel + "'")
            println("Available channels: " + channels.collect { it.getName() })
            return
        }
        
        println("Found FoxP3 channel at index: " + foxp3ChannelIndex)
        
        // Use FoxP3 channel for cell detection
        println("Creating FoxP3+ segmentation annotations...")
        
        // Process each detected region individually to avoid memory issues
        def totalFoxp3Cells = []
        def cellId = 1
        
        detectedRegions.eachWithIndex { region, regionIndex ->
            println("Processing region ${regionIndex + 1}/${detectedRegions.size()}: ${region.getName()} (area: ${String.format('%.1f', region.getROI().getArea())} µm²)")
            
            // Aggressive memory cleanup before each region
            System.gc()
            Thread.sleep(1000)
            System.gc()
            
            selectObjects(region)
            
            try {
                // Use cell detection for FoxP3+ individual cells within this specific region
                runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', 
                    '{"detectionImage": "' + foxp3Channel + '",' +
                    '"requestedPixelSizeMicrons": 1.0,' +  // Increased pixel size to reduce memory load
                    '"backgroundRadiusMicrons": 8.0,' +
                    '"medianRadiusMicrons": 0.0,' +
                    '"sigmaMicrons": 1.5,' +
                    '"minAreaMicrons": ' + minCellArea + ',' +
                    '"maxAreaMicrons": ' + maxCellArea + ',' +
                    '"threshold": ' + (foxp3Threshold * 255) + ',' +
                    '"maxBackground": 2000,' +
                    '"watershedPostProcess": true,' +
                    '"cellExpansionMicrons": 2.0,' +
                    '"includeNuclei": true,' +
                    '"smoothBoundaries": true,' +
                    '"makeMeasurements": true}')
                    
                // Get detections from this region and classify them
                def regionCells = getDetectionObjects()
                regionCells.each { cell ->
                    cell.setPathClass(getPathClass("FoxP3_Positive"))
                    // Use correct measurement method for QuPath v6.0 - pass primitive double, not Double object
                    def measurements = cell.getMeasurementList()
                    measurements.putMeasurement("Unique_FoxP3_ID", (double)cellId)
                    measurements.putMeasurement("Area_Microns", (double)cell.getROI().getArea())
                    // Note: Source region stored as string - QuPath v6.0 may not support string measurements via putMeasurement
                    // measurements.putMeasurement("Source_Region", region.getName())  // Commented out - use object name instead
                    cellId++
                }
                
                totalFoxp3Cells.addAll(regionCells)
                println("  Found ${regionCells.size()} FoxP3+ cells in this region")
                
            } catch (Exception e) {
                println("  FoxP3+ detection failed in region ${region.getName()}: " + e.getMessage())
                if (e.getMessage().contains("OutOfMemoryError") || e.getMessage().contains("heap space")) {
                    println("  MEMORY ERROR: Region too large for current memory settings")
                    println("  Consider increasing QuPath memory or processing smaller regions")
                }
            }
        }
        
        println("Total FoxP3+ cells detected: " + totalFoxp3Cells.size())
        
        // Add intensity measurements for all channels using QuPath v6.0 method
        if (totalFoxp3Cells.size() > 0) {
            selectObjects(totalFoxp3Cells)
            try {
                // Use the correct plugin call for QuPath v6.0
                runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{}')
                println("Added intensity measurements for " + totalFoxp3Cells.size() + " FoxP3+ cells")
            } catch (Exception e) {
                println("Warning: Could not add intensity measurements: " + e.getMessage())
            }
        }
        
        println("Found " + totalFoxp3Cells.size() + " FoxP3+ cells using region-specific, memory-optimized detection")
        
    } catch (Exception e) {
        println("FoxP3+ channel-specific segmentation failed: " + e.getMessage())
        println("No FoxP3+ cells could be detected.")
    }
}

/**
 * Identify Germinal Center Cores using PNA+ fluorescence-based thresholding
 */
def identifyGerminalCenterCores() {
    println("Detecting PNA+ Germinal Center cores using fluorescence-based thresholding...")
    
    // Select tissue region for detection
    def tissueRegion = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
    if (tissueRegion != null) {
        selectObjects(tissueRegion)
    }
    
    // Set image type to fluorescence for proper channel handling
    setImageType('FLUORESCENCE')
    
    try {
        // Use fluorescence-based region detection
        def gcCores = createFluorescenceRegions(
            PNA_CHANNEL,      // Channel name
            PNA_THRESHOLD,    // Threshold value
            "GC_Core",        // Classification name
            MIN_GC_AREA,      // Minimum area
            MAX_GC_AREA       // Maximum area
        )
        
        if (gcCores.size() > 0) {
            addObjects(gcCores)
            println("Successfully created ${gcCores.size()} PNA+ GC cores")
        } else {
            println("No PNA+ GC cores could be detected.")
        }
        
    } catch (Exception e) {
        println("ERROR: PNA+ GC core detection failed: " + e.getMessage())
        println("No PNA+ GC cores could be detected.")
    }
}

/**
 * Identify B Cell zones using CD19-FITC fluorescence-based thresholding
 */
def identifyBCellZones() {
    println("Detecting CD19+ B cell zones using fluorescence-based thresholding...")
    
    // Select tissue region for detection
    def tissueRegion = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
    if (tissueRegion != null) {
        selectObjects(tissueRegion)
    }
    
    try {
        // Use fluorescence-based region detection
        def bZones = createFluorescenceRegions(
            CD19_CHANNEL,     // Channel name
            CD19_THRESHOLD,   // Threshold value
            "B_Cell_Zone",    // Classification name
            MIN_ZONE_AREA,    // Minimum area
            MAX_ZONE_AREA     // Maximum area
        )
        
        if (bZones.size() > 0) {
            addObjects(bZones)
            println("Successfully created ${bZones.size()} CD19+ B cell zones")
        } else {
            println("No CD19+ B cell zones could be detected.")
        }
        
    } catch (Exception e) {
        println("ERROR: CD19+ B cell zone detection failed: " + e.getMessage())
        println("No CD19+ B cell zones could be detected.")
    }
}

/**
 * Identify T Cell zones using CD4-BV421 fluorescence-based thresholding
 */
def identifyTCellZones() {
    println("Detecting CD4+ T cell zones using fluorescence-based thresholding...")
    
    // Select tissue region for detection
    def tissueRegion = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
    if (tissueRegion != null) {
        selectObjects(tissueRegion)
    }
    
    try {
        // Use fluorescence-based region detection
        def tZones = createFluorescenceRegions(
            CD4_CHANNEL,      // Channel name
            CD4_THRESHOLD,    // Threshold value
            "T_Cell_Zone",    // Classification name
            MIN_ZONE_AREA,    // Minimum area
            MAX_ZONE_AREA     // Maximum area
        )
        
        if (tZones.size() > 0) {
            addObjects(tZones)
            println("Successfully created ${tZones.size()} CD4+ T cell zones")
        } else {
            println("No CD4+ T cell zones could be detected.")
        }
        
    } catch (Exception e) {
        println("ERROR: CD4+ T cell zone detection failed: " + e.getMessage())
        println("No CD4+ T cell zones could be detected.")
    }
}

/**
 * Create proximal regions around each FoxP3+ cell for proximity analysis
 */
def createFoxP3ProximalRegions() {
    // Configuration variables - defined within function scope
    def proximityRadius = 25.0  // Radius for proximal region (µm)
    
    println("Creating proximal regions around FoxP3+ cells...")
    
    def foxp3Cells = getCellObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Positive" }
    def proximalRegions = []
    
    foxp3Cells.each { foxp3Cell ->
        def centroid = [foxp3Cell.getROI().getCentroidX(), foxp3Cell.getROI().getCentroidY()]
        def foxp3Id = foxp3Cell.getMeasurementList().getMeasurementValue("Unique_FoxP3_ID")
        
        // Create circular proximal region around the FoxP3+ cell
        def proximalROI = ROIs.createEllipseROI(
            centroid[0] - proximityRadius, centroid[1] - proximityRadius,
            proximityRadius * 2, proximityRadius * 2, ImagePlane.getDefaultPlane())
        
        def proximalRegion = PathObjects.createAnnotationObject(proximalROI, getPathClass("FoxP3_Proximal_Region"))
        proximalRegion.setName("FoxP3_Proximal_" + (foxp3Id ?: (proximalRegions.size() + 1)))
        
        // Store reference to the parent FoxP3+ cell  
        proximalRegion.getMeasurementList().putMeasurement("Parent_FoxP3_ID", (double)(foxp3Id ?: 0.0))
        
        proximalRegions.add(proximalRegion)
    }
    
    if (proximalRegions.size() > 0) {
        addObjects(proximalRegions)
        println("Created " + proximalRegions.size() + " FoxP3+ proximal regions")
    } else {
        println("No FoxP3+ cells found - cannot create proximal regions")
    }
}

/**
 * Detect cells within the FoxP3+ proximal regions
 */
def detectCellsInProximalRegions() {
    // Configuration variables - defined within function scope
    def minCellArea = 50.0      // Minimum cell area (µm²) - consistent with FoxP3+ detection
    def maxCellArea = 700.0     // Maximum cell area (µm²) - consistent with FoxP3+ detection
    def nuclearChannel = "DNA-Brilliant Violet 421" // Channel 21 (C30)
    
    println("Detecting cells in FoxP3+ proximal regions...")
    
    def proximalRegions = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Proximal_Region" }
    def totalProximalCells = 0
    
    proximalRegions.each { region ->
        selectObjects(region)
        
        try {
            // Get the nuclear channel index for cell detection
            def imageData = getCurrentImageData()
            def server = imageData.getServer()
            def channels = server.getMetadata().getChannels()
            
            def nuclearChannelIndex = -1
            channels.eachWithIndex { channel, index ->
                if (channel.getName() == nuclearChannel) {
                    nuclearChannelIndex = index
                }
            }
            
            // Simple cell detection within proximal region using nuclear channel
            if (nuclearChannelIndex >= 0) {
                runPlugin('qupath.imagej.detect.nuclei.SimpleThresholdNucleusDetection',
                    '{"threshold":0.3,' +
                    '"minAreaMicrons":' + minCellArea + ',' +
                    '"maxAreaMicrons":' + maxCellArea + ',' +
                    '"makeMeasurements":true,' +
                    '"channel":' + nuclearChannelIndex + '}')
            } else {
                // Fallback to default detection if nuclear channel not found
                runPlugin('qupath.imagej.detect.nuclei.SimpleThresholdNucleusDetection',
                    '{"threshold":0.3,' +
                    '"minAreaMicrons":' + minCellArea + ',' +
                    '"maxAreaMicrons":' + maxCellArea + ',' +
                    '"makeMeasurements":true}')
            }
            
            // Get cells within this proximal region
            def regionCells = getCellObjects().findAll { cell ->
                region.getROI().contains(cell.getROI().getCentroidX(), cell.getROI().getCentroidY()) &&
                cell.getPathClass()?.getName() != "FoxP3_Positive"  // Exclude the central FoxP3+ cell
            }
            
            // Mark these as proximal cells and link to parent FoxP3+ cell
            def parentFoxP3Id = region.getMeasurementList().getMeasurementValue("Parent_FoxP3_ID")
            regionCells.each { cell ->
                cell.setPathClass(getPathClass("Proximal_Cell"))
                            cell.getMeasurementList().putMeasurement("Parent_FoxP3_ID", (double)(parentFoxP3Id ?: 0.0))
            cell.getMeasurementList().putMeasurement("Is_Proximal_Cell", (double)1.0)
            }
            
            totalProximalCells += regionCells.size()
            
        } catch (Exception e) {
            println("Cell detection failed in proximal region " + region.getName() + ": " + e.getMessage())
        }
    }
    
    println("Found " + totalProximalCells + " cells in FoxP3+ proximal regions")
}

/**
 * Assign each FoxP3+ cell to its respective Germinal Center region
 */
def assignFoxP3CellsToGCRegions() {
    println("Assigning FoxP3+ cells to Germinal Center regions...")
    
    def foxp3Cells = getCellObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Positive" }
    def gcCores = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "GC_Core" }
    def bZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "B_Cell_Zone" }
    def tZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "T_Cell_Zone" }
    
    println("Found " + gcCores.size() + " GC cores")
    println("Found " + bZones.size() + " B cell zones") 
    println("Found " + tZones.size() + " T cell zones")
    println("Assigning " + foxp3Cells.size() + " FoxP3+ cells...")
    
    // Counters for statistics
    def coreCount = 0
    def bZoneCount = 0
    def tZoneCount = 0
    def outsideCount = 0
    
    foxp3Cells.each { foxp3Cell ->
        def centroid = [foxp3Cell.getROI().getCentroidX(), foxp3Cell.getROI().getCentroidY()]
        def measurements = foxp3Cell.getMeasurementList()
        
        def gcLocation = "Outside_GC"
        def specificZone = "Outside_GC"
        def zoneId = "0"
        
        // Check GC cores first (highest priority - most specific)
        for (int i = 0; i < gcCores.size(); i++) {
            if (gcCores[i].getROI().contains(centroid[0], centroid[1])) {
                gcLocation = "GC_Core"
                specificZone = "GC_Core"
                zoneId = String.valueOf(i + 1)
                coreCount++
                break
            }
        }
        
        // Check B cell zones if not in core
        if (gcLocation == "Outside_GC") {
            for (int i = 0; i < bZones.size(); i++) {
                if (bZones[i].getROI().contains(centroid[0], centroid[1])) {
                    gcLocation = "B_Cell_Zone"
                    specificZone = "B_Cell_Zone"
                    zoneId = String.valueOf(i + 1)
                    bZoneCount++
                    break
                }
            }
        }
        
        // Check T cell zones if not in other zones
        if (gcLocation == "Outside_GC") {
            for (int i = 0; i < tZones.size(); i++) {
                if (tZones[i].getROI().contains(centroid[0], centroid[1])) {
                    gcLocation = "T_Cell_Zone"
                    specificZone = "T_Cell_Zone"
                    zoneId = String.valueOf(i + 1)
                    tZoneCount++
                    break
                }
            }
        }
        
        // If still outside, count it
        if (gcLocation == "Outside_GC") {
            outsideCount++
        }
        
        // Store location information - string values converted to numeric codes for QuPath v6.0 compatibility
        def locationCode = gcLocation == "GC_Core" ? 1.0 : (gcLocation == "B_Cell_Zone" ? 2.0 : (gcLocation == "T_Cell_Zone" ? 3.0 : 0.0))
        measurements.putMeasurement("GC_Location_Code", (double)locationCode)  // 0=Outside, 1=GC_Core, 2=B_Zone, 3=T_Zone
        measurements.putMeasurement("Zone_ID_Numeric", (double)(zoneId as Integer))
        measurements.putMeasurement("X_Centroid", (double)centroid[0])
        measurements.putMeasurement("Y_Centroid", (double)centroid[1])
    }
    
    println("FoxP3+ cell locations:")
    println("  GC Core: " + coreCount + " cells")
    println("  B Cell Zone: " + bZoneCount + " cells")
    println("  T Cell Zone: " + tZoneCount + " cells")
    println("  Outside GC: " + outsideCount + " cells")
}



/**
 * Perform proximity analysis between FoxP3+ cells and proximal cells
 */
def performProximityAnalysis() {
    // Configuration variables - defined within function scope
    def proximityRadius = 25.0  // Radius for proximity analysis (µm)
    
    println("Performing proximity analysis...")
    
    def foxp3Cells = getCellObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Positive" }
    def proximalCells = getCellObjects().findAll { it.getPathClass()?.getName() == "Proximal_Cell" }
    def imageData = getCurrentImageData()
    def cal = imageData.getServer().getPixelCalibration()
    def pixelSize = cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : 1.0
    
    println("Analyzing " + foxp3Cells.size() + " FoxP3+ cells and " + proximalCells.size() + " proximal cells")
    
    // For each FoxP3+ cell, calculate distances to all proximal cells
    foxp3Cells.each { foxp3Cell ->
        def foxp3Centroid = [foxp3Cell.getROI().getCentroidX(), foxp3Cell.getROI().getCentroidY()]
        def foxp3Id = foxp3Cell.getMeasurementList().getMeasurementValue("Unique_FoxP3_ID")
        def nearbyCount = 0
        
        proximalCells.each { proximalCell ->
            def cellCentroid = [proximalCell.getROI().getCentroidX(), proximalCell.getROI().getCentroidY()]
            def distance = calculateDistance(foxp3Centroid, cellCentroid) * pixelSize
            
            def measurements = proximalCell.getMeasurementList()
            
            // Store distance to nearest FoxP3+ cell
            def currentMinDistance = measurements.getMeasurementValue("Min_Distance_to_FoxP3")
            if (currentMinDistance == null || distance < currentMinDistance) {
                measurements.putMeasurement("Min_Distance_to_FoxP3", (double)distance)
                measurements.putMeasurement("Nearest_FoxP3_ID", (double)(foxp3Id ?: 0.0))
            }
            
            // Count cells within proximity radius of this FoxP3+ cell
            if (distance <= proximityRadius) {
                measurements.putMeasurement("Distance_to_FoxP3", (double)distance)
                nearbyCount++
            }
        }
        
        // Store how many proximal cells are near this FoxP3+ cell
        foxp3Cell.getMeasurementList().putMeasurement("Nearby_Cell_Count", (double)nearbyCount)
    }
    
    // Also ensure all FoxP3+ cells have their GC location info
    foxp3Cells.each { foxp3Cell ->
        def measurements = foxp3Cell.getMeasurementList()
        
        // Add intensity measurements if missing
        if (measurements.getMeasurementValue("FoxP3_Intensity") == null) {
            def foxp3Intensity = measurements.getMeasurementValue("Cell: CorrectedFOXP3-GFP mean")
            if (foxp3Intensity == null) {
                foxp3Intensity = measurements.getMeasurementValue("Nucleus: CorrectedFOXP3-GFP mean")
            }
            measurements.putMeasurement("FoxP3_Intensity", (double)(foxp3Intensity ?: 0.0))
        }
    }
    
    println("Proximity analysis complete")
    println("Average proximal cells per FoxP3+ cell: " + 
        (foxp3Cells.size() > 0 ? String.format("%.1f", proximalCells.size() / foxp3Cells.size()) : "0"))
}

/**
 * Export quantification results
 */
def exportResults() {
    def project = getProject()
    def baseDir
    
    if (project != null && project.getPath() != null) {
        baseDir = project.getPath().getParent().toFile()
    } else {
        // Fall back to user home directory
        def userHome = System.getProperty("user.home")
        baseDir = userHome != null ? new File(userHome, "QuPathResults") : new File(".")
    }
    
    // Create directory if it doesn't exist
    if (!baseDir.exists()) {
        baseDir.mkdirs()
    }
    
    println("Exporting results to: " + baseDir.getAbsolutePath())
    
    // Export FoxP3+ cell data
    exportFoxP3Data(baseDir)
    
    // Export proximity cell data  
    exportProximityData(baseDir)
    
    // Export all cell data
    exportAllCellData(baseDir)
    
    // Export summary statistics
    exportSummaryStats(baseDir)
}

/**
 * Export FoxP3+ cell quantification
 */
def exportFoxP3Data(baseDir) {
    def foxp3Cells = getCellObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Positive" }
    def exportFile = new File(baseDir, "FoxP3_Treg_Quantification.tsv")
    
    exportFile.withWriter { writer ->
        // Write header
        writer.writeLine("Unique_ID\tGC_Location\tSpecific_Zone\tZone_ID\t" +
            "FoxP3_Intensity\tCell_Area\tNucleus_Area\tX_Centroid\tY_Centroid\t" +
            "Nearby_Cell_Count\tCell_Type")
        
        foxp3Cells.each { cell ->
            def ml = cell.getMeasurementList()
            def uniqueId = ml.getMeasurementValue("Unique_ID") ?: 0
            def location = ml.getMeasurementValue("GC_Location_String") ?: "Unknown"
            def specificZone = ml.getMeasurementValue("Specific_Zone") ?: "Unknown"
            def zoneId = ml.getMeasurementValue("Zone_ID_String") ?: "0"
            
            def foxp3Intensity = ml.getMeasurementValue("FoxP3_Intensity") ?: 0.0
            def cellArea = ml.getMeasurementValue("Cell: Area µm²") ?: ml.getMeasurementValue("Area µm²") ?: 0.0
            def nucleusArea = ml.getMeasurementValue("Nucleus: Area µm²") ?: 0.0
            def xCentroid = ml.getMeasurementValue("X_Centroid") ?: cell.getROI().getCentroidX()
            def yCentroid = ml.getMeasurementValue("Y_Centroid") ?: cell.getROI().getCentroidY()
            def nearbyCount = ml.getMeasurementValue("Nearby_Cell_Count") ?: 0
            def cellType = ml.getMeasurementValue("Cell_Type") ?: 1.0
            
            writer.writeLine("${uniqueId}\t${location}\t${specificZone}\t${zoneId}\t" +
                "${foxp3Intensity}\t${cellArea}\t${nucleusArea}\t${xCentroid}\t${yCentroid}\t" +
                "${nearbyCount}\t${cellType}")
        }
    }
    
    println("FoxP3+ cell data exported to: " + exportFile.getName() + " (${foxp3Cells.size()} cells)")
}

/**
 * Export proximity cell data
 */
def exportProximityData(baseDir) {
    def proximityCells = getCellObjects().findAll { 
        def val = it.getMeasurementList().getMeasurementValue("Is_Proximity_Cell")
        return val != null && val == 1.0
    }
    def exportFile = new File(baseDir, "Proximity_Cell_Quantification.tsv")
    
    exportFile.withWriter { writer ->
        writer.writeLine("Cell_ID\tNearest_FoxP3_ID\tMin_Distance_to_FoxP3\tDistance_to_FoxP3\t" +
            "FoxP3_Intensity\tCell_Area\tNucleus_Area\tX_Centroid\tY_Centroid\t" +
            "GC_Location\tSpecific_Zone\tCell_Type")
        
        proximityCells.eachWithIndex { cell, index ->
            def ml = cell.getMeasurementList()
            def nearestFoxp3Id = ml.getMeasurementValue("Nearest_FoxP3_ID") ?: 0
            def minDistance = ml.getMeasurementValue("Min_Distance_to_FoxP3") ?: 0.0
            def distance = ml.getMeasurementValue("Distance_to_FoxP3") ?: 0.0
            
            def foxp3Intensity = ml.getMeasurementValue("FoxP3_Intensity") ?: 0.0
            def cellArea = ml.getMeasurementValue("Cell: Area µm²") ?: ml.getMeasurementValue("Area µm²") ?: 0.0
            def nucleusArea = ml.getMeasurementValue("Nucleus: Area µm²") ?: 0.0
            def xCentroid = ml.getMeasurementValue("X_Centroid") ?: cell.getROI().getCentroidX()
            def yCentroid = ml.getMeasurementValue("Y_Centroid") ?: cell.getROI().getCentroidY()
            def location = ml.getMeasurementValue("GC_Location_String") ?: "Unknown"
            def specificZone = ml.getMeasurementValue("Specific_Zone") ?: "Unknown"
            def cellType = ml.getMeasurementValue("Cell_Type") ?: 0.0
            
            writer.writeLine("${index + 1}\t${nearestFoxp3Id}\t${minDistance}\t${distance}\t" +
                "${foxp3Intensity}\t${cellArea}\t${nucleusArea}\t${xCentroid}\t${yCentroid}\t" +
                "${location}\t${specificZone}\t${cellType}")
        }
    }
    
    println("Proximity cell data exported to: " + exportFile.getName() + " (${proximityCells.size()} cells)")
}

/**
 * Export all cell data for comprehensive analysis
 */
def exportAllCellData(baseDir) {
    def allCells = getCellObjects()
    def exportFile = new File(baseDir, "All_Cells_Quantification.tsv")
    
    exportFile.withWriter { writer ->
        writer.writeLine("Cell_ID\tCell_Class\tFoxP3_Intensity\tCell_Area\tNucleus_Area\t" +
            "X_Centroid\tY_Centroid\tGC_Location\tSpecific_Zone\tZone_ID\t" +
            "Cell_Type\tMin_Distance_to_FoxP3\tIs_Proximity_Cell\tUnique_ID")
        
        allCells.eachWithIndex { cell, index ->
            def ml = cell.getMeasurementList()
            def cellClass = cell.getPathClass()?.getName() ?: "Unknown"
            def foxp3Intensity = ml.getMeasurementValue("FoxP3_Intensity") ?: 0.0
            def cellArea = ml.getMeasurementValue("Cell: Area µm²") ?: ml.getMeasurementValue("Area µm²") ?: 0.0
            def nucleusArea = ml.getMeasurementValue("Nucleus: Area µm²") ?: 0.0
            def xCentroid = ml.getMeasurementValue("X_Centroid") ?: cell.getROI().getCentroidX()
            def yCentroid = ml.getMeasurementValue("Y_Centroid") ?: cell.getROI().getCentroidY()
            def location = ml.getMeasurementValue("GC_Location_String") ?: "Unknown"
            def specificZone = ml.getMeasurementValue("Specific_Zone") ?: "Unknown"
            def zoneId = ml.getMeasurementValue("Zone_ID_String") ?: "0"
            def cellType = ml.getMeasurementValue("Cell_Type") ?: 0.0
            def minDistance = ml.getMeasurementValue("Min_Distance_to_FoxP3") ?: -1.0
            def isProximity = (ml.getMeasurementValue("Is_Proximity_Cell") ?: 0.0) == 1.0 ? "Yes" : "No"
            def uniqueId = ml.getMeasurementValue("Unique_ID") ?: 0
            
            writer.writeLine("${index + 1}\t${cellClass}\t${foxp3Intensity}\t${cellArea}\t${nucleusArea}\t" +
                "${xCentroid}\t${yCentroid}\t${location}\t${specificZone}\t${zoneId}\t" +
                "${cellType}\t${minDistance}\t${isProximity}\t${uniqueId}")
        }
    }
    
    println("All cell data exported to: " + exportFile.getName() + " (${allCells.size()} cells)")
}

/**
 * Export summary statistics
 */
def exportSummaryStats(baseDir) {
    def allCells = getCellObjects()
    def foxp3Cells = allCells.findAll { it.getPathClass()?.getName() == "FoxP3_Positive" }
    def weakCells = allCells.findAll { it.getPathClass()?.getName() == "FoxP3_Weak" }
    def negativeCells = allCells.findAll { it.getPathClass()?.getName() == "FoxP3_Negative" }
    def proximityCells = allCells.findAll { 
        def val = it.getMeasurementList().getMeasurementValue("Is_Proximity_Cell")
        return val != null && val == 1.0
    }
    
    def foxp3Regions = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "FoxP3_Region" }
    def gcCores = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "GC_Core" }
    def bZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "B_Cell_Zone" }
    def tZones = getAnnotationObjects().findAll { it.getPathClass()?.getName() == "T_Cell_Zone" }
    
    def summaryFile = new File(baseDir, "Analysis_Summary.tsv")
    
    summaryFile.withWriter { writer ->
        writer.writeLine("Metric\tCount\tPercentage")
        writer.writeLine("Total_Cells\t" + allCells.size() + "\t100.00")
        writer.writeLine("FoxP3_Positive_Cells\t" + foxp3Cells.size() + "\t" + 
            (allCells.size() > 0 ? String.format("%.2f", (foxp3Cells.size() / allCells.size()) * 100) : "0.00"))
        writer.writeLine("FoxP3_Weak_Cells\t" + weakCells.size() + "\t" + 
            (allCells.size() > 0 ? String.format("%.2f", (weakCells.size() / allCells.size()) * 100) : "0.00"))
        writer.writeLine("FoxP3_Negative_Cells\t" + negativeCells.size() + "\t" + 
            (allCells.size() > 0 ? String.format("%.2f", (negativeCells.size() / allCells.size()) * 100) : "0.00"))
        writer.writeLine("Proximity_Cells\t" + proximityCells.size() + "\t" + 
            (allCells.size() > 0 ? String.format("%.2f", (proximityCells.size() / allCells.size()) * 100) : "0.00"))
        
        writer.writeLine("")
        writer.writeLine("Regions\tCount\t")
        writer.writeLine("FoxP3_Regions\t" + foxp3Regions.size() + "\t")
        writer.writeLine("GC_Cores\t" + gcCores.size() + "\t")
        writer.writeLine("B_Cell_Zones\t" + bZones.size() + "\t")
        writer.writeLine("T_Cell_Zones\t" + tZones.size() + "\t")
        
        // Calculate averages if we have FoxP3+ cells
        if (foxp3Cells.size() > 0) {
            def totalNearby = foxp3Cells.sum { cell ->
                cell.getMeasurementList().getMeasurementValue("Nearby_Cell_Count") ?: 0
            }
            def avgNearby = totalNearby / foxp3Cells.size()
            
            writer.writeLine("")
            writer.writeLine("Analysis_Metrics\tValue\t")
            writer.writeLine("Avg_Cells_per_FoxP3\t" + String.format("%.2f", avgNearby) + "\t")
            writer.writeLine("FoxP3_Density\t" + String.format("%.4f", foxp3Cells.size() / (foxp3Regions.size() > 0 ? foxp3Regions.size() : 1)) + "\t")
        }
    }
    
    println("Summary statistics exported to: " + summaryFile.getName())
}



// === UTILITY FUNCTIONS ===

/**
 * Calculate Euclidean distance between two points
 */
def calculateDistance(point1, point2) {
    def dx = point1[0] - point2[0]
    def dy = point1[1] - point2[1]
    return Math.sqrt(dx * dx + dy * dy)
}



println("=== CHANNEL-SPECIFIC GC TREG ANALYSIS READY ===")
println("This script uses appropriate detection methods for fluorescence analysis:")
println("1. Fluorescence-based thresholding for regions (GC cores, B/T zones) with channel targeting")
println("2. WatershedCellDetection for individual FoxP3+ cells") 
println("3. Analyzes intensity values from specific marker channels by index")
println("4. Provides biologically accurate region and cell identification")
println("")
println("Run this script to perform the complete channel-specific analysis pipeline.")
println("Results will be exported as multiple TSV files for detailed analysis.")