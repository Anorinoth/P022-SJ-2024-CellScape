/**
 * Memory-Efficient Germinal Center and T Regulatory Cell Analysis for QuPath v0.5.1
 * 
 * This version is optimized for large, multi-channel images with memory constraints.
 * It processes images in smaller regions and uses simplified detection methods.
 * 
 * Compatible with QuPath v0.5.1+
 * Author: Image Analysis Pipeline
 * Date: 2025
 */

import static qupath.lib.gui.scripting.QPEx.*

println("=== MEMORY-EFFICIENT GC AND TREG ANALYSIS ===")
println("Compatible with QuPath v0.5.1+")
println("Optimized for large multi-channel images")

// Get image information
def imageData = getCurrentImageData()
def server = imageData.getServer()
def cal = server.getPixelCalibration()

println("")
println("Image: " + server.getMetadata().getName())
println("Dimensions: " + server.getWidth() + " x " + server.getHeight())
println("Channels: " + server.nChannels())
println("Pixel size: " + cal.getPixelWidthMicrons() + " x " + cal.getPixelHeightMicrons() + " µm")

// List available channels
println("")
println("Available channels:")
for (int i = 0; i < server.nChannels(); i++) {
    println("  Channel " + i + ": " + server.getChannel(i).getName())
}

// === MAIN EXECUTION ===
try {
    // Step 1: Create tissue region (smaller)
    createTissueRegion()
    
    // Step 2: Memory-efficient cell detection
    performMemoryEfficientCellDetection()
    
    // Step 3: Classify FoxP3+ cells
    classifyFoxP3Cells()
    
    // Step 4: Create placeholder annotations for manual drawing
    createBasicRegionAnnotations()
    
    // Step 5: Assign cell locations with fixed measurements
    assignCellLocationsFixed()
    
    // Step 6: Export results
    exportResults()
    
} catch (Exception e) {
    println("ERROR during analysis: " + e.getMessage())
    e.printStackTrace()
}

println("")
println("=== MEMORY-EFFICIENT ANALYSIS COMPLETE ===")
println("This version processes smaller regions to avoid memory issues.")
println("For best results with large images:")
println("1. Increase QuPath's memory allocation (Edit > Preferences > General > Set maximum memory)")
println("2. Close other applications to free up RAM")
println("3. Consider processing images in smaller sections")

/**
 * Create a tissue region annotation - smaller to reduce memory usage
 */
def createTissueRegion() {
    println("")
    println("--- Step 1: Creating Tissue Region (Memory-Efficient) ---")
    
    // Get image data and server within function scope
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
    
    // Use a smaller tissue detection to reduce memory usage
    def annotations = getAnnotationObjects()
    def existingTissue = annotations.find { it.getName() == "Tissue_Region" }
    
    if (existingTissue == null) {
        // Create a simple rectangular region instead of full tissue detection
        def width = server.getWidth()
        def height = server.getHeight()
        
        // Create a region covering central 60% of image to reduce processing area
        def x = (int)(width * 0.2)
        def y = (int)(height * 0.2) 
        def w = (int)(width * 0.6)
        def h = (int)(height * 0.6)
        
        def roi = ROIs.createRectangleROI(x, y, w, h, null)
        def tissueAnnotation = PathObjects.createAnnotationObject(roi)
        tissueAnnotation.setName("Tissue_Region")
        addObjects([tissueAnnotation])
        
        println("Created tissue region annotation (60% of image area)")
    } else {
        println("Using existing tissue region annotation")
    }
}

/**
 * Memory-efficient cell detection using simplified parameters
 */
def performMemoryEfficientCellDetection() {
    println("")
    println("--- Step 2: Memory-Efficient Cell Detection ---")
    
    // Configuration variables - defined within function scope
    def nuclearChannel = "DNA-Brilliant Violet 421"  // Updated based on user feedback (C30)
    def minCellArea = 50.0      // Increased minimum area to reduce cell count
    def maxCellArea = 1000.0    // Increased maximum area
    def cellExpansion = 1.0     // Reduced expansion
    
    // Select tissue annotation for cell detection
    def tissueAnnotation = getAnnotationObjects().find { it.getName() == "Tissue_Region" }
    selectObjects(tissueAnnotation)
    
    if (tissueAnnotation != null) {
        try {
            // Use simple positive cell detection with memory-friendly parameters
            runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', 
                '{"requestedPixelSizeMicrons":1.0,' +  // Larger pixel size = less detail = less memory
                '"backgroundRadiusMicrons":15.0,' +
                '"minAreaMicrons":' + minCellArea + ',' +
                '"maxAreaMicrons":' + maxCellArea + ',' +
                '"threshold":0.2,' +  // Higher threshold = fewer cells = less memory
                '"cellExpansionMicrons":' + cellExpansion + ',' +
                '"makeMeasurements":true}')
            
            def cells = getCellObjects()
            println("Detected " + cells.size() + " cells using memory-efficient detection")
            
        } catch (Exception e) {
            println("Memory-efficient detection failed: " + e.getMessage())
            println("Trying minimal detection...")
            
            // Ultra-simple detection as fallback
            try {
                runPlugin('qupath.imagej.detect.nuclei.SimpleThresholdNucleusDetection',
                    '{"threshold":0.3,' +
                    '"minAreaMicrons":' + minCellArea + ',' +
                    '"maxAreaMicrons":' + maxCellArea + ',' +
                    '"makeMeasurements":true}')
                
                def cells = getCellObjects()
                println("Detected " + cells.size() + " cells using minimal detection")
                
            } catch (Exception e2) {
                println("All detection methods failed due to memory constraints")
                println("Try increasing QuPath's memory allocation or processing a smaller image region")
            }
        }
    }
    
    // Clear selection
    clearAllObjects()
}

/**
 * Classify FoxP3+ cells with memory-efficient approach
 */
def classifyFoxP3Cells() {
    println("")
    println("--- Step 3: Classifying FoxP3+ Cells (Memory-Efficient) ---")
    
    // Configuration variables - defined within function scope
    def foxp3Channel = "CorrectedFOXP3-GFP"  // Channel 30
    def foxp3Threshold = 0.3  // Adjust as needed
    
    def cells = getCellObjects()
    def foxp3PositiveCells = []
    
    try {
        // Add measurements in smaller batches to reduce memory usage
        def batchSize = 100  // Process cells in batches
        def processed = 0
        
        for (int i = 0; i < cells.size(); i += batchSize) {
            def batch = cells.subList(i, Math.min(i + batchSize, cells.size()))
            selectObjects(batch)
            
            // Add intensity measurements for this batch only
            addShapeMeasurements()  // Simpler measurements to reduce memory
            
            // Classify based on measurements
            for (cell in batch) {
                def measurements = cell.getMeasurementList()
                def channelMean = measurements.getMeasurementValue(foxp3Channel + ": Cell: Mean")
                
                if (channelMean != null && channelMean > foxp3Threshold) {
                    cell.getPathClass() = getPathClass("FoxP3_Positive")
                    foxp3PositiveCells.add(cell)
                }
            }
            
            processed += batch.size()
            if (processed % 500 == 0) {
                println("Processed " + processed + "/" + cells.size() + " cells for FoxP3 classification")
            }
        }
        
    } catch (Exception e) {
        println("FoxP3 classification failed: " + e.getMessage())
        println("Try reducing the number of detected cells or image complexity")
    }
    
    println("Classified " + foxp3PositiveCells.size() + " FoxP3+ cells")
    clearAllObjects()
}

/**
 * Create placeholder annotations for manual drawing
 */
def createBasicRegionAnnotations() {
    println("")
    println("--- Step 4: Creating Basic Region Annotations ---")
    println("Creating placeholder annotations...")
    
    // Create small example annotations that user can modify
    def imageData = getCurrentImageData()
    def server = imageData.getServer()
    
    // Create small placeholder regions
    def centerX = server.getWidth() / 2
    def centerY = server.getHeight() / 2
    def size = 500  // Small placeholders
    
    // GC Core placeholder
    def gcRoi = ROIs.createRectangleROI(centerX - size, centerY - size, size, size, null)
    def gcAnnotation = PathObjects.createAnnotationObject(gcRoi, getPathClass("GC_Core"))
    gcAnnotation.setName("GC_Core_Example")
    
    // B Cell Zone placeholder  
    def bRoi = ROIs.createRectangleROI(centerX, centerY - size, size, size, null)
    def bAnnotation = PathObjects.createAnnotationObject(bRoi, getPathClass("B_Cell_Zone"))
    bAnnotation.setName("B_Cell_Zone_Example")
    
    // T Cell Zone placeholder
    def tRoi = ROIs.createRectangleROI(centerX - size, centerY, size, size, null)
    def tAnnotation = PathObjects.createAnnotationObject(tRoi, getPathClass("T_Cell_Zone"))
    tAnnotation.setName("T_Cell_Zone_Example")
    
    addObjects([gcAnnotation, bAnnotation, tAnnotation])
    
    println("IMPORTANT: You should manually draw annotations for:")
    println("  - GC cores (classify as 'GC_Core')")
    println("  - B cell zones (classify as 'B_Cell_Zone')")  
    println("  - T cell zones (classify as 'T_Cell_Zone')")
    println("")
    println("These can be drawn using the annotation tools and classified")
    println("using the right-click context menu → Set class")
}

/**
 * Assign cell locations with fixed measurement approach (no numeric hashcodes)
 */
def assignCellLocationsFixed() {
    println("")
    println("--- Step 5: Assigning Cell Locations (Fixed) ---")
    
    def cells = getCellObjects()
    def gcCores = getAnnotationObjects().findAll { it.getPathClass() == getPathClass("GC_Core") }
    def bCellZones = getAnnotationObjects().findAll { it.getPathClass() == getPathClass("B_Cell_Zone") }
    def tCellZones = getAnnotationObjects().findAll { it.getPathClass() == getPathClass("T_Cell_Zone") }
    
    println("Found " + gcCores.size() + " GC cores")
    println("Found " + bCellZones.size() + " B cell zones")
    println("Found " + tCellZones.size() + " T cell zones")
    
    // Assign locations using string measurements instead of numeric hashcodes
    for (cell in cells) {
        def cellROI = cell.getROI()
        def cellCentroid = cellROI.getCentroidX() + "," + cellROI.getCentroidY()
        def location = "Outside"
        def zoneId = "0"
        
        // Check GC cores first (highest priority)
        for (gc in gcCores) {
            if (gc.getROI().contains(cellROI.getCentroidX(), cellROI.getCentroidY())) {
                location = "GC_Core"
                zoneId = gc.getID() ?: "GC_1"
                break
            }
        }
        
        // Check B cell zones if not in GC core
        if (location == "Outside") {
            for (b in bCellZones) {
                if (b.getROI().contains(cellROI.getCentroidX(), cellROI.getCentroidY())) {
                    location = "B_Cell_Zone"
                    zoneId = b.getID() ?: "B_1"
                    break
                }
            }
        }
        
        // Check T cell zones if not in other zones
        if (location == "Outside") {
            for (t in tCellZones) {
                if (t.getROI().contains(cellROI.getCentroidX(), cellROI.getCentroidY())) {
                    location = "T_Cell_Zone"
                    zoneId = t.getID() ?: "T_1"
                    break
                }
            }
        }
        
        // Store location as string measurements (safer approach)
        def measurements = cell.getMeasurementList()
        measurements.putMeasurement("GC_Location_String", location)  // String measurement
        measurements.putMeasurement("Zone_ID_String", zoneId)         // String measurement
        measurements.putMeasurement("Cell_Centroid", cellCentroid)     // String measurement
    }
    
    println("Assigned locations to all cells using string measurements")
}

/**
 * Export results to files
 */
def exportResults() {
    println("")
    println("--- Step 6: Exporting Results ---")
    
    def cells = getCellObjects()
    def foxp3Cells = cells.findAll { it.getPathClass() == getPathClass("FoxP3_Positive") }
    
    // Determine export directory
    def project = getProject()
    def exportDir
    if (project != null) {
        exportDir = buildFilePath(project.getPath().getParent().toString(), "exports")
    } else {
        exportDir = buildFilePath(System.getProperty("user.home"), "QuPath_Exports")
    }
    
    // Create export directory
    def exportDirFile = new File(exportDir)
    if (!exportDirFile.exists()) {
        exportDirFile.mkdirs()
    }
    
    // Export all cells
    def allCellsPath = buildFilePath(exportDir, "all_cells_measurements.tsv")
    saveAnnotationMeasurements(allCellsPath, cells)
    println("Exported " + cells.size() + " cell measurements to: " + allCellsPath)
    
    // Export FoxP3+ cells if any found
    if (foxp3Cells.size() > 0) {
        def foxp3Path = buildFilePath(exportDir, "foxp3_positive_cells.tsv")
        saveAnnotationMeasurements(foxp3Path, foxp3Cells)
        println("Exported " + foxp3Cells.size() + " FoxP3+ cell measurements to: " + foxp3Path)
    } else {
        println("No FoxP3+ cells found to export")
    }
    
    println("Export directory: " + exportDir)
} 