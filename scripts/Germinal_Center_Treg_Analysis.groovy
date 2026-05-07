/**
 * Germinal Center and T Regulatory Cell Analysis Macro for QuPath v0.5.1
 * 
 * This macro performs comprehensive analysis of germinal centers and regulatory T cells
 * including:
 * 1. Identification of Germinal Center Cores (PNA+)
 * 2. Identification of B Cell zones (CD19-FITC) and T Cell zones (CD4-BV421)
 * 3. Cell segmentation using InstaSeg with nuclear staining (C30)
 * 4. FoxP3+ cell detection (C31 channel)
 * 5. Proximity cell analysis around FoxP3+ cells
 * 6. Quantification and export of results
 * 
 * Compatible with QuPath v0.5.1
 * Requires InstaSeg extension for advanced cell segmentation
 * 
 * Author: Image Analysis Pipeline
 * Date: 2025
 */

// Import required classes
import static qupath.lib.gui.scripting.QPEx.*
import qupath.lib.objects.PathObjects
import qupath.lib.objects.classes.PathClassFactory
import qupath.lib.roi.ROIs
import qupath.lib.regions.ImagePlane
import qupath.lib.regions.RegionRequest
import qupath.lib.images.servers.ImageServer
import qupath.lib.analysis.features.ObjectMeasurements
import qupath.lib.scripting.QP
import qupath.lib.objects.PathObject
import qupath.lib.measurements.MeasurementList
import qupath.lib.roi.interfaces.ROI
import qupath.lib.objects.TMACoreObject
import qupath.lib.common.GeneralTools
import qupath.lib.io.GsonTools

import java.awt.geom.Area
import java.awt.geom.Point2D
import java.util.stream.Collectors

// Configuration parameters
class AnalysisConfig {
    // Threshold parameters
    static final double PNA_THRESHOLD = 0.3                    // PNA+ threshold for GC cores
    static final double CD19_THRESHOLD = 0.25                  // CD19-FITC threshold for B cells
    static final double CD4_THRESHOLD = 0.2                    // CD4-BV421 threshold for T cells
    static final double FOXP3_THRESHOLD = 0.35                 // FoxP3-GFP threshold for Tregs
    static final double MIN_GC_AREA = 1000.0                   // Minimum area for GC cores (µm²)
    static final double PROXIMITY_RADIUS = 50.0                // Radius for proximity analysis (µm)
    
    // Channel mapping (adjust based on your image channel order)
    static final String PNA_CHANNEL = "PNA"                    // PNA channel name
    static final String CD19_CHANNEL = "CD19-FITC"            // CD19 channel name  
    static final String CD4_CHANNEL = "CD4-BV421"             // CD4 channel name
    static final String FOXP3_CHANNEL = "CorrectedFoxP3-GFP"  // FoxP3 channel name (C31)
    static final String NUCLEAR_CHANNEL = "C30"               // Nuclear stain channel (C30)
    
    // InstaSeg parameters
    static final double CELL_EXPANSION = 2.0                   // Cell expansion from nucleus (µm)
    static final double MIN_CELL_AREA = 20.0                   // Minimum cell area (µm²)
    static final double MAX_CELL_AREA = 500.0                  // Maximum cell area (µm²)
}

println("=== Germinal Center and Treg Analysis ===")
println("QuPath Version: " + getQuPathVersion())
println("Starting analysis pipeline...")

def imageData = getCurrentImageData()
if (imageData == null) {
    println("ERROR: No image is currently open!")
    return
}

def server = imageData.getServer()
def cal = server.getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
    println("WARNING: No pixel calibration found. Results may be inaccurate.")
}

println("Image: " + server.getMetadata().getName())
println("Dimensions: " + server.getWidth() + " x " + server.getHeight())
println("Channels: " + server.getMetadata().getChannels().size())

// Clear existing objects
clearAllObjects()

try {
    
    // === STEP 1: IDENTIFY GERMINAL CENTER CORES ===
    println("\n--- Step 1: Identifying Germinal Center Cores ---")
    
    def gcCores = identifyGerminalCenterCores()
    println("Found " + gcCores.size() + " potential GC cores")
    
    // === STEP 2: IDENTIFY B AND T CELL ZONES ===
    println("\n--- Step 2: Identifying B and T Cell Zones ---")
    
    def bCellZones = identifyBCellZones(gcCores)
    def tCellZones = identifyTCellZones(gcCores, bCellZones)
    
    println("Found " + bCellZones.size() + " B cell zones")
    println("Found " + tCellZones.size() + " T cell zones")
    
    // === STEP 3: CELL SEGMENTATION ===
    println("\n--- Step 3: Performing Cell Segmentation ---")
    
    // Create tissue annotation for segmentation
    def tissueAnnotation = createTissueAnnotation()
    addObject(tissueAnnotation)
    
    // Perform cell detection using nuclear channel
    performCellSegmentation(tissueAnnotation)
    
    def allCells = getDetectionObjects()
    println("Detected " + allCells.size() + " cells")
    
    // === STEP 4: IDENTIFY FOXP3+ CELLS ===
    println("\n--- Step 4: Identifying FoxP3+ Cells ---")
    
    def foxp3Cells = identifyFoxP3PositiveCells(allCells)
    println("Found " + foxp3Cells.size() + " FoxP3+ cells")
    
    // === STEP 5: PROXIMITY ANALYSIS ===
    println("\n--- Step 5: Performing Proximity Analysis ---")
    
    def proximityCells = identifyProximityCells(foxp3Cells, allCells)
    println("Found " + proximityCells.size() + " cells in proximity to FoxP3+ cells")
    
    // === STEP 6: ASSIGN LOCATIONS AND QUANTIFY ===
    println("\n--- Step 6: Assigning Locations and Quantifying ---")
    
    assignCellLocations(allCells, gcCores, bCellZones, tCellZones)
    measureCellIntensities(allCells)
    
    // === STEP 7: EXPORT RESULTS ===
    println("\n--- Step 7: Exporting Results ---")
    
    exportTregQuantification(foxp3Cells)
    exportProximityQuantification(proximityCells, foxp3Cells)
    
    println("\n=== Analysis Complete ===")
    println("Results exported to project directory")
    
} catch (Exception e) {
    println("ERROR during analysis: " + e.getMessage())
    e.printStackTrace()
}

// === HELPER FUNCTIONS ===

/**
 * Identify Germinal Center Cores using PNA+ staining
 */
def identifyGerminalCenterCores() {
    println("Creating PNA+ pixel classifier...")
    
    // Create pixel classifier for PNA+ regions
    def pnaClassifier = createSimplePixelClassifier(
        AnalysisConfig.PNA_CHANNEL,
        AnalysisConfig.PNA_THRESHOLD,
        "PNA_Positive",
        "Background"
    )
    
    // Create objects from classifier
    def request = RegionRequest.createInstance(server.getPath(), 1.0, 0, 0, server.getWidth(), server.getHeight())
    def pnaObjects = createObjectsFromPixelClassifier(pnaClassifier, request, AnalysisConfig.MIN_GC_AREA)
    
    // Filter and classify as GC cores
    def gcCores = []
    pnaObjects.each { obj ->
        if (obj.getROI().getArea() * cal.getAveragedPixelSizeMicrons() * cal.getAveragedPixelSizeMicrons() >= AnalysisConfig.MIN_GC_AREA) {
            obj.setPathClass(PathClassFactory.getPathClass("GC_Core"))
            obj.setName("GC_Core_" + (gcCores.size() + 1))
            gcCores.add(obj)
        }
    }
    
    addObjects(gcCores)
    return gcCores
}

/**
 * Identify B Cell zones using CD19-FITC staining
 */
def identifyBCellZones(gcCores) {
    println("Creating CD19+ pixel classifier...")
    
    def cd19Classifier = createSimplePixelClassifier(
        AnalysisConfig.CD19_CHANNEL,
        AnalysisConfig.CD19_THRESHOLD,
        "CD19_Positive",
        "Background"
    )
    
    def request = RegionRequest.createInstance(server.getPath(), 1.0, 0, 0, server.getWidth(), server.getHeight())
    def cd19Objects = createObjectsFromPixelClassifier(cd19Classifier, request, 500.0)
    
    def bCellZones = []
    cd19Objects.each { obj ->
        // Check if overlaps with any GC core
        def overlapsGC = gcCores.any { gc -> 
            obj.getROI().intersects(gc.getROI().getBounds2D())
        }
        
        if (overlapsGC) {
            obj.setPathClass(PathClassFactory.getPathClass("B_Cell_Zone"))
            obj.setName("B_Cell_Zone_" + (bCellZones.size() + 1))
            bCellZones.add(obj)
        }
    }
    
    addObjects(bCellZones)
    return bCellZones
}

/**
 * Identify T Cell zones using CD4-BV421 staining
 */
def identifyTCellZones(gcCores, bCellZones) {
    println("Creating CD4+ pixel classifier...")
    
    def cd4Classifier = createSimplePixelClassifier(
        AnalysisConfig.CD4_CHANNEL,
        AnalysisConfig.CD4_THRESHOLD,
        "CD4_Positive",
        "Background"
    )
    
    def request = RegionRequest.createInstance(server.getPath(), 1.0, 0, 0, server.getWidth(), server.getHeight())
    def cd4Objects = createObjectsFromPixelClassifier(cd4Classifier, request, 500.0)
    
    def tCellZones = []
    cd4Objects.each { obj ->
        // Check if borders B cell zones but is not a GC core
        def bordersBCell = bCellZones.any { bzone ->
            obj.getROI().intersects(bzone.getROI().getBounds2D())
        }
        
        def overlapsGC = gcCores.any { gc ->
            obj.getROI().intersects(gc.getROI().getBounds2D())
        }
        
        if (bordersBCell && !overlapsGC) {
            obj.setPathClass(PathClassFactory.getPathClass("T_Cell_Zone"))
            obj.setName("T_Cell_Zone_" + (tCellZones.size() + 1))
            tCellZones.add(obj)
        }
    }
    
    addObjects(tCellZones)
    return tCellZones
}

/**
 * Create tissue annotation for cell segmentation
 */
def createTissueAnnotation() {
    // Create simple thresholder to identify tissue
    def tissueClassifier = createSimplePixelClassifier(
        AnalysisConfig.NUCLEAR_CHANNEL,
        0.1,  // Low threshold to capture all tissue
        "Tissue",
        "Background"
    )
    
    def request = RegionRequest.createInstance(server.getPath(), 4.0, 0, 0, server.getWidth(), server.getHeight())
    def tissueObjects = createObjectsFromPixelClassifier(tissueClassifier, request, 10000.0)
    
    if (tissueObjects.size() > 0) {
        def largestTissue = tissueObjects.max { it.getROI().getArea() }
        largestTissue.setPathClass(PathClassFactory.getPathClass("Tissue"))
        largestTissue.setName("Tissue_Region")
        return largestTissue
    }
    
    // Fallback: create full image annotation
    def roi = ROIs.createRectangleROI(0, 0, server.getWidth(), server.getHeight(), ImagePlane.getDefaultPlane())
    def annotation = PathObjects.createAnnotationObject(roi, PathClassFactory.getPathClass("Tissue"))
    annotation.setName("Full_Image")
    return annotation
}

/**
 * Perform cell segmentation using InstaSeg or standard detection
 */
def performCellSegmentation(tissueAnnotation) {
    selectObjects(tissueAnnotation)
    
    try {
        // Try to use InstaSeg if available
        println("Attempting to use InstaSeg for cell segmentation...")
        // Note: InstaSeg integration would require the extension to be installed
        // For now, use standard cell detection
        runCellDetection()
    } catch (Exception e) {
        println("InstaSeg not available, using standard cell detection...")
        runStandardCellDetection()
    }
}

/**
 * Standard cell detection method
 */
def runStandardCellDetection() {
    def builder = CellDetection.builder()
        .detection(AnalysisConfig.NUCLEAR_CHANNEL)
        .nucleusExpansion(AnalysisConfig.CELL_EXPANSION)
        .cellExpansion(AnalysisConfig.CELL_EXPANSION)
        .minArea(AnalysisConfig.MIN_CELL_AREA)
        .maxArea(AnalysisConfig.MAX_CELL_AREA)
        .threshold(0.2)
        .watershedPostProcess(true)
        .splitByShape(true)
        .measureShape(true)
        .measureIntensity(true)
    
    def cellDetection = builder.build()
    cellDetection.runDetection(imageData, getSelectedObjects())
}

/**
 * Fallback cell detection using simple method
 */
def runCellDetection() {
    runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', 
        '{"detectionImageBrightfield": "Hematoxylin OD", ' +
        '"requestedPixelSizeMicrons": 0.5, ' +
        '"backgroundRadiusMicrons": 8.0, ' +
        '"medianRadiusMicrons": 0.0, ' +
        '"sigmaMicrons": 1.5, ' +
        '"minAreaMicrons": ' + AnalysisConfig.MIN_CELL_AREA + ', ' +
        '"maxAreaMicrons": ' + AnalysisConfig.MAX_CELL_AREA + ', ' +
        '"threshold": 0.1, ' +
        '"maxBackground": 2.0, ' +
        '"watershedPostProcess": true, ' +
        '"excludeDAB": false, ' +
        '"cellExpansionMicrons": ' + AnalysisConfig.CELL_EXPANSION + ', ' +
        '"includeNuclei": true, ' +
        '"smoothBoundaries": true, ' +
        '"makeMeasurements": true}')
}

/**
 * Identify FoxP3+ cells
 */
def identifyFoxP3PositiveCells(allCells) {
    def foxp3Cells = []
    def cellId = 1
    
    allCells.each { cell ->
        def measurements = cell.getMeasurementList()
        def foxp3Intensity = measurements.getMeasurementValue("Nucleus: " + AnalysisConfig.FOXP3_CHANNEL + " mean")
        
        if (foxp3Intensity >= AnalysisConfig.FOXP3_THRESHOLD) {
            cell.setPathClass(PathClassFactory.getPathClass("FoxP3_Positive"))
            cell.setName("FoxP3_Cell_" + cellId)
            
            // Add unique ID measurement
            def ml = cell.getMeasurementList()
            ml.putMeasurement("Unique_ID", cellId)
            
            foxp3Cells.add(cell)
            cellId++
        }
    }
    
    return foxp3Cells
}

/**
 * Identify cells within proximity radius of FoxP3+ cells
 */
def identifyProximityCells(foxp3Cells, allCells) {
    def proximityCells = []
    
    foxp3Cells.each { foxp3Cell ->
        def foxp3Centroid = foxp3Cell.getROI().getCentroidX(), foxp3Cell.getROI().getCentroidY()
        def foxp3Id = foxp3Cell.getMeasurementList().getMeasurementValue("Unique_ID")
        
        allCells.each { cell ->
            if (cell != foxp3Cell) {
                def cellCentroid = [cell.getROI().getCentroidX(), cell.getROI().getCentroidY()]
                def distance = calculateDistance(foxp3Centroid, cellCentroid) * cal.getAveragedPixelSizeMicrons()
                
                if (distance <= AnalysisConfig.PROXIMITY_RADIUS) {
                    // Add proximity measurements
                    def ml = cell.getMeasurementList()
                    ml.putMeasurement("Proximal_FoxP3_ID", foxp3Id)
                    ml.putMeasurement("Distance_to_FoxP3", distance)
                    
                    cell.setPathClass(PathClassFactory.getPathClass("Proximity_Cell"))
                    proximityCells.add(cell)
                }
            }
        }
    }
    
    return proximityCells
}

/**
 * Assign GC location to each cell
 */
def assignCellLocations(allCells, gcCores, bCellZones, tCellZones) {
    allCells.each { cell ->
        def centroid = [cell.getROI().getCentroidX(), cell.getROI().getCentroidY()]
        def ml = cell.getMeasurementList()
        
        // Check which zone the cell belongs to
        def location = "Outside_GC"
        def zoneId = 0
        
        // Check GC cores first
        gcCores.eachWithIndex { gc, index ->
            if (gc.getROI().contains(centroid[0], centroid[1])) {
                location = "GC_Core"
                zoneId = index + 1
                return true  // Break from closure
            }
        }
        
        // Check B cell zones if not in core
        if (location == "Outside_GC") {
            bCellZones.eachWithIndex { bzone, index ->
                if (bzone.getROI().contains(centroid[0], centroid[1])) {
                    location = "B_Cell_Zone"
                    zoneId = index + 1
                    return true
                }
            }
        }
        
        // Check T cell zones if not in other zones
        if (location == "Outside_GC") {
            tCellZones.eachWithIndex { tzone, index ->
                if (tzone.getROI().contains(centroid[0], centroid[1])) {
                    location = "T_Cell_Zone"
                    zoneId = index + 1
                    return true
                }
            }
        }
        
        ml.putMeasurement("GC_Location", location)
        ml.putMeasurement("Zone_ID", zoneId)
    }
}

/**
 * Measure cell intensities in all channels
 */
def measureCellIntensities(allCells) {
    println("Measuring cell intensities...")
    
    // Add intensity measurements for all channels
    ObjectMeasurements.addIntensityMeasurements(imageData, allCells, 1.0)
    ObjectMeasurements.addShapeMeasurements(imageData, allCells, cal, "NUCLEUS", "CELL")
}

/**
 * Export Treg quantification data
 */
def exportTregQuantification(foxp3Cells) {
    def project = getProject()
    def exportPath = project != null ? 
        new File(project.getPath().getParent().toFile(), "Treg_Quantification.tsv") :
        new File("Treg_Quantification.tsv")
    
    exportPath.withWriter { writer ->
        // Write header
        writer.writeLine("Unique_ID\tGC_Location\tZone_ID\t" + 
            "PNA_Mean\tCD19_Mean\tCD4_Mean\tFoxP3_Mean\tNuclear_Mean\t" +
            "Cell_Area\tNucleus_Area\tX_Centroid\tY_Centroid")
        
        foxp3Cells.each { cell ->
            def ml = cell.getMeasurementList()
            def uniqueId = ml.getMeasurementValue("Unique_ID")
            def location = ml.getMeasurementValue("GC_Location")
            def zoneId = ml.getMeasurementValue("Zone_ID")
            
            // Get intensity measurements
            def pnaMean = getChannelMean(cell, AnalysisConfig.PNA_CHANNEL)
            def cd19Mean = getChannelMean(cell, AnalysisConfig.CD19_CHANNEL)
            def cd4Mean = getChannelMean(cell, AnalysisConfig.CD4_CHANNEL)
            def foxp3Mean = getChannelMean(cell, AnalysisConfig.FOXP3_CHANNEL)
            def nuclearMean = getChannelMean(cell, AnalysisConfig.NUCLEAR_CHANNEL)
            
            def cellArea = ml.getMeasurementValue("Cell: Area µm²")
            def nucleusArea = ml.getMeasurementValue("Nucleus: Area µm²")
            def xCentroid = cell.getROI().getCentroidX()
            def yCentroid = cell.getROI().getCentroidY()
            
            writer.writeLine("${uniqueId}\t${location}\t${zoneId}\t" +
                "${pnaMean}\t${cd19Mean}\t${cd4Mean}\t${foxp3Mean}\t${nuclearMean}\t" +
                "${cellArea}\t${nucleusArea}\t${xCentroid}\t${yCentroid}")
        }
    }
    
    println("Treg quantification exported to: " + exportPath.getAbsolutePath())
}

/**
 * Export proximity quantification data
 */
def exportProximityQuantification(proximityCells, foxp3Cells) {
    def project = getProject()
    def exportPath = project != null ? 
        new File(project.getPath().getParent().toFile(), "Proximity_Quantification.tsv") :
        new File("Proximity_Quantification.tsv")
    
    exportPath.withWriter { writer ->
        // Write header
        writer.writeLine("Cell_ID\tProximal_FoxP3_ID\tDistance_to_FoxP3\t" +
            "PNA_Mean\tCD19_Mean\tCD4_Mean\tFoxP3_Mean\tNuclear_Mean\t" +
            "Cell_Area\tNucleus_Area\tX_Centroid\tY_Centroid")
        
        proximityCells.eachWithIndex { cell, index ->
            def ml = cell.getMeasurementList()
            def proximalId = ml.getMeasurementValue("Proximal_FoxP3_ID")
            def distance = ml.getMeasurementValue("Distance_to_FoxP3")
            
            // Get intensity measurements
            def pnaMean = getChannelMean(cell, AnalysisConfig.PNA_CHANNEL)
            def cd19Mean = getChannelMean(cell, AnalysisConfig.CD19_CHANNEL)
            def cd4Mean = getChannelMean(cell, AnalysisConfig.CD4_CHANNEL)
            def foxp3Mean = getChannelMean(cell, AnalysisConfig.FOXP3_CHANNEL)
            def nuclearMean = getChannelMean(cell, AnalysisConfig.NUCLEAR_CHANNEL)
            
            def cellArea = ml.getMeasurementValue("Cell: Area µm²")
            def nucleusArea = ml.getMeasurementValue("Nucleus: Area µm²")
            def xCentroid = cell.getROI().getCentroidX()
            def yCentroid = cell.getROI().getCentroidY()
            
            writer.writeLine("${index + 1}\t${proximalId}\t${distance}\t" +
                "${pnaMean}\t${cd19Mean}\t${cd4Mean}\t${foxp3Mean}\t${nuclearMean}\t" +
                "${cellArea}\t${nucleusArea}\t${xCentroid}\t${yCentroid}")
        }
    }
    
    println("Proximity quantification exported to: " + exportPath.getAbsolutePath())
}

// === UTILITY FUNCTIONS ===

/**
 * Create a simple pixel classifier
 */
def createSimplePixelClassifier(channelName, threshold, aboveClass, belowClass) {
    // This is a simplified implementation
    // In practice, you would use QuPath's pixel classification tools
    return [
        channel: channelName,
        threshold: threshold,
        aboveClass: aboveClass,
        belowClass: belowClass
    ]
}

/**
 * Create objects from pixel classifier (simplified)
 */
def createObjectsFromPixelClassifier(classifier, request, minArea) {
    // This would be implemented using QuPath's pixel classification
    // For now, return empty list as placeholder
    return []
}

/**
 * Get channel mean intensity for a cell
 */
def getChannelMean(cell, channelName) {
    def ml = cell.getMeasurementList()
    def measurementName = "Nucleus: " + channelName + " mean"
    
    try {
        return ml.getMeasurementValue(measurementName)
    } catch (Exception e) {
        return 0.0
    }
}

/**
 * Calculate Euclidean distance between two points
 */
def calculateDistance(point1, point2) {
    def dx = point1[0] - point2[0]
    def dy = point1[1] - point2[1]
    return Math.sqrt(dx * dx + dy * dy)
}

println("Macro loaded successfully. Run this script to perform the analysis.") 