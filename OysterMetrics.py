
#import modules
import arcpy, os, sys
arcpy.CheckOutExtension("Spatial")
from arcpy.sa import *
from arcpy import env
import arcpy.mp

import pandas as pd
import numpy as np

#set relative paths
scriptPath = sys.argv[0]
scriptWS = os.path.dirname(scriptPath) #gives you script folder
rootWS = os.path.dirname(scriptWS) #goes up one folder into project folder
dataWS = os.path.join(rootWS, "data") #go into data folder from project folder'
resultsWS = os.path.join(rootWS, "results") #go into results folder from project folder

mxd = arcpy.mp.ArcGISProject("CURRENT")
mxd.relativePaths = True

#set environmental variables
arcpy.env.workspace = dataWS
arcpy.env.overwriteOutput = True

## USER INPUTS
input_scratchWS = arcpy.GetParameterAsText(0) # input scratch folder
arcpy.env.scratchWorkspace = input_scratchWS #set the environmental variable for the scratch folder after knowing the scratch folder
#input files
input_dsm = arcpy.GetParameterAsText(1) #input oyster dsm
input_reefPolygon = arcpy.GetParameterAsText(2) # line along reef edge
input_reefName = arcpy.GetParameterAsText(3) # enter the reef name, no spaces
input_sampleNum = arcpy.GetParameterAsText(4) # enter the number of samples you want to pick

# Determine Parent Arc Folder (assuming scratch nested within - change this to more adaptable to other users later)
ArcPar = os.path.dirname(input_scratchWS) # get arc parent folder
arcpy.AddMessage(ArcPar)

# Create subfolder based on reef name
## allows you to save files for each reef rather than overwriting each time
scratchWS = input_scratchWS + "\\" + input_reefName # set scratch workspace to subfolder with reef name
if not os.path.exists(scratchWS):
    os.makedirs(scratchWS)
arcpy.AddMessage("folder for this reef created at: " + scratchWS)

# Create path to gdb
gdb = ArcPar + "\\OysterMetrics.gdb" # get path of geodatabase

# 1. Normalize DSM plane

# a. Buffer input reef polygon (outside)
outBufferReefEdge = gdb + "\\outerBufferReef" # save within geodatabase
arcpy.Buffer_analysis(input_reefPolygon, outBufferReefEdge, "0.25 Meters")
arcpy.AddMessage("(1a) outer reef edge buffer created")

# b. Convert Outer Buffer Polygon to line
outBufferReefLine = gdb + "\\outerBufferLine"
arcpy.PolygonToLine_management(outBufferReefEdge, outBufferReefLine) # create outer buffer line
arcpy.AddMessage("(1b) outer reef buffer converted to line")

# c. Generate Samples along line
normalizeSamplePoints = gdb + "\\normalizeSamplePts"
arcpy.GeneratePointsAlongLines_management(outBufferReefLine, normalizeSamplePoints,
    'PERCENTAGE', Percentage = 10)
arcpy.AddMessage("(1c) normalizing samples created")

# d. Copy Features to prepare points with object IDs
normalizeSamplePoints_CF = gdb + "\\normalizeSamplePtsCF"
arcpy.CopyFeatures_management(normalizeSamplePoints, normalizeSamplePoints_CF)

# e. Extract values to points
normalizeSamplePts_DSMAttached = gdb + "\\normalizeSamplePts_DSMAttached"
ExtractValuesToPoints(normalizeSamplePoints_CF, input_dsm, normalizeSamplePts_DSMAttached)
arcpy.AddMessage("(1d) normalizing samples extracted DSM values")

# f. Create Plane from Sample Point DSM values)
normalPlane = scratchWS + "\\normalPlane.tif" # delete previous version
if os.path.exists(normalPlane):
    os.remove(normalPlane)

normalPlane = scratchWS + "\\normalPlane.tif"
outTrend = Trend(normalizeSamplePts_DSMAttached, "RASTERVALU", order = 1)
outTrend.save(normalPlane)
arcpy.AddMessage("(1e) normalized plane generated and saved:" + normalPlane)

# g. Subtract normalized plane from DSM to create normalized dsm
normalDSM = scratchWS + "\\normalDSM.tif" # delete previous version
if os.path.exists(normalDSM):
    os.remove(normalDSM)

normalDSM = scratchWS + "\\normalDSM.tif"
outMinus = Raster(input_dsm) - Raster(normalPlane)
outMinus.save(normalDSM)
arcpy.AddMessage("(1f) DSM normalized to height and saved:" + normalDSM)

### DSM NORMALIZED ####

# 2. Downsample normalized DSM to 2xGSD
normalDSMx2 = scratchWS + "\\normalDSMx2.tif"
gsdx = arcpy.GetRasterProperties_management(normalDSM, property_type = "CELLSIZEX") # original gsdx
gsdy = arcpy.GetRasterProperties_management(normalDSM, property_type = "CELLSIZEY") # original gsdy
gsd = float(gsdx[0])*2 # mulitply GSD (x) by 2 for GSDx2
arcpy.Resample_management(normalDSM, normalDSMx2, gsd, "BILINEAR") # use gsd as an input for cell size
arcpy.AddMessage("(2) normalized DSM resampled to 2xGSD: " + str(gsd))

### DSM DOWNSAMPLED ###


## 3. GENERATE SAMPLES
# a. Crest

# 1. Buffer further inside the reef edge for crest (to avoid overlap with reef edge)
BufferReefCrest = gdb + "\\crestBufferReef"
arcpy.Buffer_analysis(input_reefPolygon, BufferReefCrest, "-0.7 Meters") # 2x the diagonal
arcpy.AddMessage("(3a, 1) crest reef buffer created")

# 2. Convert Crest Buffer Polygon to line
BufferReefCrestLine = gdb + "\\crestBufferLine"
arcpy.PolygonToLine_management(BufferReefCrest, BufferReefCrestLine) # create crest buffer line
arcpy.AddMessage("(3a, 2) crest reef buffer converted to line")

# 3. Generate transects across reef edge inner buffer (Step 8)
crestTransects = gdb + "\\crestTransects"
arcpy.GenerateTransectsAlongLines_management(BufferReefCrestLine, crestTransects,
    '0.1 Meters', '100 Meters') # create transects at each 0.1m interval 100m in length (big enough to cross any oyster reef)
arcpy.AddMessage("(3a, 3) reef crest transects created")

# 4. Clip transects to area of reef edge
crestTransectsClip = gdb + "\\crestTransectsClip"
arcpy.Clip_analysis(crestTransects, BufferReefCrest, crestTransectsClip) # clip transects to reef edge boundary
arcpy.AddMessage("(3a, 4) reef crest transects clipped")

# 5. Copy features to attach shape length information for export
crestTransectsClipCF = gdb + "\\crestTransectsClipCF"
arcpy.CopyFeatures_management(crestTransectsClip, crestTransectsClipCF)
arcpy.AddMessage("(3a, 5) transects copied")

# 6. Export Attribute Table to find longest line
transectLines = scratchWS + "\\transectLines.csv"
arcpy.CopyRows_management(crestTransectsClipCF, transectLines)
arcpy.AddMessage("(3a, 6) reef crest transect table exported")

# 7. Select Max Length of longest line
transectLines_df  = pd.read_csv(transectLines) # read transect length csv
maxlength = pd.to_numeric(transectLines_df['Shape_Length']).max() # grab max length
arcpy.AddMessage("(3a, 7) the longest line is " + maxlength.astype(str)) # print longest length

# 8. Select Object ID of longest Line
ismax = pd.to_numeric(transectLines_df['Shape_Length']) == maxlength  # make True/False column for longest line
maxID = (transectLines_df[ismax]).iloc[0]['OBJECTID'] # filter for True rows, select OID of the first row
arcpy.AddMessage("(3a, 8) the longest line has Object ID:" + maxID.astype(str))

# 9. Select crest transect by OID and copy to new feature class
maxCrestTransectFile = gdb + "\\maxCrestTransect"
query = 'OBJECTID =' + maxID.astype(str)  # build the query
arcpy.AddMessage("the query is: " + query) # print the query
arcpy.management.MakeFeatureLayer(crestTransectsClipCF, 'maxCrestTransect', query, None) # make layer with longest transect
arcpy.CopyFeatures_management('maxCrestTransect', maxCrestTransectFile) # save that layer to a shapefile
arcpy.AddMessage("(3a, 9) copy features")

# 10. Generate Samples along reef crest
reefCrestSamplePoints = gdb + "\\reefCrestSamplePoints"

maxSamples = int(maxlength/0.4) # extract max # of samples possible given longest crest line
arcpy.AddMessage(maxSamples)
numSamplesChosen = int(input_sampleNum) # convert input sample # to integer
numSamplesUsed_Crest = 0 # create dummy variable to fill in for the number of samples actually used

if(numSamplesChosen > maxSamples): # if the # of samples chosen exceeds the # possible
    percentageSep = 100/maxSamples # how far the samples should be spaced %
    arcpy.management.GeneratePointsAlongLines(maxCrestTransectFile, reefCrestSamplePoints,
        'PERCENTAGE', None, percentageSep, None) # take crest samples spaced by percentageSep
    numSamplesUsed_Crest = maxSamples # update num samples used
    arcpy.AddWarning("Your reef is not big enough for" + str(numSamplesChosen) + " samples. Using " +
                str(maxSamples) + " samples")
    arcpy.AddMessage("(3a, 10) " + str(maxSamples) + " crest samples created")

if (numSamplesChosen <= maxSamples): # if the # of samples chosen is less than or equal to # possible
    percentageSep = 100/numSamplesChosen # how far the samples should be spaced by %
    arcpy.management.GeneratePointsAlongLines(maxCrestTransectFile, reefCrestSamplePoints,
        'PERCENTAGE', None, percentageSep, None) # take edge samples spaced by percentageSep
    numSamplesUsed_Crest = numSamplesChosen # update num samples used
    arcpy.AddMessage("(3a, 10) " + str(numSamplesChosen) + " crest samples created")

# 11. Buffer reef crest Samples
reefCrestSamplesBuffered = gdb + "\\reefCrestSamplesBuffered"
arcpy.Buffer_analysis(reefCrestSamplePoints, reefCrestSamplesBuffered, "0.125 Meters") # radius of 0.125 to get diameter (and square side) of 0.25m
arcpy.AddMessage("(3a, 11) crest sample points buffered")

# 12. Make Buffer Plots Square
reefCrestSamplesSquareBuf = gdb + "\\reefCrestSamplesSquareBuf"
arcpy.MinimumBoundingGeometry_management(reefCrestSamplesBuffered, reefCrestSamplesSquareBuf,
    "ENVELOPE")
arcpy.AddMessage("(3a, 12) crest sample points converted to square buffer")

# 13. Copy out features from gdb to spatially reference plots later
reefCrestSamplePlots = scratchWS + "//reefCrestSamplePlots.shp"
arcpy.CopyFeatures_management(reefCrestSamplesSquareBuf, reefCrestSamplePlots)

## b. Edge Samples

# 1. Buffer input reef polygon (inside)
innerBufferReefEdge = gdb + "\\innerBufferReef"
arcpy.Buffer_analysis(input_reefPolygon, innerBufferReefEdge, "-0.35 Meters") # buffer in by diagonal length
arcpy.AddMessage("(3b, 1) reef edge buffer created")

# 2. Convert Inner Buffer Polygon to line
innerBufferReefLine = gdb + "\\innerBufferLine"
arcpy.PolygonToLine_management(innerBufferReefEdge, innerBufferReefLine) # create inner buffer line
arcpy.AddMessage("(3b, 2) reef edge buffer converted to line")

# 3. Generate Samples along inner reef edge
numSamplesUsed_Edge = numSamplesUsed_Crest # set the number of samples to use as the same as crest
percentageSep_Edge = 100/numSamplesUsed_Crest # spacing of samples in %
reefEdgeSamplePoints = gdb + "\\reefEdgeSamplePoints"
arcpy.management.GeneratePointsAlongLines(innerBufferReefLine, reefEdgeSamplePoints,
    'PERCENTAGE', None, percentageSep_Edge, None)
arcpy.AddMessage("(3b, 3) " + str(numSamplesUsed_Edge) + " edge samples created")

# 4. Buffer reef edge Samples
reefEdgeSamplesBuffered = gdb + "\\reefEdgeSamplesBuffered"
arcpy.Buffer_analysis(reefEdgeSamplePoints, reefEdgeSamplesBuffered, "0.25 Meters")
arcpy.AddMessage("(3b, 4) edge sample points buffered")

# 5. Make Buffer Plots Square
reefEdgeSamplesSquareBuf = gdb + "\\reefEdgeSamplesSquareBuf"
arcpy.MinimumBoundingGeometry_management(reefEdgeSamplesBuffered, reefEdgeSamplesSquareBuf,
    "ENVELOPE")
arcpy.AddMessage("(3b, 5) edge sample points converted to square buffer")

# 6. Copy out features from gdb to spatially reference plots later
reefEdgeSamplePlots = scratchWS + "//reefEdgeSamplePlots.shp"
arcpy.CopyFeatures_management(reefEdgeSamplesSquareBuf, reefEdgeSamplePlots)

## METRICS
arcpy.AddMessage("----Begin Metrics Calculations-----")
tables = scratchWS + "\\tables" # output location for metric tables
if not os.path.exists(tables):
    os.makedirs(tables)

# 1. Elevation
# a. Edge
elevationEdgeDBF = tables + "\\elevationEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, "OBJECTID", normalDSMx2, elevationEdgeDBF,
    "DATA", "MEAN")

elevationEdge = tables + "\\elevationEdge.csv"
arcpy.conversion.TableToTable(elevationEdgeDBF, tables, "elevationEdge.csv") # convert to csv format
arcpy.AddMessage("(1a) zonal statistics: elevation at reef edge exported")

# b. Crest
elevationCrestDBF = tables + "\\elevationCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, "OBJECTID", normalDSMx2, elevationCrestDBF,
    "DATA", "MEAN")

elevationCrest = tables + "\\elevationCrest.csv"
arcpy.conversion.TableToTable(elevationCrestDBF, tables, "elevationCrest.csv") # convert to csv format
arcpy.AddMessage("(1b) zonal statistics: elevation at reef crest exported")

# 2. Mean Slope

# Create Slope Raster
slopeDSM = tables + "//slopeDSM.tif"
outSlopeDSM = arcpy.sa.Slope(normalDSMx2, "DEGREE", 1, "PLANAR", "METER")
outSlopeDSM.save(slopeDSM) # save the slope raster
arcpy.AddMessage("(2) slope raster generated")

# a. Edge
slopeEdgeDBF = tables + "\\slopeEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeEdgeDBF,
    "DATA", "MEAN")

slopeEdge = tables + "\\slopeEdge.csv"
arcpy.conversion.TableToTable(slopeEdgeDBF, tables, "slopeEdge.csv") # convert to csv format
arcpy.AddMessage("(2a) zonal statistics: slope at reef edge exported")

# b. Slope
slopeCrestDBF = tables + "\\slopeCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeCrestDBF,
    "DATA", "MEAN")

slopeCrest = tables + "\\slopeCrest.csv"
arcpy.conversion.TableToTable(slopeCrestDBF, tables, "slopeCrest.csv") # convert to csv format
arcpy.AddMessage("(2b) zonal statistics: slope at reef crest exported")

# 3. Standard Deviation of Slope
# a. Edge
slopeStdEdgeDBF = tables + "\\slopeStdEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeStdEdgeDBF,
    "DATA", "STD")
slopeStdEdge = tables + "\\slopeStdEdge.csv"
arcpy.conversion.TableToTable(slopeStdEdgeDBF, tables, "slopeStdEdge.csv") # convert to csv format
arcpy.AddMessage("(3a) zonal statistics: slope standard deviation at reef edge exported")

# b. Crest
slopeStdCrestDBF = tables + "\\slopeStdCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeStdCrestDBF,
    "DATA", "STD")
slopeStdCrest = tables + "\\slopeStdCrest.csv"
arcpy.conversion.TableToTable(slopeStdCrestDBF, tables, "slopeStdCrest.csv") # convert to csv format
arcpy.AddMessage("(3b) zonal statistics: slope standard deviation at reef crest exported")

# 4. Rugosity

# a. Edge
rugosityEdge = tables + "\\rugosityEdge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, normalDSMx2,
    "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, rugosityEdge, None)
arcpy.AddMessage("(4a) rugosity at reef edge exported")

# b. Crest
rugosityCrest = tables + "\\rugosityCrest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, normalDSMx2,
    "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, rugosityCrest, None)
arcpy.AddMessage("(4b) rugosity at reef crest exported")

# 5. Fractal Dimension

desc = arcpy.Describe(normalDSM) # grab the information for the original raster
rasExtent = desc.extent # grab the extent object to set extent for resampling
arcpy.env.extent = rasExtent

# resample to 2x gsd
gsd2 = gsd*2

# resample to 4x gsd
gsd4 =  gsd*4
dsmx4 = scratchWS + "\\gsdx4.tif"
arcpy.Resample_management(normalDSM, dsmx4, gsd4, "BILINEAR")

# resample to 6x gsd
gsd6 =  gsd*6
dsmx6 = scratchWS + "\\gsdx6.tif"
arcpy.Resample_management(normalDSM, dsmx6, gsd6, "BILINEAR")

# resample to 8x gsd
gsd8 =  gsd*8
dsmx8 = scratchWS + "\\gsdx8.tif"
arcpy.Resample_management(normalDSM, dsmx8, gsd8, "BILINEAR")

# resample to 10x gsd
gsd10 =  gsd*10
dsmx10 = scratchWS + "\\gsdx10.tif"
arcpy.Resample_management(normalDSM, dsmx10, gsd10, "BILINEAR")

# resample to 12x gsd
gsd12 =  gsd*12
dsmx12 = scratchWS + "\\gsdx12.tif"
arcpy.Resample_management(normalDSM, dsmx12, gsd12, "BILINEAR")

arcpy.AddMessage("(5) rasters resampled for Fractal Dimension")

# a. Edge

# 2 x GSD
fractal2_edge = tables + "\\fractal2_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, normalDSMx2, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal2_edge, None)

# 4 x GSD
fractal4_edge = tables + "\\fractal4_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, dsmx4, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal4_edge, None)

# 6 x GSD
fractal6_edge = tables + "\\fractal6_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, dsmx6, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal6_edge, None)

# 8 x GSD
fractal8_edge = tables + "\\fractal8_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, dsmx8, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal8_edge, None)

# 10 x GSD
fractal10_edge = tables + "\\fractal10_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, dsmx10, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal10_edge, None)

# 12 x GSD
fractal12_edge = tables + "\\fractal12_edge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, dsmx12, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, fractal12_edge, None)

arcpy.AddMessage("(5a) Surface area for fractal D. extracted at edge plots")

# b. Crest

# 2 x GSD
fractal2_crest = tables + "\\fractal2_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, normalDSMx2, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal2_crest, None)

# 4 x GSD
fractal4_crest = tables + "\\fractal4_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, dsmx4, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal4_crest, None)

# 6 x GSD
fractal6_crest = tables + "\\fractal6_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, dsmx6, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal6_crest, None)

# 8 x GSD
fractal8_crest = tables + "\\fractal8_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, dsmx8, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal8_crest, None)

# 10 x GSD
fractal10_crest = tables + "\\fractal10_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, dsmx10, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal10_crest, None)

# 12 x GSD
fractal12_crest = tables + "\\fractal12_crest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, dsmx12, "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, fractal12_crest, None)

arcpy.AddMessage("(5b) Surface area for fractal D. extracted at crest plots")

# 6. Curvature

# generate curvature rasters
profileCurv = scratchWS + "\\profileCurv.tif"
planformCurv = scratchWS + "\\planformCurv.tif"
genCurv = scratchWS + "\\generalCurv.tif"
arcpy.ddd.Curvature(normalDSMx2, genCurv, 1, profileCurv, planformCurv)
arcpy.AddMessage("curvature rasters generated")

# a. Edge
# profile
profileCurvEdgeDBF = tables + "\\profileCurvEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', profileCurv,
    profileCurvEdgeDBF, "DATA", "MEAN")
profileCurvEdge = tables + "\\profileCurvEdge.csv"
arcpy.conversion.TableToTable(profileCurvEdgeDBF, tables, "profileCurvEdge.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: profile curvature at edge exported")

# planform
planformCurvEdgeDBF = tables + "\\planformCurvEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', planformCurv,
    planformCurvEdgeDBF, "DATA", "MEAN")
planformCurvEdge = tables + "\\planformCurvEdge.csv"
arcpy.conversion.TableToTable(planformCurvEdgeDBF, tables, "planformCurvEdge.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: planform curvature at edge exported")

# b. Crest
# profile
profileCurvCrestDBF = tables + "\\profileCurvCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', profileCurv,
    profileCurvCrestDBF, "DATA", "MEAN")
profileCurvCrest = tables + "\\profileCurvCrest.csv"
arcpy.conversion.TableToTable(profileCurvCrestDBF, tables, "profileCurvCrest.csv")
arcpy.AddMessage("zonal statistics: profile curvature at crest exported")

# planform
planformCurvCrestDBF = tables + "\\planformCurvCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', planformCurv,
    planformCurvCrestDBF, "DATA", "MEAN")
planformCurvCrest = tables + "\\planformCurvCrest.csv"
arcpy.conversion.TableToTable(planformCurvCrestDBF, tables, "planformCurvCrest.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: planform curvature at crest exported")

### End Metric Extraction ###

# Data Analysis

# 1. Elevation
# a. Edge
elevation_edge = pd.read_csv(elevationEdge)
elevation_edge = elevation_edge[["OBJECTID", "MEAN"]]
elevation_edge = elevation_edge.assign(gsd = gsd2, metric = "elevation", area = "edge")
elevation_edge.columns = ['value' if x == 'MEAN' else x for x in elevation_edge.columns]

# b. Crest
elevation_crest = pd.read_csv(elevationCrest)
elevation_crest = elevation_crest[["OBJECTID", "MEAN"]]
elevation_crest = elevation_crest.assign(gsd = gsd2, metric = "elevation", area = "crest")
elevation_crest.columns = ['value' if x == 'MEAN' else x for x in elevation_crest.columns]

# 2. Mean Slope

# a. Mean Slope, Edge
slope_edge = pd.read_csv(slopeEdge)
slope_edge = slope_edge[["OBJECTID", "MEAN"]]
slope_edge = slope_edge.assign(gsd = gsd2, metric = "slope", area = "edge")
slope_edge.columns =  ['value' if x == 'MEAN' else x for x in slope_edge.columns]

# b. Mean Slope, Crest
slope_crest = pd.read_csv(slopeCrest)
slope_crest = slope_crest[["OBJECTID", "MEAN"]]
slope_crest = slope_crest.assign(gsd = gsd2, metric = "slope", area = "crest")
slope_crest.columns =  ['value' if x == 'MEAN' else x for x in slope_crest.columns]

# 3. Slope Standard Deviation

# a. Standard Deviation of Slope, Edge
slopeStd_edge = pd.read_csv(slopeStdEdge)
slopeStd_edge = slopeStd_edge[["OBJECTID", "STD"]]
slopeStd_edge = slopeStd_edge.assign(gsd = gsd2, metric = "slopestd", area = "edge")
slopeStd_edge.columns =  ['value' if x == 'STD' else x for x in slopeStd_edge.columns]

# b. Standard Deviation of Slope, Crest
slopeStd_crest = pd.read_csv(slopeStdCrest)
slopeStd_crest = slopeStd_crest[["OBJECTID", "STD"]]
slopeStd_crest = slopeStd_crest.assign(gsd = gsd2, metric = "slopestd", area = "crest")
slopeStd_crest.columns =  ['value' if x == 'STD' else x for x in slopeStd_crest.columns]

# 4. Rugosity

# a. Edge
rugosity_edge = pd.read_csv(rugosityEdge)
rugosity_edge = rugosity_edge[["OBJECTID", "SArea"]]
rugosity_edge = rugosity_edge.assign(gsd = gsd2, metric = "rugosity", area = "edge")
rugosity_edge["SArea"] = rugosity_edge["SArea"]/0.25 # make rugosity ratio (3D SA (SArea) / 2D SA (0.25m))
rugosity_edge.columns = ['value' if x == 'SArea' else x for x in rugosity_edge.columns]

# b. Crest
rugosity_crest = pd.read_csv(rugosityCrest)
rugosity_crest = rugosity_crest[["OBJECTID", "SArea"]]
rugosity_crest = rugosity_crest.assign(gsd = gsd2, metric = "rugosity", area = "crest")
rugosity_crest["SArea"] = rugosity_crest["SArea"]/0.25 # make rugosity ratio (3D SA (SArea) / 2D SA (0.25m))
rugosity_crest.columns = ['value' if x == 'SArea' else x for x in rugosity_crest.columns]

# 5. Fractal Dimension
# a. edge

# 2 x GSD
fractald2_edge = pd.read_csv(fractal2_edge)
fractald2_edge = fractald2_edge[["OBJECTID", "SArea"]]
fractald2_edge = fractald2_edge.assign(gsd = gsd2, area = "edge")

# 4 x GSD
fractald4_edge = pd.read_csv(fractal4_edge)
fractald4_edge = fractald4_edge[["OBJECTID", "SArea"]]
fractald4_edge = fractald4_edge.assign(gsd = gsd4, metric = "fractaldimension", area = "edge")

# 6 x GSD
fractald6_edge = pd.read_csv(fractal6_edge)
fractald6_edge = fractald6_edge[["OBJECTID", "SArea"]]
fractald6_edge = fractald6_edge.assign(gsd = gsd6, metric = "fractaldimension", area = "edge")

# 8 x GSD
fractald8_edge = pd.read_csv(fractal8_edge)
fractald8_edge = fractald8_edge[["OBJECTID", "SArea"]]
fractald8_edge = fractald8_edge.assign(gsd = gsd8, metric = "fractaldimension", area = "edge")

# 10 x GSD
fractald10_edge = pd.read_csv(fractal10_edge)
fractald10_edge = fractald10_edge[["OBJECTID", "SArea"]]
fractald10_edge = fractald10_edge.assign(gsd = gsd10, metric = "fractaldimension", area = "edge")

# 12 x GSD
fractald12_edge = pd.read_csv(fractal12_edge)
fractald12_edge = fractald12_edge[["OBJECTID", "SArea"]]
fractald12_edge = fractald12_edge.assign(gsd = gsd12, metric = "fractaldimension", area = "edge")

# b. crest

# 2 x GSD
fractald2_crest = pd.read_csv(fractal2_crest)
fractald2_crest = fractald2_crest[["OBJECTID", "SArea"]]
fractald2_crest = fractald2_crest.assign(gsd = gsd2, metric = "fractaldimension", area = "crest")

# 4 x GSD
fractald4_crest = pd.read_csv(fractal4_crest)
fractald4_crest = fractald4_crest[["OBJECTID", "SArea"]]
fractald4_crest = fractald4_crest.assign(gsd = gsd4, metric = "fractaldimension", area = "crest")

# 6 x GSD
fractald6_crest = pd.read_csv(fractal6_crest)
fractald6_crest = fractald6_crest[["OBJECTID", "SArea"]]
fractald6_crest = fractald6_crest.assign(gsd = gsd6, metric = "fractaldimension", area = "crest")

# 8 x GSD
fractald8_crest = pd.read_csv(fractal8_crest)
fractald8_crest = fractald8_crest[["OBJECTID", "SArea"]]
fractald8_crest = fractald8_crest.assign(gsd = gsd8, metric = "fractaldimension", area = "crest")

# 10 x GSD
fractald10_crest = pd.read_csv(fractal10_crest)
fractald10_crest = fractald10_crest[["OBJECTID", "SArea"]]
fractald10_crest = fractald10_crest.assign(gsd = gsd10, metric = "fractaldimension", area = "crest")

# 12 x GSD
fractald12_crest = pd.read_csv(fractal12_crest)
fractald12_crest = fractald12_crest[["OBJECTID", "SArea"]]
fractald12_crest = fractald12_crest.assign(gsd = gsd12, metric = "fractaldimension", area = "crest")

list_fractals = [fractald2_edge, fractald4_edge, fractald6_edge, fractald8_edge, fractald10_edge,
    fractald12_edge, fractald2_crest, fractald4_crest, fractald6_crest, fractald8_crest,
    fractald10_crest, fractald12_crest] # dfs

df_fractal = pd.concat(list_fractals) # merge

# format master fractal table
  # CURRENTLY TABLE IS JUST OF SURFACE AREAS
df_fractal.columns = ['value' if x == 'SArea' else x for x in df_fractal.columns]
fractalWrite = tables + "\\fractalmerge.csv"
df_fractal.to_csv(fractalWrite)


# 6. Curvature

# a. Edge
# planform
planformcurvature_edge = pd.read_csv(planformCurvEdge)
planformcurvature_edge = planformcurvature_edge[["OBJECTID", "MEAN"]]
planformcurvature_edge = planformcurvature_edge.assign(gsd = gsd2, metric = "planformcurvature", area = "edge")
planformcurvature_edge.columns =  ['value' if x == 'MEAN' else x for x in planformcurvature_edge.columns]

# profile
profilecurvature_edge = pd.read_csv(profileCurvEdge)
profilecurvature_edge = profilecurvature_edge[["OBJECTID", "MEAN"]]
profilecurvature_edge = profilecurvature_edge.assign(gsd = gsd2, metric = "profilecurvature", area = "edge")
profilecurvature_edge.columns =  ['value' if x == 'MEAN' else x for x in profilecurvature_edge.columns]

# b. Crest
planformcurvature_crest = pd.read_csv(planformCurvCrest)
planformcurvature_crest = planformcurvature_crest[["OBJECTID", "MEAN"]]
planformcurvature_crest = planformcurvature_crest.assign(gsd = gsd2, metric = "planformcurvature", area = "crest")
planformcurvature_crest.columns =  ['value' if x == 'MEAN' else x for x in planformcurvature_crest.columns]

# profile
profilecurvature_crest = pd.read_csv(profileCurvCrest)
profilecurvature_crest = profilecurvature_crest[["OBJECTID", "MEAN"]]
profilecurvature_crest = profilecurvature_crest.assign(gsd = gsd2, metric = "profilecurvature", area = "crest")
profilecurvature_crest.columns =  ['value' if x == 'MEAN' else x for x in profilecurvature_crest.columns]

# 1. Merge tables
list_dfs = [elevation_edge, elevation_crest, slope_edge, slope_crest,
    slopeStd_edge, slopeStd_crest, rugosity_edge, rugosity_crest,
    profilecurvature_edge, profilecurvature_crest, planformcurvature_edge,
    planformcurvature_crest, df_fractal] # dfs

list_dfnames =  ["elevation_edge", "elevation_crest", "slope_edge", "slope_crest",
    "slopeStd_edge", "slopeStd_crest", "rugosity_edge", "rugosity_crest",
    "planformcurvature_edge", "planformcurvature_crest", "profilecurvature_crest",
    "profilecurvature_edge", "df_fractal"] # df names

df = pd.concat(list_dfs, keys = list_dfnames, sort = False) # merge

resultFolder = scratchWS + "\\results" # make a folder for the result dataframe
if not os.path.exists(resultFolder):
    os.makedirs(resultFolder)

outputDf = resultFolder + "\\metrics_joined.csv" # write out the result dataframe
df.to_csv(outputDf)
arcpy.AddMessage("data joined and exported to: " + outputDf)
