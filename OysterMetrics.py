
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
gsd = float(gsdx[0])*2 # target GSD; 2x original GSD
arcpy.Resample_management(normalDSM, normalDSMx2, gsd, "BILINEAR") # use gsd as an input for cell size
arcpy.AddMessage("(2) normalized DSM resampled to 2xGSD: " + str(gsd))

### DSM DOWNSAMPLED ###


## GENERATE SAMPLES
# a. Crest

# 12a. Buffer further inside the reef edge for crest (to avoid overlap with reef edge)
BufferReefCrest = gdb + "\\crestBufferReef"
arcpy.Buffer_analysis(input_reefPolygon, BufferReefCrest, "-0.5 Meters")
arcpy.AddMessage("crest reef buffer created")

# 12b. Convert Crest Buffer Polygon to line
BufferReefCrestLine = gdb + "\\crestBufferLine"
arcpy.PolygonToLine_management(BufferReefCrest, BufferReefCrestLine) # create crest buffer line
arcpy.AddMessage("crest reef buffer converted to line")

# 12c. Generate transects across reef edge inner buffer (Step 8)
crestTransects = gdb + "\\crestTransects"
arcpy.GenerateTransectsAlongLines_management(BufferReefCrestLine, crestTransects,
    '0.1 Meters', '100 Meters') # create transects at each 0.1m interval 100m in length (big enough to cross any oyster reef)
arcpy.AddMessage("reef crest transects created")

# 13. Clip transects to area of reef edge
crestTransectsClip = gdb + "\\crestTransectsClip"
arcpy.Clip_analysis(crestTransects, BufferReefCrest, crestTransectsClip) # clip transects to reef edge boundary
arcpy.AddMessage("reef crest transects clipped")

# 14. Copy features to attach shape length information for export
crestTransectsClipCF = gdb + "\\crestTransectsClipCF"
arcpy.CopyFeatures_management(crestTransectsClip, crestTransectsClipCF)
arcpy.AddMessage("transects copied")

# 15. Export Attribute Table to find longest line
transectLines = scratchWS + "\\transectLines.csv"
arcpy.CopyRows_management(crestTransectsClipCF, transectLines)
arcpy.AddMessage("reef crest transect table exported")

# 16. Select OBJECT ID of longest line
transectLines_df  = pd.read_csv(transectLines) # read transect length csv
maxlength = pd.to_numeric(transectLines_df['Shape_Length']).max() # grab max length
arcpy.AddMessage("the longest line is " + maxlength.astype(str)) # print longest length

ismax = pd.to_numeric(transectLines_df['Shape_Length']) == maxlength  # make True/False column for longest line
maxID = (transectLines_df[ismax]).iloc[0]['OBJECTID'] # filter for True rows, select OID of the first row
arcpy.AddMessage("the longest line has Object ID:" + maxID.astype(str))

# 17. Select crest transect by OID and copy to new feature class
maxCrestTransectFile = gdb + "\\maxCrestTransect"
query = 'OBJECTID =' + maxID.astype(str)  # build the query
arcpy.AddMessage("the query is: " + query) # print the query
arcpy.management.MakeFeatureLayer(crestTransectsClipCF, 'maxCrestTransect', query, None) # make layer with longest transect
arcpy.CopyFeatures_management('maxCrestTransect', maxCrestTransectFile) # save that layer to a shapefile

# 18. Generate Samples along reef crest
reefCrestSamplePoints = gdb + "\\reefCrestSamplePoints"
arcpy.GeneratePointsAlongLines_management(maxCrestTransectFile, reefCrestSamplePoints,
    'PERCENTAGE', Percentage = 10)
arcpy.AddMessage("crest sample points created")

# 19. Buffer reef crest Samples
reefCrestSamplesBuffered = gdb + "\\reefCrestSamplesBuffered"
arcpy.Buffer_analysis(reefCrestSamplePoints, reefCrestSamplesBuffered, "0.25 Meters")
arcpy.AddMessage("crest sample points buffered")

# 20. Make Buffer Plots Square
reefCrestSamplesSquareBuf = gdb + "\\reefCrestSamplesSquareBuf"
arcpy.MinimumBoundingGeometry_management(reefCrestSamplesBuffered, reefCrestSamplesSquareBuf,
    "ENVELOPE")
arcpy.AddMessage("crest sample points converted to square buffer")





## Edge Samples
# a. Buffer input reef polygon (inside)
innerBufferReefEdge = gdb + "\\innerBufferReef"
arcpy.Buffer_analysis(input_reefPolygon, innerBufferReefEdge, "-0.35 Meters") # buffer in by diagonal length
arcpy.AddMessage("(3a) inner reef edge buffer created")

# b. Convert Inner Buffer Polygon to line
innerBufferReefLine = gdb + "\\innerBufferLine"
arcpy.PolygonToLine_management(innerBufferReefEdge, innerBufferReefLine) # create inner buffer line
arcpy.AddMessage("(3b) inner reef buffer converted to line")

# c. Generate Samples along inner reef edge
reefEdgeSamplePoints = gdb + "\\reefEdgeSamplePoints"

if (int(input_sampleNum) == 4):
    arcpy.AddMessage("Default, 4 Samples chosen")
    arcpy.management.GeneratePointsAlongLines(innerBufferReefLine, reefEdgeSamplePoints,
        'PERCENTAGE', None, 25, None) # take four edge samples
else :
    arcpy.AddMessage("The default was not chosen...")
    arcpy.management.GeneratePointsAlongLines(innerBufferReefLine, reefEdgeSamplePoints,
        'PERCENTAGE', None, 25, None) # take the chosen # of edge samples

#arcpy.management.GeneratePointsAlongLines(innerBufferReefLine, reefEdgeSamplePoints,
    #'DISTANCE', "0.4 Meters") # separate each by the diagonal length (0.35) + buffer for separation
#arcpy.AddMessage("(3c) edge sample points created")

# 10. Buffer reef edge Samples
reefEdgeSamplesBuffered = gdb + "\\reefEdgeSamplesBuffered"
arcpy.Buffer_analysis(reefEdgeSamplePoints, reefEdgeSamplesBuffered, "0.25 Meters")
arcpy.AddMessage("edge sample points buffered")

# 11. Make Buffer Plots Square
reefEdgeSamplesSquareBuf = gdb + "\\reefEdgeSamplesSquareBuf"
arcpy.MinimumBoundingGeometry_management(reefEdgeSamplesBuffered, reefEdgeSamplesSquareBuf,
    "ENVELOPE")
arcpy.AddMessage("edge sample points converted to square buffer")






## METRICS
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
arcpy.AddMessage("zonal statistics: elevation at reef edge exported")

# b. Crest
elevationCrestDBF = tables + "\\elevationCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, "OBJECTID", normalDSMx2, elevationCrestDBF,
    "DATA", "MEAN")

elevationCrest = tables + "\\elevationCrest.csv"
arcpy.conversion.TableToTable(elevationCrestDBF, tables, "elevationCrest.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: elevation at reef crest exported")

# 2. Mean Slope

# Create Slope Raster
slopeDSM = tables + "//slopeDSM.tif"
outSlopeDSM = arcpy.sa.Slope(normalDSMx2, "DEGREE", 1, "PLANAR", "METER")
outSlopeDSM.save(slopeDSM) # save the slope raster
arcpy.AddMessage("slope raster generated")

# a. Edge
slopeEdgeDBF = tables + "\\slopeEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeEdgeDBF,
    "DATA", "MEAN")

slopeEdge = tables + "\\slopeEdge.csv"
arcpy.conversion.TableToTable(slopeEdgeDBF, tables, "slopeEdge.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: slope at reef edge exported")

# b. Slope
slopeCrestDBF = tables + "\\slopeCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeCrestDBF,
    "DATA", "MEAN")

slopeCrest = tables + "\\slopeCrest.csv"
arcpy.conversion.TableToTable(slopeCrestDBF, tables, "slopeCrest.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: slope at reef crest exported")

# 3. Standard Deviation of Slope
# a. Edge
slopeStdEdgeDBF = tables + "\\slopeStdEdge.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefEdgeSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeStdEdgeDBF,
    "DATA", "STD")
slopeStdEdge = tables + "\\slopeStdEdge.csv"
arcpy.conversion.TableToTable(slopeStdEdgeDBF, tables, "slopeStdEdge.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: slope standard deviation at reef edge exported")

# b. Crest
slopeStdCrestDBF = tables + "\\slopeStdCrest.dbf"
arcpy.sa.ZonalStatisticsAsTable(reefCrestSamplesSquareBuf, 'OBJECTID', slopeDSM, slopeStdCrestDBF,
    "DATA", "STD")
slopeStdCrest = tables + "\\slopeStdCrest.csv"
arcpy.conversion.TableToTable(slopeStdCrestDBF, tables, "slopeStdCrest.csv") # convert to csv format
arcpy.AddMessage("zonal statistics: slope standard deviation at reef crest exported")

# 24. Rugosity

# a. Edge
rugosityEdge = tables + "\\rugosityEdge.csv"
arcpy.ddd.AddSurfaceInformation(reefEdgeSamplesSquareBuf, normalDSMx2,
    "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefEdgeSamplesSquareBuf, rugosityEdge, None)
arcpy.AddMessage("rugosity at reef edge exported")

# b. Crest
rugosityCrest = tables + "\\rugosityCrest.csv"
arcpy.ddd.AddSurfaceInformation(reefCrestSamplesSquareBuf, normalDSMx2,
    "SURFACE_AREA", "BILINEAR", None, 1, 0, None)
arcpy.management.CopyRows(reefCrestSamplesSquareBuf, rugosityCrest, None)
arcpy.AddMessage("rugosity at reef crest exported")

# 4. Curvature

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
elevation_edge = elevation_edge.assign(metric = "elevation", area = "edge")
elevation_edge.columns = ['value' if x == 'MEAN' else x for x in elevation_edge.columns]
arcpy.AddMessage(elevation_edge)

# b. Crest
elevation_crest = pd.read_csv(elevationCrest)
elevation_crest = elevation_crest[["OBJECTID", "MEAN"]]
elevation_crest = elevation_crest.assign(metric = "elevation", area = "crest")
elevation_crest.columns = ['value' if x == 'MEAN' else x for x in elevation_crest.columns]
arcpy.AddMessage(elevation_crest)

# 2. Mean Slope

# a. Mean Slope, Edge
slope_edge = pd.read_csv(slopeEdge)
slope_edge = slope_edge[["OBJECTID", "MEAN"]]
slope_edge = slope_edge.assign(metric = "slope", area = "edge")
slope_edge.columns =  ['value' if x == 'MEAN' else x for x in slope_edge.columns]

# b. Mean Slope, Crest
slope_crest = pd.read_csv(slopeCrest)
slope_crest = slope_crest[["OBJECTID", "MEAN"]]
slope_crest = slope_crest.assign(metric = "slope", area = "crest")
slope_crest.columns =  ['value' if x == 'MEAN' else x for x in slope_crest.columns]

# 3. Slope Standard Deviation

# a. Standard Deviation of Slope, Edge
slopeStd_edge = pd.read_csv(slopeStdEdge)
slopeStd_edge = slopeStd_edge[["OBJECTID", "STD"]]
slopeStd_edge = slopeStd_edge.assign(metric = "slopestd", area = "edge")
slopeStd_edge.columns =  ['value' if x == 'STD' else x for x in slopeStd_edge.columns]

# b. Standard Deviation of Slope, Crest
slopeStd_crest = pd.read_csv(slopeStdCrest)
slopeStd_crest = slopeStd_crest[["OBJECTID", "STD"]]
slopeStd_crest = slopeStd_crest.assign(metric = "slopestd", area = "crest")
slopeStd_crest.columns =  ['value' if x == 'STD' else x for x in slopeStd_crest.columns]

# 4. Rugosity

# a. Edge
rugosity_edge = pd.read_csv(rugosityEdge)
rugosity_edge = rugosity_edge[["OBJECTID", "SArea"]]
rugosity_edge = rugosity_edge.assign(metric = "rugosity", area = "edge")
rugosity_edge["SArea"] = rugosity_edge["SArea"]/0.25 # make rugosity ratio (3D SA (SArea) / 2D SA (0.25m))
rugosity_edge.columns = ['value' if x == 'SArea' else x for x in rugosity_edge.columns]

# b. Crest
rugosity_crest = pd.read_csv(rugosityCrest)
rugosity_crest = rugosity_crest[["OBJECTID", "SArea"]]
rugosity_crest = rugosity_crest.assign(metric = "rugosity", area = "crest")
rugosity_crest["SArea"] = rugosity_crest["SArea"]/0.25 # make rugosity ratio (3D SA (SArea) / 2D SA (0.25m))
rugosity_crest.columns = ['value' if x == 'SArea' else x for x in rugosity_crest.columns]

# 5. Curvature

# a. Edge
# planform
planformcurvature_edge = pd.read_csv(planformCurvEdge)
planformcurvature_edge = planformcurvature_edge[["OBJECTID", "MEAN"]]
planformcurvature_edge = planformcurvature_edge.assign(metric = "planformcurvature", area = "edge")
planformcurvature_edge.columns =  ['value' if x == 'MEAN' else x for x in planformcurvature_edge.columns]

# profile
profilecurvature_edge = pd.read_csv(profileCurvEdge)
profilecurvature_edge = profilecurvature_edge[["OBJECTID", "MEAN"]]
profilecurvature_edge = profilecurvature_edge.assign(metric = "profilecurvature", area = "edge")
profilecurvature_edge.columns =  ['value' if x == 'MEAN' else x for x in profilecurvature_edge.columns]

# b. Crest
planformcurvature_crest = pd.read_csv(planformCurvCrest)
planformcurvature_crest = planformcurvature_crest[["OBJECTID", "MEAN"]]
planformcurvature_crest = planformcurvature_crest.assign(metric = "planformcurvature", area = "crest")
planformcurvature_crest.columns =  ['value' if x == 'MEAN' else x for x in planformcurvature_crest.columns]

# profile
profilecurvature_crest = pd.read_csv(profileCurvCrest)
profilecurvature_crest = profilecurvature_crest[["OBJECTID", "MEAN"]]
profilecurvature_crest = profilecurvature_crest.assign(metric = "profilecurvature", area = "crest")
profilecurvature_crest.columns =  ['value' if x == 'MEAN' else x for x in profilecurvature_crest.columns]

# 1. Merge tables
list_dfs = [elevation_edge, elevation_crest, slope_edge, slope_crest,
    slopeStd_edge, slopeStd_crest, rugosity_edge, rugosity_crest,
    profilecurvature_edge, profilecurvature_crest, planformcurvature_edge,
    planformcurvature_crest] # dfs

list_dfnames =  ["elevation_edge", "elevation_crest", "slope_edge", "slope_crest",
    "slopeStd_edge", "slopeStd_crest", "rugosity_edge", "rugosity_crest",
    "planformcurvature_edge", "planformcurvature_crest", "profilecurvature_crest",
    "profilecurvature_edge"] # df names

df = pd.concat(list_dfs, keys = list_dfnames) # merge

resultFolder = scratchWS + "\\results" # make a folder for the result dataframe
if not os.path.exists(resultFolder):
    os.makedirs(resultFolder)

outputDf = resultFolder + "\\metrics_joined.csv" # write out the result dataframe
df.to_csv(outputDf)
arcpy.AddMessage("data joined and exported to: " + outputDf)
