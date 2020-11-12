
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
scratchWS = arcpy.GetParameterAsText(0) #scratch folder
arcpy.env.scratchWorkspace = scratchWS #set the environmental variable for the scratch folder after knowing the scratch folder
#input files
input_dsm = arcpy.GetParameterAsText(1) #input oyster dsm
input_orthomosaic = arcpy.GetParameterAsText(2) # input oyster ortho
input_reefPolygon = arcpy.GetParameterAsText(3) # line along reef edge


# 0. Get Parent Arc Folder (assuming scratch nested within - change this to
# more adaptable to other users later)
ArcPar = os.path.dirname(scratchWS) # get arc parent folder
arcpy.AddMessage(ArcPar)

# 1. Buffer input reef polygon (outside)
gdb = ArcPar + "\\OysterMetrics.gdb" # get path of geodatabase

outBufferReefEdge = gdb + "\\outerBufferReef" # save within geodatabase
arcpy.Buffer_analysis(input_reefPolygon, outBufferReefEdge, "0.25 Meters")
arcpy.AddMessage("outer reef edge buffer created")

# 2. Convert Outer Buffer Polygon to line
outBufferReefLine = gdb + "\\outerBufferLine"
arcpy.PolygonToLine_management(outBufferReefEdge, outBufferReefLine) # create outer buffer line
arcpy.AddMessage("outer reef buffer converted to line")

# 3. Generate Samples along line
normalizeSamplePoints = gdb + "\\normalizeSamplePts"
arcpy.GeneratePointsAlongLines_management(outBufferReefLine, normalizeSamplePoints,
    'PERCENTAGE', Percentage = 10)
arcpy.AddMessage("normalizing samples created")

# 4. Extract DSM Values at Sample Points
normalizeSamplePts_DSMAttached = gdb + "\\normalizeSamplePts_DSMAttached"
ExtractValuesToPoints(normalizeSamplePoints, input_dsm, normalizeSamplePts_DSMAttached)
arcpy.AddMessage("normalizing samples extracted DSM values")

# 5. Create Plane from Sample Point DSM values)
normalPlane = scratchWS + "\\normalPlane.tif" # delete previous version
if os.path.exists(normalPlane):
    os.remove(normalPlane)

normalPlane = scratchWS + "\\normalPlane.tif"
outTrend = Trend(normalizeSamplePts_DSMAttached, "RASTERVALU", order = 1)
outTrend.save(normalPlane)
arcpy.AddMessage("normalized plane generated and saved:" + normalPlane)

# 6. Subtract normalized plane from DSM to create normalized dsm
normalDSM = scratchWS + "\\normalDSM.tif" # delete previous version
if os.path.exists(normalDSM):
    os.remove(normalDSM)

normalDSM = scratchWS + "\\normalDSM.tif"
outMinus = Raster(input_dsm) - Raster(normalPlane)
outMinus.save(normalDSM)
arcpy.AddMessage("DSM normalized to height and saved:" + normalDSM)


## Take Reef Samples from Normalized Reef DSM

# 7. Buffer input reef polygon (inside)
innerBufferReefEdge = gdb + "\\innerBufferReef"
arcpy.Buffer_analysis(input_reefPolygon, innerBufferReefEdge, "-0.25 Meters")
arcpy.AddMessage("inner reef edge buffer created")

# 8. Convert Inner Buffer Polygon to line
innerBufferReefLine = gdb + "\\innerBufferLine"
arcpy.PolygonToLine_management(innerBufferReefEdge, innerBufferReefLine) # create inner buffer line
arcpy.AddMessage("inner reef buffer converted to line")

# 9. Generate Samples along inner reef edge
reefEdgeSamplePoints = gdb + "\\reefEdgeSamplePoints"
arcpy.GeneratePointsAlongLines_management(innerBufferReefLine, reefEdgeSamplePoints,
    'PERCENTAGE', Percentage = 10)
arcpy.AddMessage("edge sample points created")

# 10. Buffer reef edge Samples
reefEdgeSamplesBuffered = gdb + "\\reefEdgeSamplesBuffered"
arcpy.Buffer_analysis(reefEdgeSamplePoints, reefEdgeSamplesBuffered, "0.25 Meters")
arcpy.AddMessage("edge sample points buffered")

# 11. Make Buffer Plots Square
reefEdgeSamplesSquareBuf = gdb + "\\reefEdgeSamplesSquareBuf"
arcpy.MinimumBoundingGeometry_management(reefEdgeSamplesBuffered, reefEdgeSamplesSquareBuf,
    "ENVELOPE")
arcpy.AddMessage("edge sample points converted to square buffer")

## Extract Crest points

# 12. Convert reef edge polygon to line
reefEdgeLine = gdb + "\\reefEdgeLine"
arcpy.PolygonToLine_management(input_reefPolygon, reefEdgeLine) # create reef edge line
arcpy.AddMessage("reef edge converted to line")

# 13. Generate transects along reef reefEdgeLine
crestTransects = gdb + "\\crestTransects"
arcpy.GenerateTransectsAlongLines_management(reefEdgeLine, crestTransects,
    '0.1 Meters', '100 Meters') # create transects at each 0.1m interval 100m in length (big enough to cross any oyster reef)
arcpy.AddMessage("reef crest transects created")

# 14. Clip transects to area of reef edge
crestTransectsClip = gdb + "\\crestTransectsClip"
arcpy.Clip_analysis(crestTransects, input_reefPolygon, crestTransectsClip) # clip transects to reef edge boundary
arcpy.AddMessage("reef crest transects clipped")

# 15. Copy features to attach shape length information for export
crestTransectsClipCF = gdb + "\\crestTransectsClipCF"
arcpy.CopyFeatures_management(crestTransectsClip, crestTransectsClipCF)
arcpy.AddMessage("transects copied")

# 16. Export Attribute Table to find longest line
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
