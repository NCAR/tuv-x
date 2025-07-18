#!/usr/bin/env python3

###############################################################
# tuvUtils.py
#
# Created:  October 2023
# Author:   Carl Drews
#
# This module contains TUV-x utility functions.
#
# (C) Copyright 2023 by Atmospheric Chemistry Observations & Modeling,
# University Corporation for Atmospheric Research (UCAR).
# 3450 Mitchell Lane
# Boulder, Colorado USA  80301
#
###############################################################

import sys

import urllib
import urllib.request
import http.client

import netCDF4
import numpy
import pandas



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Find the index corresponding to the requested height.
# Method is by proximity fuze - proceed until distance starts increasing
# altitude = array of heights above sea level
# specHeight = look for closest altitude to this height
# return array index closest to specHeight
def getHeightIndex(altitudes, specHeight):
   #progress("altitudes = {}".format(altitudes))

   distance =  abs(altitudes[0] - specHeight) + 1
   gotGreater = False		# found a greater value within the array
   for levelIndex in range(len(altitudes)):
      newDist = abs(altitudes[levelIndex] - specHeight)
      if (newDist >= distance):
         gotGreater = True
         break
      distance = newDist

   if (gotGreater):
      levelIndex -= 1
   return(levelIndex)



# Use linear interpolation to extract value from array.
# refArray and valArray are the same size.
# refValue = looking for value at this height
# refArray = vertical profile of height
# valArray = extract temperature or density (or height) from here
# return the value calculated between reference levels
def linearInterpolate(refValue, refArray, valArray):
   #progress("refValue = {}".format(refValue))
   #progress("refArray = {}".format(refArray))
   #progress("valArray = {}".format(valArray))

   # take care of the trivial edge cases
   if (refValue <= refArray[0]):
      return(valArray[0])
   if (refValue >= refArray[-1]):
      return(valArray[-1])

   # straddle the reference value
   rfi = 0
   while (refArray[rfi+1] <= refValue):
      rfi += 1

   # Calculate the reference fractional change
   # and apply fraction to the value delta.
   #progress("rfi = {}".format(rfi))
   #progress("straddle {} and {}".format(refArray[rfi], refArray[rfi+1]))
   refDelta = refArray[rfi+1] - refArray[rfi]
   valDelta = valArray[rfi+1] - valArray[rfi]
   fraction = (refValue - refArray[rfi]) / refDelta
   value = valArray[rfi] + fraction * valDelta
   #progress("value = {}".format(value))

   return(value)



# Retrieve one web page from the Internet.
# return page contents or None if cannot load
def downloadPage(pageURL):

   exception = None

   try:
      # handle unicode chars in the URL
      urlFile = urllib.request.urlopen(pageURL)
      #pageText = str(urlFile.read())
      pageText = urlFile.read()
      actualURL = urlFile.geturl()
      urlFile.close()

   except (urllib.error.URLError, http.client.InvalidURL,
      http.client.RemoteDisconnected, ConnectionResetError) as oops:
      progress("{}: {}".format(oops, pageURL))
      pageText = None
      actualURL = None
      exception = oops

   return(pageText)



# Retrieve some JSON TUV-x profile from configuration file.
# profileURL = CGI call to retrieve list of numerical values
# return list of floats
def getProfile(profileURL):
   allFloatStr = downloadPage(profileURL).decode("utf-8")
   #progress("allFloatStr = {}".format(allFloatStr))

   strList = allFloatStr.split('\n')[:-1]
   #progress("strList = {}".format(strList))
   values = list(map(float, strList))
   #progress("values = {}".format(values))

   return(values)



# Get values of temperature and density at specified height,
# using linear interpolation between layers if necessary.
# atHeight = measurement altitude
# return tuple of air (temperature, density)
def getStandardProfile(atHeight):

   # use the Javascript-intended web service for looking up standard profile
   heights = getProfile("https://www.acom.ucar.edu/cgi-bin/acom/tuvLookup.py?profile=height")
   #progress("heights = {}".format(heights))
   temperatures = getProfile("https://www.acom.ucar.edu/cgi-bin/acom/tuvLookup.py?profile=temperature")
   #progress("temperatures = {}".format(temperatures))
   densities = getProfile("https://www.acom.ucar.edu/cgi-bin/acom/tuvLookup.py?profile=density")
   #progress("densities = {}".format(densities))

   temperature = linearInterpolate(atHeight, heights, temperatures)
   density = linearInterpolate(atHeight, heights, densities)

   return(temperature, density)



# Modify aerosol optical depth by inserting layers at certain heights,
# dividing the existing layer instead of interpolating across it.
# heights = the vertical profile
# myAOD = from the JSON list of grids
# heightsToInsert = list of new levels
# return the modified grid
def insertOpticalDepths(heights, myAOD, heightsToInsert):
   #progress("heights = {} {} length {}".format(heights, type(heights), len(heights)))
   #progress("myAOD = {} {} length {}".format(myAOD, type(myAOD), len(myAOD)))
   #progress("insert = {} length {}".format(heightsToInsert, len(heightsToInsert)))

   # extend the AOD list to match size of height profile - 1
   myAOD.extend(numpy.zeros(len(heights) - len(myAOD) - 1))
   #progress("Extended myAOD = {} {} length {}".format(myAOD, type(myAOD), len(myAOD)))

   for heightToInsert in heightsToInsert:

      needInsert = True
      for atIndex in range(len(heights) + 1):
         if (atIndex >= len(heights)):
            # insertion is after end of height
            break

         # check for exact match with existing level
         if (heightToInsert == heights[atIndex]):
            needInsert = False
            break

         # check for height inside current layer
         if (heightToInsert < heights[atIndex]):
            break
 
      if (not needInsert):
         continue

      #progress("heightToInsert = {}   atIndex = {}".format(heightToInsert, atIndex))

      if (atIndex >= len(heights)):
         # insertion is after end of height
         heights.insert(atIndex, heightToInsert)
         myAOD.insert(atIndex, myAOD[atIndex - 1])
         continue

      # divide the AOD proportionally between the lower layer and the new layer
      span = heights[atIndex] - heights[atIndex - 1]
      lowerSpan = heightToInsert - heights[atIndex - 1]
      upperSpan = heights[atIndex] - heightToInsert
      #progress("spans = {} + {} = {}".format(lowerSpan, upperSpan, span))

      # add the new level
      heights.insert(atIndex, heightToInsert)

      # set lower layer and new layer by allocation, not interpolation
      valueToSplit = myAOD[atIndex - 1]
      myAOD[atIndex - 1] = valueToSplit * lowerSpan / span
      #progress("myAOD02 = {} {}".format(myAOD, type(myAOD)))
      myAOD.insert(atIndex, valueToSplit * upperSpan / span)
      #progress("myAOD03 = {}".format(myAOD))

   return(heights)



# Load vertical grid from JSON configuration grid at equal intervals.
# myGrid = from the JSON list of grids
# return the created grid as list of altitudes
def intervalToValues(myGrid):
   # verify that the type is what we know how to modify
   if (myGrid["type"].lower() != "equal interval"):
      progress("Warning: Unknown grid type {}<br>".format(myGrid["type"]))
      return(myGrid)

   # collect interval parameters
   begins = myGrid["begins at"]
   ends = myGrid["ends at"]
   delta = myGrid["cell delta"]

   # create the numeric grid
   heights = numpy.arange(begins, ends + delta, delta)

   myGrid["values"] = heights.tolist()
   #progress("Created heights are now {}".format(myGrid["values"]))

   myGrid["type"] = "from config file"
   del myGrid["begins at"]
   del myGrid["ends at"]
   del myGrid["cell delta"]

   return(myGrid)



# Set ground elevation of grid by compressing it upwards.
# Raise the first grid element from zero to groundHeight,
# and preserve the highest level in the last grid element.
# Everything in between is scaled linear upwards.
# heights = array of increasing altitues
# groundHeight = make this value the lowest level in the grid
#       Example: Leadville, Colorado is at elevation 3.66 km.
# return the modified heights
def raiseGround(heights, groundHeight):
   if (groundHeight is None):
      return(heights)
   if (groundHeight == heights[0]):
      return(heights)

   maxHeight = heights[-1]
   squeezeFactor = ((maxHeight - groundHeight)
      / (maxHeight - heights[0]))
   #progress("Ground elevation interval squeezeFactor = {}".format(squeezeFactor))

   for hi in range(len(heights)):
      diff = maxHeight - heights[hi]
      heights[hi] = maxHeight - diff * squeezeFactor

   return(heights)


# common JSON tags
filePathAttr = "file path"
title01 = "Geometric Altitude (km)"
title02 = "Number density n (cm-3)"           # could also be temperature

# Load data profile from a CSV file.
# myProfile = from the JSON list of profiles
# tuvBuildDir = the data file is read from here
# return the loaded DataFrame of Altitude and (Density or Temperature)
def fileToValues(myProfile, tuvBuildDir):

   # load contents of existing CSV file
   csvName = myProfile[filePathAttr]
   #progress("csvName = {}".format(csvName))
   csvPath = "{}{}".format(tuvBuildDir, csvName)
   progress("csvPath = {}".format(csvPath))

   csvContent = pandas.read_csv(csvPath, skiprows=3,
      names=[title01, title02],
      sep=" ", skipinitialspace=True,
      dtype={title01: float})

   #progress("Loaded csvContent {}".format(csvContent))
   return(csvContent)



# Modify profile by interpolating onto height grid.
# myProfile = already loaded, JSON changed to enumerated values
# myDataFrame = altitude and data that we loaded earlier
# heights = list of altitudes in raised and inserted grid
# customHeight = probably the measurement height
# customValue = set this air density or temperature
# return the modified profile
def interpolateProfile(myProfile, myDataFrame, heights, customHeight, customValue=None):
   #progress("interpolateProfile() myProfile = {}".format(myProfile))
   #progress("interpolateProfile() heights = {}".format(heights))

   # extract profile columns as arrays
   proHeights = numpy.array(myDataFrame[title01])
   proValues = numpy.array(myDataFrame[title02])
   #progress("interpolateProfile() proHeights = {}".format(proHeights))
   #progress("interpolateProfile() proValues = {}".format(proValues))

   # calculate the profile values at the supplied heights
   newValues = numpy.zeros(len(heights))
   for vi in range(len(newValues)):
      newValues[vi] = linearInterpolate(heights[vi],
         proHeights, proValues)
   #progress("interpolateProfile() newValues = {}".format(newValues)) 

   # set custom value at the custom height
   if (customValue is not None):
      heightIndex = getHeightIndex(heights, customHeight)
      newValues[heightIndex] = customValue

   # change profile JSON to enumerated values
   # remove the JSON file path
   del myProfile[filePathAttr]
   myProfile["values"] = newValues.tolist()
   myProfile["type"] = "from config file"

   # if height grid is not there, add it???
   if (not "grid" in myProfile):
      #progress("Adding grid to profile {}".format(myProfile["name"]))
      myProfile["grid"] = {"name": "height", "units": "km"}

   return(myProfile)



# Create NetCDF file containing derived variables.
# filename = with .nc extension
# dimensions = dictionary of tuples [name, size] for all variables
# nameList = list of variable names to create
# dimList = dimensions for each item in varList
# valueList = values for each item in varList
# unitList = units for each variable in varList
def writeFile(filename, dimensions, nameList,
   dimList, valueList, unitList):
   file = netCDF4.Dataset(filename, mode='w')

   for dim in dimensions:
      file.createDimension(dim[0], size=dim[1])

   # loop through and create variables, assign values
   for ni, di, vi, ui in zip(nameList, dimList, valueList, unitList):
      #progress("{}: {}".format(ni, di))
      newVar = file.createVariable(ni, 'f', di)
      newVar[:] = vi[:]
      newVar.units = ui

   file.close()
   return



# Write NetCDF file containing clouds radiator.
# tuvParams = TUV numerical parameters in a tuvCore object
# numLevels = size over dimension for vertical levels
# numWavelengths = size of dimension for wavelength
# heightGrid = the vertical levels
# toFilename = path and filename to create
def writeClouds(tuvParams, numLevels, numWavelengths,
   heightGrid, toFilename):

   # set up the dimensions
   cloudDims = [
      ["vertical_level", numLevels],
      ["wavelength", numWavelengths]
   ]

   # set up the variables
   varNames = ["optical_depth", "single_scattering_albedo", "asymmetry_factor"]
   varDims = [
      (cloudDims[1][0], cloudDims[0][0]),
      (cloudDims[1][0], cloudDims[0][0]),
      (cloudDims[1][0], cloudDims[0][0])
   ]

   # create zero arrays the proper size
   od = numpy.zeros((numWavelengths, numLevels))
   ssa = numpy.zeros((numWavelengths, numLevels))
   af = numpy.zeros((numWavelengths, numLevels))
   varValues = [ od, ssa, af ]

   # assign numerical values strictly *between* base and top
   baseIndex = getHeightIndex(heightGrid, tuvParams.cloudBase)
   topIndex = getHeightIndex(heightGrid, tuvParams.cloudTop)
   #progress("baseIndex = {}   topIndex = {}".format(baseIndex, topIndex))

   #progress("Fill in clouds from height index {} up to {}".format(baseIndex, topIndex))
   ssa[:, baseIndex:topIndex] = tuvParams.cloudSingScatAlbedo
   af[:, baseIndex:topIndex] = tuvParams.cloudAsymmetry

   totalThickness = heightGrid[topIndex] - heightGrid[baseIndex]
   #progress("Allocate total cloud AOD over cloud thickness = {}".format(totalThickness))
   for cloudIndex in range(baseIndex, topIndex):
      od[:, cloudIndex] = (tuvParams.cloudOptDepth
         * (heightGrid[cloudIndex + 1] - heightGrid[cloudIndex]) / totalThickness)
   #progress("Cloud AOD[7] = {}".format(od[7]))

   varUnits = ["", "", ""]

   # create the NetCDF file
   writeFile(toFilename, cloudDims, varNames,
      varDims, varValues, varUnits)

   return



# Calculated the normalized column per layer
# from the standard profile of optical depths.
# stdOptDepths = standard optical depths, probably loaded from standard config.json
# levelHeights = the altitudes; should be 1 larger than stdOptDepths
def normColumnLayer(stdOptDepths, levelHeights):

   if (False):
      progress("in normColumnLayer() stdOptDepths = {} {}"
         .format(stdOptDepths, len(stdOptDepths)))
      progress("in normColumnLayer() levelHeights = {} {}"
         .format(levelHeights, len(levelHeights)))

   # begin with the TUV level density (AOD/km)
   numDepths = len(stdOptDepths)
   normalized = numpy.zeros(numDepths)
   normalized[:] = stdOptDepths

   #progress("numDepths = {}".format(numDepths))

   # calculate the TUV layer AOT 340 nm
   layerAOT = numpy.zeros(numDepths)
   for li in range(numDepths - 1):
      layerAOT[li] = (((normalized[li] + normalized[li+1]) / 2)
         * (levelHeights[li+1] - levelHeights[li]))
   #progress("layerAOT = {} {}".format(layerAOT, len(layerAOT)))

   # normalize so that sum = 1.0
   sumLayerAOT = numpy.sum(layerAOT)
   #progress("sum layerAOT = {}".format(sumLayerAOT))
   normalized = layerAOT / sumLayerAOT

   return(normalized)



# Write NetCDF file containing aerosols radiator.
# tuvParams = TUV numerical parameters in a tuvCore object
# numLevels = size over dimension for vertical levels
# numWavelengths = size of dimension for wavelength
# heightGrid = the altitudes
# wavelengthGrid = the wavelengths
# toFilename = path and filename to create
def writeAerosols(tuvParams, numLevels, numWavelengths,
   heightGrid, wavelengthGrid, toFilename):

   # set up the dimensions
   cloudDims = [
      ["vertical_level", numLevels],
      ["wavelength", numWavelengths]
   ]

   # set up the variables
   varNames = ["optical_depth", "single_scattering_albedo", "asymmetry_factor"]
   varDims = [
      (cloudDims[1][0], cloudDims[0][0]),
      (cloudDims[1][0], cloudDims[0][0]),
      (cloudDims[1][0], cloudDims[0][0])
   ]

   # create zero arrays the proper size
   od = numpy.zeros((numWavelengths, numLevels))
   ssa = numpy.zeros((numWavelengths, numLevels))
   af = numpy.zeros((numWavelengths, numLevels))
   varValues = [ od, ssa, af ]

   # assign constant numerical values
   ssa[:] = tuvParams.aerosolSingScatAlbedo
   af[:] = tuvParams.aerosolAsymmetry

   # retrieve and modify the standard profile of optical depths
   normalizedColumn = normColumnLayer(tuvParams.stdOpticalDepths,
      heightGrid)
   if (False):
      progress("normalizedColumn = {} len {} sum {}"
         .format(normalizedColumn, len(normalizedColumn),
         numpy.sum(normalizedColumn)))
   normalizedColumn *= tuvParams.aerosolOptDepth

   # calculate the aerosol optical depth by wavelength
   od[:] = tuvParams.aerosolOptDepth
   for wi in range(len(wavelengthGrid)):
      wavelength = wavelengthGrid[wi]

      for li in range(numLevels):
         # AOD(lambda) formula by Kirk Ullmann - November 2023
         od[wi, li] = (normalizedColumn[li]
            * (tuvParams.aerosolLambda0 / wavelength) ** tuvParams.aerosolAlpha)
         #od[wi, 100] = wi / 100.0		# bogus, for ncview plotting
         #progress("{} = {} * {} / {}".format(od[wi,li], normalizedColumn[li],
         #   tuvParams.aerosolLambda0, wavelength))

   varUnits = ["", "", ""]

   # create the NetCDF file
   writeFile(toFilename, cloudDims, varNames,
      varDims, varValues, varUnits)

   return

