#!/usr/bin/env python3

###############################################################
# tuv-x-remote.py
#
# Created:  February 2023
# Author:   Carl Drews
#
# This script carries out a TUV-x calculation.
#
# Usage:
#   This script is called from the POST method of an input form.
#
# (C) Copyright 2023 by Atmospheric Chemistry Observations & Modeling,
# University Corporation for Atmospheric Research (UCAR).
# 3450 Mitchell Lane
# Boulder, Colorado USA  80301
#
###############################################################

import sys
import datetime
import tuvCore
import utilsLite

import json
import os
import subprocess
import contextlib

import netCDF4
import numpy
import pandas

import tuvUtils



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Replace all key values in a dictionary at all levels.
# Recursively traverse all levels of the dictionary.
# myObject = Python object, probably starting with JSON dictionary
#	Could also be a list.
# myKey = look for this key
# myValue = replace value of the key with this
def replaceAllKeys(myObject, myKey, myValue):

   if (isinstance(myObject, dict)):
      #progress("Dictionary")
      for key, item in myObject.items():
         if (key == myKey):
            #progress("dict key = {}: {}".format(key, item))
            myObject[key] = myValue	# replace the key's value
         replaceAllKeys(item, myKey, myValue)

   elif (isinstance(myObject, list)):
      #progress("List")
      for item in myObject:
         #progress("list item = {}".format(item))
         replaceAllKeys(item, myKey, myValue)

   return



# Replace all specific keys in the specified dictionary
# myObject = Python object, probably starting with JSON dictionary
#	Could also be a list.
# dictName = search for dictionaries with this name
# myKey = search for this key
# myValue = key values becomes this value (replaced or new)
def setAllDictKeys(myObject, dictName, myKey, myValue):

   if (isinstance(myObject, dict)):
      #progress("Dictionary")

      for key, item in myObject.items():
         #progress("\tdict key = {}".format(key))
         if (key == dictName):
            #progress("\t\tsetting key {} to {}".format(myKey, myValue))
            item[myKey] = myValue	# add or replace the key's value
         setAllDictKeys(item, dictName, myKey, myValue)

   elif (isinstance(myObject, list)):
      #progress("List")
      for item in myObject:
         setAllDictKeys(item, dictName, myKey, myValue)

   return



CLOUDS_FILENAME = "clouds.nc"
AEROSOLS_FILENAME = "aerosols.nc"

# Substitute new values into an existing JSON structure.
# jStruct = some values will get replaced
# tuvValues = replacements will come from this structure
# tuvBuildDir = TUVX build directory containing config files below
def replaceJson(jStruct, tuvValues, tuvBuildDir):
   # turn all diagnostics off
   replaceAllKeys(jStruct, "enable diagnostics", False)

   # collect the new vertical levels to insert
   heightsToInsert = [tuvValues.height]
   if (tuvValues.cloudOptDepth > 0.0):
      heightsToInsert.extend([tuvValues.cloudBase, tuvValues.cloudTop])

   # make JSON substitions to profiles
   densityProfile = densityDataFrame = None
   temperatureProfile = temperatureDataFrame = None
   for profile in jStruct["profiles"]:
      #progress("profile = {}".format(profile))
      if (profile["name"].lower() == "solar zenith angle"):
         #progress("profile name {}".format(profile["name"]))
         if (tuvValues.inputOption == 1):
            profile["latitude"] = tuvValues.latitude
            profile["longitude"] = tuvValues.longitude

         if (tuvValues.inputOption == 2):
            profile["type"] = "from config file"
            profile["values"] = [tuvValues.zenithAngle]
            profile["grid"] = {
               "name": "time",
               "units": "hours"
            }

            del profile["latitude"]
            del profile["longitude"]
            del profile["year"]
            del profile["month"]
            del profile["day"]

      if ((profile["name"].lower() == "solar zenith angle"
            and tuvValues.inputOption == 1)
         or profile["name"].lower() == "earth-sun distance"):
         #progress("profile name {}".format(profile["name"]))
         myDateStr = tuvValues.date
         year = myDateStr[0:4]
         month = myDateStr[4:6]
         day = myDateStr[6:8]
         profile["year"] = year
         profile["month"] = month
         profile["day"] = day

      if (profile["name"].lower() == "air"):
         densityProfile = profile
         densityDataFrame = tuvUtils.fileToValues(profile, tuvBuildDir)
         #progress("Loaded air density profile {}.<br>".format(densityProfile))

      if (profile["name"].lower() == "temperature"):
         temperatureProfile = profile
         temperatureDataFrame = tuvUtils.fileToValues(profile, tuvBuildDir)
         #progress("Loaded temperature profile {}.<br>".format(temperatureProfile))

      if (profile["name"].lower() == "o3"):
         profile["reference column"] = tuvValues.ozone
      if (profile["name"].lower() == "surface albedo"):
         profile["uniform value"] = tuvValues.albedo

   # make JSON substitions to grids
   heightGrid = None
   for grid in jStruct["grids"]:
      #progress("grid = {}".format(grid))
      if (grid["name"].lower() == "time"):
         myTime = tuvValues.time.split(":")
         hour = int(myTime[0])
         minute = int(myTime[1])
         second = int(myTime[2])
         hourFloat = hour + minute / 60.0 + second / 3600.0
         grid["values"][0] = hourFloat
         if (False):
            hourFloat += 3.0	# What are the multiple time values? What would they mean?
            grid["values"][1] = hourFloat
         else:
            # reduce template to one time value
            del grid["values"][1]

      if (grid["name"].lower() == "height"):
         # load the standard height grid
         heightGrid = tuvUtils.intervalToValues(grid)
         #progress("Loaded height profile {}<br>.".format(heightGrid))

         # adjust vertical profile to ground elevation
         compressedHeights = tuvUtils.raiseGround(heightGrid["values"],
            tuvValues.gElevation)
         heightGrid["values"] = compressedHeights
         #progress("Raised height profile {}<br>.".format(heightGrid))

      if (grid["name"].lower() == "wavelength"
         and tuvValues.customWavelength):
         # convert wavelength grid to equal interval
         grid["type"] = "equal interval"
         del grid["file path"]
         grid["begins at"] = tuvValues.waveStart
         grid["ends at"] = tuvValues.waveEnd
         grid["cell delta"] = ((tuvValues.waveEnd - tuvValues.waveStart)
            / tuvValues.waveIntervals)

   # Actinic flux
   radXfer = jStruct["radiative transfer"]
   #progress("radXfer = {}<br>".format(radXfer))
   radXfer["__output"] = tuvValues.actinicFlux
   #progress("radXfer2 = {}<br>".format(radXfer))

   radiators = radXfer["radiators"]

   # clouds
   if (tuvValues.cloudOptDepth > 0.0):
      clouds = {
         "name": "clouds",
         "type": "from netcdf file",
         "file path": CLOUDS_FILENAME
      }
      radiators.append(clouds)

   # aerosols
   radiatorToDelete = None
   aerosols = None
   for radi, radiator in enumerate(radiators):
      #progress("radiator = {}".format(radiator))
      if radiator["name"] == "aerosols":
         # retrieve the standard optical depth profile
         tuvValues.stdOpticalDepths = radiator["optical depths"]
         if (False):
            progress("Standard opticalDepthGrid = {} length {} sum {}"
               .format(tuvValues.stdOpticalDepths, len(tuvValues.stdOpticalDepths),
               numpy.sum(tuvValues.stdOpticalDepths)))

         if (False):
            # modify the JSON template
            radiator["single scattering albedo"] = tuvValues.aerosolSingScatAlbedo
            radiator["asymmetry factor"] = tuvValues.aerosolAsymmetry
            radiator["550 nm optical depth"] = tuvValues.aerosolOptDepth
         else:
            # remove JSON and replace with NetCDF file
            radiatorToDelete = radi

            aerosols = {
               "name": "aerosols",
               "type": "from netcdf file",
               "file path": AEROSOLS_FILENAME
            }

   # Insert the new layers of optical depth;
   # these two profiles have to be inserted together.
   insertedHeights = tuvUtils.insertOpticalDepths(heightGrid["values"],
      tuvValues.stdOpticalDepths, heightsToInsert)
   heightGrid["values"] = insertedHeights
   if (False):
      progress("Inserted height profile {} length {}<br>."
         .format(heightGrid["values"], len(heightGrid["values"])))
      progress("Inserted opticalDepthGrid = {} length {} sum {}"
         .format(tuvValues.stdOpticalDepths, len(tuvValues.stdOpticalDepths),
         numpy.sum(tuvValues.stdOpticalDepths)))

   # map air density and temperature onto the new height grid
   tuvUtils.interpolateProfile(densityProfile, densityDataFrame,
      heightGrid["values"], tuvValues.height, tuvValues.density);
   tuvUtils.interpolateProfile(temperatureProfile, temperatureDataFrame,
      heightGrid["values"], tuvValues.height, tuvValues.temperature);

   if (radiatorToDelete is not None):
      del radiators[radiatorToDelete]
   if (aerosols is not None):
      radiators.append(aerosols)

   # solver
   if (tuvValues.numStreams != 2):
      solver = radXfer["solver"]
      solver["type"] = "discrete ordinate"
      solver["number of streams"] = tuvValues.numStreams

   # Cross sections
   setAllDictKeys(jStruct["photolysis"], "cross section", "__output", tuvValues.crossSections)

   # Quantum yields
   setAllDictKeys(jStruct["photolysis"], "quantum yield", "__output", tuvValues.quantumYields)

   return



# Safely move to new directory and guarantee return.
@contextlib.contextmanager
def pushd(new_dir):
   previous_dir = os.getcwd()
   os.chdir(new_dir)
   try:
      yield
   finally:
      os.chdir(previous_dir)



# Determine and return dimensions for clouds
# myConfig = the already-modified JSON configuration
# myStage = stage directory for this config
# return tuple of (numLevels, numWavelengths, heightLevels, wavelengths)
def getCloudDimSizes(myConfig, myStage):
   # set up some defaults in case grids not found
   numLevels = 234
   heightValues = None
   numWavelengths = 345
   wavelengthFile = None
   wavelengthValues = None

   # find the height grid
   # Is there a better way to find the named grids than looping?
   #    Carl Drews - October 15, 2023
   for grid in myConfig["grids"]:
      if (grid["name"] == "height"):
         heightValues = grid["values"]

      if (grid["name"] == "wavelength"):
         if (grid["type"] == "from csv file"):
            wavelengthFile = grid["file path"]

         if (grid["type"] == "equal interval"):
            beginWave = grid["begins at"]
            endWave = grid["ends at"]
            delta = grid["cell delta"]
            numWavelengths = round((endWave - beginWave) / delta) + 1
            progress("numWavelengths = {}".format(numWavelengths))
            wavelengthValues = numpy.linspace(beginWave, endWave, numWavelengths)

   if (heightValues is not None):
      numLevels = len(heightValues)

   if (wavelengthFile is not None):
      progress("wavelength file = {}".format(myStage + wavelengthFile))
      # first line of this file contains grid size
      fpWave = open(myStage + wavelengthFile, 'r')
      firstLine = fpWave.readline()
      numWavelengths = int(firstLine)

      wavelengthValues = numpy.zeros(numWavelengths)
      for wi, oneLine in enumerate(fpWave):
         wavelengthValues[wi] = float(oneLine.strip())
      fpWave.close()
      #progress("wavelengthValues = {}".format(wavelengthValues))

   # both dimensions refer to the layers *between* levels
   numLevels -= 1
   numWavelengths -= 1

   if (False):
      # calculate the heights between the original levels
      betweenHeights = numpy.zeros(numLevels)
      for hi in range(numLevels):
         betweenHeights[hi] = (heightValues[hi] + heightValues[hi + 1]) / 2
      #progress("betweenHeights = {}".format(betweenHeights))

   # calculate the wavelengths between the original grid
   betweenWavelengths = numpy.zeros(numWavelengths)
   for wi in range(numWavelengths):
      betweenWavelengths[wi] = (wavelengthValues[wi] + wavelengthValues[wi + 1]) / 2

   return((numLevels, numWavelengths, heightValues, betweenWavelengths))



# Extract photolysis numerical output from NetCDF dataset.
# fpDataset = created by opening a NetCDF file
# myHeight = height to display, or None to display surface
# return block of HTML text, and a pandas DataFrame for CSV
def buildPhotoOutput(fpDataset, myHeight):
   dataText = ""
   dataText += "<pre>"

   # set up for pandas DataFrame
   csvNames = []
   csvValues = []

   levelIndex = 0
   waveIndex = 0
   timeIndex = 0

   if (myHeight is not None):
      altValues = fpDataset.variables["altitude"][:]
      levelIndex = tuvUtils.getHeightIndex(altValues, myHeight)

   #progress("levelIndex = {}".format(levelIndex))

   # loop through all the variables
   varNames = fpDataset.variables
   for vari, varName in enumerate(varNames):
      #progress("varName[{}] = {}<br>".format(vari, varName))
      dataText += "{:>4} ".format(vari)
      dataText += "{:<39} ".format(varName)

      # extract a useful single cell for display
      myVar = fpDataset.variables[varName]
      #progress("myVar dimensions = {}<br>".format(myVar.dimensions))

      # loop through the dimensions and extract 3 indexes
      indexes = []
      for dimName in myVar.dimensions:
         #progress(" {}".format(dimName), endString="")
         if ("level" in dimName):
            indexes.append(levelIndex)
         if ("wavelength" in dimName):
            indexes.append(waveIndex)
         if ("time" in dimName):
            indexes.append(timeIndex)
      #progress(" ")
      #progress("indexes = {}".format(indexes))

      # retrieve that cell by indexes
      varValue = -99999.9
      varValues = fpDataset.variables[varName][:]
      if (len(varValues.shape) == 1):
         varValue = varValues[indexes[0]]
      elif (len(varValues.shape) == 2):
         varValue = varValues[indexes[0], indexes[1]]
      elif (len(varValues.shape) == 3):
         varValue = varValues[indexes[0], indexes[1], indexes[2]]

      # format floating-point outpu
      dataText += "{:1.5e}<br>".format(varValue)

      csvNames.append(varName)
      csvValues.append(varValue)

   dataText += "</pre>"

   # assemble the pandas DataFrame
   csvFrame = pandas.DataFrame(
      data= {
         "Variable": csvNames,
         "Value": csvValues
      }
   )

   return(dataText, csvFrame)



# index constants for readability
#DIRECT = 0
#SUNDOWN = 1
#SUNUP = 2

# Extract flux numerical output from NetCDF dataset.
# fpDataset = created by opening a NetCDF file
# myHeight = height to display, or None to display surface
# sunlight[] = array of multipliers for direct, down, and up
# return pandas DataFrame
def buildFluxOutput(fpDataset, myHeight, sunlight):
   levelIndex = 0
   timeIndex = 0

   if (myHeight is not None):
      altValues = fpDataset.variables["altitude"][:]
      levelIndex = tuvUtils.getHeightIndex(altValues, myHeight)

   if (False):
      fluxFrame = pandas.DataFrame(
         data={"Hello": [1.2, 3.4, 5.6],
            "Flux": [7,8,9],
            "Output": [100,200,300]}
      )

   # start with the wavelengths
   colTitle = "wavelength"
   fluxFrame = pandas.DataFrame(
      data={colTitle: fpDataset.variables[colTitle][:]}
   )

   # append the radiations
   radVars = ["direct radiation", "downward radiation", "upward radiation"]
   for colTitle in radVars:
      fluxFrame.insert(len(fluxFrame.columns),
         colTitle, fpDataset.variables[colTitle][:, levelIndex, 0])

   # calculate and append the total
   total = None
   for cti, colTitle in enumerate(radVars):
      if (total is None):
         total = fpDataset.variables[colTitle][:, levelIndex, 0] * sunlight[cti]
      else:
         total += fpDataset.variables[colTitle][:, levelIndex, 0] * sunlight[cti]

   fluxFrame.insert(len(fluxFrame.columns), "total", total)

   return(fluxFrame)



# Extract photolysis numerical output from NetCDF dataset,
# retrieving all variables that have a certain prefix.
# fpDataset = created by opening a NetCDF file
# myHeight = height to display, or None to display surface
# varPrefix = probably "cross section" or "quantum yield"
# return pandas DataFrame
def buildPrefixOutput(fpDataset, myHeight, varPrefix):
   levelIndex = 0
   timeIndex = 0

   if (myHeight is not None):
      altValues = fpDataset.variables["altitude"][:]
      levelIndex = tuvUtils.getHeightIndex(altValues, myHeight)

   #progress("levelIndex = {}".format(levelIndex))

   # start with the wavelengths
   colTitle = "wavelength"
   varsFrame = pandas.DataFrame(
      data={colTitle: fpDataset.variables[colTitle][:]}
   )

   # append the variables with the specified prefix
   varNames = fpDataset.variables
   for vari, varName in enumerate(varNames):
      #progress("prefix = {}   varName = {}".format(varPrefix, varName))
      if (not varName.startswith(varPrefix)):
         continue
      varsFrame.insert(len(varsFrame.columns),
         varName, fpDataset.variables[varName][:, levelIndex, 0])

   return(varsFrame)



# Create URL to web directory where user can download.
# myStageBase = base directory in the local file system
# dateTimeStamp = unique ID for directory within staging area
# return URL of indexed directory suitable for downloading
def createStageLink(myStageBase, dateTimeStamp):
   # index the stage directory
   cmd = "python3 webIndex.py {}tuvxDir{}.txt".format(myStageBase, dateTimeStamp)
   #progress("cmd = {}".format(cmd))
   proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

   # create http link to the directory
   stageURL = "https://www.acom.ucar.edu/tuv-x/STAGE/{}/".format(dateTimeStamp)
   stageLink = "<a href=\"{}\">{}</a>".format(stageURL, stageURL)

   return(stageLink)



TUVX_BUILD_DIR = "tuvx-build/tuv-x/build/"
STAGE_BASEDIR = "/data/TUV-x/stage/"

# Main program begins here.
def main():
   mainValue = 0

   # provide information for the log file
   progress("<h3>Hello, TUV-x Remote World!</h3>")

   # collect the command-line arguments
   args = tuvCore.tuvParams()
   for argPair in sys.argv:
      #progress("argPair = {}<br>".format(argPair))
      pairValue = argPair.split('=')
      #progress("pairValue = {}<br>".format(pairValue))
      if (len(pairValue) < 2):
         continue

      if (pairValue[0].lower() == "customwavelength"):
         args.customWavelength = (pairValue[1].lower() == "true")
      if (pairValue[0].lower() == "wavestart"):
         args.waveStart = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "waveend"):
         args.waveEnd = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "waveintervals"):
         args.waveIntervals = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "latitude"):
         args.latitude = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "longitude"):
         args.longitude = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "date"):
         args.date = utilsLite.sanitize(pairValue[1])
      if (pairValue[0].lower() == "time"):
         args.time = utilsLite.sanitize(pairValue[1])

      if (pairValue[0].lower() == "inputoption"):
         args.inputOption = int(pairValue[1])
      if (pairValue[0].lower() == "zenithangle"):
         args.zenithAngle = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "photolysis"):
         args.photolysis = (pairValue[1].lower() == "true")
      if (pairValue[0].lower() == "actinic"):
         args.actinicFlux = (pairValue[1].lower() == "true")
      if (pairValue[0].lower() == "crosssections"):
         args.crossSections = (pairValue[1].lower() == "true")
      if (pairValue[0].lower() == "quantumyields"):
         args.quantumYields = (pairValue[1].lower() == "true")

      if (pairValue[0].lower() == "ozone"):
         args.ozone = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "albedo"):
         args.albedo = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "gelevation"):
         args.gElevation = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0].lower() == "height"):
         args.height = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "customprofile"):
         args.customProfile = (pairValue[1].lower() == "true")
      if (pairValue[0].lower() == "temperature"):
         args.temperature = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0].lower() == "density"):
         args.density = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0] == "cloudOptDepth"):
         args.cloudOptDepth = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "cloudSingScatAlbedo"):
         args.cloudSingScatAlbedo = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "cloudAsymmetry"):
         args.cloudAsymmetry = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "cloudTop"):
         args.cloudTop = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "cloudBase"):
         args.cloudBase = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0] == "aerosolOptDepth"):
         args.aerosolOptDepth = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "aerosolSingScatAlbedo"):
         args.aerosolSingScatAlbedo = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "aerosolAsymmetry"):
         args.aerosolAsymmetry = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "aerosolAlpha"):
         args.aerosolAlpha = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "aerosolLambda0"):
         args.aerosolLambda0 = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0] == "directSun"):
         args.directSunlight = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "diffDown"):
         args.diffuseDown = utilsLite.safeFloat(pairValue[1])
      if (pairValue[0] == "diffUp"):
         args.diffuseUp = utilsLite.safeFloat(pairValue[1])

      if (pairValue[0] == "numStreams"):
         args.numStreams = int(pairValue[1])

   progress("TUV-x-remote args = {}<br>".format(args.toString()))

   # Create a unique ID to stage this request in a directory.
   # The date-time string should provide a unique filename.
   nowDate = datetime.datetime.now()
   uniqueID = nowDate.strftime("%Y%m%d-%H%M%S--%f")        # year through microseconds
   progress("uniqueID = {}<br>".format(uniqueID))

   if (not args.customProfile):
      # load standard profiles at the specified measurement height
      progress("Loading standard profile at altitude . . .")
      args.temperature, args.density = tuvUtils.getStandardProfile(args.height)
      progress("height = {}   temperature = {}   density = {}"
         .format(args.height, args.temperature, args.density))
      progress(". . . done.")

   # Read template json configuration file.
   jsonConfig = None
   jsonFile = TUVX_BUILD_DIR + "examples/full_config.json"
   progress("jsonFile = {}".format(jsonFile))
   try:
      with open(jsonFile, 'r') as jFile:
         jsonConfig = json.load(jFile)
   except FileNotFoundError as oops:
      progress("{}<br>".format(oops))
      return(1)

   #progress("jsonConfig = {}<br>"
   #   .format(json.dumps(jsonConfig, sort_keys=True, indent=3)))

   # Substitute entered values.
   try:
      replaceJson(jsonConfig, args, TUVX_BUILD_DIR)
   except KeyError as oops:
      progress("KeyError in replaceJson(): {}<br>".format(oops))
      return(1)
   #progress("modified json[] = {}<br>"
   #   .format(json.dumps(jsonConfig["profiles"], indent=3)))

   #  Write a new json file to the staging area.
   stageDir = STAGE_BASEDIR + "{}/".format(uniqueID)
   try:
      os.mkdir(stageDir)
   except FileNotFoundError as oops:
      progress("{}<br>".format(oops))
      return(1)

   jsonFilename = "tuvx_config.json"
   jsonFilePath = "{}{}".format(stageDir, jsonFilename)
   fpJson = open(jsonFilePath, "w")
   json.dump(jsonConfig, fpJson, indent=3)
   fpJson.close()

   # create symbolic link in stage to data directory under TUV-x build directory
   workingDir = os.getcwd() + "/"
   progress("workingDir = {}<br>".format(workingDir))
   os.symlink(workingDir + TUVX_BUILD_DIR + "data", stageDir + "data")

   numLevels, numWaves, myHeightLevels, myWavelengths = getCloudDimSizes(
      jsonConfig, stageDir)
   if (args.cloudOptDepth > 0.0):
      # create clouds NetCDF file containing cloud radiator
      cloudFilename = "{}{}".format(stageDir, CLOUDS_FILENAME)
      progress("Cloud file {} will contain {} levels and {} wavelengths."
         .format(cloudFilename, numLevels, numWaves))
      tuvUtils.writeClouds(args, numLevels, numWaves, myHeightLevels, cloudFilename)

   # create aerosols NetCDF file containing aerosol radiator
   aerosolFilename = "{}{}".format(stageDir, AEROSOLS_FILENAME)
   progress("Aerosol file {} will contain {} levels and {} wavelengths."
      .format(aerosolFilename, numLevels, numWaves))
   tuvUtils.writeAerosols(args, numLevels, numWaves,
      myHeightLevels, myWavelengths, aerosolFilename)

   # move to the stage directory for relative file paths
   with pushd(stageDir):
      progress("Moving to stage directory {}<br>".format(stageDir))

      # Call Fortran executable tuv-x.
      cmd = "export LD_LIBRARY_PATH=\"/opt/local/lib64:/opt/local/lib:$LD_LIBRARY_PATH\""
      cmd += "; " + workingDir + TUVX_BUILD_DIR + "tuv-x" + " " + jsonFilename
      progress("cmd = {}<br>".format(cmd))
      proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
      progress("stdout = {}<br>".format(proc.stdout.decode("utf-8")))
      progress("stderr = {}<br>".format(proc.stderr.decode("utf-8")))

      # Read Fortran output from the staging area.
      fpPhoto = netCDF4.Dataset("photolysis_rate_constants.nc", 'r')

      # Extract results from NetCDF and pass them back to web3 for display.
      photoText, photoData = buildPhotoOutput(fpPhoto, args.height)
      #progress("photoData = {}".format(photoData))

      fluxOutput = None
      if (args.actinicFlux):
         # extract just the actinic flux
         fluxOutput = buildFluxOutput(fpPhoto, args.height,
            [args.directSunlight, args.diffuseDown, args.diffuseUp])

         # save flux to CSV file
         fluxOutput.to_csv("actinicFlux.csv", index_label="index")

         # save flux to text file for easier browsing
         with open("actinicFlux.txt", 'w') as fpText:
            fluxString = fluxOutput.to_string()
            fpText.write(fluxString)

      prefixes = ["cross section", "quantum yield"]
      baseFilenames = ["crossSections", "quantumYields"]
      for pfi, prefix in enumerate(prefixes):
         if (prefix == "cross section" and not args.crossSections):
            continue
         if (prefix == "quantum yield" and not args.quantumYields):
            continue

         prefixOutput = buildPrefixOutput(fpPhoto, args.height, prefix)
         #progress("{} Output = {}".format(prefix, prefixOutput))

         # save flux to CSV file
         baseFilename = baseFilenames[pfi]
         prefixOutput.to_csv("{}.csv".format(baseFilename), index_label="index")

         # save flux to text file for easier browsing
         with open("{}.txt".format(baseFilename), 'w') as fpText:
            prefixString = prefixOutput.to_string()
            fpText.write(prefixString)

      # close NetCDF file
      fpPhoto.close()

      # save j-values to CSV file
      photoData.to_csv("photoRateConstants.csv", index_label="index")

      # save j-values to text file for easier browsing
      with open("photoRateConstants.txt", 'w') as fpText:
         photoString = photoData.to_string()
         fpText.write(photoString)

      if (False):
         # Read Fortran output from the staging area.
         fpDose = netCDF4.Dataset("dose_rates.nc", 'r')

         # What are we going to retrieve from this file? Anything?

         # close NetCDF file
         fpDose.close()

   # create directory file for webIndex.py
   dirFile = STAGE_BASEDIR + "tuvxDir{}.txt".format(uniqueID)
   fp = open(dirFile, 'w')
   fp.write("{}\n".format(stageDir))
   fp.close()

   # remove symbolic link so download does not follow it
   os.unlink(stageDir + "data")

   # pass text summary back to the browser on web3
   progress("Date = {}   time = {}<br>".format(args.date, args.time))
   progress(photoText)

   # create URL indexed directory for user to download files
   stageURL = createStageLink(STAGE_BASEDIR, uniqueID)
   progress("Download TUV-x output at: {}<br>".format(stageURL))

   return(mainValue)	# 0 = no error



# call the main
progress("{}<br>".format(__file__))
progress("Start time remote: {} Local<br>".format(datetime.datetime.now()))

retValue = main()

progress("End time remote: {} Local<br>".format(datetime.datetime.now()))
sys.exit(retValue)

