#!/usr/bin/env python3

# tuvCore.py
# Python module to provide basic support for TUV-x calculations.
#
# Author: Carl Drews, Atmospheric Chemistry Observations & Modeling
# Created: February 2023
# Copyright (C) 2023-2025 University Corporation for Atmospheric Research



import sys



# Display progress message on the console.
# endString = set this to '' for no return
def progress(message, endString='\n'):
   if (True):   # disable here in production
      print(message, end=endString)
      sys.stdout.flush()



# Encapsulate the parameters passed to TUV-x.
class tuvParams:
   # Create default set of parameters.
   def __init__(self, lat=0.0, lon=0.0,
      myDate = "20230308", myTime = "12:34:56"):

      # custom wavelength
      self.customWavelength = False
      self.waveStart = 100
      self.waveEnd = 500
      self.waveIntervals = 80

      # geolocation
      self.latitude = lat
      self.longitude = lon

      # date and time for this calculation
      self.date = myDate
      self.time = myTime

      # Input Option and solar zenith angle
      self.inputOption = 1
      self.zenithAngle = None

      # output
      self.photolysis = False
      self.actinicFlux = False
      self.crossSections = False
      self.quantumYields = False

      # other input parameters
      self.ozone = 280.0
      self.albedo = 0.16
      self.gElevation = 1.655		# km above sea level

      # single values at a certain level
      self.customProfile = False
      self.height = None
      self.temperature = None
      self.density = None

      # clouds
      self.cloudOptDepth = 0.245
      self.cloudSingScatAlbedo = 0.945
      self.cloudAsymmetry = 0.845
      self.cloudTop = 5.9
      self.cloudBase = 4.1

      # aerosols
      self.aerosolOptDepth = 0.24
      self.aerosolSingScatAlbedo = 0.9
      self.aerosolAsymmetry = 0.6
      self.aerosolAlpha = 1.0
      self.aerosolLambda0 = 550
      self.stdOpticalDepths = None	# radiators | optical depths[] in the JSON config

      # sunlight
      self.directSunlight = 1.0
      self.diffuseDown = 1.0
      self.diffuseUp = 1.0

      # solver
      self.numStreams = 2
      return

   def toString(self):
      myString = ""

      myString += (" customWavelength:{}".format(self.customWavelength))
      myString += (" waveStart:{}".format(self.waveStart))
      myString += (" waveEnd:{}".format(self.waveEnd))
      myString += (" waveIntervals:{}".format(self.waveIntervals))

      myString += (" latitude:{:06f} longitude:{:06f}"
         .format(self.latitude, self.longitude))
      myString += (" date:{} time:{}"
         .format(self.date, self.time))

      myString += (" inputOption:{}".format(self.inputOption))
      myString += (" zenithAngle:{}".format(self.zenithAngle))

      myString += (" photolysis:{}".format(self.photolysis))
      myString += (" actinic flux:{}".format(self.actinicFlux))
      myString += (" cross sections:{}".format(self.crossSections))
      myString += (" quantum yields:{}".format(self.quantumYields))

      myString += (" ozone:{}".format(self.ozone))
      myString += (" albedo:{}".format(self.albedo))
      myString += (" gElevation:{}".format(self.gElevation))

      myString += (" customProfile:{}".format(self.customProfile))
      myString += (" height:{}".format(self.height))
      myString += (" temperature:{}".format(self.temperature))
      myString += (" density:{}".format(self.density))

      myString += (" cloudOptDepth:{}".format(self.cloudOptDepth))
      myString += (" cloudSingScatAlbedo:{}".format(self.cloudSingScatAlbedo))
      myString += (" cloudAsymmetry:{}".format(self.cloudAsymmetry))
      myString += (" cloudTop:{}".format(self.cloudTop))
      myString += (" cloudBase:{}".format(self.cloudBase))

      myString += (" aerosolOptDepth:{}".format(self.aerosolOptDepth))
      myString += (" aerosolSingScatAlbedo:{}".format(self.aerosolSingScatAlbedo))
      myString += (" aerosolAsymmetry:{}".format(self.aerosolAsymmetry))
      myString += (" aerosolAlpha:{}".format(self.aerosolAlpha))
      myString += (" aerosolLambda0:{}".format(self.aerosolLambda0))

      myString += (" directSunlight:{}".format(self.directSunlight))
      myString += (" diffuseDown:{}".format(self.diffuseDown))
      myString += (" diffuseUp:{}".format(self.diffuseUp))

      myString += (" numStreams:{}".format(self.numStreams))

      return(myString)

