#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

import numpy as np
import sys
import json
from xsqy_subs import xform_to_netCDF

#-----------------------------------------------------
#  json config file is only possible argument
#-----------------------------------------------------
if( len(sys.argv) > 2 ):
  print(f'\n{sys.argv[0]}: requires one or no arguments')
  sys.exit( -1)
elif( len(sys.argv) == 2 ):
  filespec = sys.argv[1]
else:
  filespec = 'photo.config.json'

  
#-----------------------------------------------------
#  open json photo config file
#-----------------------------------------------------
#ilespec = 'photo.config.tst.json'
try:
  fp = open(filespec,'r')
except:
  print(f"Failed to open {filespec}")
  sys.exit(-1)

#-----------------------------------------------------
#  transfer config file into dictionary
#-----------------------------------------------------
try:
  photDict = json.load(fp)
except:
  print(f"Failed to load json file {filespec}")
  sys.exit(-1)

#-----------------------------------------------------
#  done with json input file
#-----------------------------------------------------
fp.close()

#-----------------------------------------------------
#  loop through photo reactions in dictionary
#-----------------------------------------------------
list = photDict['photoreactions']
molecule = ""
for rxt in list:
#-----------------------------------------------------
#  call xform_to_netCDF
#-----------------------------------------------------
  if( molecule != rxt['molecule']):
    nFile = 1
    molecule = rxt['molecule']
  else:
    nFile += 1
  xform_to_netCDF(nFile,rxt,'./')

print('\n')
print(f'\nThere are {len(list)} photoreactions in {filespec}')
