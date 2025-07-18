#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

import numpy as np
from numpy import unravel_index
import sys
from scipy.io import FortranFile as FF

"""
The main program
"""
while True:
#--------------------------
# get the variable name
#--------------------------
   prompt = '\nEnter diagnostic variable: '
   try:
      variable = input(prompt)
   except EOFError:
      print(f"Error with input {variable} diagnostic")
      sys.exit(-1)
   if( variable == 'quit' ): break
#--------------------------
# get TUV variable
#--------------------------
   filespec = 'OUTPUTS/' + variable
   try:
      fileold = FF(filespec,'r')
   except:
      print(f"Error opening file {filespec}")
      sys.exit(-1)

   try:
      varold = fileold.read_reals(dtype=np.float64)
   except:
      print(f"Error reading file {filespec}")
      sys.exit(-1)
   
   sizeold = varold.size

   nz = 120
   nbins = 156
   if( (varold.size % nz) != 0 ):
      nz = 121
      if( (varold.size % nz) != 0 ):
         print(f"Error {variable} spatial dim is neither 120 or 121")
         sys.exit(-1)
   print(f"\nnz = {nz}")

   is2d = variable == 'dtrl' or variable == 'radField'

   print(f"\n{variable}")
   print(varold)

