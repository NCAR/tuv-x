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
   filespec = 'OUTPUTS/' + variable + '.old'
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
#--------------------------
# get photo-decomp variable
#--------------------------
   filespec = 'OUTPUTS/' + variable + '.new'
   try:
      filenew = FF(filespec,'r')
   except:
      print(f"Error opening file {filespec}")
      sys.exit(-1)

   try:
      varnew = filenew.read_reals(dtype=np.float64)
   except:
      print(f"Error reading file {filespec}")
      sys.exit(-1)

   sizenew = varnew.size
   if( sizeold != sizenew ):
     print(f"size {variable}.old,{variable}.new are not the same") 
     print(f"size {variable}.old = {sizeold}")
     print(f"size {variable}.new = {sizenew}")

   if( variable != 'radField' ):
      nz = 120
   else:
      nz = 121
   nbins = 156
#  if( (varnew.size % nz) != 0 ):
#     nz = 121
#     if( (varnew.size % nz) != 0 ):
#        print(f"Error {variable} spatial dim is neither 120 or 121")
#        sys.exit(-1)
#  print(f"\nnz = {nz}")

#  is2d = variable == 'dtrl' or variable == 'radField' or variable == 'dtaer' or variable == 'dto3'
   is2d = varnew.size/nz > 1

#--------------------------
# zero value check
#--------------------------
   if( np.any(varnew == 0.) ):
      if( np.all(varnew == 0.) ):
         print(f"{variable}.new == 0")
      else:
         print(f"some {variable}.new == 0")

   diff = np.zeros(varnew.size)

   for n in range(varnew.size):
      if( varnew[n] == 0. ):
         if( varold[n] == 0. ):
            diff[n] = 0.
         else:
            diff[n] = -100.
      else:
         diff[n] = (varold[n]/varnew[n] - 1.)*100.

   diff = abs(diff)

   if( is2d ):
     varnew = np.reshape( varnew,(nz,nbins))
     varold = np.reshape( varold,(nz,nbins))
     diff = np.reshape( diff,(nz,nbins))
     if( variable == 'radField' ):
        print(f"\nradField.new  @ ground")
        print(f"{varnew[0,:]}")
        print(f"\nradField.old  @ ground")
        print(f"{varold[0,:]}")
        print(f"\nradField.new @ top atm")
        print(f"{varnew[120,:]}")
        print(f"\nradField.old @ top atm")
        print(f"{varold[120,:]}")
        print(f"\nradField.new @ nlambda = 1")
        print(f"{varnew[:,0]}")
        print(f"\nradField.old @ nlambda = 1")
        print(f"{varold[:,0]}")

   mindiff = np.amin(diff)
   maxdiff = np.amax(diff)
   print(f"\n{variable} min,max diff = {mindiff}%, {maxdiff}%")
   print(f"\n{variable} RMS-diff = {np.sqrt(np.mean(diff*diff))}%")
   maxind = unravel_index(diff.argmax(),diff.shape)
   print(f"{variable} diff max index   = {diff.argmax()}")
   if( is2d ):
      print(f"maxind = {maxind[0]}, {maxind[1]}")
      print(f"{variable}.old @ max = {varold[maxind[0],maxind[1]]}")
      print(f"{variable}.new @ max = {varnew[maxind[0],maxind[1]]}")
      print(f"{variable} max diff = {diff[maxind[0],maxind[1]]}")
   else:
      print(f"{variable}.old,new @ max diff = {varold[maxind]}, {varnew[maxind]}")
