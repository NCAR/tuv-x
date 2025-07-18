import numpy as np
import sys
from netCDF4 import Dataset
import netCDF4 as ncd
from datetime import datetime as dt

"""
Function to read the data file(s)
"""
def read_data_file(data_dictionary,dataTray):

  InpFileSpec = data_dictionary['filespec']
  try:
    InpFile = open(InpFileSpec,'r')
  except:
    print(f'Failed to open data file {InpFileSpec}')
    sys.exit(-3)

  print(f'Opened data file {InpFileSpec}')
  nLines  = len(InpFile.readlines())
  InpFile.seek(0)
  nskipHdr = data_dictionary['nPreSkip'] if 'nPreSkip' in data_dictionary else 0
  nRead    = data_dictionary['nRead'] if 'nRead' in data_dictionary else 0

  header = ''
# if header lines exist then read them
  if( nskipHdr > 0 ):
    for ndx in range(nskipHdr):
      header += InpFile.readline()
    InpFile.seek(0)

  nskipHdr = abs(nskipHdr)
  nskipEnd = nLines - (nskipHdr + nRead)
  print(f'nLines,nskipHdr,nRead,nskipEnd = {nLines},{nskipHdr},{nRead},{nskipEnd}')
  try:
    data = np.genfromtxt(InpFile,dtype='float64',skip_header=nskipHdr,skip_footer=nskipEnd,comments=None)
    print(f'Read cross section file {InpFileSpec}')
  except:
    print(f'Failed to read data file {InpFileSpec}')
    sys.exit(-2)

  try:
    dataTray.append(data)
  except:
    print('Failed to append data array to dataTray')
    sys.exit(-2)

  InpFile.close()
  print(f'Closed data file {InpFileSpec}')
  return(header)

"""
Function to write the netCDF file
"""
def stuff_netCDF_file(ncFile,interpolationTemps,dataTray,hasLambdaGrid,InpFileSpecs,var_typ,Headers):
  
  ndataVars = len(dataTray)
  print(f'ndataVars = {ndataVars}')

  for dataVarNdx in range(ndataVars):
    nparameterRow,nparameterCol = np.shape(dataTray[dataVarNdx])
    if( hasLambdaGrid ):
      nparameterCol -= 1
    ntemperature = min( len(interpolationTemps),nparameterCol )
    print(f'data array is ({nparameterRow},{nparameterCol})')
    DataTag = var_typ + "_parameters"
# define dimensions
    RowDimName = 'bins'
    ColDimName = 'parameters'
    TempDimName = 'temperatures'
    print(f'Variable type = {var_typ}')
    print(f'RowDimName,ColDimName = {RowDimName},{ColDimName}')
    ncFile.createDimension(RowDimName,nparameterRow)
    ncFile.createDimension(ColDimName,nparameterCol)
    ncFile.createDimension(TempDimName,ntemperature)
# create wavelength grid
    if( hasLambdaGrid ):
      Var = ncFile.createVariable('wavelength',np.float64,(RowDimName))
      Var.units = 'nm'
# write wavelength grid
      Var[:] = dataTray[dataVarNdx][:,0]
# create interpolation temperature array
    if( len(interpolationTemps) > 0 ):
      Var = ncFile.createVariable('temperature',np.float64,(TempDimName))
      Var.units = 'K'
# write interpolation temperature array
      Var[:] = interpolationTemps[:ntemperature]
# create cross section or quantum yield data array
    Var = ncFile.createVariable(DataTag,np.float64,(ColDimName,RowDimName))
    Var.hdr = Headers[dataVarNdx]
# write data array
    if( hasLambdaGrid ):
      Var[:,:] = np.transpose(dataTray[dataVarNdx][:,1:])
    else:
      Var[:,:] = np.transpose(dataTray[dataVarNdx])
    if( var_typ == 'cross_section' ):
      if( hasLambdaGrid ):
        Var.units = 'cm^2'
      else:
        Var.units = 'see source code'
    else:
      Var.units = 'fraction'

# global attributes
  version = '1.0'
  ncFile.Author = 'TUV Data Xformer ' + version
  now = dt.now()
  ncFile.created = now.strftime("%Y-%m-%d %H:%M:%S")
  if( var_typ == 'cross_section' ):
    ncFile.title = 'Cross section parameters'
  else:
    ncFile.title = 'Quantum yield parameters'
  ncFile.file = InpFileSpecs

"""
Transform ascii data file(s) to netCDF counterpart
"""
def xform_to_netCDF(nFile,phtDictionary,ncd_path):

# cross section
  if( 'cross-sections' in phtDictionary ):
# form netCDF file for cross sections
    molecule = phtDictionary['molecule']

    print(f'\nProcessing {molecule} cross section')
    xsects = phtDictionary['cross-sections']
    nxsects = len(xsects)
    print(f'  There are {nxsects} cross section files')
    dataTray = []
    interpolationTemps = []
    Headers = []

# loop over ascii input data files
    for xsect in xsects:
      ncdFilespec = ncd_path + '/' + molecule + '.nc'
# create the netcdf dataset
      print(f'\nCreating netCDF file {ncdFilespec}')

      try:
        ncFile = Dataset(ncdFilespec,mode='w',format='NETCDF4_CLASSIC')
      except:
        print(f'Failed to create netCDF4 dataset {ncdFilespec}')
        sys.exit(-1)

# 1st data column wavelength grid?
      if( 'has lambda grid' in xsect ):
        hasLambdaGrid = xsect['has lambda grid']
      else:
        hasLambdaGrid = True

# interpolation temperatures?
      if( 'interpolation temperature' in xsect ):
        interpolationTemps = xsect['interpolation temperature']
        print("\nInterpolation temperatures:")
        print(interpolationTemps)

      header = ''
      header = read_data_file(xsect,dataTray)
      Headers.append(header)

      InpFileSpecs = xsect['filespec']

      print(f'\nThere are {len(dataTray)} arrays in dataTray')
      print(f'\nShape data array is {dataTray[0].shape}')
      stuff_netCDF_file(ncFile,interpolationTemps,dataTray,hasLambdaGrid,InpFileSpecs,'cross_section',Headers)
# close netcdf file
      ncFile.close()

  print('\n')

# quantum yield
  if( 'quantum-yields' in phtDictionary ):
# form netCDF file for cross sections
    molecule = phtDictionary['molecule']

    print(f'\nProcessing {molecule} quantum yield')
    qylds = phtDictionary['quantum-yields']
    nqylds = len(qylds)
    print(f'  There are {nqylds} quantum yield files')
    dataTray = []
    interpolationTemps = []
    Headers = []

# loop over ascii input data files
    for qyld in qylds:
      ncdFilespec = ncd_path + '/' + molecule + '_quantum_yield_' + str(nFile) + '.nc'
# create the netcdf dataset
      print(f'\nCreating netCDF file {ncdFilespec}')

      try:
        ncFile = Dataset(ncdFilespec,mode='w',format='NETCDF4_CLASSIC')
      except:
        print(f'Failed to create netCDF4 dataset {ncdFilespec}')
        sys.exit(-1)

# 1st data column wavelength grid?
      if( 'has lambda grid' in qyld ):
        hasLambdaGrid = qyld['has lambda grid']
      else:
        hasLambdaGrid = True

# interpolation temperatures?
      if( 'interpolation temperature' in qyld ):
        interpolationTemps = qyld['interpolation temperature']
        print("\nInterpolation temperatures:")
        print(interpolationTemps)

      header = ''
      header = read_data_file(qyld,dataTray)
      Headers.append(header)

      InpFileSpecs = qyld['filespec']

      print(f'\nThere are {len(dataTray)} arrays in dataTray')
      print(f'\nShape data array is {dataTray[0].shape}')
      stuff_netCDF_file(ncFile,interpolationTemps,dataTray,hasLambdaGrid,InpFileSpecs,'quantum_yield',Headers)
# close netcdf file
      ncFile.close()
