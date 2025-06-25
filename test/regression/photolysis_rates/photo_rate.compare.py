import sys
import numpy as np
from netCDF4 import Dataset

def compare_netcdf(file1, file2, tol=1.5e-1):
  ds1 = Dataset(file1, 'r')
  ds2 = Dataset(file2, 'r')

  # Create normalized variable name dictionaries
  ds1_vars = {k.replace(' ', '').replace('+hv','').lower(): v for k, v in ds1.variables.items()}
  ds2_vars = {k.replace(' ', '').replace('+hv','').lower(): v for k, v in ds2.variables.items()}

  vars1 = set(ds1_vars.keys())
  vars2 = set(ds2_vars.keys())

  # Remove non-dose rates from tuv-x
  vars2.discard('temperature')
  vars2.discard('Earth-Sun distance')
  vars2.discard('solar zenith angle')
  vars2.discard('altitude')
  vars2.discard('time')
  vars2.discard('wavelength')
  
  for var in vars1:
    if var not in vars2:
      continue
    arr1 = ds1_vars[var][:121, :]
    arr2 = ds2_vars[var][:]
    if arr1.shape != arr2.shape:
      print(f"Shape mismatch for variable '{var}': {arr1.shape} vs {arr2.shape}")
      sys.exit(1)
    if not np.allclose(arr1, arr2, rtol=tol, equal_nan=True):
      print(f"Values differ for variable '{var}' beyond tolerance {tol}")
      print("Values from file1 and file2:")
      for row1, row2 in zip(arr1, arr2):
        print(np.column_stack((row1, row2)))
      sys.exit(1)

  print("All variables match within tolerance.")

  ds1.close()
  ds2.close()

if __name__ == "__main__":
  if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <file1.nc> <file2.nc>")
    sys.exit(1)
  compare_netcdf(sys.argv[1], sys.argv[2])