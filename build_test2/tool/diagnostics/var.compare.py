import numpy
from numpy import unravel_index
import sys
from scipy.io import FortranFile
import os
import json

new_data_path = os.path.relpath("output")
old_data_path = os.path.relpath("odat/OUTPUTS")
n_vertical_bins = 120
n_wavelength_bins = 156
arrays_2D = [ "dtrl", "radField" ]

# Reads and returns a real array from a Fortran binary file
def get_variable_from_file(file_path, var_type):
    try:
        ffile = FortranFile(file_path, 'r')
    except:
        print(f"Error opening file {file_path}")
        sys.exit(1)
    try:
        var = ffile.read_reals(dtype=var_type)
    except:
        print(f"Error reading file {file_path}")
        sys.exit(1)
    return var


# Compares values of real arrays from Fortran binary files
# using provided tolerances
def compare_var(var_name):
    global data_path
    global n_vertical_bins
    global n_wavelength_bins
    global arrays_2D

    var_old = get_variable_from_file(os.path.join(old_data_path,var_name + ".old"), numpy.float64)
    var_new = get_variable_from_file(os.path.join(new_data_path,var_name + ".new"), numpy.float64)

    # compare array sizes
    if var_old.size != var_new.size :
        print(f"size of {var_name}.old and {var_name}.new are not the same")
        print(f"size {var_name}.old = {var_old.size}")
        print(f"size {var_name}.new = {var_new.size}")
        sys.exit(-1)

    # check for valid vertical dimensions
    if var_new.size % n_vertical_bins != 0 :
        n_vertical_bins += 1
        if var_new.size % n_vertical_bins != 0 :
            print(f"Error: {var_name} invalid spatial dimension")
            sys.exit(-1)

    # get percent difference
    diff = numpy.zeros( var_old.size )
    max_val = max( numpy.max(numpy.abs(var_old)),numpy.max(numpy.abs(var_new)) )
    threshold = 1.0e-20*max_val
    for n in range( var_old.size ):
        if( max( abs(var_old[n]),abs(var_new[n]) ) > threshold ):
            diff[n] = abs(var_old[n] / var_new[n] - 1.0) * 100.0

#   if numpy.any( var_new == 0.0 ) :
#       diff = var_old
#   else:
#       diff = ( var_old / var_new - 1.0 ) * 100.0

    # reshape 2D arrays
    if var_name in arrays_2D :
        var_new = numpy.reshape( var_new, ( n_vertical_bins, n_wavelength_bins ) )
        var_old = numpy.reshape( var_old, ( n_vertical_bins, n_wavelength_bins ) )
        diff    = numpy.reshape(    diff, ( n_vertical_bins, n_wavelength_bins ) )

    # get comparison stats
    results = {}
    results["minimum difference"] = numpy.amin( numpy.abs( diff ) )
    results["maximum difference"] = numpy.amax( numpy.abs( diff ) )
    results["RMS difference"]     = numpy.sqrt( numpy.mean( diff * diff ) )
    return results


# Perform an analysis of TUV output based on given metrics
def analyze_output(config):
    for var_name, options in config.items() :
        results = compare_var(var_name)
        for metric, tolerance in options.items() :
            if not metric in results.keys() :
                print(f"Error: invalid comparison metric for {var_name}: {metric}")
                sys.exit(-1)
            if tolerance < results[ metric ] :
                print(f"Error: comparison failure for {var_name}:")
                print(f"{var_name} {metric} = {results[metric]} > {tolerance}")
                print(f"       Max diff = {results['maximum difference']}")
                sys.exit(-1)
            else:
                print(f"{var_name} {metric} within tolerance: {results[metric]} <= {tolerance}")


"""
TUV output comparison
"""
if len(sys.argv) == 0 :
    print(f"Usage: python3 var.compare.py config.json [config2.json ...]")
    print(f"where config.json contains:")
    print(f"{{")
    print(f"  \"variable name\" : {{")
    print(f"    \"minimum difference\" : 12.3,")
    print(f"    \"maximum difference\" : 32.4,")
    print(f"    \"RMS difference\" : 13.7")
    print(f"  }}")
    print(f"}}")
    print(f"")
    print(f"Muliple variables per file are allowed and each difference tolerance is optional.")
    sys.exit(0)

# process each specified configuration file
for config_file_path in sys.argv[1:] :
    try:
        with open(config_file_path) as config_file :
            analyze_output(json.load(config_file))
    except IOError:
        print(f"Error opening file {config_file_path}")
        sys.exit(-1)
    except ValueError as e:
        print(f"Error parsing json file {config_file_path}: {e}")
        sys.exit(-1)


