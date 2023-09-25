#!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3

import os
import numpy as np
from scipy.io import FortranFile as FF
import json
import sys

#return index of annotated label
def get_list_index( list, match  ):
   index = -1
   for n in np.arange(len(list)):
      if( list[n].strip() == match ):
         index = n
         break
   return index

# Compares values of real arrays from Fortran binary files
# using provided tolerances
def compare_var(var_name,tolerance,var_old,var_new):

    # compare array sizes
    if var_old.size != var_new.size :
        raise RuntimeError('Mismatched arrays: {var_name}.old size = {var_old.size}; {var_name}.new size = {var_new.size}')

    # get percent difference
    diff = np.zeros( var_old.size )
    max_val = max( np.max(np.abs(var_old)),np.max(np.abs(var_new)) )
    threshold = 1.0e-20*max_val
    for n in range( var_old.size ):
        if( max( abs(var_old[n]),abs(var_new[n]) ) > threshold ):
            diff[n] = abs(var_old[n] / var_new[n] - 1.0) * 100.0

    diff = np.abs(diff)
    diff_2d = np.reshape( diff,[156,121] )
    print(f"Shape diff_2d = {np.shape(diff_2d)}")
    print(f"Max diff @ (row,col) = {np.unravel_index(np.argmax(diff_2d),diff_2d.shape)}")
    print(f"Max diff @  {np.argmax( diff_2d )}")
    # get comparison stats
    results = {}
    results["minimum difference"] = np.amin( diff )
    results["maximum difference"] = np.amax( diff )
    results["maximum difference index"] = np.argmax( diff )
    results["maximum difference 2dindex"] = np.unravel_index(np.argmax(diff_2d),diff_2d.shape)
    results["RMS difference"]     = np.sqrt( np.mean( diff * diff ) )
    results["fail count"]         = np.count_nonzero( diff > tolerance )
    return results


# Returns paths to the script folder and the output folder
def get_paths():
    argc = len(sys.argv)
    if( argc == 4 ):
        script_path = str(sys.argv[1])
        old_output_path = str(sys.argv[2])
        new_output_path = str(sys.argv[3])
    else:
        print(f"Usage: python xsqy.compare.py path/to/scripts path/to/old/output path/to/new/output")
        sys.exit(2)
    return script_path, old_output_path, new_output_path


# Returns labels for new and old photolysis reactions
def get_labels(script_path):
    with open(os.path.join(script_path, "annotatedjlabels.old"), 'r') as f :
        label_old = f.readlines()
    with open(os.path.join(script_path, "annotatedjlabels.new"), 'r') as f :
        label_new = f.readlines()
    return label_new, label_old


# Compares new output file against old output file
def compare_output(fxsqy_new_path, fxsqy_old_path, labels_new, labels_old, config):
    fxsqy_old = FF( fxsqy_old_path, 'r' )
    fxsqy_new = FF( fxsqy_new_path, 'r' )
    xsqy_old   = fxsqy_old.read_reals(dtype=np.float64)
    xsqy_new   = fxsqy_new.read_reals(dtype=np.float64)

    nlabels_new = len( labels_new )
    nlabels_old = len( labels_old )

    print(f"\nxsqy.new type  = {xsqy_new.dtype}")
    ndata = int(xsqy_new.shape[0]/nlabels_new)
    print(f"xsqy.new size = {ndata}")

    print( f"\nThere are {nlabels_new} new arrays")
    print( f"There are {nlabels_old} old arrays\n")

    XSQY_new = np.reshape( xsqy_new,[nlabels_new,ndata] )
    print(f"\nXSQY_new type  = {XSQY_new.dtype}")
    print(f"XSQY_new shape = {XSQY_new.shape}")
    print( XSQY_new[0,:] )
    print("")

    # check reshape
    maxind = np.argmax( xsqy_new )
    print(f"\nMax val xsqy_new @ {maxind}")
    print(f" {maxind-4} <= n <= {maxind+4}")
    print("xsqy_new near Max")
    print( xsqy_new[maxind-4:maxind-1] )
    print( xsqy_new[maxind] )
    print( xsqy_new[maxind+1:maxind+4] )
    maxind = np.unravel_index( np.argmax(XSQY_new),XSQY_new.shape )
    print(f"\nMax val XSQY_new @ {maxind}")
    print("XSQY_new near Max")
    print( XSQY_new[maxind[0],maxind[1]-4:maxind[1]-1] )
    print( XSQY_new[maxind[0],maxind[1]] )
    print( XSQY_new[maxind[0],maxind[1]+1:maxind[1]+4] )
    print("")

    XSQY_old = np.reshape( xsqy_old,[nlabels_old,ndata] )
    print( XSQY_old.dtype )
    print( XSQY_old.shape )
    print( XSQY_old[0,:] )
    print("")

    success = True

    # make sure requested field is in both datasets
    for match, options in config.items():
        indatasets = False
        ndx_new = get_list_index( labels_new,match )
        indatasets = ndx_new > -1
        if( not indatasets ):
            print(f"\nNo match for {match} in new dataset")
            continue
        ndx_old = get_list_index( labels_old,match )
        indatasets = ndx_old > -1
        if( not indatasets ):
            print(f"\nNo match for {match} in old dataset")
            continue
        print(f"\n{match} in both datasets; (old,new) = {ndx_old},{ndx_new}")
        # compare datasets; old first
        print("old dataset")
        print(f"Min = {np.amin(XSQY_old[ndx_old,:])}")
        print(f"Max = {np.amax(XSQY_old[ndx_old,:])}")
        print(f"Non-zero count = {np.count_nonzero(XSQY_old[ndx_old,:])}")
        # new last
        print("\nnew dataset")
        print(f"Min = {np.amin(XSQY_new[ndx_new,:])}")
        print(f"Max = {np.amax(XSQY_new[ndx_new,:])}")
        print(f"Non-zero count = {np.count_nonzero(XSQY_new[ndx_new,:])}\n")
        results = compare_var( match, options["maximum difference"], XSQY_old[ndx_old,:], XSQY_new[ndx_new,:] )
        for metric, tolerance in options.items():
            if not metric in results.keys():
                print(f"Error: invalid comparison metric for {match}: {metric}")
                sys.exit(1)
            if tolerance < results[metric]:
                XSQY_old_2D = np.reshape(XSQY_old[ndx_old,:],[156,121])
                XSQY_new_2D = np.reshape(XSQY_new[ndx_new,:],[156,121])
                print(f"Error: comparison FAILURE for {match}:")
                print(f"shape XSQY_2D = {np.shape(XSQY_old_2D)}")
                print(f"{match} {metric} = {results[metric]} > {tolerance}")
                print(f"       Max diff = {results['maximum difference']}")
                print(f"   Max diff ndx = {results['maximum difference index']}")
                print(f"Max diff 2d ndx = {results['maximum difference 2dindex']}")
                maxind = results['maximum difference index']
                print(f"   old,new vals = {XSQY_old[ndx_old,maxind],XSQY_new[ndx_new,maxind]}")
                print(f"   Max diff @ (row,col) = {np.unravel_index(results['maximum difference index'],[121,156])}")
                print(f"       Fail cnt = {results['fail count']}\n")
                success = False
                continue
            else:
                print(f"{match} {metric} within tolerance: {results[metric]}% <= {tolerance}%")

    # close open files
    fxsqy_old.close()
    fxsqy_new.close()

    return success


# main
script_path, old_output_path, new_output_path = get_paths()
labels_new, labels_old = get_labels(script_path)

print("\njlabels.new\n-----------")
for label in labels_new:
    print(label.strip())

with open(os.path.join(script_path, f"xsqy.compare.json"),"r") as f :
    config=json.load(f)

file_indices = [ '01', '02', '03', '04', '05' ]
for file_index in file_indices:
    success = True
    success = success and compare_output(os.path.join(new_output_path, f"xsqy.{file_index}.new"), \
                                         os.path.join(old_output_path, f"xsqy.{file_index}.old"), \
                                         labels_new, labels_old, config)

if( not success ):
    raise RuntimeError('Test failed')
