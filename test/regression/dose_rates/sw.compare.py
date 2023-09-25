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
    diff_2d = np.reshape( diff,[156,1] )
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
        print(f"Usage: python sw.compare.py path/to/scripts path/to/old/output path/to/new/output")
        sys.exit(2)
    return script_path, old_output_path, new_output_path


# Returns labels for new and old photolysis reactions
def get_labels(script_path):
    with open(os.path.join(script_path, "annotatedslabels.old"), 'r') as f :
        label_old = f.readlines()
    with open(os.path.join(script_path, "annotatedslabels.new"), 'r') as f :
        label_new = f.readlines()
    return label_new, label_old


# Compares new output file against old output file
def compare_output(fsw_new_path, fsw_old_path, labels_new, labels_old, config):
    fsw_old = FF( fsw_old_path, 'r' )
    fsw_new = FF( fsw_new_path, 'r' )
    sw_old   = fsw_old.read_reals(dtype=np.float64)
    sw_new   = fsw_new.read_reals(dtype=np.float64)

    nlabels_new = len( labels_new )
    nlabels_old = len( labels_old )

    print(f"\nsw.new type  = {sw_new.dtype}")
    ndata = int(sw_new.shape[0]/nlabels_new)
    print(f"sw.new size = {ndata}")

    print( f"\nThere are {nlabels_new} new arrays")
    print( f"There are {nlabels_old} old arrays\n")

    SW_new = np.reshape( sw_new,[nlabels_new,ndata] )
    print(f"\nSW_new type  = {SW_new.dtype}")
    print(f"SW_new shape = {SW_new.shape}")
    print( SW_new[0,:] )
    print("")

    # check reshape
    maxind = np.argmax( sw_new )
    print(f"\nMax val sw_new @ {maxind}")
    print(f" {maxind-4} <= n <= {maxind+4}")
    print("sw_new near Max")
    print( sw_new[maxind-4:maxind-1] )
    print( sw_new[maxind] )
    print( sw_new[maxind+1:maxind+4] )
    maxind = np.unravel_index( np.argmax(SW_new),SW_new.shape )
    print(f"\nMax val SW_new @ {maxind}")
    print("SW_new near Max")
    print( SW_new[maxind[0],maxind[1]-4:maxind[1]-1] )
    print( SW_new[maxind[0],maxind[1]] )
    print( SW_new[maxind[0],maxind[1]+1:maxind[1]+4] )
    print("")

    SW_old = np.reshape( sw_old,[nlabels_old,ndata] )
    print( SW_old.dtype )
    print( SW_old.shape )
    print( SW_old[0,:] )
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
        print(f"Min = {np.amin(SW_old[ndx_old,:])}")
        print(f"Max = {np.amax(SW_old[ndx_old,:])}")
        print(f"Non-zero count = {np.count_nonzero(SW_old[ndx_old,:])}")
        # new last
        print("\nnew dataset")
        print(f"Min = {np.amin(SW_new[ndx_new,:])}")
        print(f"Max = {np.amax(SW_new[ndx_new,:])}")
        print(f"Non-zero count = {np.count_nonzero(SW_new[ndx_new,:])}\n")
        results = compare_var( match, options["maximum difference"], SW_old[ndx_old,:], SW_new[ndx_new,:] )
        for metric, tolerance in options.items():
            if not metric in results.keys():
                print(f"Error: invalid comparison metric for {match}: {metric}")
                sys.exit(1)
            if tolerance < results[metric]:
                SW_old_2D = np.reshape(SW_old[ndx_old,:],[156,1])
                SW_new_2D = np.reshape(SW_new[ndx_new,:],[156,1])
                print(f"Error: comparison FAILURE for {match}:")
                print(f"shape XSQY_2D = {np.shape(SW_old_2D)}")
                print(f"{match} {metric} = {results[metric]} > {tolerance}")
                print(f"       Max diff = {results['maximum difference']}")
                print(f"   Max diff ndx = {results['maximum difference index']}")
                print(f"Max diff 2d ndx = {results['maximum difference 2dindex']}")
                maxind = results['maximum difference index']
                print(f"   old,new vals = {SW_old[ndx_old,maxind],SW_new[ndx_new,maxind]}")
                print(f"   Max diff @ (row,col) = {np.unravel_index(results['maximum difference index'],[1,156])}")
                print(f"       Fail cnt = {results['fail count']}\n")
                success = False
                continue
            else:
                print(f"{match} {metric} within tolerance: {results[metric]}% <= {tolerance}%")

    # close open files
    fsw_old.close()
    fsw_new.close()

    return success


# main
script_path, old_output_path, new_output_path  = get_paths()
labels_new, labels_old = get_labels(script_path)

print("\nslabels.new\n-----------")
for label in labels_new:
    print(label.strip())

with open(os.path.join(script_path, f"sw.compare.json"),"r") as f :
    config=json.load(f)

file_indices = [ '01', '02', '03', '04', '05' ]
for file_index in file_indices:
    success = True
    success = success and compare_output(os.path.join(new_output_path, f"sw.{file_index}.new"), \
                                         os.path.join(old_output_path, f"sw.{file_index}.old"), \
                                         labels_new, labels_old, config)

if( not success ):
    raise RuntimeError('Test failed')
