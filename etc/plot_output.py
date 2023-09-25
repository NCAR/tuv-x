#!/usr/bin/env python3

from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse


class NetCdfData:
    """NetCDF dataset to plot"""
    def __init__(self, file_path):
        self.dataset = Dataset(f"{file_path}", mode='r')


def get_args():
    """Gets all TUV-x plotter command-line options"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, type=str, help="NetCDF file containing data to plot")
    parser.add_argument("-p", "--photo-rate", nargs='+', type=str, help="Plot photolysis rate [s-1]")
    parser.add_argument("-sza", "--solar-zenith-angle", type=float, help="Solar zenith angle [degrees]. " \
                        "If not specified, the first solar zenith angle in the file will be used.")
    parser.add_argument("-xs", "--cross-section", nargs='+', type=str, help="Plot cross section(s)")
    parser.add_argument("-rf", "--radiation-field", action="store_true", help="Plot the radiation field")
    parser.add_argument("-qy", "--quantum-yield", nargs='+', type=str, help="Plot quantum yield(s)")
    parser.add_argument("-ht", "--height", nargs='+', type=float, help="Plot specific vertical layers [km]")
    parser.add_argument("-wl", "--wavelength", nargs='+', type=float, help="Plot specific wavelengths [nm]")
    parser.add_argument("-t", "--temperature", nargs='+', type=float, help="Plot specific temperatures [K]")
    parser.add_argument("-lht", "--list-heights", action="store_true", help="List all heights [km]")
    parser.add_argument("-lwl", "--list-wavelengths", action="store_true", help="List all wavelengths [nm]")
    parser.add_argument("-lt", "--list-temperatures", action="store_true", help="List all temperatures [K]")
    parser.add_argument("-lsza", "--list-solar-zenith-angles", action="store_true", help="List all solar " \
                        "zenith angles [degrees]")
    parser.add_argument("-log", "--log", action="store_true", help="Plot with logarithmic axis")
    parser.add_argument("-wllb", "--wavelength-lower-bound", type=float, help="Set lower wavelength bound")
    parser.add_argument("-wlub", "--wavelength-upper-bound", type=float, help="Set upper wavelength bound")
    parser.add_argument("-htlb", "--height-lower-bound", type=float, help="Set lower height bound")
    parser.add_argument("-htub", "--height-upper-bound", type=float, help="Set upper height bound")
    args = parser.parse_args( )
    if args.cross_section and args.quantum_yield:
        parser.error("--cross-section and --quantum-yield cannot be used in combination")
    if args.photo_rate and args.height:
        parser.error("--photo-rate cannot be used with --height specified")
    if args.photo_rate and args.wavelength:
        parser.error("--photo-rate cannot be used with --wavelength specified")
    return args


def find_nearest_indices(array, values):
    results = np.empty((len(values),), dtype=int)
    for i, value in enumerate(values):
        results[i] = (np.abs(array - value)).argmin()
    return results


def plot_by_wavelength(args, title, labels, plot_data, wavelengths, heights, height_ids, \
                       height_units, szas, sza_id):
    """Plots a set of profiles by wavelength"""
    plt.title(f"{title}\nSZA = {szas[sza_id]}")
    plt.ylabel(title)
    plt.xlabel('wavelength [nm]')
    if args.wavelength_lower_bound:
        l_wl = find_nearest_indices(wavelengths[:], [args.wavelength_lower_bound])[0]
    else:
        l_wl = 0
    if args.wavelength_upper_bound:
        u_wl = find_nearest_indices(wavelengths[:], [args.wavelength_upper_bound])[0]
    else:
        u_wl = len(wavelengths[:])-1
    for label in labels:
        for height_id in height_ids:
            trace_label = f"{label} {heights[height_id]} {height_units}"
            plt.plot(wavelengths[l_wl:u_wl], plot_data[label][l_wl:u_wl,height_id,sza_id], \
                     label=trace_label)
    plt.legend(loc = 'upper right')
    if args.log:
        plt.yscale('log')
    plt.show()


def plot_by_height(args, title, labels, plot_data, wavelengths, heights, wavelength_ids, \
                   wavelength_units, szas, sza_id):
    """Plots a set of profiles by height"""
    plt.title(f"{title}\nSZA = {szas[sza_id]}")
    plt.ylabel('height [km]')
    plt.xlabel(title)
    if args.height_lower_bound:
        l_ht = find_nearest_indices(heights[:], [args.height_lower_bound])[0]
    else:
        l_ht = 0
    if args.height_upper_bound:
        u_ht = find_nearest_indices(heights[:], [args.height_upper_bound])[0]
    else:
        u_ht = len(heights[:])-1
    for label in labels:
        if len(wavelengths[:]) == 0:
            plt.plot(plot_data[label][l_ht:u_ht,1], heights[l_ht:u_ht], label=label)
        for wavelength_id in wavelength_ids:
            trace_label = f"{label} {wavelengths[wavelength_id]} {wavelength_units}"
            plt.plot(plot_data[label][wavelength_id,l_ht:u_ht,sza_id], heights[l_ht:u_ht], label=trace_label)
    plt.legend(loc = 'lower right')
    if args.log:
        plt.xscale('log')
    plt.show()


"""
TUV-x output plotter
"""

args = get_args()
tuvx_data = NetCdfData(args.file)
heights = tuvx_data.dataset.variables["altitude"]
wavelengths = tuvx_data.dataset.variables["wavelength"]
szas = tuvx_data.dataset.variables["solar zenith angle"]
temperatures = tuvx_data.dataset.variables["temperature"]
sza_id = 0

if args.list_heights:
    print(f"\nheights [km]: {heights[:]}")

if args.list_wavelengths:
    print(f"\nwavelengths [nm]: {wavelengths[:]}")

if args.list_solar_zenith_angles:
    print(f"\nsolar zenith angles [degrees]: {szas[:]}")

if args.list_temperatures:
    print(f"\ntemperatures [K]: {temperatures[:]}")

if args.cross_section:
    title = "cross section [cm2 molecule-1]"
    labels = args.cross_section
    plot_data = {}
    for label in labels:
        plot_data[label] = tuvx_data.dataset.variables[f"cross section {label}"]

if args.quantum_yield:
    title = "quantum yield [unitless]"
    labels = args.quantum_yield
    plot_data = {}
    for label in labels:
        plot_data[label] = tuvx_data.dataset.variables[f"quantum yield {label}"]

if args.radiation_field:
    title = "radation field [GET_UNITS]"
    labels = { "direct radiation", "upward radiation", "downward radiation" }
    plot_data = {}
    for label in labels:
        plot_data[label] = tuvx_data.dataset.variables[label]

if args.solar_zenith_angle:
    sza_id = find_nearest_indices(szas[:], [args.solar_zenith_angle])

if args.height:
    height_ids = find_nearest_indices(heights[:], args.height)
    plot_by_wavelength(args, title, labels, plot_data, wavelengths, heights, height_ids, \
                       "km", szas, sza_id)

if args.temperature:
    temperature_ids = find_nearest_indices(temperatures[:,sza_id], args.temperature)
    plot_by_wavelength(args, title, labels, plot_data, wavelengths, temperatures[:,sza_id], \
                       temperature_ids, "K", szas, sza_id)

if args.wavelength:
    wavelength_ids = find_nearest_indices(wavelengths[:], args.wavelength)
    plot_by_height(args, title, labels, plot_data, wavelengths, heights, wavelength_ids, \
                   "nm", szas, sza_id)

if args.photo_rate:
    title = "photolysis rate constant [s-1]"
    labels = args.photo_rate
    plot_data = {}
    for label in labels:
        plot_data[label] = tuvx_data.dataset.variables[label]
    plot_by_height(args, title, labels, plot_data, [], heights, [], "km", szas, sza_id)
