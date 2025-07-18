<h1 align="center">
<img src="../../docs/source/_static/logo.svg" width="300">
</h1><br>

# Python driver for TUV-x

Copyright (C) 2024-2025 University Corporation for Atmospheric Research

TUV-x is a Fortran application that calculates photolysis rates.
ACOM has developed a [Quick TUV-x Calculator](https://www.acom.ucar.edu/Models/TUV/Interactive_TUV/tuv-x.shtml) that accepts user input and runs a calculation.
The web interface for the TUV-x calculator runs on ACOM machine web3, while the Fortran executable runs on machine modeling2.
A python script **tuv-x-remote.py** collects user input from the calculator form, formats the configuration into a json file,
stores the cloud and aerosol radiators as NetCDF, and then launches the Fortran executable on modeling2.

When the calculation is complete, tuv-x-remote.py retrieves the results from a NetCDF file.
The script then writes the results into a set of CSV, text, and NetCDF files, and packages them for user download in temporary directory.

Under this software architecture, the python script tuv-x-remote.py can also be used as a stand-alone driver for the Fortran engine.
The batch script **goTuv-x-remote** provides several examples of calling tuv-x-remote.py from the Unix command line.

```
python3 tuv-x-remote.py latitude=40.0 longitude=-105.0 date=20240320 time=19:00:00 photolysis=False actinic=False crossSections=False quantumYields=False customWavelength=False waveStart=280.0 waveEnd=700.0 waveIntervals=420.0 inputOption=1 zenithAngle=0.0 ozone=300.0 albedo=0.1 gElevation=2.0 height=4.888 customProfile=False temperature=3.14159 density=2.71828 cloudOptDepth=0.25 cloudSingScatAlbedo=0.9999 cloudAsymmetry=0.85 cloudTop=5.555 cloudBase=4.444 aerosolOptDepth=0.235 aerosolSingScatAlbedo=0.99 aerosolAsymmetry=0.62 aerosolAlpha=1.0 aerosolLambda0=550.0 directSun=1.0 diffDown=1.0 diffUp=1.0 numStreams=2
```

# License

- [Apache 2.0](/LICENSE)
- Copyright (C) 2024-2025 University Corporation for Atmospheric Research

