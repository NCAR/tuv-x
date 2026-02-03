<h1 align="center">
<img src="docs/source/_static/logo.svg" width="300">
</h1><br>

Tropospheric ultraviolet-extended (TUV-x): A photolysis rate calculator

[![License](https://img.shields.io/github/license/NCAR/tuv-x.svg)](https://github.com/NCAR/tuv-x/blob/main/LICENSE)
[![Ubuntu](https://github.com/NCAR/tuv-x/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/ubuntu.yml)
[![Mac](https://github.com/NCAR/tuv-x/actions/workflows/mac.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/mac.yml)
[![Windows](https://github.com/NCAR/tuv-x/actions/workflows/windows.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/windows.yml)
[![Docker](https://github.com/NCAR/tuv-x/actions/workflows/docker.yml/badge.svg)](https://github.com/NCAR/tuv-x/actions/workflows/docker.yml)
[![codecov](https://codecov.io/gh/NCAR/tuv-x/branch/main/graph/badge.svg?token=H46AAEAQF9)](https://codecov.io/gh/NCAR/tuv-x)
[![DOI](https://zenodo.org/badge/396946468.svg)](https://zenodo.org/badge/latestdoi/396946468)

Copyright (C) 2020-2026 University Corporation for Atmospheric Research

# Building and installing
To build and install TUV-x locally, you must have the following libraries installed:

- [yaml-cpp](https://github.com/jbeder/yaml-cpp/)
- [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (both C and Fortran libraries)

You must also have CMake installed on your machine. 

To install TUV-x locally,
open a terminal window, navigate to a folder where you would like the TUV-x files to exist,
and run the following commands:

## Build and run (Docker version)

To build and run the stand-alone version of TUV-x, you must have [Docker Desktop](https://www.docker.com/get-started) installed and running. With Docker Desktop running, open a terminal window and run the following command to start the TUV-x container:

```
docker run -it ghcr.io/ncar/tuv-x:release bash
```

Inside the container, you can run the TUV-x tests from the `/build/` folder:

```
cd build/
# to run the tests
make test
# to use the standalone tool
./tuv-x examples/tuv_5_4.json
# or 
./tuv-x examples/ts1_tsmlt.json
```

### Sharing data between your computer and the docker container

To easily retrive the output data from tuv-x, or to use your own data in the configuration,
you will need to use [docker bind mounts](https://docs.docker.com/storage/bind-mounts/).

Start docker with a bind mount volume. **You must use a full system path.** 
This command will create a link between the folder you provided into and out of the running container.

For example, this command below will start tuv-x and map the directory `/output` in the running container
to your downloads directory.
```
docker run -v /Users/$USER/Downloads:/output -it tuvx
```

To run the full example config file with the tuv-x tool from this directory, you'll need to copy over the data files.
The full example uses data files from the `data` directory, which is why you need to do this.

```
cd /output
cp -r /build/data .
# to use the standalone tool
./tuv-x examples/tuv_5_4.json
# or 
./tuv-x examples/ts1_tsmlt.json
```

Now, in your downloads folder, you should have to nc files, `photolysis_rate_constants.nc` and `dose_rates.nc`.

If you have your own data files you'd like to use, you can copy them into the Downloads directory, or whichever directory you mapped in, and use them in the container.


## Build and run (local build version)

```
git clone https://github.com/NCAR/tuv-x.git
cd tuv-x
mkdir build
cd build
ccmake ..
make -j 8
```

You will now have a runnable exectubable for `tuv-x` and the tests in the build directory.

```
# to use the standalone tool
./tuv-x examples/tuv_5_4.json
# or 
./tuv-x examples/ts1_tsmlt.json
```
Inspect the output file `photolysis_rate_constants.nc` to see the results!

## Install

After completing the previous step run `sudo make install`.
This wll install the tuvx static library, the tuv-x configuration
and runtime data, as well as the standalone `tuv-x` exectuable, which can be
added to your system path to make the executable useable from any directory.

If you would later lake to uninstall tuv-x, you can run
`sudo make uninstall` from the build directory. 


# Citation

The following bibtex can be used to cite the work which originally developed
this tool.

A recommended citation is 

> Madronich, Sasha, and Siri Flocke (1999), The role of solar radiation in atmospheric chemistry, in Handbook of Environmental Chemistry, edited by P. Boule, pp. 1-26, Springer-Verlag, Heidelberg.

However, you are encouraged to use the format matching whichever style
you prefer.

```
@incollection{madronich_role_1999,
	address = {Berlin, Heidelberg},
	series = {The {Handbook} of {Environmental} {Chemistry}},
	title = {The {Role} of {Solar} {Radiation} in {Atmospheric} {Chemistry}},
	isbn = {978-3-540-69044-3},
	url = {https://doi.org/10.1007/978-3-540-69044-3_1},
	language = {en},
	booktitle = {Environmental {Photochemistry}},
	publisher = {Springer},
	author = {Madronich, Sasha and Flocke, Siri},
	editor = {Boule, Pierre},
	year = {1999},
	doi = {10.1007/978-3-540-69044-3_1},
	keywords = {Earth-Sun geometry., photolysis rate coefficients, radiative transfer, solar radiation, spectral actinic flux},
	pages = {1--26},
}
```

# Community and contributions
We welcome contributions and feedback from anyone, everything from updating
the content or appearance of the documentation to new and
cutting edge science.

- [Contact](https://github.com/NCAR/tuv-x/discussions)
  - If you'd like to get in touch with us, feel free to start a conversation
on our [discussion board](https://github.com/NCAR/tuv-x/discussions) 
or email us at musica-info@ucar.edu. 

- [Collaboration](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf)
  - Anyone interested in scientific collaboration
which would add new software functionality should read the [MUSICA software development plan](https://github.com/NCAR/musica/blob/main/docs/Software%20Development%20Plan.pdf).

- [Contributor's guide](https://ncar.github.io/tuv-x/contributing/contributors_guide.html)
  - Before submiitting a PR, please thouroughly read this to you understand our expectations. We reserve the right to reject any PR not meeting our guidelines.


# Documentation
Please see the [TUV-x documentation](https://ncar.github.io/tuv-x/) for detailed
installation and usage instructions.

# License

- [Apache 2.0](/LICENSE)
- Copyright (C) 2020-2026 University Corporation for Atmospheric Research
