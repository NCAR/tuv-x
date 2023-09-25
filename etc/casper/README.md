# Building TUV-x on CASPER

## Get the source code

- Copy the build script you wish to use from this folder to GLADE.
  There are scripts for GNU and Intel compilers, with and without MPI support.

- Log in to CASPER and start an interactive session (BASH shell)

- Create a directory to build TUV-x in:

```
mkdir my-tuvx-build
```

- Create an environment variable named `TUVX_HOME` pointing to the absolute path of your build directory:

```
export TUVX_HOME=/path/to/my-tuvx-build
```

## Build TUV-x

Replace `/path/to/build_tuvx_casper_X.sh` with the path to the build script you copied to GLADE, in the following:

```
cd $TUVX_HOME
. /path/to/build_tuvx_casper_X.sh
```

## Run TUV-x
- Whenever you go to run TUV-x after it has been built, make sure you have the correct environment modules loaded.
  Look at the top of the build script you used to find the modules required.

- Run a test. The tests use the python numpy package, and we access this through the conda environment module on CASPER:

```
module load conda/latest
conda activate npl
cd $TUVX_HOME/tuv-x/build
make test
```

- Run a configuration (this does not require any Python packages)

```
cd $TUVX_HOME/tuv-x/build
./tuv-x examples/full_config.json
```

