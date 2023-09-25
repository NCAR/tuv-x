# Building TUV-x on modeling2

1. Log in to modeling2 (BASH shell)
1. Create a directory where TUV-x will be built.
    `mkdir tuvx-build`
1. Set the `TUVX_HOME` shell variable to the location of the build directory. 
    For example, if you made the above directory in your user's home `export TUVX_HOME=/home/$USER/tuvx-build`
1. Go into this directory `cd $TUVX_HOME`
1. Download the build script from ['here'](https://github.com/NCAR/tuv-x/blob/main/etc/modeling2/build_tuvx_modeling2_gnu.sh). 
  `curl -o build_tuvx_modeling2_gnu.sh https://raw.githubusercontent.com/NCAR/tuv-x/main/etc/modeling2/build_tuvx_modeling2_gnu.sh`
1. Make the script executable `chmod +x build_tuvx_modeling2_gnu.sh`
1. Set library paths
```
export PKG_CONFIG_PATH="/opt/local/lib/pkgconfig:$PKG_CONFIG_PATH"
export LD_LIBRARY_PATH="/opt/local/lib64:/opt/local/lib:$LD_LIBRARY_PATH"
```
8. Execute the script to build TUV-x. `./build_tuvx_modeling2_gnu.sh`
1. Run the tests
```
cd $TUVX_HOME/tuv-x/build
make test ARGS=-j4
```
10. Run a configuration
```
cd $TUVX_HOME/tuv-x/build
./tuv-x examples/full_config.json
```

