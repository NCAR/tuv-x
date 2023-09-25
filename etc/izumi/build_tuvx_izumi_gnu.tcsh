# Downloads and builds TUV-x and its dependencies on CASPER using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script


module purge
module load compiler/gnu/9.3.0
module load tool/netcdf/4.7.4/gnu/9.3.0

if ("${TUVX_HOME}" == "") then
  echo "You must set the TUVX_HOME environment variable to the directory where TUV-x should be build."
  exit
endif

if (! -d "${TUVX_HOME}" ) then
  echo "TUVX_HOME must point to an existing directory"
  exit
endif

echo "Building TUV-x"

# get the source code
cd ${TUVX_HOME}
curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.3.0.tar.gz
git clone --recurse-submodules https://github.com/NCAR/tuv-x.git

# extract
cd ${TUVX_HOME}
tar -zxf 8.3.0.tar.gz

set INSTALL_ROOT=$TUVX_HOME/install
mkdir -p $INSTALL_ROOT

# json-fortran
set JSON_FORTRAN_ROOT=$TUVX_HOME/json-fortran-8.3.0
set JSON_FORTRAN_HOME=$INSTALL_ROOT/jsonfortran-gnu-8.3.0
cd $JSON_FORTRAN_ROOT
sed -i 's/\-C $<CONFIG>//' CMakeLists.txt
sed -i -e '$a\ ' src/json_initialize_dummy_arguments.inc 
sed -i -e '$a\ ' src/json_initialize_arguments.inc 
sed -i -e '/#include "json_initialize_dummy_arguments.inc"/r src/json_initialize_dummy_arguments.inc' -e '/#include "json_initialize_dummy_arguments.inc"/d' src/json_value_module.F90
sed -i -e '/#include "json_initialize_dummy_arguments.inc"/r src/json_initialize_dummy_arguments.inc' -e '/#include "json_initialize_dummy_arguments.inc"/d' src/json_file_module.F90
sed -i -e '/#include "json_initialize_arguments.inc"/r src/json_initialize_arguments.inc' -e '/#include "json_initialize_arguments.inc"/d' src/json_value_module.F90
sed -i -e '/#include "json_initialize_arguments.inc"/r src/json_initialize_arguments.inc' -e '/#include "json_initialize_arguments.inc"/d' src/json_file_module.F90
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran \
      -D SKIP_DOC_GEN:BOOL=TRUE \
      -D CMAKE_INSTALL_PREFIX=$INSTALL_ROOT \
      ..
make install
mkdir -p $JSON_FORTRAN_HOME/lib/shared
mv $JSON_FORTRAN_HOME/lib/*.so* $JSON_FORTRAN_HOME/lib/shared

# TUV-x
set TUVX_ROOT=$TUVX_HOME/tuv-x
cd $TUVX_ROOT
git checkout release
git submodule update
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran \
      -D CMAKE_BUILD_TYPE=release \
      -D JSON_INCLUDE_DIR=$JSON_FORTRAN_HOME/lib \
      -D JSON_LIB=$JSON_FORTRAN_HOME/lib/libjsonfortran.a \
      -D NETCDF_INCLUDE_DIR=$NETCDF_PATH/include \
      -D NETCDF_C_LIB=$NETCDF_PATH/lib/libnetcdf.so \
      -D NETCDF_FORTRAN_LIB=$NETCDF_PATH/lib/libnetcdff.so \
      -D ENABLE_COVERAGE=OFF \
      ..
make
