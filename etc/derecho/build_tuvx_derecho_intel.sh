# Downloads and builds TUV-x and its dependencies on Derecho using Intel compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script

module purge
module load ncarenv/23.09
module load intel/2023.2.1
module load netcdf/4.9.2
module load mkl/2023.2.0
module load ncarcompilers/1.0.0
module load cmake/3.26.3

if [[ -z "${TUVX_HOME}" ]]; then
  echo "You must set the TUVX_HOME environment variable to the directory where TUV-x should be build."
  return
fi

if [[ ! -d "${TUVX_HOME}" ]]; then
  echo "TUVX_HOME must point to an existing directory"
  return
fi

echo "Building JSON Fortran"

# get & build the source code of JSON Fortran

cd ${TUVX_HOME}
curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.3.0.tar.gz
tar -zxvf 8.3.0.tar.gz
cd json-fortran-8.3.0
mkdir build
cd build
INSTALL_DIR=$TUVX_HOME/json-fortran-8.3.0
cmake -D SKIP_DOC_GEN:BOOL=TRUE -D CMAKE_INSTALL_PREFIX=$INSTALL_DIR ..
make install

echo "Building TUV-x"

# get & build the source code of TUV-x

cd ${TUVX_HOME}
git clone git@github.com:NCAR/tuv-x.git
cd tuv-x
mkdir build
cd build
export JSON_FORTRAN_HOME=$INSTALL_DIR/jsonfortran-intel-8.3.0
cmake -D CMAKE_BUILD_TYPE=release -D ENABLE_MEMCHECK=OFF -DBLAS_LIBRARIES=-lmkl_rt -DLAPACK_LIBRARIES=-lmkl_rt ..
make -j 8
