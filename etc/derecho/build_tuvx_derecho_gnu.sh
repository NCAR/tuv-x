# Downloads and builds TUV-x and its dependencies on Derecho using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script

module purge
module load ncarenv/23.09
module load craype/2.7.20
module load gcc/12.2.0
module load cray-libsci/23.02.1.1
module load netcdf/4.9.2
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

# download and build TUV-X
echo "Downloading and Building TUV-x"
cd ${TUVX_HOME}
git clone git@github.com:NCAR/tuv-x.git
cd tuv-x
mkdir build
cd build
cmake -D CMAKE_BUILD_TYPE=release -D TUVX_ENABLE_MEMCHECK=OFF -D LAPACK_LIBRARIES=-lsci_gnu -D BLAS_LIBRARIES=-lsci_gnu ..   
make -j 8
