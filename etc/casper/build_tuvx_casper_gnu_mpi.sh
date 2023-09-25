# Downloads and builds TUV-x and its dependencies on CASPER using GNU compilers
#
# The TUVX_HOME environment variable must be set to the directory to build TUV-x
# in prior to calling this script


module purge
module load gnu/11.2.0
module load ncarenv/1.3
module load ncarcompilers/0.5.0
module load cmake/3.22.0
module load netcdf/4.8.1
module load openmpi/4.1.1

if [[ -z "${TUVX_HOME}" ]]; then
  echo "You must set the TUVX_HOME environment variable to the directory where TUV-x should be build."
  return
fi

if [[ ! -d "${TUVX_HOME}" ]]; then
  echo "TUVX_HOME must point to an existing directory"
  return
fi

echo "Building TUV-x"

# get the source code
cd ${TUVX_HOME}
curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.1.tar.gz
git clone https://github.com/NCAR/musica-core.git
git clone --recurse-submodules https://github.com/NCAR/tuv-x.git

# extract
cd ${TUVX_HOME}
tar -zxf 8.2.1.tar.gz

INSTALL_ROOT=$TUVX_HOME/install
mkdir -p $INSTALL_ROOT

# json-fortran
JSON_FORTRAN_ROOT=$TUVX_HOME/json-fortran-8.2.1
export JSON_FORTRAN_HOME=$INSTALL_ROOT/jsonfortran-gnu-8.2.1
cd $JSON_FORTRAN_ROOT
sed -i 's/\-C $<CONFIG>//' CMakeLists.txt
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=gfortran \
      -D SKIP_DOC_GEN:BOOL=TRUE \
      -D CMAKE_INSTALL_PREFIX=$INSTALL_ROOT \
      ..
make install
mkdir -p $JSON_FORTRAN_HOME/lib/shared
mv $JSON_FORTRAN_HOME/lib/*.so* $JSON_FORTRAN_HOME/lib/shared

# musica-core
MUSICA_CORE_ROOT=${TUVX_HOME}/musica-core
export MUSICA_CORE_HOME=${INSTALL_ROOT}/musica-core-0.1.0
export MUSICA_CORE_PACKAGE=${INSTALL_ROOT}/musicacore-0.1.0/cmake/musicacore-0.1.0
cd ${MUSICA_CORE_ROOT}
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=mpifort \
      -D CMAKE_BUILD_TYPE=release \
      -D CMAKE_INSTALL_PREFIX=${INSTALL_ROOT} \
      -D ENABLE_MPI=ON \
      ..
make install

# TUV-x
TUVX_ROOT=$TUVX_HOME/tuv-x
cd $TUVX_ROOT
git checkout develop-output-cs-qy
git submodule update
mkdir -p build
cd build
cmake -D CMAKE_Fortran_COMPILER=mpifort \
      -D CMAKE_BUILD_TYPE=release \
      -D musicacore_DIR=${MUSICA_CORE_PACKAGE} \
      -D ENABLE_MPI=ON \
      -D ENABLE_MEMCHECK=OFF \
      -D ENABLE_COVERAGE=OFF \
      -D CMAKE_INSTALL_PREFIX=${TUVX_HOME} \
      ..
make install
