FROM fedora:latest

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc \
        gcc-c++ \
        gcc-fortran \
        gdb \
        git \
        lapack-devel \
        lcov \
        make \
        netcdf-fortran-devel \
        python3 \
        python3-pip \
        valgrind \
        yaml-cpp-devel \
    && dnf clean all

RUN pip3 install numpy scipy netCDF4

ENV LD_LIBRARY_PATH=/usr/local/lib64

# build the tuv-x tool
COPY . /tuv-x
RUN cd /tuv-x \
      && mkdir build \
      && cd build \
      && cmake -D CMAKE_BUILD_TYPE=release \
               -D TUVX_INSTALL_INCLUDE_DIR=/usr/local/include \
               -D TUVX_INSTALL_MOD_DIR=/usr/local/include \
               .. \
      && make install -j 8

WORKDIR /tuv-x/build
