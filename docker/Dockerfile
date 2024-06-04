FROM fedora:37

RUN dnf -y update \
    && dnf -y install \
        gcc-fortran \
        gcc-c++ \
        gcc \
        gdb \
        git \
        netcdf-fortran-devel \
        cmake \
        make \
        lcov \
        valgrind \
        python3 \
        python3-pip \
        lapack-devel \
        yaml-cpp-devel \
    && dnf clean all

RUN pip3 install numpy scipy

ENV LD_LIBRARY_PATH=/usr/local/lib64

# build the tuv-x tool
COPY . /tuv-x/
RUN mkdir /build \
      && cd /build \
      && cmake -D CMAKE_BUILD_TYPE=release \
               -D TUVX_INSTALL_INCLUDE_DIR=/usr/local/include \
               -D TUVX_INSTALL_MOD_DIR=/usr/local/include \
               /tuv-x \
      && make install -j 8

WORKDIR /build
