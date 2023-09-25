# https://www.gitpod.io/docs/configure/workspaces/workspace-image#use-a-custom-dockerfile
# configured following some instructions from above address

FROM fedora:35

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
    && dnf clean all

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && sudo make install -j 4

ENV JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0"
ENV LD_LIBRARY_PATH=/usr/local/jsonfortran-gnu-8.2.0/lib/

RUN git clone https://github.com/NCAR/tuv-x \
    && cd tuv-x \
    && mkdir build \
    && cd build \
    && cmake -D CMAKE_BUILD_TYPE=release -D ENABLE_MEMCHECK=OFF -D ENABLE_TESTS=OFF .. \
    && make -j 4 \
    && sudo make install

# Create the gitpod user. UID must be 33333.
# https://github.com/gitpod-io/gitpod/issues/7071
RUN useradd -l -u 33333 -G wheel -md /home/gitpod -s /bin/bash -p gitpod gitpod

USER gitpod
