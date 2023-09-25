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
        lapack-devel \
    && dnf clean all

RUN pip3 install numpy scipy

# install json-fortran
RUN curl -LO https://github.com/jacobwilliams/json-fortran/archive/8.2.0.tar.gz \
    && tar -zxvf 8.2.0.tar.gz \
    && cd json-fortran-8.2.0 \
    && export FC=gfortran \
    && mkdir build \
    && cd build \
    && cmake -D SKIP_DOC_GEN:BOOL=TRUE .. \
    && sudo make install

# add a symlink
# Create symlinks in the Docker container
RUN ln -s /usr/local/jsonfortran-gnu-8.2.0/lib/libjsonfortran.a /usr/local/lib64/libjsonfortran.a && \
    ln -s /usr/local/jsonfortran-gnu-8.2.0/lib/libjsonfortran.so.8.2 /usr/local/lib64/libjsonfortran.so.8.2 && \
    ln -s /usr/local/jsonfortran-gnu-8.2.0/lib/libjsonfortran.so /usr/local/lib64/libjsonfortran.so && \
    ln -s /usr/local/jsonfortran-gnu-8.2.0/lib/libjsonfortran.so.8.2.0 /usr/local/lib64/libjsonfortran.so.8.2.0

ENV LD_LIBRARY_PATH=/usr/local/lib64

# build the tuv-x tool
COPY . /tuv-x/
RUN mkdir /build \
      && cd /build \
      && export JSON_FORTRAN_HOME="/usr/local/jsonfortran-gnu-8.2.0" \
      && cmake -D CMAKE_BUILD_TYPE=release \
               -D ENABLE_MEMCHECK=OFF \
               /tuv-x \
      && make install -j 8

WORKDIR /build
