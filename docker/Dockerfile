FROM fedora:latest

RUN dnf -y update \
    && dnf -y install \
        cmake \
        gcc \
        gcc-c++ \
        gdb \
        git \
        lapack-devel \
        lcov \
        make \
        valgrind \
    && dnf clean all

# build the tuv-x library
COPY . /tuv-x
RUN cd /tuv-x \
      && mkdir build \
      && cd build \
      && cmake -D CMAKE_BUILD_TYPE=release \
               .. \
      && make -j 8

WORKDIR /tuv-x/build
