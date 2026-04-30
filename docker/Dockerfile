FROM ubuntu:24.04

RUN apt-get update && apt-get install -y \
    cmake \
    g++ \
    git \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

COPY . /tuv-x/
WORKDIR /tuv-x/build

RUN cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DTUVX_ENABLE_TESTS=ON \
    && make -j$(nproc)
