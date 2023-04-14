FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    g++-9 \
    wget \
    make \
    libssl-dev

ENV CC=/bin/gcc-9
ENV CXX=/bin/g++-9

RUN ln -s /bin/gcc-9 /bin/gcc
RUN ln -s /bin/g++-9 /bin/g++

# cmake
RUN wget https://cmake.org/files/v3.20/cmake-3.20.0.tar.gz \
   && tar -xzvf cmake-3.20.0.tar.gz \
   && cd cmake-3.20.0 \
   && ./bootstrap \
   && make -j16 \
   && make install \
   && cd ../ \
   && rm -rf cmake-3.20.0

# boost
RUN wget https://sourceforge.net/projects/boost/files/boost/1.73.0/boost_1_73_0.tar.gz \
  && tar -xzvf boost_1_73_0.tar.gz \
  && cd boost_1_73_0 \
  && ./bootstrap.sh --prefix=/usr/local/ --with-libraries=program_options \
  && ./b2 install \
  && cd ../ \
  && rm -rf boost_1_73_0


RUN $CC  --version
RUN $CXX --version

