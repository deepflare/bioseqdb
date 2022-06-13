FROM postgres:14-bullseye

LABEL maintainer="bioseqdb - https://github.com/covid-genomics/bioseqdb"


RUN set -ex \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
      autoconf \
      automake \
      autotools-dev \
	  build-essential \
      ca-certificates \
      cmake \
      g++ \
      git \
      make \
      postgresql-server-dev-14 \
	  libbz2-dev \
	  zlib1g-dev \
	  liblzma-dev

ENV LIB_FLAGS "-g -Wall -O2 -fPIC -march=native -DSIMDE_ENABLE_NATIVE_ALIASES"
ENV BWA_TAG master
ENV HTS_TAG 1.5

ADD . /bioseqdb

WORKDIR /bioseqdb

RUN cd hello_world && mkdir build && cd build && cmake .. && make