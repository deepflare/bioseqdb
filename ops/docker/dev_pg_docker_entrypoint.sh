#!/bin/bash

cp -Rf /source/. /bioseqdb

cd /bioseqdb \
&& rm -rf /build > /dev/null 2> /dev/null \
&& mkdir build \
&& cd build

cmake ..
make 2>&1 | tee /logs/docker_dev_build.log


