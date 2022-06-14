#!/bin/bash

cd /bioseqdb \
    && rm -rf /build > /dev/null 2> /dev/null \
    && mkdir build \
    && cd build

while :
do
    echo "[i] Build the source code"
    cd /bioseqdb/build && cmake ..
    if [[ "$?" == "0" ]]; then
        make 2>&1 | tee /logs/docker_dev_build.log
    fi
    if [[ "$?" == "0" ]]; then
        make install 2>&1 | tee /logs/docker_dev_build.log
    fi
    echo "[i] Run postgres now"
    pkill postgres > /dev/null 2> /dev/null
    export POSTGRES_PASSWORD=password
    export POSTGRES_USER=bioseqdb
    export POSTGRES_DB=bioseqdb
    bash /usr/local/bin/docker-entrypoint.sh postgres
    echo "[i] Wait for any change in the source code..."
    fswatch -1 /bioseqdb/*
done


