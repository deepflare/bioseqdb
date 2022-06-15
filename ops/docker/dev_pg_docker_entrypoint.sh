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
    su - postgres -c 'cd /bioseqdb/packages/db_tests && poetry install && poetry run pytest'
    fswatch -1 /bioseqdb/*
done


