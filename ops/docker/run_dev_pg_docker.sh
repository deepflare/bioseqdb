#!/bin/bash

docker build . -t bioseqdb-dev -f ops/docker/Dockerfile.dev.pg \
&& docker run -it \
    -v $(pwd):/bioseqdb \
    -v $(pwd)/_logs:/logs \
    bioseqdb-dev