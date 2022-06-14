#!/bin/bash

docker build . -t bioseqdb-dev -f ops/docker/Dockerfile.dev.pg \
&& docker run -it \
    -v $(pwd):/source \
    -v $(pwd)/_logs:/logs \
    bioseqdb-dev