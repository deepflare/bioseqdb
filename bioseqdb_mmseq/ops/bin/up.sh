#!/bin/bash

PROJECT_DIR=$(dirname $(dirname "$0"))

cd "$PROJECT_DIR""/docker" && BUILDKIT_PROGRESS=plain docker-compose -p bioseqdb up --build --remove-orphans
