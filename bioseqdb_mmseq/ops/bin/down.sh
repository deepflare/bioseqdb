#!/bin/bash

PROJECT_DIR=$(dirname $(dirname "$0"))

cd "$PROJECT_DIR""/docker" && docker-compose -p bioseqdb down --remove-orphans
