#!/bin/bash


echo "Generating sources for extension ${POSTGRES_EXTENSION_NAME} (version: ${POSTGRES_EXTENSION_VERSION})"

mkdir -p gen > /dev/null 2> /dev/null \
&& poetry run yasha \
    -o gen/extension.cpp \
    extension/extension.cpp.jinja \
&& poetry run yasha \
    -o gen/${POSTGRES_EXTENSION_NAME}--${POSTGRES_EXTENSION_VERSION}.sql \
    extension/extension.sql.jinja \
&& poetry run yasha \
    -o gen/${POSTGRES_EXTENSION_NAME}.control \
    extension/extension.control.jinja \
&& poetry run yasha \
    -o gen/config.h \
    extension/config.h.jinja
