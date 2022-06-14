#pragma once

extern "C" {
#include <postgres.h>
#include <fmgr.h>
#include <funcapi.h>
#include <miscadmin.h>
#include <executor/spi.h>
#include <catalog/pg_type.h>
}

#include <algorithm>
#include <array>
#include <chrono>
#include <charconv>
#include <string>
#include <string_view>
#include <optional>
#include <stdint.h>
#include <cstdlib>

#define raise_pg_error(code, msg) ereport(ERROR, (errcode(code)), msg);

#define POSTGRES_FUNCTION(name) \
    extern "C" { \
    PG_FUNCTION_INFO_V1(name); \
    Datum name(PG_FUNCTION_ARGS)

#define POSTGRES_FUNCTION_END() \
}