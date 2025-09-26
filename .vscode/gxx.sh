#!/usr/bin/env bash
# Make GCC 14’s runtime visible to cpptools
export LD_LIBRARY_PATH="/opt/bwhpc/common/compiler/gnu/14.1.0/lib64:${LD_LIBRARY_PATH}"
exec /gpfs/bwfor/software/common/compiler/gnu/14.1.0/bin/g++ "$@"