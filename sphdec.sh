#!/bin/sh

LD_LIBRARY_PATH=$(pwd)/src/fplll/lib/:$(pwd)/src/fplll/include/
export LD_LIBRARY_PATH
./sphdec "$@"
