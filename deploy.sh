#!/usr/bin/env bash

source bootstrap.sh
mkdir build
cd build
../configure
make
