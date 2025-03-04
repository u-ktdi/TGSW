#!/bin/bash

mkdir -p build
cd build

cmake ..
make

./test_tgsw