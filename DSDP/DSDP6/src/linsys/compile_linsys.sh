#!/bin/bash

rm -rf build
mkdir build
cd build
cmake ..
make
mv *.dylib /Users/gaowenzhi/Desktop/gwz/cmake/src/linsys/
