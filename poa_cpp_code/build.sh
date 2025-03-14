#!/bin/bash

rm -rf build
mkdir build
/usr/bin/cmake -DCMAKE_BUILD_TYPE:STRING=Release -S./ -B./build -G "Unix Makefiles"
/usr/bin/cmake --build build/ --config Release --target all -j40 --
/usr/bin/cmake --build build/ --config Release --target all --
