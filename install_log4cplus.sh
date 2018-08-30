#!/bin/sh
set -ex
mkdir log4cplus && cd log4cplus
git clone --recursive https://github.com/log4cplus/log4cplus.git --branch REL_2_0_2 .
mkdir build && cd build && cmake .. && make && sudo make install
