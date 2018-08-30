#!/bin/sh
set -ex
wget https://github.com/log4cplus/log4cplus/archive/REL_2_0_1.zip
unzip REL_2_0_1.zip && cd log4cplus-REL_2_0_1
mkdir build && cd build && cmake .. && make && sudo make install
