#!/bin/bash
# see https://vaneyckt.io/posts/safer_bash_scripts_with_set_euxo_pipefail/
set -euo pipefail

if [ "$TRAVIS_OS_NAME" == "linux" ]; then
  OS=Linux
elif [ "$TRAVIS_OS_NAME" == "macos" ]; then
  OS=Darwin
fi

URL=https://github.com/Kitware/CMake/releases/download/v{$CMAKE_VERSION}/cmake-${CMAKE_VERSION}-${OS}-x86_64.tar.gz

mkdir -p $HOME/cmake
curl -fsSL $URL | tar xz --strip-components=1 -C $HOME/cmake 
