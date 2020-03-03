# Installation

## Build Prequesites

The following commands can be used to install the required and optional
dependencies on your system.
The c++ dependencies are vendored using git submodules and are included in the
release tar balls.

### Ubuntu / Debian

```
$ apt install g++ \
  cmake \
  doxygen
```

### Arch Linux

```
$ pacman -S g++ \
  cmake \
  doxygen
```

### RHEL / Centos

Note: RHEL provides both cmake version 2.x and 3.x, the version 3.x executable is
called `cmake3`.
PROPOSAL needs cmake3.

```
$ yum install g++ cmake3 doxygen
```

### MacOS

```
$ xcode-select --install
$ brew install cmake \
  doxygen
```

## Building PROPOSAL

1. Download a release tarball from <https://github.com/tudo-astroparticlephysics/PROPOSAL/releases>
  and extract it or use git (the recursive is needed because we use submodules for dependencies):
  ```
  $ git clone --recursive https://github.com/tudo-astroparticlephysics/PROPOSAL
  ```

  To use a specific version of PROPOSAL, use `git checkout vX.Y.Z` after cloning.

  In case you want to perform the unit tests, you need to download the TestFiles.
  ```
  $ curl https://proposal.app.tu-dortmund.de/resources/2020-02-26/TestFiles.tar.gz --output tests/TestFiles.tar.gz
  ```

1.  Move to the build directory and generate the Makefile with cmake:

  ```
  mkdir build
  cmake ..
  ```

  If you don't want to compile the pybindings, call cmake with

  ```
  cmake .. -DADD_PYTHON=OFF
  ```

  To specify an installation prefix other than `/usr/local` use

  ```
  cmake .. -DCMAKE_INSTALL_PREFIX=/custom/prefix
  ```

  The prefix is used to install PROPOSAL on your system (See
  item 6).  
  To show further installation options use `cmake ..` and/or
  visit the [documentation](https://cmake.org/documentation/).
  Also have a look at the additional cmake options down below.

  #### Note

  * The option `CMAKE_INSTALL_PREFIX` adds the given path also to the
  include directories. So if you have installed PROPOSAL with
  `CMAKE_INSTALL_PREFIX` and are modifying the header files, you will have to 
  to uninstall PROPOSAL before the next build otherwise your local
  changes won't be used.

  * To ensure, that cmake finds the right python paths use these
  cmake options and adjust the version number to your python version:
  ```
  -DPYTHON_LIBRARY=$(python-config --prefix)/lib/libpython2.7.dylib
  -DPYTHON_INCLUDE_DIR=$(python-config --prefix)/include/python2.7
  ```

1.  Compile the project:

    make

  or

    make -j#

  with `#` being the number of processors you can choose to compile
  the project on multiple processes simultaneously.

1.  Install the library into `-DCMAKE_INSTALL_PREFIX`, by default `/usr/share`

  ```
  make install
  ```

# Build types

CMake uses `CMAKE_BUILD_TYPE` when building with make, the default
is set to `Debug` when in a git checkout and to `Release` otherwise.
Two other options exist: `RelWithDebInfo` and `MinSizeRel`.
See https://cmake.org/cmake/help/v3.17/variable/CMAKE_BUILD_TYPE.html

# Additional Cmake options #

| Option | Default value | Description |
| --- | --- | --- |
| `ADD_PYTHON` | ON | Compile the python wrapper |
| `ADD_PERFORMANCE_TEST` | OFF | Compile the performace test source |
| `ADD_ROOT` | ON | Compile PROPOSAL with ROOT support |
| `ADD_TESTS` | OFF | Compile unit tests. This downloads [googletest](https://github.com/google/googletest). |

**Examples**

  cmake -DADD_PYTHON=OFF <further options>

# Compiling #

After installation PROPOSAL can be used as a C++ library and easily included in any `*.cxx` file with the command

    #include "PROPOSAL/PROPOSAL.h"

Assuming `PROPOSAL.h` has been included in a file with the name `example.cxx`, this file can be compiled with

    g++ example.cxx -std=c++11 -lPROPOSAL <further options>
 
 or
 
    gcc -lstdc++ example.cxx -std=c++11 -lPROPOSAL <further options>


# Uninstalling #

It is also possible to uninstall PROPOSAL with

  make uninstall

This will remove all files listed in `install_mainfest.txt` which should
have been created in your build directory after the installation.
