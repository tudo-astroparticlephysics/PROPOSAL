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
  doxygen \
  curl
```

### Arch Linux

```
$ pacman -S g++ \
  cmake \
  doxygen \
  curl
```

### RHEL / Centos

Note: RHEL provides both cmake version 2.x and 3.x, the version 3.x executable is
called `cmake3`.
PROPOSAL needs cmake3.

```
$ yum install g++ cmake3 doxygen curl
```

### MacOS

```
$ xcode-select --install
$ brew install cmake \
  doxygen 
```


## Installing from PyPI

If you only want to use PROPOSAL from python, go ahead and 
```
$ pip install [--user] proposal
```

## Building PROPOSAL

1. Download a release tarball from <https://github.com/tudo-astroparticlephysics/PROPOSAL/releases>
    and extract it or use git (the recursive is needed because we use submodules for dependencies):
    ```sh
    $ git clone --recursive https://github.com/tudo-astroparticlephysics/PROPOSAL
    ```

    To use a specific version of PROPOSAL, use `git checkout vX.Y.Z` after cloning.


1. Create a build directory and generate the build configuration:
    ```sh
    $ mkdir build
    $ cmake ..
    ```
    
    If you don't want to compile the pybindings, call cmake with
    ```sh
    $ cmake .. -DADD_PYTHON=OFF
    ```
    
    To specify an installation prefix other than `/usr/local` use
    ```sh
    $ cmake .. -DCMAKE_INSTALL_PREFIX=/custom/prefix
    ```
    
    The prefix is used to install PROPOSAL on your system.  
    To show further installation options use `cmake ..` and/or
    visit the [documentation](https://cmake.org/documentation/).
    Also have a look at the additional cmake options down below.
    
    **Note**
    
    * The option `CMAKE_INSTALL_PREFIX` adds the given path also to the
    include directories. So if you have installed PROPOSAL with
    `CMAKE_INSTALL_PREFIX` and are modifying the header files, you will have to 
    to uninstall PROPOSAL before the next build otherwise your local
    changes won't be used.
    
    * To ensure, that cmake finds the right python paths use these
    cmake options and adjust the version number to your python version:
    ```sh
    -DPYTHON_LIBRARY=/path/to/python/shared/library
    -DPYTHON_INCLUDE_DIR=/path/to/python/include/directory
    ```


1. Compile the project:
    ```sh
    $ cmake --build . [-j CORES]
    ```
    with `#` being the number of processors you can choose to compile
    the project on multiple processes simultaneously.

1.  Install the library into `-DCMAKE_INSTALL_PREFIX`, by default `/usr/share`
    ```sh
    $ cmake --build . --target install
    ```

# Build types

CMake uses `CMAKE_BUILD_TYPE` when building with make, the default
is set to `Debug` when in a git checkout and to `Release` otherwise.
Two other options exist: `RelWithDebInfo` and `MinSizeRel`.
See https://cmake.org/cmake/help/v3.17/variable/CMAKE_BUILD_TYPE.html

# Additional cmake options

| Option | Default value | Description |
| --- | --- | --- |
| `ADD_PYTHON` | ON | Compile the python wrapper |
| `ADD_PERFORMANCE_TEST` | OFF | Compile the performance test source |
| `ADD_ROOT` | ON | Compile PROPOSAL with ROOT support |


# Compiling your executables using PROPOSAL


## Simple executable

After installation PROPOSAL can be used as a C++ library and easily included in any `*.cxx` file with the command

```cpp
#include "PROPOSAL/PROPOSAL.h"
```

Assuming `PROPOSAL.h` has been included in a file with the name `example.cxx`, this file can be compiled with

```sh
$ g++ example.cxx -std=c++11 -lPROPOSAL <further options>
```

##  CMake Project

PROPOSAL exports a cmake config file, suitable for use with `find_package`.
To use PROPOSAL as dependency in your own CMake project, use:

```cmake
add_executable(example ...)
find_package(PROPOSAL REQUIRED)
target_link_libraries(example PRIVATE PROPOSAL::PROPOSAL)  # or PUBLIC
```

In case you did install PROPOSAL in a custom prefix, use `PROPOSAL_DIR` to tell
cmake where to find PROPOSAL:

```
$ PROPOSAL_DIR=/path/to/proposal/prefix cmake ...
```

# Tests

We are using `gtest`, which is vendored via a git submodule and which is also
included in the release tarballs.
Tests are automatically compiled when building PROPOSAL itself directly rather
than as dependency in your own CMake project.
To deactivate the compilation of tests pass `-DBUILD_TESTING=off` to the cmake call.

To run the tests, go into the build directory and call
```sh
$ `ctest [-V]` 
```

There are several very long running regression tests, disable them by
invoking
```sh
$ ctest -V -E '(Brems|Photo|Epair|Mupair)'
```

Test resources are not included in the git repository for size reasons and
are downloaded by cmake from our server <proposal.app.tu-dortmund.de/resources>.
 

# Uninstalling

It is also possible to uninstall PROPOSAL with

```
$ make uninstall
```

This will remove all files listed in `install_mainfest.txt` which should
have been created in your build directory after the installation.
