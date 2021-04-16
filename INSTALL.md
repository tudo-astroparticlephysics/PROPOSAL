# Installation

The installation process is based on **CMake** (version `>3.9` required).
To handle dependencies and for simple usage in other projects, the package manager **conan** (version `>1.33` required) can be used additionally.
If you are only interested to use PROPOSAL in python, you can install PROPOSAL using **pip**.

All different installation approaches are going to be explained in the following.

For more detailed information about the specific building tools, see the listed documentations: 

- [conan documentation](https://docs.conan.io/en/latest/)
- [CMake documentation](https://cmake.org/cmake/help/latest/)
- [pip documentation](https://pip.pypa.io/en/stable/)


## Building using conan (recommended for C++ users)

For this installation approach, all dependencies will be fetched by conan, meaning that you don't have to install them by yourself. If you have not installed conan yet, you can do so, for example:

```sh
$ pip install conan
``` 

Clone the repository and create a build directory

```sh
$ git clone https://github.com/tudo-astroparticlephysics/PROPOSAL.git
$ cd PROPOSAL && mkdir build && cd build      
```

Use conan to prepare all dependencies. You can pass additional options to conan.

```sh             
$ conan install .. -o with_python=True	# other optional dependencies
```

The following options can be passed to `conan install`:

| Option.              | default | Description                                   |
| -------------------- | ------- | --------------------------------------------- |
| `with_python`        | False   | Build and install python interface.           |
| `with_testing`       | False   | Build TestFiles for Python.                   |
| `with_documentation` | False   | Build doxygen documentation of C++ code (WIP) |

Build and install PROPOSAL. You may require root privileges when installing, depending on the installation location:

```sh
$ conan build ..
# cmake --install .
```

*Note:* As an alternative, you may create a local conan package and use it in your project. See the [conan documentation](https://docs.conan.io/en/latest/) for more information.

## Building using pip (recommended for Python users)

If you only want to use PROPOSAL in Python, the easiest way is to use **pip**:

```sh
$ pip install proposal
```

This will install the most recent version of PROPOSAL. 
See the [pip documentation](https://pip.pypa.io/en/stable/) for more information.

## Building using CMake (recommended for advanced users)

If you want to provide all dependencies by your own, you can skip conan and only use CMake.
In this case, you need to provide the following dependencies for CMake:

- [spdlog](https://github.com/gabime/spdlog)
- [CubicInterpolation](https://github.com/MaxSac/cubic_interpolation)
- [pybind11](https://github.com/pybind/pybind11.git) (only if you want to install python bindings)
- [gtest](https://github.com/google/googletest) (only if you want to use the UnitTests)
- [doxygen](https://github.com/doxygen/doxygen) (only if you want to build the documentation)

Clone the repository and create a build directory

```sh
$ git clone https://github.com/tudo-astroparticlephysics/PROPOSAL.git
$ cd PROPOSAL && mkdir build && cd build      
```

Use CMake and make to build and install PROPOSAL. We recommend building PROPOSAL in Release since performance will be better by several orders of magnitude:

```sh
$ cmake .. -DCMAKE_BUILD_TYPE=Release	# or other CMake options
$ make -j
# sudo make install
```

There are several CMake options you may pass:

| Option.               | default | Description                                   |
| --------------------- | ------- | --------------------------------------------- |
| `BUILD_PYTHON`        | OFF     | Build and install python interface.           |
| `BUILD_TESTING`       | OFF     | Build TestFiles for Python.                   |
| `BUILD_DOCUMENTATION` | OFF     | Build doxygen documentation of C++ code (WIP) |


# Minimal working example

If you have installed PROPOSAL at a relocatable place, linking against it
should be straight forward.

Create a simple Cmake file called `CMakeLists.txt`,
which searches for the target PROPOSAL and link your executable against it.

```cmake
CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

PROJECT(PROPOSAL_test)

find_package(PROPOSAL REQUIRED)

add_executable(test test.cpp)
target_link_libraries(test PROPOSAL::PROPOSAL)
```

Run CMake and build the executable.

```sh
$ cmake -B build .
$ make -C build 
$ ./build/bin/test
```

Hopefully everything has been built correctly and your program will be executed.
Otherwise, create an issue with the corresponding error code.

Enjoy!

# FAQ:

This section provides help for common problems during installation. 

##### I am experiencing linker errors when importing `proposal` in python or when trying to link PROPOSAL.
If you have installed PROPOSAL using conan, check your conan profile for the `compiler.libcxx` setting, for example with

```sh
$ cat .conan/profiles/default`
``` 

if you are using the default conan profile. This should be set to `libstdc++11`.
If not, you can change this setting with

```sh
$ conan profile update settings.compiler.libcxx=libstdc++11 default
```

before installing for your default profile.

