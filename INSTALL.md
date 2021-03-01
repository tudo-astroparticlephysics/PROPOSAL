# Installation

The installation is based on CMake and Conan. 
Therefore a cmake version `>3.9` and a conan version `>1.33` is required.
More information about how the programs are operated can be read in the 
listed documentations. 

- [conan documentation](https://docs.conan.io/en/latest/)
- [cmake documentation](https://cmake.org/cmake/help/latest/)

In fact there is no good reason why a simple user must interact with cmake more
than a simple installation of a version which is capable with the requirements. 

There are different targets which can be build by the building process. 
The targets can be activated or deactivated by the installation process

- `build_python`: python interface to interact with the lib
- `build_documentation`: doxygen cpp documentation
- `build_testing`: testing of the source code


## Build Prequesites

Requiremenst will be handled by conan.
Simply pass the optional targets behind the conan install process 
and let it do the work for you.

```sh
$ mkdir build && cd build                   
$ conan install .. -o build_python=True     # other optional dependencies
```

## Building PROPOSAL

Building could be achieved by using conan or cmake. 
The more straight forward approach is, that conan do the work for you.

```sh
$ conan build ..
```

That will populate the build directory with the set targets.

## Packageing and installation at relocatable place

It's very common to locate package at place where the compiler and linker can
find them. 
For that choose a place (`/usr/local/` for example is often a good choice) and 
create the package there. 
Dependend where you want to store the package, the necessary privileques are
required. 

```sh
# conan package .. -pf /usr/local
```

## Minimal working example

If you had installed PROPOSAL at a relocatable place place linking against it
should be straight forward.

Create a simple cmake file called `CMakeLists.txt`,
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

Hopefully everything had be build correctly and your programm would be executed.
Otherwise create an Issue with the corresponding error code.

Enjoy.
