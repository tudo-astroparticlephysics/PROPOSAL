
# Installation #

## Install dependencies ##

The following commands can be used to install the required and optional
dependencies on your system.

### Ubuntu 16.10 ###

	apt install cmake \
		doxygen \
		liblog4cplus-dev
		libboost-dev \

### Arch Linux ###

	pacman -S cmake \
		doxygen \
		log4cplus \
		gtest \
		boost-libs

### Mac OS X ###

	brew install cmake \
		doxygen \
		log4cplus \
		boost \

## Install PROPOSAL ##


1. 	Make a directory where the whole project will be, e.g.:

		mkdir PROPOSAL

2.	Create a build and src directory, e.g.:

		mkdir PROPOSAL/src PROPOSAL/build

3. 	Extract the sources from the
	[hompage](http://app.tu-dortmund.de/cms/de/Projekte/PROPOSAL/) or
	gitlab to the folder, e.g.:

		unzip PROPOSAL.zip PROPOSAL/

	or

		git clone https://github.com/tudo-astroparticlephysics/PROPOSAL.git PROPOSAL/src

4.	Move to the build directory and generate the Makefile with cmake:

		cd PROPOSAL/build
		cmake ../src

	If you don't want to compile the pybindings, call cmake with

		cmake ../src -DADD_PYTHON=OFF

	To specify an installation prefix other than `/usr/local` use

		cmake ../src -DCMAKE_INSTALL_PREFIX=/custom/prefix

	The prefix is used to install PROPOSAL on your system (See
	item 6).<br>
	To show further installation options use `ccmake ../src` and/or
	visit the [documentation](https://cmake.org/documentation/).
	Also have a look at the additional cmake option down below.

	#### **Note** ####

	* The option `CMAKE_INSTALL_PREFIX` adds the given path also to the
	  include directories. So if you have installed PROPOSAL with
	  `CMAKE_INSTALL_PREFIX` and are modifying the header files, you will make
	  sure to uninstall PROPOSAL before the next build otherwise your local
	  changes won't be used.

	* To ensure, that cmake finds the right python paths use these
	  cmake options and adjust the version number to your python version:

			-DPYTHON_LIBRARY=$(python-config --prefix)/lib/libpython2.7.dylib
			-DPYTHON_INCLUDE_DIR=$(python-config --prefix)/include/python2.7

6.  Compile the project:

		make

	or

		make -j#

	with `#` being the number of processors you can choose to compile
	the project on multiple processors simultaneously.

7.	Install the library on your system

		make install

	or e.g.

		make install DESTDIR=$HOME

	The latter command will install PROPOSAL in `$HOME/<prefix>`, where
	the `prefix` was defined by `CMAKE_INSTALL_PREFIX` which defaults
	to `usr/local`.

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
