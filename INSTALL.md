
# Installation #

## Install dependencies ##

The following commands can be used to install the required and optional
dependencies on your system.

### Ubuntu 16.10 ###

	apt install cmake \
		doxygen \
		liblog4cplus-dev
		libgtest-dev \
		libboost-dev \
		libboost-python-dev \

#### **Note** ####

The package `libgtest-dev` only installs the source files.
To create the gtest libraries you have compile these source files.
Therefore you can use:

	cd /usr/src/gtest
	sudo cmake CMakeLists.txt
	sudo make

	# copy or symlink libgtest.a and libgtest_main.a to your /usr/lib folder
	sudo cp *.a /usr/lib

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
		boost-python

#### **Note** ####

For Mac OS X `GTest` must be installed from source.
Therefore clone the repo from
[googletest](https://github.com/google/googletest)
into a local build directory and install googletest:

	cd && mkdir build
	cd build && mkdir gtest
	cd gtest && mkdir src build
	git clone https://github.com/google/googletest src
	cd build
	cmake ../src
	make
	make install

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

		git clone ... PROPOSAL/src

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

	The option `CMAKE_INSTALL_PREFIX` adds the given path also to the
	include directories. So if you have installed PROPOSAL with
	`CMAKE_INSTALL_PREFIX` and are modifying the header files, you will make
	sure to uninstall PROPOSAL before the next build otherwise your local
	changes won't be used.

6.  Compile the project:

		make

	or

		make -j#

	with `#` being the number of processors you can choose to compile
	the project on multiple processors simultaneously.

7.	Install the library on your system

		make install

	or for e.g.

		make install DESTDIR=$HOME

	The latter command will install PROPOSAL in `$HOME/<prefix>`, where
	the `prefix` was defined by `CMAKE_INSTALL_PREFIX` which defaults
	to `usr/local`.

# Additional Cmake options #

| Option | Default value | Description |
| --- | --- | --- |
| `ADD_PYTHON` | ON | Choose to compile the python wrapper |
| `ADD_PERFORMANCE_TEST` | OFF | Choose to compile the performace test source |
| `ADD_ROOT` | ON | Choose to compile PROPOSAL with ROOT support |

The test file generator is useful ensure not break the basic
functionality, if use decide to modify the library.

**Examples**

	cmake -DADD_PYTHON=OFF <further options>

# Uninstalling #

It is also possible to uninstall PROPOSAL with

	make uninstall

This will remove all files listed in `install_mainfest.txt` which should
have been created in your build directory after the installation.
