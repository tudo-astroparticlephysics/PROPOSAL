
# Installation #

1. 	Make a directory where the whole project will be, e.g:

		mkdir PROPOSAL

2.	Create a build and src directory, e.g.:

		mkdir PROPOSAL/src PROPOSAL/build

3. 	Extract the sources from the homepage/gitlab to the folder, e.g.:

		unzip PROPOSAL.zip PROPOSAL/

	or

		git clone "..." PROPOSAL/src

4.	Move to the build directory and generate the MakeFile with cmake:

		cd PROPOSAL/build
		cmake ../src

	If you don't want to compile the pybindings, call cmake with

		cmake ../src -DADD_PYTHON=OFF

	To specify an installation location other than `/usr/local` use

		cmake ../src -DCMAKE_INSTALL_PREFIX=/path/to/install

	To show further installation options use `ccmake ../src`

6.  Compile the project:

		make

	or

		make -j#

	with `#` being the number of processors you can choose to compile
	the project on multiple processors simultaneously.

7.	Install the library on your system

		make install
