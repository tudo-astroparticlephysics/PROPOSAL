
# Install #

1. 	Make a directory where the whole project will be, e.g:
	```sh
	mkdir PROPOSAL
	```
2.	Create a build and src directory, e.g.:
	```sh
	mkdir PROPOSAL/src PROPOSAL/build
	```

3. 	Extract the sources from the homepage/gitlab to the folder, e.g.:
	```sh
	unzip PROPOSAL.zip PROPOSAL/
	```
	or
	```sh
	git clone "..." PROPOSAL/src
	```
4.	Move to the build directory and generate the MakeFile with cmake:
	```sh
	cd PROPOSAL/build
	cmake ../src
	```
	If you don't want to compile the pybindings, call cmake with
	```sh
	cmake ../src -DADD_PYTHON=OFF
	```
6.  Compile the project:
	```sh
	make
	```
	or
	```sh
	make -j#
	```
	with `#` being the number of processors you can choose to compile
	the project on multiple processors simultaneously.
