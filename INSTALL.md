	eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
	e                                                                                              e
	e      K#5z#zy   W#5z9#X,      zE9y    K#Xz#zW      XEEX      WGeGy      9yXK    EXG           e
	e      Geeeeeeee eeeeeeeee  Weeeeeeee  eeeeeeeeD  eeeeeeee  ,eeeueeee   Weeee    eee           e
	e      Xee   eee 9ee   Eee  eee    eee Dee   eee eee    eee Deee        ee eee   eee           e
	e      Xeeeeeee  EeeeeeeX   eee    eee Deeeeeee  eee    eee   Weeeeez  eee  eey  eee           e
	e      Xee,      9ee  eee,  eee   Weee Dee       eee    eee eee   eee  eeeeeeee  eee           e
	e      eeee      eee9  eeeE  eeeeeeee  eeeD       eeeeeeee  KeeeeeeeW eee KuXeee eeeeeeee      e
	e       uu        ,u     uz     uu      uK           uu        uWK    X      u z eeeeeeee      e
	e                                                                                              e
	e                                                                                              e
	e                                         5   Ku                  yKKy#K                       e
	e                    ee ee                ee  ee                eeeeeee#                       e
	e                   9ee eee               e9  ee                   eD                          e
	e                   eey                   e9  ee                   eK                          e
	e                    eeeee5               eeeEDee                  eW 9                        e
	e                     ,ED                 e                        eee#                        e
	e                                         e                                                    e
	e                                                                                              e
	eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeze


# Install #

1. 	Make a directory where the whole project will be, e.g:
	```bash
	mkdir PROPOSAL
	```
2.	Create a build and src directory, e.g.:
	```bash
	mkdir PROPOSAL/src PROPOSAL/build
	```

3. 	Extract the sources from the homepage/gitlab to the folder, e.g.:
	```bash
	unzip PROPOSAL.zip PROPOSAL/
	```
	or
	```bash
	git clone "..." PROPOSAL/src
	```
4.	Move to the build directory and generate the MakeFile with cmake:
	```bash
	cd PROPOSAL/build
	cmake ../src
	```
	If you don't want to compile the pybindings, call cmake with
	```bash
	cmake ../src -DADD_PYTHON=OFF
	```
6.  Compile the project:
	```bash
	make
	```
	or
	```bash
	make -j#
	```
	with `#` being the number of processors you can choose to compile
	the project on multiple processors simultaneously.
