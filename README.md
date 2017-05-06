# PROPOSAL #

This is the restructured PROPOSAL Version 2, which should run in the IceCube Simulation.

PROPOSAL was tested on Mac OS X V. 10.7.5, Ubuntu 12.04, SUSE Enterprise 10 and PCLinuxos. Since
all these OS are UNIX based it should be fine to run and compile PROPOSAL on a UNIX based OS.

## Requirements ##

- Boost Library 1.48 or higher
- log4cplus 1.1.0 or higher
- CMake 2.8 or higher

## Recommended ##

- Doxygen 	(For pdf and html documentation of the code)
- GTEST 	(To be sure not to break program code if you change something)
- ROOT		(This is highly recommended since there are lots of example plots and you are able to save the output in an root file.)

## Installation ##

Before installation of icesim, one should copy some files from the
resources folder to other IceCube Simulation projects to provide compatibility:

```
cp src/PROPOSAL/resources/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/
cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/
cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/
cp src/PROPOSAL/resources/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/
```
---

The three following files have not been changed yet:


* `cp src/PROPOSAL/resources/MuonGun/shower_and_propagate.py src/MuonGun/resources/scripts/`

* `cp src/PROPOSAL/resources/MuonGun/utils.py src/MuonGun/resources/scripts/`

Further install instruction are found in [install](INSTALL.txt).


## Usage ##

You can execute the PROPOSAL main from from the build directory. (e.g.
`./bin/PROPOSALtest`) The main part of the configuration of the propagator
routine are the configuration files which you can find in resources. The file
is fully documented and should guide you through your configuration.  Even if
you haven't installed root you should find some interesting code in the
`root_examples`.  All particle coordinates take the detector as the origin of the
coordinate system.

## Issues ##

When you encounter any errors or misunderstandings don't hesitate and write a mail to
Tomasz.Fuchs@tu-dortmund.de
Jan-Hendrik.Koehne@tu-dortmund.de

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

Boost_INCLUDE_DIR (ADVANCED)
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc

It seems that you have not installed boost or didn't set the include directory for it properly.
To install boost following
e.g. http://www.linuxfromscratch.org/blfs/view/svn/general/boost.html	(02.Jan.2014)

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

LOG4CPLUS_INCLUDE_DIR (ADVANCED)
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
LOG4CPLUS_LIBRARY (ADVANCED)
	linked by target "PROPOSAL" in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc

It seems that you have not installed log4cplus or didn't set the directorys for it properly.
To install log4cplus use
e.g. http://sourceforge.net/projects/log4cplus/files/log4cplus-stable/1.1.1/	(02.Jan.2014)

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

In general when some variables are not found by cmake you can have a look into the FindModules
files of cmake to see what is actually happening.
In Ubuntu they are located at /usr/share/cmake-2.8/Modules

Some modules are found in the main source directory or in the resource folder.

## Authors ##

*Mario Dunsch*, *Jan Soedingrekso*

## Former Developer and Maintainer ##

*Jan-Hendrick Koehne*, *Tomasz Fuchs*
