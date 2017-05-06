```
###############################################################################
#                                                                             #
#            _____  _____   ____  _____   ____   _____         _              #
#           |  __ \|  __ \ / __ \|  __ \ / __ \ / ____|  /\   | |             #
#           | |__) | |__) | |  | | |__) | |  | | (___   /  \  | |             #
#           |  ___/|  _  /| |  | |  ___/| |  | |\___ \ / /\ \ | |             #
#           | |    | | \ \| |__| | |    | |__| |____) / ____ \| |____         #
#           |_|    |_|  \_\\____/|_|     \____/|_____/_/    \_\______|        #
#                                                                             #
#                            _   _    ___   ___                               #
#                           | | | |  / _ \ (   )                              #
#                           | |_| | |  __/  | |                               #
#                           | ._,_|  \___|   \_)                              #
#                           | |                                               #
#                           |_|                                               #
#                                                                             #
###############################################################################
```


# PROPOSAL #

This is the restructured PROPOSAL Version 2, which should run in the IceCube Simulation.

PROPOSAL was tested on Mac OS X V. 10.7.5, Ubuntu 12.04, SUSE Enterprise 10 and PCLinuxos. Since
all these OS are UNIX based it should be fine to run and compile PROPOSAL on a UNIX based OS.

## Requirements ##

- Boost Library 1.48 or higher
- [log4cplus](https://github.com/log4cplus/log4cplus) 1.1.0 or higher
- CMake 2.8 or higher

## Recommended ##

- Doxygen (For pdf and html documentation of the code)
- [googletest](https://github.com/google/googletest)
  (To be sure not to break program code if you change something)

- [ROOT](https://root.cern.ch/)
  (This is highly recommended since there are lots of example plots and you are able to save the output in an root file.)

## Installation ##

#### Standalone ####

Install instruction for the standalone installation
are found in [install](INSTALL.md).

---

#### IceSim ####

Before installation of icesim, one should copy some files from the
resources folder to other IceCube Simulation projects to provide compatibility:

```bash
cp src/PROPOSAL/resources/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/
cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/
cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/
cp src/PROPOSAL/resources/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/
```
---

The three following files have not been changed yet:


* `cp src/PROPOSAL/resources/MuonGun/shower_and_propagate.py src/MuonGun/resources/scripts/`

* `cp src/PROPOSAL/resources/MuonGun/utils.py src/MuonGun/resources/scripts/`


## Usage ##

### Console ###

You can execute the PROPOSAL main from from the build directory. (e.g.
`./bin/PROPOSALtest`) The main part of the configuration of the propagator
routine are the configuration files which you can find in resources. The file
is fully documented and should guide you through your configuration.  Even if
you haven't installed root you should find some interesting code in the
`root_examples`.  All particle coordinates take the detector as the origin of the
coordinate system.

### Python ###

How to use PROPOSAL within Python is demonstrated with some example
scripts you can find in
[resources/examples/standalone](resources/examples/standalone).

For a short demonstration the following snippet will create data you can use to
show the distribution of muon ranges and the number of interactions in ice.

```python
import pyPROPOSAL

ptype = pyPROPOSAL.ParticleType.MuMinus

med = pyPROPOSAL.Medium("ice")
E = pyPROPOSAL.EnergyCutSettings()
prop = pyPROPOSAL.Propagator(med, E, ptype, "resources/tables")

mu_length = list()
n_daughters = list()

for i in range(1000):
    prop.reset_particle()
    prop.particle.energy = 1e8  # Unit in MeV
    d = p.propagate()

    mu_length.append(prop.particle.propagated_distance / 100)
    n_daughters.append(len(d))
```


## Issues ##

When you encounter any errors or misunderstandings don't hesitate and write a mail to
Tomasz.Fuchs@tu-dortmund.de
Jan-Hendrik.Koehne@tu-dortmund.de

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

```
Boost_INCLUDE_DIR (ADVANCED)
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
```

It seems that you have not installed boost or didn't set the include directory for it properly.
To install boost following
e.g. http://www.linuxfromscratch.org/blfs/view/svn/general/boost.html	(02.Jan.2014)

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

```
LOG4CPLUS_INCLUDE_DIR (ADVANCED)
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
   used as include directory in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
LOG4CPLUS_LIBRARY (ADVANCED)
	linked by target "PROPOSAL" in directory /home/tomasz/PROPOSALUpgrade/PROPOSALsrc
```

It seems that you have not installed log4cplus or didn't set the directorys for it properly.
To install log4cplus use
e.g. http://sourceforge.net/projects/log4cplus/files/log4cplus-stable/1.1.1/	(02.Jan.2014)

---
CMake Error: The following variables are used in this project, but they are set to NOTFOUND.

In general when some variables are not found by cmake you can have a look into the FindModules
files of cmake to see what is actually happening.
In Ubuntu they are located at /usr/share/cmake-2.8/Modules

Some modules are found in the main source directory or in the resource folder.

## License ##

[License](LICENSE.md)

## Authors ##

*Mario Dunsch*, *Jan Soedingrekso*

## Former Developer and Maintainer ##

*Jan-Hendrick Koehne*, *Tomasz Fuchs*
