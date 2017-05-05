# PROPOSAL

This is the restructured PROPOSAL Version 2, which should run in the IceCube Simulation.

PROPOSAL was tested on Mac OS X V. 10.7.5, Ubuntu 12.04, SUSE Enterprise 10 and PCLinuxos. Since
all these OS are UNIX based it should be fine to run and compile PROPOSAL on a UNIX based OS.

## Requirements:
	- Boost Library 1.48 or higher
	- log4cplus 1.1.0 or higher
	- CMake 2.8 or higher

## Recommended:
	- Doxygen 	(For pdf and html documentation of the code)
	- GTEST 	(To be sure not to break program code if you change something)
	- ROOT		(This is highly recommended since there are lots of example plots and you are able to save the output in an root file.)

## Installation

Before installation of icesim, one should copy some files from the
resources folder to other IceCube Simulation projects to provide compatibility:

* `cp src/PROPOSAL/resources/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/`

* `cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/`

* `cp src/PROPOSAL/resources/simprod-scripts/proposal_alternate_cross_sections.py src/simprod-scripts/resources/examples/`

* `cp src/PROPOSAL/resources/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/`

---

The three following files have not been changed yet:


* `cp src/PROPOSAL/resources/MuonGun/shower_and_propagate.py src/MuonGun/resources/scripts/`

* `cp src/PROPOSAL/resources/MuonGun/utils.py src/MuonGun/resources/scripts/`

Further install instruction are found in [install](INSTALL.txt).

## Authors:

*Mario Dunsch*, *Jan Soedingrekso*

## Former Developer and Maintainer:

*Jan-Hendrick Koehne*, *Tomasz Fuchs*
