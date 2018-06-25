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

[![Build Status](https://travis-ci.org/tudo-astroparticlephysics/PROPOSAL.svg?branch=master)](https://travis-ci.org/tudo-astroparticlephysics/PROPOSAL)

# PROPOSAL #

PROPOSAL (Propagator with optimal precision and optimized speed for all
leptons) is presented as a public tool for muon propagation through transparent
media. Up-to-date cross sections for ionization, bremsstrahlung, photonuclear
interactions, electron pair production, Landau–Pomeranchuk–Migdal and
Ter-Mikaelian effects, muon and tau decay, as well as Molière scattering are
implemented for different parametrizations.
The full Paper can be found
[here](http://www.sciencedirect.com/science/article/pii/S0010465513001355)

PROPOSAL was tested on Mac OS X V. 10.7.5, Ubuntu 12.04, SUSE Enterprise 10 and PCLinuxos. Since
all these OS are UNIX based it should be fine to run and compile PROPOSAL on a UNIX based OS.

## Requirements ##

- Boost Library 1.48 or higher
- [log4cplus](https://github.com/log4cplus/log4cplus) 1.1.0 or higher
- CMake 2.8 or higher

## Recommended ##

- Doxygen (For pdf and html documentation of the code)
- [Boost.Python](http://www.boost.org/doc/libs/1_64_0/libs/python/doc/html/index.html)
  (To build the python wrapper)

## Installation ##

### Standalone ###

Install instruction for the standalone installation
are found in [install](INSTALL.md).

---

### IceSim ###

Before installation of icesim, one should copy some files from the
resources folder to other IceCube Simulation projects to provide compatibility:

```sh
cp src/PROPOSAL/resources/icesim/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/
cp src/PROPOSAL/resources/icesim/MuonGun/shower_and_propagate.py src/MuonGun/resources/scripts/
cp src/PROPOSAL/resources/icesim/MuonGun/utils.py src/MuonGun/resources/scripts/
cp src/PROPOSAL/resources/MuonGun/propagate_muons.py src/MuonGun/resources/scripts/
cp src/PROPOSAL/resources/icesim/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/
cp src/PROPOSAL/resources/icesim/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/
cp src/PROPOSAL/resources/icesim/dataclasses/I3Particle.cxx src/dataclasses/private/dataclasses/physics/
cp src/PROPOSAL/resources/icesim/dataclasses/I3Particle.h src/dataclasses/public/dataclasses/physics/
```
---

## Usage ##

### Deployment ###

PROPOSAL is build as library. So you can include this project in your own
c++ project by including the header files. The following snippet uses the
[configuration](resources/config.json) to propagte muons and
store the muon ranges for further proceeds.
The parameters of the configuration file are described
[here](resources/config_docu.md).

```c++
#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

Propagator prop(MuMinusDef::Get(), "resources/config.json");
Particle& mu = prop.GetParticle();
Particle mu_backup(mu);

mu_backup.SetEnergy(9e6);
mu_backup.SetDirection(Vector3D(0, 0, -1));

std::vector<double> ranges;

for (int i = 0; i < 10; i++)
{
  mu.InjectState(mu_backup);

  prop.Propagate();

  ranges.push_back(mu.GetPropagatedDistance());
}

// ... Do stuff with ranges, e.g. plot histogram

```

Supposing this snippet is the content of `foo.cxx` within the
following minimal code structure

    my_program
    ├── CMakeLists.txt
    ├── resources
    │   ├── configuration
    │   └── tables
    └── source
        └── foo.cpp

the `CMakeLists.txt` could look like

```
cmake_minimum_required(VERSION 2.6)

add_executable(foo source/foo.cpp)

find_library(PROPOSAL_LIBRARIES REQUIRED NAMES PROPOSAL)

if (PROPOSAL_LIBRARIES)
  include_directories(${PROPOAL_INCLUDE_DIRS})
  target_link_libraries (foo ${PROPOSAL_LIBRARIES})
endif ()
```

### Python ###

How to use PROPOSAL within Python is demonstrated with some example
scripts you can find in
[resources/examples/standalone](resources/examples/standalone).

For a short demonstration the following snippet will create data you can use to
show the distribution of muon ranges and the number of interactions in ice.
The parameters of the given configuration file are described
[here](resources/config_docu.md).

```python
import pyPROPOSAL as pp

prop = pp.Propagator(pp.MuMinusDef.get(), "path/to/config.json")
mu = prop.particle
mu_backup = pp.Particle(mu)

mu_backup.energy = 9e6
mu_backup.direction = pp.Vector3D(0, 0, -1)

mu_length = []
mu_secondaries = []

for i in range(1000):
    mu.inject_state(mu_backup)
    secondaries = prop.propagate()

    mu_length.append(prop.particle.propagated_distance / 100)
    mu_secondaries.append(len(secondaries))
```

## Documentation ##

The C++ API can be build using

	make doc

A documentation of the configuration file can be found
[here](resources/config_docu.md).

## Issues ##

When you encounter any errors or misunderstandings don't hesitate and write a mail to
[Jan Soedingrekso](mailto:jan.soedingrekso@tu-dortmund.de),
[Mario Dunsch](mailto:mario.dunsch@tu-dortmund.de),
[Tomasz Fuchs](mailto:Tomasz.Fuchs@tu-dortmund.de) or
[Jan-Hendrik Koehne](mailto:Jan-Hendrik.Koehne@tu-dortmund.de).

## License ##

This software may be modified and distributed under the terms of
a modified LGPL License. See the LICENSE for details of the LGPL License.

Modifcations of the LGPL [License](LICENSE.md):

1. The user shall acknowledge the use of PROPOSAL by citing the following reference:

	> J.H. Koehne et al.
	> Comput.Phys.Commun. 184 (2013) 2070-2090
	> DOI: 10.1016/j.cpc.2013.04.001

2. The user should report any bugs/errors or improvments to the current maintainer of PROPOSAL.

## Developers and Maintainers ##

*Mario Dunsch*, *Jan Soedingrekso*

## Former Developers and Maintainers ##

*Jan-Hendrick Koehne*, *Tomasz Fuchs*

## Acknowledgement ##

![SFB876](https://raw.githubusercontent.com/wiki/tudo-astroparticlephysics/Cor-PlusPlus/images/sfb876.png)
This work was created as part of the project [C3](http://sfb876.tu-dortmund.de/SPP/sfb876-c3.html) of the [SFB876](http://sfb876.tu-dortmund.de/index.html).
