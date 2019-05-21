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
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3066511.svg)](https://doi.org/10.5281/zenodo.3066511)

# PROPOSAL #

PROPOSAL (Propagator with optimal precision and optimized speed for all
leptons) is presented as a public tool for muon propagation through transparent
media. Up-to-date cross sections for ionization, bremsstrahlung, photonuclear
interactions, electron pair production, Landau–Pomeranchuk–Migdal and
Ter-Mikaelian effects, muon and tau decay, as well as Molière scattering are
implemented for different parametrizations.
The full Paper can be found
[here](https://doi.org/10.1016/j.cpc.2013.04.001).
Recent improvements are documented [here](https://doi.org/10.1016/j.cpc.2019.03.021).

PROPOSAL was tested on Mac OS X V. 10.13.6, Ubuntu 12.04, SUSE Enterprise 10 and PCLinuxos. Since
all these OS are UNIX based it should be fine to run and compile PROPOSAL on a UNIX based OS.

PROPOSAL is now a C++11 library using pybind11 Python bindings!


## How to cite PROPOSAL?

If you use PROPOSAL, please cite the PROPOSAL paper

```bibtex
@article{koehne2013proposal,
  title     ={PROPOSAL: A tool for propagation of charged leptons},
  author    = {Koehne, Jan-Hendrik and
               Frantzen, Katharina and
               Schmitz, Martin and
               Fuchs, Tomasz and
               Rhode, Wolfgang and
               Chirkin, Dmitry and
               Tjus, J Becker},
  journal   = {Computer Physics Communications},
  volume    = {184},
  number    = {9},
  pages     = {2070--2090},
  year      = {2013},
  publisher = {Elsevier},
  doi       = {10.1016/j.cpc.2013.04.001}
}
```
and our zenodo entry of the version you use
```bibtex
@misc{mario_dunsch_2019_2598747,
  author = {Dunsch, Mario and
            Soedingrekso, Jan and
            Koehne, Jan-Hendrik and
            Fuchs, Tomasz and
            Alameddine, Jean-Marco and
            Sackel, Maximilian and
            van Santen, Jacob and
            Kopper, Claudio and
            Krings, Kai and
            Olivas, Alex},
  title  = {tudo-astroparticlephysics/PROPOSAL: Zenodo},
  month  = may,
  year   = 2019,
  doi    = {10.5281/zenodo.3066511},
  url    = {https://doi.org/10.5281/zenodo.3066511}
}
```
and if you want to cite the latest improvements
```bibtex
@online{dunsch_2018_proposal_improvements,
  title       = {Recent Improvements for the Lepton Propagator PROPOSAL},
  author      = {Dunsch, Mario and
                 Soedingrekso, Jan and
                 Sandrock, Alexander and
                 Meier, Max and
                 Menne, Thorben and
                 Rhode, Wolfgang},
  year        = {2018},
  eprint      = {1809.07740},
  eprinttype  = {arxiv},
  eprintclass = {hep-ph},
  journal     = {Computer Physics Communications},
  doi         = {10.1016/j.cpc.2019.03.021}
}
```


## Requirements ##

- Boost Library 1.48 or higher
- [log4cplus](https://github.com/log4cplus/log4cplus) 2.0.2
- CMake 2.8 or higher

## Recommended ##

- Doxygen (For pdf and html documentation of the code)
- [pybind11](https://github.com/pybind/pybind11)
  (To build the python wrapper)
  If you decide to build the python wrapper and pybind11 is not
  provided on your system, pybind11 will be cloned to the project
  source folder.

## Installation ##

Install and compiling instructions for the standalone installation
are found in [install](INSTALL.md).


## Usage ##

### Deployment ###

PROPOSAL is built as library. So you can include this project in your own
c++ project by including the header files. The following snippet uses the
[configuration](resources/config.json) to propagate muons and
store the muon ranges for further proceeds.
The parameters of the configuration file are described
[here](resources/config_docu.md).

```c++
#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main(){
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

}
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
set (CMAKE_CXX_STANDARD 11)

add_executable(foo source/foo.cpp)

find_library(PROPOSAL_LIBRARIES REQUIRED NAMES PROPOSAL)

if (PROPOSAL_LIBRARIES)
  target_link_libraries (foo ${PROPOSAL_LIBRARIES})
endif ()
```

The file can then be compiled with
    
    cmake . 
    
and

    make

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

prop = pp.Propagator(
	particle_def=pp.particle.MuMinusDef.get(),
	config_file="path/to/config.json"
)

mu = prop.particle
mu_backup = pp.particle.Particle(mu)

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

The C++ API can be built using

	make doc

A documentation of the configuration file can be found
[here](resources/config_docu.md).

## Issues ##

When you encounter any errors or misunderstandings don't hesitate and write a mail to
[Jan Soedingrekso](mailto:jan.soedingrekso@tu-dortmund.de),
[Alexander Sandrock](mailto:alexander.sandrock@tu-dortmund.de).

## License ##

This software may be modified and distributed under the terms of
a modified LGPL License. See the LICENSE for details of the LGPL License.

Modifcations of the LGPL [License](LICENSE.md):

1. The user shall acknowledge the use of PROPOSAL by citing the following reference:

	> J.H. Koehne et al.
	> Comput.Phys.Commun. 184 (2013) 2070-2090
	> DOI: 10.1016/j.cpc.2013.04.001

2. The user should report any bugs/errors or improvements to the current maintainer of PROPOSAL.

## Developers and Maintainers ##

*Jan Soedingrekso*, *Alexander Sandrock*, *Jean-Marco Alameddine*, *Maximilian Sackel*

## Former Developers and Maintainers ##

*Jan-Hendrik Koehne*, *Tomasz Fuchs*, *Mario Dunsch*

## Acknowledgement ##

![SFB876](https://raw.githubusercontent.com/wiki/tudo-astroparticlephysics/Cor-PlusPlus/images/sfb876.png)
This work was created as part of the project [C3](http://sfb876.tu-dortmund.de/SPP/sfb876-c3.html) of the [SFB876](http://sfb876.tu-dortmund.de/index.html).
