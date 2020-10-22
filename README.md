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


# PROPOSAL [![Build Status](https://travis-ci.org/tudo-astroparticlephysics/PROPOSAL.svg?branch=master)](https://travis-ci.org/tudo-astroparticlephysics/PROPOSAL) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1484180.svg)](https://doi.org/10.5281/zenodo.1484180)

PROPOSAL (Propagator with optimal precision and optimized speed for all
leptons) is presented as a public tool for propagating leptons and gamma rays
through media.
Up-to-date cross sections for ionization, bremsstrahlung, photonuclear
interactions, electron pair production, Landau–Pomeranchuk–Migdal and
Ter-Mikaelian effects, muon and tau decay, as well as Molière scattering are
implemented for different parametrizations.
The full Paper can be found
[here](https://doi.org/10.1016/j.cpc.2013.04.001).
Recent improvements are documented [here](https://doi.org/10.1016/j.cpc.2019.03.021).

PROPOSAL is developed and tested on macOS and linux. 
Continuous integration is setup on travis and tests several version of gcc and clang.

PROPOSAL is now a C++11 library using pybind11 Python bindings!

This version is the latest not using C++14 or more modern methods.
The next release will make use of C++14 methods.


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
  doi       = {10.1016/j.cpc.2013.04.001}
}
```
and our zenodo entry of the version you use
```bibtex
@misc{dunsch_2020_1484180,
  author = {Dunsch, Mario and
            Soedingrekso, Jan and
            Koehne, Jan-Hendrik and
            Fuchs, Tomasz and
            Alameddine, Jean-Marco and
            Sackel, Maximilian and
            Noethe, Maximilian and
            van Santen, Jacob and
            Menne, Thorben and
            Sandrock, Alexander and
            Kuo, Chung-Yun and
            Kopper, Claudio and
            Krings, Kai and
            Olivas, Alex},
  title  = {tudo-astroparticlephysics/PROPOSAL: Zenodo},
  month  = mar,
  year   = 2020,
  doi    = {10.5281/zenodo.1484180},
  url    = {https://doi.org/10.5281/zenodo.1484180}
}
```
and if you want to cite the latest improvements
```bibtex
@article{dunsch_2018_proposal_improvements,
  title       = {Recent Improvements for the Lepton Propagator PROPOSAL},
  author      = {Dunsch, Mario and
                 Soedingrekso, Jan and
                 Sandrock, Alexander and
                 Meier, Max and
                 Menne, Thorben and
                 Rhode, Wolfgang},
  journal     = {Computer Physics Communications},
  volume      = {242},
  pages       = {132--144},
  year        = {2019},
  eprint      = {1809.07740},
  eprinttype  = {arXiv},
  eprintclass = {hep-ph},
  doi         = {10.1016/j.cpc.2019.03.021}
}
```


## Requirements

- CMake 3.8 or higher
- C++11 compatible compiler

## Recommended

- Doxygen (For pdf and html documentation of the code)
- [pybind11](https://github.com/pybind/pybind11)
  (To build the python wrapper)
  If you decide to build the python wrapper and pybind11 is not
  provided on your system, pybind11 will be cloned to the project
  source folder.

## Installation

Install and compiling instructions for the standalone installation
are found in [install](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/INSTALL.md).


## Usage

### Deployment

PROPOSAL is built as library. So you can include this project in your own
c++ project by including the header files. The following snippet uses the
[configuration](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/resources/config.json) to propagate muons and
store the muon ranges for further proceeds.
The parameters of the configuration file are described
**[here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/resources/config_docu.md)**.

```c++
#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main(){
    ParticleDef mu_def = MuMinusDef::Get();
    Propagator prop(mu_def, "resources/configuration/config.json");
    DynamicData mu(mu_def.particle_type);

    mu.SetEnergy(9e6);
    mu.SetPosition(Vector3D(0, 0, 0))
    mu.SetDirection(Vector3D(0, 0, -1));

    std::vector<double> ranges;

    for (int i = 0; i < 10; i++)
    {
        Secondaries sec = prop.Propagate(mu);

        ranges.push_back(sec.GetPosition().back().magnitude());
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
        └── foo.cxx

the `CMakeLists.txt` could look like

```cmake
cmake_minimum_required(VERSION 3.8)

add_executable(foo source/foo.cxx)

find_package(PROPOSAL REQUIRED)
target_link_libraries(foo PRIVATE PROPOSAL::PROPOSAL)  # or PUBLIC
```

In case you did install PROPOSAL in a custom prefix, use `PROPOSAL_DIR` to tell
cmake where to find PROPOSAL:

```
$ mkdir build && cd build
$ PROPOSAL_DIR=/path/to/proposal/prefix cmake .. [CMAKE options]
$ cmake --build .
```


### Python ###

How to use PROPOSAL within Python is demonstrated with some example
scripts you can find in
[resources/examples/standalone](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/resources/examples/standalone).

For a short demonstration the following snippet will create data you can use to
show the distribution of muon ranges and the number of interactions in ice.
The parameters of the given configuration file are described
[here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/resources/config_docu.md).

```python
import proposal as pp

mu_def = pp.particle.MuMinusDef()
prop = pp.Propagator(
	  particle_def=mu_def,
	  config_file="path/to/config.json"
)

mu = pp.particle.DynamicData(mu_def.particle_type)

mu.energy = 9e6
mu.direction = pp.Vector3D(0, 0, -1)

mu_length = []
mu_secondaries = []

for i in range(1000):
    sec = prop.propagate(mu)

    mu_length.append(sec.position[-1].magnitude() / 100)
    mu_secondaries.append(sec.number_of_particles)
```

## Documentation ##

The C++ API can be built using

	make doc

A documentation of the configuration file can be found
[here](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/resources/config_docu.md).

## Issues ##

When you encounter any errors or misunderstandings don't hesitate and write a mail to
[Jan Soedingrekso](mailto:jan.soedingrekso@tu-dortmund.de),
[Alexander Sandrock](mailto:alexander.sandrock@tu-dortmund.de).

## License ##

This software may be modified and distributed under the terms of
a modified LGPL License. See the LICENSE for details of the LGPL License.

Modifications of the LGPL [License](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/LICENSE.md):

1. The user shall acknowledge the use of PROPOSAL by citing the following reference:

	> J.H. Koehne et al.
	> Comput.Phys.Commun. 184 (2013) 2070-2090
	> DOI: 10.1016/j.cpc.2013.04.001

2. The user should report any bugs/errors or improvements to the current maintainer of PROPOSAL.

## Developers and Maintainers ##

*Jan Soedingrekso*, *Alexander Sandrock*, *Jean-Marco Alameddine*, *Maximilian Sackel*

## Former Developers and Maintainers ##

*Jan-Hendrik Koehne*, *Tomasz Fuchs*, *Mario Dunsch*

## Acknowledgment ##

![SFB876](https://raw.githubusercontent.com/wiki/tudo-astroparticlephysics/Cor-PlusPlus/images/sfb876.png)
This work was created as part of the project [C3](http://sfb876.tu-dortmund.de/SPP/sfb876-c3.html) of the [SFB876](http://sfb876.tu-dortmund.de/index.html).
