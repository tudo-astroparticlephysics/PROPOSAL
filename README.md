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

PROPOSAL is now a C++14 library using pybind11 Python bindings!
In the next major release, C++17 methods may also be used in the core library, which is currently only used when building the tests.


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

- CMake 3.9 or higher (to build the tests CMake 3.10 is required)
- C++14 compatible compiler

Furthermore, you either need the package manager

- conan

which will provide all dependencies that are necessary for PROPOSAL or you need to provide these dependencies by your own. For more information, see [here](INSTALL.md).

## Installation

Install and compiling instructions are found in [install](INSTALL.md).

## Usage

### Usage as a C++ library

PROPOSAL is built as library. So you can include this project in your own
C++ project by including the header files. The following snippet uses the
[configuration](examples/config_minimal.json) to propagate muons and
store the muon ranges for further proceeds.
The parameters of the configuration file are described
**[here](docs/config_docu.md)**.

```c++
#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

int main(){
    auto mu_def = MuMinusDef();
    Propagator prop(mu_def, "path/to/config.json");

    Cartesian3D position(0, 0, 0);
    Cartesian3D direction(0, 0, 1);
    auto energy = 1e8; // MeV
    auto init_state = ParticleState(position, direction, energy, 0., 0.);

    std::vector<double> ranges;

    for (int i = 0; i < 10; i++)
    {
        auto track = prop.Propagate(init_state, 50000); // distance to propagate in cm

        ranges.push_back(track.back().propagated_distance);
    }

// ... Do stuff with ranges, e.g. plot histogram

}
```

To see an example on how to run this script with PROPOSAL using CMake, see [here](INSTALL.md).

### Usage in Python

How to use PROPOSAL within Python is demonstrated with some example
jupyter notebooks you can find in the
[examples](examples) folder.

For a short demonstration the following snippet will create data you can use to
show the distribution of muon ranges and the number of interactions in ice.
The parameters of the given configuration file are described
[here](docs/config_docu.md).

```python
import proposal as pp

mu_def = pp.particle.MuMinusDef()
prop = pp.Propagator(
	  particle_def=mu_def,
	  config_file="path/to/config.json"
)

init_state = pp.particle.ParticleState()
init_state.energy = 1e9 # initial energy in MeV
init_state.position = pp.Cartesian3D(0, 0, 0)
init_state.direction = pp.Cartesian3D(0, 0, 1)

mu_length = []

for i in range(1000):
    track = prop.propagate(init_state)

    mu_length.append(track.track_propagated_distances()[-1] / 100)
```


## Issues ##

When you encounter any errors or misunderstandings, you can always create an issue here on GitHub.
Furthermore, you may always contact us with your questions via
[Jean-Marco Alameddine](mailto:jean-marco.alameddine@tu-dortmund.de),
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
