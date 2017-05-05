# PROPOSAL

This is the restructured PROPOSAL Version 2, which should run in the IceCube Simulation.


## Installation

Before installation of icesim, one should copy some files from the
resources folder to other IceCube Simulation projects to provide compatibility:

* `cp src/PROPOSAL/resources/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/`

* `cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments/`

* `cp src/PROPOSAL/resources/simprod-scripts/proposal_alternate_cross_sections.py src/simprod-scripts/resources/examples/`

* `cp src/PROPOSAL/resources/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/`

---

the three following files ave not been changed yet

* `cp src/PROPOSAL/resources/MuonGun/propagate_muons.py src/MuonGun/resources/scripts/`

* `cp src/PROPOSAL/resources/MuonGun/shower_and_propagate.py src/MuonGun/resources/scripts/`

* `cp src/PROPOSAL/resources/MuonGun/utils.py src/MuonGun/resources/scripts/`

Further install instruction are found in [install](INSTALL.txt).

## Authors:
*Mario Dunsch* and *Jan Soedingrekso*

Former Developer and Maintainer:

*Jan-Hendrick Koehne* and *Tomasz Fuchs*