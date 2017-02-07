This is the restructured PROPOSAL Version 2, which should run in the IceCube Simulation.


Before installation of icesim, one should copy some files from the resources folder to other IceCube Simulation projects to provide compatibility:

cp src/PROPOSAL/resources/MuonGun/MuonPropagator.* src/MuonGun/private/MuonGun/
cp src/PROPOSAL/resources/simprod-scripts/PropagateMuons.py src/simprod-scripts/python/segments
cp src/PROPOSAL/resources/simprod-scripts/proposal_alternative_cross_sections.py src/simprod-scripts/resources/examples/
cp src/PROPOSAL/resources/clsim/PropagateMuons.py src/clsim/resources/scripts/photonPaths/

Authors:
Mario Dunsch
and
Jan Soedingrekso

Former Developer and Maintainer:
Jan-Hendrick Koehne
and
Tomasz Fuchs
