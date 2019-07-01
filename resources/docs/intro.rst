
Introduction
------------

Leptons in PROPOSAL can have different interactions with the media they are in. These stochastic energy losses are:

- Ionization
- electron-Pairproduction
- Bremsstrahlung
- Photonuclear Interaction
- "Decay"

The decay is not a real stochastic loss, but is an interaction which the charged lepton can do. Below 1TeV the energy loss of muons is dominated by Ionization and above 1TeV by Pairproduction.

Not all interactions are simulated since this would require a lot of computing time. In PROPOSAL you can set an minimum energy for an stochastic loss which will then be calculated. This can be an absolute value or a relative value of the current particle energy. All energy losses below this limit are averaged over the propagated distance.
