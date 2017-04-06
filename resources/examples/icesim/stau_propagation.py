#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

import numpy as np

import os
import sys

os.chdir(sys.path[0])

# ptype = dc.I3Particle.STauMinus
ptype = dc.I3Particle.MuMinus
propagator = PROPOSAL.I3PropagatorServicePROPOSAL(
    type=ptype,
    mediadef="../configuration_IceOnly"
    # type=ptype,
    # particleMass=1000
)

stau = dc.I3Particle()
stau.type = ptype
stau.pos = dc.I3Position(0,0,0)
stau.dir = dc.I3Direction(0,0)
stau.energy = 100 * I3Units.TeV
stau.time = 0 * I3Units.ns
stau.location_type = dc.I3Particle.InIce

mu_length = list()
n_daughters = list()

primary_energy_epair = list()
secondary_energy_epair = list()

for i in range(10000):
    stau.length = NaN
    # returns None instead of an I3MMCTrack
    daughters = propagator.Propagate(stau)
    # length of daughters is always 1
    mu_length.append(stau.length)
    n_daughters.append(len(daughters))

    for particle in daughters:
        # print("-------------------------------")
        # print("particle type: ", particle.type)
        # print("primary energy", particle.speed)
        # print("secondary energy: ", particle.energy)

        if particle.type == dc.I3Particle.PairProd:
            primary_energy_epair.append(np.log10(particle.speed))
            secondary_energy_epair.append(np.log10(particle.energy))
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("Mu Lengths")
    pylab.hist(mu_length, histtype="step", log=True, bins=100)
    pylab.xlabel(r'$l_{\mu}(\rm{m})$')
    pylab.savefig('MuonLenghts.png')

    pylab.figure()
    pylab.title("N Daughters")
    pylab.hist(n_daughters, histtype="step", log=True, bins=100)
    pylab.xlabel('N')
    pylab.savefig('Daughters.png')

    from matplotlib.colors import LogNorm

    pylab.figure()
    pylab.title("Energyloss")
    pylab.hist2d(
        primary_energy_epair,
        secondary_energy_epair,
        bins=60,
        norm=LogNorm()
    )
    pylab.colorbar()
    pylab.savefig('Energyloss.png')
except ImportError :
    print("pylab not installed.  no plots for you.")


