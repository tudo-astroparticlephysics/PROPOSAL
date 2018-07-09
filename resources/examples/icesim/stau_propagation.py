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

ptype = dc.I3Particle.STauMinus

stau = dc.I3Particle()
stau.type = ptype
stau.pos = dc.I3Position(0, 0, 0)
stau.dir = dc.I3Direction(0, 0)
stau.energy = 100 * I3Units.TeV
stau.time = 0 * I3Units.ns
stau.location_type = dc.I3Particle.InIce

propagator = PROPOSAL.I3PropagatorServicePROPOSAL(
    config_file="$I3_BUILD/PROPOSAL/resources/config_icesim.json"
)

propagator.register_particletype(ptype)

mu_length = list()
n_daughters = list()

primary_energy_epair = list()
secondary_energy_epair = list()

for i in range(1000):
    print(i)

    stau.length = 0
    stau.pos = dc.I3Position(0, 0, 0)
    stau.dir = dc.I3Direction(0, 0)
    stau.energy = 100 * I3Units.TeV
    stau.time = 0 * I3Units.ns

    daughters = propagator.Propagate(stau)
    n_daughters.append(len(daughters))

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("N Daughters (Stau)")
    pylab.hist(n_daughters, histtype="step", log=True, bins=50)
    pylab.xlabel('N')
    pylab.savefig('stau_daughters.pdf')

except ImportError:
    print("pylab not installed.  no plots for you.")
