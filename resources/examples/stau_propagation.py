#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

ptype = dc.I3Particle.STauMinus
propagator = PROPOSAL.I3PropagatorServicePROPOSAL(
    type=ptype
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
for i in range(10000):
    stau.length = NaN
    # returns None instead of an I3MMCTrack
    daughters = propagator.Propagate(stau)
    # length of daughters is always 1
    mu_length.append(stau.length)
    n_daughters.append(len(daughters))
try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title("Stau Lengths")
    pylab.hist(mu_length, histtype = "step", log = True, bins = 100)
    pylab.xlabel(r'$l_{s\tau}(\rm{m})$')
    pylab.savefig('sTau_Lenghts.png')

    pylab.figure()
    pylab.title("N Daughters of sTau")
    pylab.hist(n_daughters, histtype = "step", log = True, bins = 100)
    pylab.xlabel('N')
    pylab.savefig('sTau_Daughters.png')
except ImportError :
    print("pylab not installed.  no plots for you.")


