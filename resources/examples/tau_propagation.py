#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import PROPOSAL
from icecube import simclasses

propagator = PROPOSAL.I3PropagatorServicePROPOSAL("-tau")

tau = dc.I3Particle()
tau.type = dc.I3Particle.TauMinus
tau.pos = dc.I3Position(0,0,0)
tau.dir = dc.I3Direction(0,0)
tau.energy = 100 * I3Units.TeV
tau.time = 0 * I3Units.ns
tau.location_type = dc.I3Particle.InIce


tau_length = list()
for i in range(100000):
    tau.length = NaN
    daughters = dc.I3ParticleVect()
    # returns None instead of an I3MMCTrack
    mmctrack = propagator.propagate(tau,daughters)
    # length of daughters is always 1
    tau_length.append(tau.length)

try:
    import pylab
    pylab.title("Tau Lengths")
    pylab.hist(tau_length, histtype = "step", log = True, bins = 100)
    pylab.xlabel(r'$l_{\tau}(\rm{m})$')
    pylab.show()
except ImportError :
    print "pylab not installed.  no plots for you."
    

