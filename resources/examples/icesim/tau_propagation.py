#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

ptype = dc.I3Particle.TauMinus
propagator = PROPOSAL.I3PropagatorServicePROPOSAL(type = ptype)

tau = dc.I3Particle()
tau.type = ptype
tau.pos = dc.I3Position(0,0,0)
tau.dir = dc.I3Direction(0,0)
tau.energy = 1000 * I3Units.TeV
tau.time = 0 * I3Units.ns
tau.location_type = dc.I3Particle.InIce

tau_length = list()
n_daughters = list()

for i in range(10000):
    tau.length = NaN
    # daughters = dc.I3ParticleVect()
    # returns None instead of an I3MMCTrack
    # frame = icetray.I3Frame()
    daughters = propagator.Propagate(tau) # ,frame)
    # mmctrack = propagator.propagate(tau,daughters)
    # length of daughters is always 1
    tau_length.append(tau.length)
    n_daughters.append(len(daughters))

try:
    import pylab
    pylab.figure()
    pylab.title("Tau Lengths")
    pylab.hist(tau_length, histtype = "step", log = True, bins = 100)
    pylab.xlabel(r'$l_{\tau}(\rm{m})$')
    pylab.savefig('tau_Lengths.png')

    pylab.figure()
    pylab.title("N Daughters of Tau")
    pylab.hist(n_daughters, histtype = "step", log = True, bins = 100)
    pylab.xlabel('N')
    pylab.savefig('tau_Daughters.png')
except ImportError :
    print("pylab not installed.  no plots for you.")
