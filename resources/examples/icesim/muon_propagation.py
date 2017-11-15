#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

ptype = dc.I3Particle.MuMinus
propagator = PROPOSAL.I3PropagatorServicePROPOSAL(type = ptype)

mu = dc.I3Particle()
mu.type = ptype
mu.pos = dc.I3Position(0,0,0)
mu.dir = dc.I3Direction(0,0)
mu.energy = 100 * I3Units.TeV
mu.time = 0 * I3Units.ns
mu.location_type = dc.I3Particle.InIce

mu_length = list()
n_daughters = list()

for i in range(10000):
    mu.length = NaN
    # returns None instead of an I3MMCTrack
    daughters = propagator.Propagate(mu)
    # length of daughters is always 1
    mu_length.append(mu.length)
    n_daughters.append(len(daughters))

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt

    plt.figure()
    plt.title(str(mu.energy) + " GeV Muon Lengths")
    plt.hist(mu_length, histtype = "step", log = True, bins = 100)
    plt.xlabel(r'$l_{\mu}(\rm{m})$')
    plt.savefig('muon_Lenghts.png')

    plt.figure()
    plt.title("N Daughters of " + str(mu.energy) + "GeV Muon")
    plt.hist(n_daughters, histtype = "step", log = True, bins = 100)
    plt.xlabel('N')
    plt.savefig('muon_Daughters.png')
except ImportError :
    print("matplotlib not installed.  no plots for you.")


