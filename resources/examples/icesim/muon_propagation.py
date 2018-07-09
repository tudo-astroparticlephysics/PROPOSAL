#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

ptype = dc.I3Particle.MuMinus
propagator = PROPOSAL.I3PropagatorServicePROPOSAL(
    config_file="$I3_BUILD/PROPOSAL/resources/config_icesim.json"
)

mu = dc.I3Particle()
mu.type = ptype
mu.pos = dc.I3Position(0, 0, 0)
mu.dir = dc.I3Direction(0, 0)
mu.energy = 100 * I3Units.TeV
mu.time = 0 * I3Units.ns
mu.location_type = dc.I3Particle.InIce

mu_length = list()
n_daughters = list()

for i in range(1000):
    mu.length = 0
    mu.pos = dc.I3Position(0, 0, 0)
    mu.dir = dc.I3Direction(0, 0)
    mu.energy = 100 * I3Units.TeV
    mu.time = 0 * I3Units.ns

    daughters = propagator.Propagate(mu)

    mu_length.append(mu.length)
    n_daughters.append(len(daughters))

try:
    import matplotlib as mpl
    mpl.use('Agg')
    import pylab

    pylab.figure()
    pylab.title(str(mu.energy) + " GeV Muon Lengths")
    pylab.hist(mu_length, histtype="step", log=True, bins=50)
    pylab.xlabel(r'$l_{\mu}(\rm{m})$')
    pylab.savefig('muon_lenghts.pdf')

    pylab.figure()
    pylab.title("N Daughters of " + str(mu.energy) + "GeV Muon")
    pylab.hist(n_daughters, histtype="step", log=True, bins=50)
    pylab.xlabel('N')
    pylab.savefig('muon_daughters.pdf')
except ImportError:
    print("pylab not installed.  no plots for you.")
