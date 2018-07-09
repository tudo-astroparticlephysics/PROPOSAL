#!/usr/bin/env python

from I3Tray import I3Units
from I3Tray import NaN
from icecube import dataclasses as dc
from icecube import icetray
from icecube import PROPOSAL
from icecube import simclasses

ptype = dc.I3Particle.TauMinus

tau = dc.I3Particle()
tau.type = ptype
tau.pos = dc.I3Position(0, 0, 0)
tau.dir = dc.I3Direction(0, 0)
tau.energy = 1000 * I3Units.TeV
tau.time = 0 * I3Units.ns
tau.location_type = dc.I3Particle.InIce

propagator = PROPOSAL.I3PropagatorServicePROPOSAL(
    config_file="$I3_BUILD/PROPOSAL/resources/config_icesim.json"
)
print(tau)

tau_length = list()
n_daughters = list()

for i in range(1000):
    tau.length = 0
    tau.pos = dc.I3Position(0, 0, 0)
    tau.dir = dc.I3Direction(0, 0)
    tau.energy = 100 * I3Units.TeV
    tau.time = 0 * I3Units.ns

    daughters = propagator.Propagate(tau)
    tau_length.append(tau.length)
    n_daughters.append(len(daughters))

try:
    import pylab
    pylab.figure()
    pylab.title("Lengths (Tau)")
    pylab.hist(tau_length, histtype="step", log=True, bins=50)
    pylab.xlabel(r'$l_{\tau}(\rm{m})$')
    pylab.savefig('tau_lengths.pdf')

    pylab.figure()
    pylab.title("N Daughters (Tau)")
    pylab.hist(n_daughters, histtype="step", log=True, bins=50)
    pylab.xlabel('N')
    pylab.savefig('tau_daughters.pdf')
except ImportError:
    print("pylab not installed.  no plots for you.")
