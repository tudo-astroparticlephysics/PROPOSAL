#!/usr/bin/env python

import os

from icecube.dataclasses import I3Particle
from icecube.PROPOSAL import I3PropagatorServicePROPOSAL

I3PropagatorServicePROPOSAL(type=I3Particle.MuMinus)
I3PropagatorServicePROPOSAL(type=I3Particle.TauMinus)

# This should have generated all data files.
# Generate a marker file to tell the build system this task
# has been successfully completed.
f = open(os.path.expandvars("$I3_BUILD/PROPOSAL/resources/tables/.tables.auto_generated"), "w")
f.close()
