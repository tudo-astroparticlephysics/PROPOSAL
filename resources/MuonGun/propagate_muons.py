#!/usr/bin/env python

import numpy
#import dashi
#import pycorsika

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
from icecube import phys_services, c2j_icetray, mmc_icetray
from icecube import MuonGun
icetray.load('corsika-reader')
from I3Tray import I3Tray

from os.path import expandvars

# See:
# Y. Becherini, A. Margiotta, M. Sioli, and M. Spurio. A parameterisation of single and multiple muons in the deep water or ice. Astroparticle Physics, 25(1):1 - 13, 2006.
# http://dx.doi.org/10.1016/j.astropartphys.2005.10.005

def MMCFactory(length=10*I3Units.m, seed=12345, mediadef=expandvars('$I3_BUILD/MuonGun/resources/iceworld-mediadef')):
	jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
	jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
	jvmOpts.append("-Xms256m")
	jvmOpts.append("-Xmx512m")
	jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")

	jvm = c2j_icetray.I3JavaVM(jvmOpts)

	# Now create the MMC propagators, but first *all* of the options must be set here. 
	# There's no special options added behind the scenes.  This is much more flexible. 
	#  Below are the standard options.  To interpret them see the MMC docs.
	mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
	mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
	mmcOpts += expandvars("-mediadef=%s " % mediadef)
	mmcOpts += "-radius=100000 "
	mmcOpts += "-length=%d " % length
	mmcOpts += "-seed=%d " % seed

	return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)

from optparse import OptionParser
parser = OptionParser(usage="%prog [OPTIONS] infiles outfile")
parser.add_option("--nevents", dest="nevents", type=int, default=int(1e6), help="Number of showers per Corsika file")

opts, args = parser.parse_args()
if len(args) < 2:
	parser.error("You must supply at least one input and one output file")
infiles, outfile = args[:-1], args[-1]
# infile, outfile = '/data/sim/IceCube/2012/generated/corsika/H/DAT000000', 'DAT000000_muons.i3'

tray = I3Tray()

tray.AddService('I3GSLRandomServiceFactory', 'rng')

if infiles[0].endswith('.i3.bz2'):
	tray.AddModule('I3Reader', 'reader', filenamelist=infiles)
else:
	tray.AddModule('I3CORSIKAReader', 'reader', filenamelist=infiles,
	    NEvents=opts.nevents, CylinderHeight=0, CylinderRadius=0,
	    ParticlesToWrite=MuonGun.I3ParticleTypeSeries([dataclasses.I3Particle.MuMinus, dataclasses.I3Particle.MuPlus]),
	)


tray.AddModule('Muonitron', 'propatron',
    Depths=list(numpy.linspace(1.0, 4.5, 8)*I3Units.km),
    # Depths=list(numpy.arange(0.1, 1.0, 0.1)*I3Units.km),
    MMC=MMCFactory(), CylinderHeight=10.*I3Units.m,
)
def printy(frame):
	#mctree = frame['I3MCTree']
	#print mctree.primaries[0]
	#print mctree.primaries[0].dir
	# Drop any frames where *no* muons made it to depth
	tracks = frame['Tracks']
	return len(tracks) > 0
tray.AddModule(printy, 'printy', Streams=[icetray.I3Frame.DAQ])

tray.AddModule('Keep', 'killcruft', Keys=['I3EventHeader', 'CorsikaWeightDict', 'CorsikaWeightMap', 'Tracks', 'MCPrimary'])

tray.AddModule('Dump', 'dumply')

tray.AddModule('I3Writer', 'writer',
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    # DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=outfile)

tray.AddModule('TrashCan', 'YesWeCan')
tray.Execute()
tray.Finish()
 
