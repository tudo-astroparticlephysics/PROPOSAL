#!/usr/bin/env python

from optparse import OptionParser
from os.path import expandvars
parser = OptionParser(usage="%prog [OPTIONS] infiles outfile")
parser.add_option("--run", dest="run", type=int, default=1)
parser.add_option("--nevents", dest="nevents", type=int, default=int(1e6), help="Number of showers per Corsika file")

parser.add_option("--emin", dest="emin", type=float, default=5e2, help="Lower bound of energy range [GeV/nucleon] ")
parser.add_option("--emax", dest="emax", type=float, default=1e10, help="Upper bound of energy range [GeV]")
parser.add_option("--eslope", dest="eslope", type=float, default=-2., help="Slope of generation power law")
parser.add_option("--primary", dest="primary", default='H', choices=('H', 'He', 'N', 'Al', 'Fe'))

parser.add_option("--atmosphere", dest="atmosphere", type=int, default=11)
parser.add_option("--model", dest="model", default="SIBYLL", help="Hadronic interaction model", choices=('QGSJET', 'EPOS', 'SIBYLL', 'QGSII', 'DPMJET'))

parser.add_option("--corsika-dir", dest="rundir", type=str, default=expandvars('$I3_BUILD/bin/corsika-run'), help='Coriska run directory')

opts, args = parser.parse_args()
if len(args) < 1:
	parser.error("You must supply an output file name")
outfile = args[-1]

opts.primary = dict(H=14, He=402, N=1407, Al=2713, Fe=5626)[opts.primary]
if opts.primary > 100:
	opts.emin *= (opts.primary % 100)

print(repr(opts))

from icecube import icetray, dataclasses, dataio
from icecube.icetray import I3Units
from icecube import phys_services, sim_services
from icecube import MuonGun
icetray.load('corsika-reader')
from I3Tray import I3Tray

def MMCFactory(length=10*I3Units.m, seed=12345, impl='mmc', mediadef=None):

	if impl == 'mmc':
		from icecube import c2j_icetray, mmc_icetray
		if mediadef is None:
			mediadef=expandvars('$I3_BUILD/MuonGun/resources/iceworld-mediadef')
		# Now create the MMC propagators, but first *all* of the options must be set here. 
		# There's no special options added behind the scenes.  This is much more flexible. 
		#  Below are the standard options.  To interpret them see the MMC docs.
		mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
		mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
		mmcOpts += expandvars("-mediadef=%s " % mediadef)
		mmcOpts += "-radius=100000 "
		mmcOpts += "-length=%d " % length
		mmcOpts += "-seed=%d " % seed
	
		jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
		jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
		jvmOpts.append("-Xms256m")
		jvmOpts.append("-Xmx512m")
		jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")
		jvmOpts.append("-XX:+TraceHPI")

		jvm = c2j_icetray.I3JavaVM(jvmOpts)
		return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)
	else:
		from icecube import PROPOSAL_icetray
		# in PROPOSAL everything can be defined in the configuration file
		if mediadef is None:
			mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/config_iceworld.json')
		return PROPOSAL_icetray.I3PropagatorServicePROPOSAL(config_file=mediadef)

CORSIKA_CONFIG =\
"""
RUNNR   %(run)d                       run number
NSHOW   %(nshower)d                   number of showers to generate
PRMPAR  %(primary)d                prim. particle (1=gamma, 14=proton, ...)
ESLOPE  %(eslope).1f                  slope of primary energy spectrum
ERANGE  %(emin).1E %(emax).1E         energy range of primary particle (GeV)
THETAP  0.  89.                       range of zenith angle (degree)
PHIP    -180.  180.                   range of azimuth angle (degree)
SEED    %(seed1)d   0   0             seed for 1. random number sequence
SEED    %(seed2)d   0   0             seed for 2. random number sequence
OBSLEV  283400.                       observation level (in cm)
FIXCHI  0.                            starting altitude (g/cm**2)
MAGNET  16.59 -52.79                  magnetic field centr. Europe
HADFLG  0  1  0  1  0  2              flags hadr.interact.&fragmentation
ECUTS   300.  300.  10000000.  100000000.          energy cuts for particles
LOCUT   T   1.58                      zenith-dependent energy threshold
MUADDI  T                             additional info for muons
MUMULT  T                             muon multiple scattering angle
ELMFLG  F   T                         em. interaction flags (NKG,EGS)
STEPFC  1.0                           mult. scattering step length fact.
RADNKG  200.E2                        outer radius for NKG lat.dens.distr.
ARRANG  -119.                         rotation of array to north
ATMOD   %(atmod)d                     atmospheric model
LONGI   T  10.  F  F                  longit.distr. & step size & fit & out
ECTMAP  1.E2                          cut on gamma factor for printout
MAXPRT  0                             max. number of printed events
DIRECT  %(outdir)s                    output directory
DATBAS  F                             write .dbase file
USER    %(user)s                      user
PAROUT T F                            suppress DAT file
DEBUG   F  6  F  1000000              debug flag and log.unit for out
EXIT
"""

class CorsikaRunner(object):
	"""
	Run a modified version of CORSIKA that writes to a named pipe rather than a file.
	"""
	import getpass
	defaults = dict(run=1, nshower=1000, primary=14, eslope=-2., emin=5e2, emax=1e10, atmod=12, outdir='./', user=getpass.getuser())
	executable = './corsika6990Linux_%s_gheisha'
	fname = None
	def __init__(self, model='SIBYLL', rundir='/scratch/jvs/metaprojects/minisim/corsika-6990/run', **kwargs):
		from subprocess import Popen, PIPE
		import fcntl, random
		config = dict(self.defaults)
		config.update(kwargs)
		config['seed1'] = random.randint(0, 1<<32 - 1)
		config['seed2'] = random.randint(0, 1<<32 - 1)
		self.fname = config['outdir'] + 'DAT%.6d' % config['run']
		os.mkfifo(self.fname)
		self.p = Popen([self.executable % model], stdin=PIPE, stdout=open('/dev/null', 'w'), cwd=rundir)		
		self.p.stdin.write(CORSIKA_CONFIG % config)
		
	def __del__(self):
		if self.fname:
			os.unlink(self.fname)

import os
corsika = CorsikaRunner(run=opts.run, nshower=opts.nevents, primary=opts.primary,
    eslope=opts.eslope, emin=opts.emin, emax=opts.emax, atmod=opts.atmosphere,
    model=opts.model, outdir=os.getcwd()+'/', rundir=opts.rundir)

tray = I3Tray()

tray.AddService('I3GSLRandomServiceFactory', 'rng')

tray.AddModule('I3CORSIKAReader', 'reader', filenamelist=[corsika.fname],
    NEvents=opts.nevents, CylinderHeight=0, CylinderRadius=0,
    ParticlesToWrite=MuonGun.I3ParticleTypeSeries([dataclasses.I3Particle.MuMinus, dataclasses.I3Particle.MuPlus]),
)

# tray.AddModule('Dump', 'dump')

tray.AddModule('Muonitron', 'propatron',
    Depths=[1e3, 1.5e3, 2.e3, 2.5e3, 3.e3, 3.5e3, 4.e3, 4.5e3],
    MMC=MMCFactory(impl='PROPOSAL'), CylinderHeight=10.*I3Units.m,
)
def droppy(frame):
	# Drop any frames where *no* muons made it to depth
	tracks = frame['Tracks']
	if len(tracks) > 0:
		frame['CorsikaWeightDict']['Atmosphere'] = opts.atmosphere
	else:
		return False
tray.AddModule(droppy, 'droppy', Streams=[icetray.I3Frame.DAQ])

tray.AddModule('Keep', 'killcruft', Keys=['I3EventHeader', 'CorsikaWeightDict', 'CorsikaWeightMap', 'Tracks', 'MCPrimary'])

tray.AddModule('Dump', 'dumply')

tray.AddModule('I3Writer', 'writer',
    Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
    # DropOrphanStreams=[icetray.I3Frame.DAQ],
    filename=outfile)


tray.Execute()

 
