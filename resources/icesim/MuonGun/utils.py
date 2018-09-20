from __future__ import print_function

from icecube.icetray import I3Units
from icecube import icetray
import random, sys

import dashi
class buffering_histogram(dashi.histogram.histogram):
	"""
	A histogram that consolidates multiple calls to fill() into more efficient blocks
	"""
	maxbuf = 65535
	def __init__(self, *args, **kwargs):
		super(buffering_histogram, self).__init__(*args, **kwargs)
		self._bh_sample = numpy.empty((self.maxbuf, self.ndim))
		self._bh_weights = numpy.empty(self.maxbuf)
		self._bh_pos = 0
	def flush(self):
		super(buffering_histogram, self).fill(self._bh_sample[:self._bh_pos], self._bh_weights[:self._bh_pos])
		self._bh_pos = 0
	def fill_single(self, values, weight):
		"""
		Add a single entry with the given weight
		"""
		self._bh_sample[self._bh_pos] = values
		self._bh_weights[self._bh_pos] = weight
		self._bh_pos += 1
		if self._bh_pos == self.maxbuf:
			self.flush()
	def fill(self, sample, weight):
		"""
		Add multiple entries, each with 1/N of the given weight
		"""
		n = sample.shape[0]
		# For large entries, kick straight to vanilla fill()
		if n > self.maxbuf:
			super(buffering_histogram, self).fill(sample, numpy.ones(n)*(weight/n))
			return 
		elif self._bh_pos + n > self.maxbuf:
			self.flush()
		self._bh_sample[self._bh_pos:self._bh_pos+n] = sample
		self._bh_weights[self._bh_pos:self._bh_pos+n] = weight/n
		self._bh_pos += n
		if self._bh_pos == self.maxbuf:
			self.flush()

def fluxsum(emin=6e2, emax=1e11, eslope=-2):
	g = eslope+1
	return (emax**g - emin**g)/g

def dcorsika_spectra(gamma=[-2.]*5, normalization=[10., 5., 3., 2., 1.], emin=6e2, emax=1e11, nevents=2.5e6):
	"""
	Calculate the generation spectra for the H, He, N, Al, and Fe components of 5-component dCORSIKA
	"""
	from icecube.weighting.weighting import PowerLaw
	masses = [1, 4, 14, 27, 56] 
	fluxsums = numpy.array([n*fluxsum(emin*m, emax, g) for m, g, n in zip(masses, gamma, normalization)])
	nshower = nevents*fluxsums/fluxsums.sum()
	print(nshower)
	return [PowerLaw(n, emin*m, emax, g) for n, m, g in zip(nshower, masses, gamma)]


class EnergyWeightCollection(object):
	def __init__(self, target_fluxes, generation_spectra):
		self.target_fluxes = target_fluxes
		self.generation_spectra = generation_spectra
	def __call__(self, ptype, E, zenith):
		return self.target_fluxes[ptype](E)/self.generation_spectra[ptype](E, zenith)

class IsotropicWeight(object):
	def __init__(self, f):
		self.f = f
	def __call__(self, E, ptype, theta):
		return self.f(E, ptype)

class VolumeCorrWeight(object):
	def __init__(self, detcfg, energy_spectrum):
		self.chi = detcfg*4./numpy.pi
		self.energy_spectrum = energy_spectrum
	def __call__(self, E, ptype, theta):
		return self.energy_spectrum(E, ptype)*(numpy.cos(theta) + self.chi*numpy.sin(theta))/(1/2. + numpy.pi*self.chi/4.)


from os.path import expandvars
def MMCFactory(radius=100000*I3Units.m, length=10*I3Units.m, seed=random.randint(0, (1<<32) - 1), impl='mmc', mediadef=None):

	if impl.lower() == 'mmc':
		from icecube import c2j_icetray, mmc_icetray

		if mediadef is None:
			mediadef=expandvars('$I3_BUILD/MuonGun/resources/iceworld-mediadef')
		# Now create the MMC propagators, but first *all* of the options must be set here. 
		# There's no special options added behind the scenes.  This is much more flexible. 
		# Below are the standard options.  To interpret them see the MMC docs.
		mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
		mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
		mmcOpts += expandvars("-mediadef=%s " % mediadef)
		mmcOpts += "-radius=%d " % radius
		mmcOpts += "-length=%d " % length
		mmcOpts += "-seed=%d " % seed

		jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
		jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
		jvmOpts.append("-Xms256m")
		jvmOpts.append("-Xmx512m")
		jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")
		
		jvm = c2j_icetray.I3JavaVM(jvmOpts)
		return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)
	else:
		from icecube import PROPOSAL_icetray
		# in PROPOSAL everything can be defined in the configuration file
		if mediadef is None:
			mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/config_iceworld.json')
		return PROPOSAL_icetray.I3PropagatorServicePROPOSAL(config_file=mediadef)

def PropagatorMMC(tray, name, seed=random.randint(0, (1<<32) - 1)):
	from icecube import c2j_icetray, mmc_icetray
	jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
	jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
	jvmOpts.append("-Xms256m")
	jvmOpts.append("-Xmx512m")
	jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")
	tray.AddService('I3JavaVMFactory', 'jvm', Options=jvmOpts)
	
	mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
	mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
	mmcOpts += expandvars("-mediadef=$I3_BUILD/mmc-icetray/resources/mediadef ")
	mmcOpts += "-radius=800 "
	mmcOpts += "-length=1600 "
	mmcOpts += "-seed=%d " % seed
	tray.AddModule('I3PropagatorMMC', name, mode=-1, opts=mmcOpts)


class CrustyPropagator(icetray.I3Module):
	def __init__(self, context):
		super(CrustyPropagator, self).__init__(context)
		self.AddOutBox("OutBox")
	
	def Configure(self):
		from icecube.MuonGun import MuonPropagator, Crust, Sphere, Cylinder

		self.crust = Crust(MuonPropagator("air", ecut=-1, vcut=5e-2, rho=0.673))
		self.crust.add_layer(Sphere(1948, 6374134), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=0.832))
		self.crust.add_layer(Sphere(1748, 6373934), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005))
		self.crust.add_layer(Cylinder(1600, 800),   MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005))
	
	def DAQ(self, frame):
		
		for track in frame['I3MCTree']:
			if not track.type in (track.MuPlus, track.MuMinus):
				continue
			if track.location_type != track.InIce:
				continue
			newtrack = self.crust.ingest(track)
			if newtrack.energy > 0:
				track.length = float('inf')
			else:
				track.length = 0
			
		self.PushFrame(frame)

def MuonitronPropagator(tray, name):
	
	def find_intersection(frame):
		from icecube.MuonGun import Cylinder
		from icecube.dataclasses import I3Constants, I3VectorDouble
		primary = frame['I3MCTree'].primaries[0]
		impact = Cylinder(1600, 800).intersection(primary.pos, primary.dir)
		d0 = I3Constants.SurfaceElev - I3Constants.OriginElev - (primary.pos.z + impact.first*primary.dir.z)
		frame['Depths'] = I3VectorDouble([d0])
	
	def patch_mctree(frame):
		tracks = frame['Tracks']
		if len(tracks) > 0:
			ntracks = len(list(frame['Tracks'].values())[0])
		else:
			ntracks = 0
		i = 0
		for track in frame['I3MCTree']:
			if not track.type in (track.MuPlus, track.MuMinus):
				continue
			if track.location_type != track.InIce:
				continue
			if i < ntracks:
				track.length = float('inf')
				i += 1
			else:
				track.length = 0
	
	tray.AddModule(find_intersection, name+'_find_intersection', Streams=[icetray.I3Frame.DAQ])
	
	from icecube.MuonGun import MuonPropagator, Crust, Sphere

	crust = Crust(MuonPropagator("air", ecut=-1, vcut=5e-2, rho=0.673))
	crust.add_layer(Sphere(1948, 6374134), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=0.832))
	crust.add_layer(Sphere(1748, 6373934), MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005))

	tray.AddModule('Muonitron', name+'propagator',
	    Depths=[0],
	    Propagator=MuonPropagator("ice", ecut=-1, vcut=5e-2, rho=1.005),
	    Crust=crust,
	)
	
	tray.AddModule(patch_mctree, name+'_patch_mctree', Streams=[icetray.I3Frame.DAQ])

import dashi, numpy, tables
from icecube import icetray, dataclasses
from icecube import sim_services, MuonGun
class Router(icetray.I3Module):
	def __init__(self, context):
		icetray.I3Module.__init__(self, context)
		self.AddOutBox("OutBox")
		for ptype in 'PPlus', 'He4Nucleus', 'N14Nucleus', 'Al27Nucleus', 'Fe56Nucleus':
			self.AddOutBox(ptype)
		
		self.AddParameter("Routes", "", dict())
	
	def Configure(self):
		
		self.routes = self.GetParameter("Routes")
		
		slots = set()
		for ptype, boxes in self.routes.items():
			for box in boxes:
				if not box in slots:
					self.AddOutBox(box)
					slots.add(box)
	def DAQ(self, frame):
		ptype = frame['MCPrimary'].type
		for box in self.routes.get(ptype, list()):
			self.PushFrame(frame, box)
		


class Filler(icetray.I3ConditionalModule):
	"""
	Tabulate tracks by depth, zenith angle, multiplicity, energy, and radius from shower axis
	"""
	def __init__(self, ctx):
		icetray.I3ConditionalModule.__init__(self, ctx)
		self.AddOutBox("OutBox")
		self.AddParameter("Weight", "", None)
		self.AddParameter("Outfile", "", None)
		self.AddParameter("MinDepth", "", 1.)
		self.AddParameter("MaxDepth", "", 5.)
		self.AddParameter("DepthSteps", "", 9)
	
	def Configure(self):
		from collections import defaultdict
		
		depthbins = numpy.linspace(self.GetParameter("MinDepth"), self.GetParameter("MaxDepth"), self.GetParameter("DepthSteps"))
		depthbins -= numpy.diff(depthbins)[0]/2.
		zenbins = numpy.arccos(numpy.linspace(1, 0, 11))
		zenbins_fine = numpy.arccos(numpy.linspace(1, 0, 101))
		multbins = numpy.array([1, 2, 3, 4, 10, 20, 40, 100], dtype=float)
		rbins = numpy.array([0, 5, 10, 15, 25, 45], dtype=float)
		
		self.primary = buffering_histogram(2, (zenbins, numpy.logspace(2, 11, 101)))
		self.multiplicity = dashi.histogram.histogram(3, (zenbins_fine, depthbins, numpy.arange(1, 100)))
		self.radius = dashi.histogram.histogram(4, (zenbins, depthbins, multbins, numpy.linspace(0, numpy.sqrt(250), 101)**2))
		self.energy = dashi.histogram.histogram(5, (zenbins, depthbins, multbins, rbins, numpy.logspace(0, 6, 101)))
		
		self.multiplicity_slices = tuple([tuple([buffering_histogram(1, (numpy.arange(1, 100),)) for j in range(len(depthbins))]) for i in range(len(zenbins_fine))])
		self.radius_slices = tuple([tuple([buffering_histogram(2, (multbins, numpy.linspace(0, numpy.sqrt(250), 101)**2)) for j in range(len(depthbins))]) for i in range(len(zenbins))])
		self.energy_slices = tuple([tuple([buffering_histogram(3, (multbins, rbins, numpy.logspace(0, 6, 101))) for j in range(len(depthbins))]) for i in range(len(zenbins))])
		
		self.depthbins = depthbins
		self.zenbins = zenbins
		self.zenbins_fine = zenbins_fine
		
		self.weighter = self.GetParameter("Weight")
		
		import os
		components = os.path.split(self.GetParameter("Outfile"))
		if os.path.splitext(components[-1])[1] in ('.hdf5', '.h5'):
			self.where = '/'
			self.outfile = os.path.join(*components)
		else:
			self.where = '/' + components[-1]
			self.outfile = os.path.join(*components[:-1])
		
		if os.path.exists(self.outfile):
			os.unlink(self.outfile)
			
		self.nevents = 0
		
	def DAQ(self, frame):
		primary = frame['MCPrimary']
		
		zenith = primary.dir.zenith
		zi = max(numpy.searchsorted(self.zenbins, zenith) - 1, 0)
		zif = max(numpy.searchsorted(self.zenbins_fine, zenith) - 1, 0)
		
		weight = self.weighter[primary.type](primary.energy, zenith)
		self.primary.fill_single((zenith, primary.energy), weight)
		
		multiplicity=self.multiplicity_slices[zif]
		radius=self.radius_slices[zi]
		energy=self.energy_slices[zi]
		
		for di, (depth, tracks) in enumerate(frame['Tracks'].items()):
			kmwe = depth/I3Units.km
			mult = len(tracks)
			values = numpy.asarray([(mult, p.radius, p.energy) for p in tracks])
			
			multiplicity[di].fill_single(mult, weight)
			radius[di].fill(values[:,:2], weight)
			energy[di].fill(values, weight)
		
		self.nevents += 1
		# if self.nevents % 10000 == 0:
		# 	print '%d events' % self.nevents
		
		# self.PushFrame(frame)
		
	def Finish(self):
		for i in range(len(self.zenbins_fine)):
			for j in range(len(self.depthbins)):
				self.multiplicity_slices[i][j].flush()
				self.multiplicity._h_bincontent[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_bincontent
				self.multiplicity._h_squaredweights[i+1,j+1,:] += self.multiplicity_slices[i][j]._h_squaredweights
		for i in range(len(self.zenbins)):
			for j in range(len(self.depthbins)):
				self.radius_slices[i][j].flush()
				self.energy_slices[i][j].flush()
				self.radius._h_bincontent[i+1,j+1,:,:] += self.radius_slices[i][j]._h_bincontent
				self.radius._h_squaredweights[i+1,j+1,:,:] += self.radius_slices[i][j]._h_squaredweights
				self.energy._h_bincontent[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_bincontent
				self.energy._h_squaredweights[i+1,j+1,:,:,:] += self.energy_slices[i][j]._h_squaredweights
		self.primary.flush()
		with tables.openFile(self.outfile, 'a') as hdf:
			if self.where != '/':
				hdf.createGroup('/', self.where[1:])
			dashi.histsave(self.primary, hdf, self.where, 'primary')
			dashi.histsave(self.multiplicity, hdf, self.where, 'multiplicity')
			dashi.histsave(self.radius, hdf, self.where, 'radius')
			dashi.histsave(self.energy, hdf, self.where, 'energy')

class FillWorker(object):
	def __init__(self, flux, binner, outfile, where='/'):
		self.flux = flux
		self.binner = binner
		self.outfile = outfile
		self.where = where
		import os
		with tables.openFile(self.outfile, 'a') as hdf:
			if self.where != '/':
				parent, node = os.path.split(self.where)
				hdf.createGroup(parent, node, createparents=True)
	def consume(self, primary, tracks, norm):
		energy = primary.energy
		zenith = primary.dir.zenith
		weight = self.flux(energy, primary.type)/norm
		self.binner.consume(tracks, energy, zenith, weight)
	def save(self, hdf, hist, where):
		# save in single precision
		hist._h_bincontent = hist._h_bincontent.astype(numpy.float32)
		hist._h_squaredweights = hist._h_squaredweights.astype(numpy.float32)
		dashi.histsave(hist, hdf, self.where, where)
	def finish(self):
		with tables.openFile(self.outfile, 'a') as hdf:
			for label in ('primary', 'multiplicity', 'radius', 'energy'):
				hist = getattr(self.binner, label).to_dashi()
				self.save(hdf, hist, label)

class MultiFiller(icetray.I3Module):
	def __init__(self, ctx):
		icetray.I3Module.__init__(self, ctx)
		self.AddOutBox("OutBox")
		self.AddParameter("Fluxes", "", None)
		self.AddParameter("GenerationSpectra", "", None)
		self.AddParameter("Outfile", "", None)
		self.AddParameter("MinDepth", "", 1.)
		self.AddParameter("MaxDepth", "", 5.)
		self.AddParameter("DepthSteps", "", 9)
	
	def Configure(self):
		
		fluxes = self.GetParameter("Fluxes")
		self.generator = self.GetParameter("GenerationSpectra")
		outfile = self.GetParameter("Outfile")
		import os
		if os.path.exists(outfile):
			os.unlink(outfile)
			
		self.nevents = 0
		
		def make_binner():
			return MuonGun.TrackBinner(self.GetParameter("MinDepth"), self.GetParameter("MaxDepth"), self.GetParameter("DepthSteps"))
		
		self.workers = dict()
		
		for label, weight in fluxes.items():
			worker = FillWorker(weight, make_binner(), outfile, '/%s' % label)
			self.workers[label] = worker
		
	def DAQ(self, frame):
		primary = frame['MCPrimary']
		tracks = frame['Tracks']
		norm = self.generator(primary.energy, primary.type, primary.dir.zenith)
		for worker in self.workers.values():
			worker.consume(primary, tracks, norm)
		
		self.nevents += 1
		if self.nevents % 10000 == 0:
			sys.stderr.write('%d events\n' % self.nevents)
	
	def Finish(self):
		for worker in self.workers.values():
			worker.finish()
		


class FastFiller(icetray.I3ConditionalModule):
	"""
	Tabulate tracks by depth, zenith angle, multiplicity, energy, and radius from shower axis
	"""
	def __init__(self, ctx):
		icetray.I3ConditionalModule.__init__(self, ctx)
		self.AddOutBox("OutBox")
		self.AddParameter("Weight", "", None)
		self.AddParameter("Outfile", "", None)
		self.AddParameter("MinDepth", "", 1.)
		self.AddParameter("MaxDepth", "", 5.)
		self.AddParameter("DepthSteps", "", 9)
	
	def Configure(self):
		
		self.binner = MuonGun.TrackBinner(self.GetParameter("MinDepth"), self.GetParameter("MaxDepth"), self.GetParameter("DepthSteps"))
		
		self.weighter = self.GetParameter("Weight")
		
		import os
		components = os.path.split(self.GetParameter("Outfile"))
		if os.path.splitext(components[-1])[1] in ('.hdf5', '.h5'):
			self.where = '/'
			self.outfile = os.path.join(*components)
		else:
			self.where = '/' + components[-1]
			self.outfile = os.path.join(*components[:-1])
		
		if os.path.exists(self.outfile):
			os.unlink(self.outfile)
			
		self.nevents = 0
		
	def DAQ(self, frame):
		# print 'woop!'
		primary = frame['MCPrimary']
		
		energy = primary.energy
		zenith = primary.dir.zenith
		weight = float(self.weighter(primary.type, energy, zenith))
		# weight = 1.
		
		self.binner.consume(frame['Tracks'], energy, zenith, weight)
		
		self.nevents += 1
		# if self.nevents % 10000 == 0:
		# 	print '%d events' % self.nevents
		
		# self.PushFrame(frame)
		
	def Finish(self):
		# return
		with tables.openFile(self.outfile, 'a') as hdf:
			if self.where != '/':
				hdf.createGroup('/', self.where[1:])
			dashi.histsave(self.binner.primary.to_dashi(), hdf, self.where, 'primary')
			dashi.histsave(self.binner.multiplicity.to_dashi(), hdf, self.where, 'multiplicity')
			dashi.histsave(self.binner.radius.to_dashi(), hdf, self.where, 'radius')
			dashi.histsave(self.binner.energy.to_dashi(), hdf, self.where, 'energy')

class CylinderWeight(object):
	def __init__(self, radius, length, depthbins, energy_spectrum):
		from numpy import pi, sqrt
		def diffarea_top(ct):
			return 2*pi**2*radius*(radius*ct);
		def diffarea_sides(ct):
			2*pi*radius*(2*sqrt(1-ct*ct))

class CylinderWeighter(object):
	"""
	In VOLUMECORR mode, CORSIKA generates showers with a zenith distribution proportional
	to the projected area of the cylindrical sampling surface. To convert counts of muons
	detected at the sampling surface back to fluxes, we have to weight each entry
	by 1/(dA/dcos(theta)), the differential projected area of the sampling surface
	into a plane perpendicular to the shower axis.
	"""
	def __init__(self, radius, height, timescale, depthbins):
		from numpy import pi, sqrt
		chi = 2*height/(pi*radius)
		# convert depths (in kmwe) to z coordinates
		z = I3Constants.SurfaceElev - I3Constants.OriginElev - (depthbins*I3Units.km)
		def diffarea(ct, r, l):
			"""
			differential projected area of an upright cylinder: dA_perp d\Omega/(dcos(theta))
			"""
			return 2*pi**2*r*(r*ct + (2*l/pi)*sqrt(1-ct**2))
		def bandarea(ct, r, l):
			return diffarea(ct, r, l) - diffarea(ct, r, 0)
		def intarea(r, l):
			"""
			projected area of an upright cylinder, integrated over the upper half-sphere
			"""
			return (pi**2*r*(r+l))
		class weighter(object):
			def __init__(self, f, timescale, r, l):
				self.f = f
				self.timescale = timescale
				self.r = r
				self.l = l
			def __call__(self, ct):
				return self.timescale*self.f(ct, self.r, self.l)
		self.weights = []
		for zhi, zlo in zip(z[:-1], z[1:]):
			print(zlo, zhi, end=' ')
			if zlo > height/2. or zhi < -height/2.:
				# outside the cylinder
				weight = None
			elif zhi > height/2. and zlo <= height/2.:
				# this layer includes the top of the cylinder
				sideheight = (height/2.-zlo)
				print(sideheight, end=' ')
				weight = weighter(diffarea, timescale, radius, sideheight)
			else:
				# a side layer has no top surface
				if zlo <= -height/2.:
					# XXX HACK: ucr's target cylinder is displaced slightly w.r.t 
					# the IceCube coordinate system. Adjust the effective area of
					# the bottom-most slice to compensate.
					sideheight = zhi+height/2. - 5.
				else:
					sideheight = zhi-zlo
				weight = weighter(bandarea, timescale, radius, sideheight)
				
			print('')
			self.weights.append(weight)
		
	def __call__(self, depthidx, zenith):
		wt = 1./self.weights[depthidx](numpy.cos(zenith))
		#print depthidx, self.weights[depthidx][0]
		return wt


class DFastFiller(FastFiller):
	def __init__(self, context):
		super(self, DFastFiller).__init__(context)
		from icecube.MuonGun import Cylinder
		self.cyl = Cylinder(1600, 800)
		
		mindepth = self.GetParameter("MinDepth")
		maxdepth = self.GetParameter("MaxDepth")
		steps = self.GetParameter("DepthSteps")
		self.depthbins = numpy.linspace(mindepth, maxdepth, steps)
		
	def DAQ(self, frame):
		from icecube.dataclasses import I3Constants
		primary = frame['MCPrimary']
		
		impact = self.cyl.intersection(primary.pos, primary.dir)
		d0 = I3Constants.SurfaceElev - I3Constants.OriginElev - (primary.pos.z + impact.first*primary.dir.z)
		di = numpy.floor(len(self.depthbins)*(d0 - mindepth)/(maxdepth-mindepth))
		
		energy = primary.energy
		zenith = primary.dir.zenith
		weight = float(self.weighter(energy, zenith))
		# weight = 1.
		
		self.binner.consume(frame['Tracks'], energy, zenith, weight)
		
		self.nevents += 1
