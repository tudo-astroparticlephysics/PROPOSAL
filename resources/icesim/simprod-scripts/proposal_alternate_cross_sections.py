#!/usr/bin/env python
######################################################################
#   
#   TraySegements and IceProdModules for production
#   NuGen with varying mmc/PROPOSAL cross-sections
#   
######################################################################
import os,sys
from os.path import expandvars
import math
import logging

# Load libraries 
from I3Tray import *
from icecube import (icetray, dataio, dataclasses, phys_services,
                     interfaces, simclasses, sim_services,
                     neutrino_generator, earthmodel_service,
                     PROPOSAL, cmc, phys_services)

from icecube.simprod import ipmodule
from icecube.simprod.util import BasicCounter


class NuGen(ipmodule.ParsingModule):
   """
   Wrapper class that runs NuGen
   
   Use like:
   python -c "from proposal_alternate_cross_sections import NuGen; NuGen().ExecuteOpts({})" --gcdfile=gcd.i3 --outputfile=nugen_out.i3 --summaryfile=summary.xml --mjd=56429 --nevents=100 --NuFlavor=NuMu --RunId=1000 --RNGSeed=12345 --RNGStream=0 --RNGNumberOfStreams=1 --bs=1 --ph=3 --bb=2 --sh=2 
   """

   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('outputfile','Output filename','')
        self.AddParameter('mjd','MJD for the GCD file','')
        self.AddParameter('nevents','Number of events','')
        self.AddParameter('elogmax','Maximum energy','9.0')
        self.AddParameter('elogmin','Minimun energy','1.0')
        self.AddParameter('gamma','Gamma index','2.0')
        self.AddParameter('NuFlavor','Neutrino Flavor','')
        self.AddParameter('UseCSS','Use CSS cross total sections',False)
        self.AddParameter("RunId","Configure run ID",0)
        self.AddParameter("RNGSeed","RNG seed",0)
        self.AddParameter("RNGStream","RNG stream number",0)
        self.AddParameter("RNGNumberOfStreams","Number of RNG streams",1)
        self.AddParameter('bs','Bremsstrahlung Parametrization',1)
        self.AddParameter('ph','Photonuclear Parametrization Family',3)
        self.AddParameter('bb','Photonuclear Parametrization',2)
        self.AddParameter('sh','Nuclear Shadowing Parametrization',2)

   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        # Instantiate a tray 
        tray = I3Tray()

        randomService = phys_services.I3SPRNGRandomService(self.rngseed, self.rngnumberofstreams, self.rngstream)
        tray.context['I3RandomService'] = randomService

        summary = dataclasses.I3MapStringDouble()
        tray.context["I3SummaryService"] = summary
         
        tray.AddSegment(NuGenTraySegment,"nugen",
             gcdfile=self.gcdfile,
             nuflavor=self.nuflavor,
             mjd=self.mjd, 
             nevents=self.nevents,
             cssflag=self.usecss, 
             runId=self.runid,
             randomService=randomService,
             elogmax=self.elogmax,
             elogmin=self.elogmin,
             gamma=self.gamma,
             stats=stats,
             bs=self.bs,
             ph=self.ph,
             bb=self.bb,
             sh=self.sh)

        tray.AddModule("I3Writer","writer")(
            ("filename",self.outputfile),
            ("streams",[icetray.I3Frame.DAQ])
        )

        

        # Execute the Tray
        tray.Execute()
        

        # Free memory
        del tray
        return 0


class CorsikaGenerator(ipmodule.ParsingModule):
   """
   Wrapper class that runs Corsika
   
   Use like:
   I3_TOPDIR=$PWD python -c "from proposal_alternate_cross_sections import CorsikaGenerator; CorsikaGenerator().ExecuteOpts({})" --procnum=0 --seed=12345 --nproc=1 --gcdfile=gcd.i3 --outputfile=out_corsika.i3 --summaryfile=summary.xml --mjd=56429 --RunId=1 --nshowers=1000 --RunCorsika --bs=1 --ph=2 --bb=3 --sh=2
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('nshowers','Number of generated CR showers',1)
        self.AddParameter('procnum','process number',0)
        self.AddParameter('seed','RNG seed',1)
        self.AddParameter('nproc','Number of processes for (RNG)',1)
        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('outputfile','Output filename','corsika.i3')
        self.AddParameter('inputfile',"Input filename (only if you are not generating file)",'')
        self.AddParameter('RunCorsika','Run CORSIKA or only generate INPUTS file',False)
        self.AddParameter('mjd','MJD to set the event times',55697)
        self.AddParameter('atmospheres','Atmospheric models',[11,12,13,14]) # mar, jul, oct, dec
        self.AddParameter('eslope','CR spectral index (only if ranpri=0)',-2.7)
        self.AddParameter('ranpri','CR spectrum: 0=individual-nuclei, 1=Wiebel-Sooth, 2=Hoerandel, 3=5-component',2)
        self.AddParameter('pnorm','5-component relative contribution H,He,N,Al,Fe',[10.,5.,3.,2.,1.])
        self.AddParameter('pgam','5-component spectral indices H,He,N,Al,Fe',[2.0,2.0,2.0,2.0,2.0])
        self.AddParameter('dCorsikaVersion','version of dCorsika to run','v6900')
        self.AddParameter("RunId","Configure run ID",0)
        self.AddParameter("locut","Enables skew angle cutfoff",1.58)

        self.AddParameter("dslope","Change in spectral index",0)
        self.AddParameter("eprimarymax",'CR max energy',1e5)
        self.AddParameter("eprimarymin","CR min energy",600)
        self.AddParameter("fluxsum","", 0.131475115)
        self.AddParameter("length","",1600)
        self.AddParameter("model","","sibyll")
        self.AddParameter("oversampling","",1)
        self.AddParameter("radius","",800)

        self.AddParameter('cthmin','Min theta of injected cosmic rays',0.0)  
        self.AddParameter('cthmax','Max theta of injected cosmic rays',89.99)  
  
        self.AddParameter('ecuts1','hadron min energy (see corsika docs)',273)  
        self.AddParameter('ecuts2','muon min energy (see corsika docs)',273)  
        self.AddParameter('ecuts3','electron min energy (see corsika docs)',0.003)  
        self.AddParameter('ecuts4','photon min energy (see corsika docs)',0.003)  
        self.AddParameter("CutoffType","Sets SPRIC=T (EnergyPerNucleon) or F (EnergyPerParticle) ","EnergyPerNucleon")
        self.AddParameter("RepoURL","URL of repository containing corsika tarballs","http://convey.icecube.wisc.edu/data/sim/sim-new/downloads")

        self.AddParameter('bs','Bremsstrahlung Parametrization',1)
        self.AddParameter('ph','Photonuclear Parametrization Family',3)
        self.AddParameter('bb','Photonuclear Parametrization',2)
        self.AddParameter('sh','Nuclear Shadowing Parametrization',2)

 
   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        from I3Tray import I3Tray
        from icecube import icetray,phys_services, dataio, dataclasses
        
   
        if self.cutofftype == "EnergyPerNucleon" : self.spric = True
        elif cutofftype == "EnergyPerParticle" : self.spric = False
        else: raise Exception, "Undefined CutoffType %s" % cutoff_typ
        print self.atmospheres
        cors = configure_corsika(self)
        cors.f2k = "F"

        print self.runcorsika, "runcorsika",self.inputfile
        if self.runcorsika and not self.inputfile:
           # Run corsika
           retval = cors.Execute(stats) 

           if retval: 
              raise Exception, "dCorsika exited with return value %s" % retval
        elif self.inputfile:
           cors.outfile = self.inputfile
           print self.runcorsika, "runcorsika",self.inputfile
        elif not self.runcorsika:
             cors.configure()
             cors.write_steering()
    
        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        summary = dataclasses.I3MapStringDouble()
        tray.context["I3SummaryService"] = summary
         
        randomService = phys_services.I3SPRNGRandomService(self.seed, self.nproc, self.procnum)#, oustatefile="rng.state")
        tray.context["I3RandomService"] = randomService

        # Configure TraySegment that actually does stuff
        tray.AddSegment(CorsikaReaderTraySegment,"corsika-reader",
              gcdfile=self.gcdfile, mjd=self.mjd, runId=self.runid, randomService=randomService,
              oversampling=self.oversampling, fluxsum=self.fluxsum, cors=cors, stats={}, 
             CylinderLength=self.length,
             CylinderRadius=self.radius,
             bs=self.bs,
             ph=self.ph,
             bb=self.bb,
             sh=self.sh)

        tray.AddModule(BasicCounter,"count_g", Streams = [icetray.I3Frame.DAQ], 
              name = "Generated Events", Stats = stats)
       
        tray.AddModule("I3Writer","writer", filename = self.outputfile, streams =[icetray.I3Frame.DAQ] )

        

        # Execute the Tray
        tray.Execute()
        

        del tray
        return 0




class MuonGunGenerator(ipmodule.ParsingModule):
   """
   Wrapper class that runs MuonGun
   
   Use like:
   python -c "from proposal_alternate_cross_sections import MuonGunGenerator; MuonGunGenerator().ExecuteOpts({})" --procnum=0 --seed=12345 --nproc=1 --gcdfile=gcd.i3 --outputfile=out_corsika.i3 --summaryfile=summary.xml --mjd=56429 --RunId=1 --nshowers=1000 --bs=1 --ph=2 --bb=3 --sh=2
   """
   def __init__(self):
        ipmodule.ParsingModule.__init__(self)

        self.AddParameter('nshowers','Number of generated CR showers',1)
        self.AddParameter('procnum','process number',0)
        self.AddParameter('seed','RNG seed',1)
        self.AddParameter('nproc','Number of processes for (RNG)',1)
        self.AddParameter('gcdfile','GeoCalibDetStatus filename','')
        self.AddParameter('outputfile','Output filename','corsika.i3')
        self.AddParameter('mjd','MJD to set the event times',55697)
        self.AddParameter('model','MuonGun Model','Hoerandel5_atmod12_SIBYLL')
        self.AddParameter("RunId","Configure run ID",0)

        self.AddParameter("gamma",'CR max energy',2)
        self.AddParameter("eprimarymax",'CR max energy',1e5)
        self.AddParameter("eprimarymin","CR min energy",600)
        self.AddParameter("length","",1600)
        self.AddParameter("radius","",800)
        self.AddParameter('cthmin','Min theta of injected cosmic rays',0.0)  
        self.AddParameter('cthmax','Max theta of injected cosmic rays',89.99)  

        self.AddParameter('bs','Bremsstrahlung Parametrization',1)
        self.AddParameter('ph','Photonuclear Parametrization Family',3)
        self.AddParameter('bb','Photonuclear Parametrization',2)
        self.AddParameter('sh','Nuclear Shadowing Parametrization',2)
 
   def Execute(self,stats):
        if not ipmodule.ParsingModule.Execute(self,stats): return 0

        import icecube.icetray
        from icecube import icetray,dataclasses, phys_services, dataio
        from icecube.MuonGun import load_model, StaticSurfaceInjector, Cylinder, OffsetPowerLaw
        from icecube.MuonGun.segments import GenerateBundles

        from I3Tray import I3Tray, I3Units
        
        modelstr = self.GetParameter('model')
        emax   = self.GetParameter("eprimarymax")
        emin   = self.GetParameter("eprimarymin")
        length = self.GetParameter("length")
        radius = self.GetParameter("radius")
        cthmin = self.GetParameter('cthmin')
        cthmax = self.GetParameter('cthmax')
        gamma  = self.GetParameter('gamma')

        seed   = self.GetParameter('seed')
        nproc  = self.GetParameter('nproc')
        procnum  = self.GetParameter('procnum')
        gcdfile = self.GetParameter('gcdfile')
    
        # Instantiate a tray 
        tray = I3Tray()

        # Configure IceTray services
        summary = dataclasses.I3MapStringDouble()
        tray.context["I3SummaryService"] = summary
         
        randomService = phys_services.I3SPRNGRandomService(seed, nproc, procnum)
        tray.context["I3RandomService"] = randomService


        # Use Hoerandel as a template for generating muons
        model = load_model(modelstr)
        # Generate only single muons, no bundles
        model.flux.max_multiplicity = 1
        # Center the sampling surface on the barycenter of IC79 strings
        #surface = Cylinder(1600*I3Units.m, 800*I3Units.m, dataclasses.I3Position(31.25, 19.64, 0))
        surface = Cylinder(length*I3Units.m, radius*I3Units.m)
        # Draw energies from an E^-2 power law broken at 1 TeV, from 10 TeV to 10 PeV
        spectrum = OffsetPowerLaw(gamma, 1*I3Units.TeV, 10*I3Units.TeV, 10*I3Units.PeV)
        # Set up the generator. This gets stored in a special frame for later reference
        generator = StaticSurfaceInjector(surface, model.flux, spectrum, model.radius)

        # Configure TraySegment that actually does stuff
        tray.AddSegment(GenerateBundles,"corsika-background", Generator=generator, 
              RunNumber=self.runid, NEvents=self.nshowers, 
              GCDFile=self.gcdfile,
              FromTime=dataclasses.I3Time(self.mjd),
              ToTime=dataclasses.I3Time(self.mjd))

        def renameMCTree(frame):
          mctree = frame["I3MCTree"]
          del frame["I3MCTree"]
          frame["I3MCTree_preMuonProp"] = mctree
        tray.AddModule(renameMCTree, "renameMCTree", Streams=[icetray.I3Frame.DAQ])
        tray.AddSegment(PropagateMuons,'propagator', RandomService=randomService, 
             CylinderLength=self.length,
             CylinderRadius=self.radius,
             bs=self.bs,
             ph=self.ph,
             bb=self.bb,
             sh=self.sh)

        tray.AddModule(BasicCounter,"count_g", Streams = [icetray.I3Frame.DAQ], 
              name = "Generated Events", Stats = stats)
       
        tray.AddModule("I3Writer","writer", filename = self.outputfile, streams =[icecube.icetray.I3Frame.DAQ] )

        

        # Execute the Tray
        tray.Execute()
        

        del tray
        return 0

def MakePropagator(
    radius=800*I3Units.m,
    length=1600*I3Units.m,
    particleType=dataclasses.I3Particle.ParticleType.MuMinus,
    impl='proposal',
    mediadef=None,
    ):
        """
        Create a muon propagator service.

        :param radius: radius of the target cylinder
        :param length: full height of the target cylinder
        :param impl: if "mmc", use MMC, otherwise use PROPOSAL
        :param mediadef: path to MMC media definition file
        """

        seed = 12345 # fixed seed, RNG will be re-set by the propagator module

        if impl.lower() == 'mmc':
            if mediadef is None:
                mediadef=expandvars('$I3_BUILD/mmc-icetray/resources/mediadef')

            from icecube import sim_services, c2j_icetray, mmc_icetray, icetray
            jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
            jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
            jvmOpts.append("-Xms256m")
            jvmOpts.append("-Xmx512m")
            jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")

            # Now create the MMC propagators, but first *all* of the options must be set here.
            # There's no special options added behind the scenes.  This is much more flexible.
            #  Below are the standard options.  To interpret them see the MMC docs.
            mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=%d -ph=%d -bb=%d -sh=%d -frho -cont "%(bs,ph,bb,sh)
            mmcOpts += expandvars("-tdir=$I3_BUILD/mmc-icetray/resources ")
            mmcOpts += expandvars("-mediadef=%s " % mediadef)
            mmcOpts += "-radius=%d " % radius
            mmcOpts += "-length=%d " % length
            mmcOpts += "-seed=%d " % seed

            if particleType in [dataclasses.I3Particle.ParticleType.MuMinus, dataclasses.I3Particle.ParticleType.MuPlus]:
                pass # no extra options necessary
            elif particleType in [dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.TauPlus]:
                mmcOpts += "tau "
            else:
                raise RuntimeError("Cannot propagate particle type!")

            jvm = c2j_icetray.I3JavaVM(jvmOpts)
            return mmc_icetray.I3PropagatorServiceMMC(jvm,mmcOpts)
        elif impl.lower() == 'proposal':
            from icecube import sim_services, PROPOSAL
            # in PROPOSAL everything can be defined in the configuration file
            if mediadef is None:
                mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/config_icesim.json')
            return PROPOSAL.I3PropagatorServicePROPOSAL(
                config_file=mediadef)
        else:
            raise RuntimeError("unknown propagator: %s" % impl)

@icetray.traysegment
def PropagateMuons(tray, name,
    RandomService = None,
    CylinderRadius = 1200.*I3Units.m,
    CylinderLength = 1700.*I3Units.m,
    SaveState = True,
    InputMCTreeName="I3MCTree_preMuonProp",
    OutputMCTreeName="I3MCTree",
    **kwargs
    ):

    from I3Tray import I3Units

    from icecube import icetray, dataclasses, phys_services, sim_services, simclasses
    from icecube import cmc

    muPropagator      = MakePropagator(radius=CylinderRadius, length=CylinderLength, particleType=dataclasses.I3Particle.ParticleType.MuMinus,**kwargs)
    tauPropagator     = MakePropagator(radius=CylinderRadius, length=CylinderLength, particleType=dataclasses.I3Particle.ParticleType.TauMinus,**kwargs)
    cascadePropagator = cmc.I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG

    # set up propagators
    propagators = sim_services.I3ParticleTypePropagatorServiceMap()
    for pt in 'MuMinus', 'MuPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = muPropagator
    for pt in 'TauMinus', 'TauPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = tauPropagator
    for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cascadePropagator

    if SaveState:
        additionalArgs = dict(RNGStateName="RNGState")
    else:
        additionalArgs = dict()

    tray.AddModule('I3PropagatorModule', name+'_propagator',
        PropagatorServices=propagators,
        RandomService=RandomService,
        InputMCTreeName=InputMCTreeName,
        OutputMCTreeName=OutputMCTreeName,
        **additionalArgs)

    # add empty MMCTrackList objects for events that have none
    def addEmptyMMCTrackList(frame):
        if "MMCTrackList" not in frame:
            frame["MMCTrackList"] = simclasses.I3MMCTrackList()
    tray.AddModule(addEmptyMMCTrackList, name+'_addEmptyMMCTrackList',
        Streams=[icetray.I3Frame.DAQ])


class DAQCounter(icetray.I3Module):
    count   = 0
    nevents = 0

    def __init__(self,ctx):
        icetray.I3Module.__init__(self,ctx)
        self.AddParameter("NEvents","name of event counter",self.nevents)
        self.AddOutBox("OutBox");

    def Configure(self):
        self.nevents = self.GetParameter("NEvents")

    def DAQ(self,frame):
        self.count += 1
        if self.count > self.nevents: 
           self.RequestSuspension()
        self.PushFrame(frame)


    
def NuGenTraySegment(tray,name,
    gcdfile, mjd, nevents, nuflavor, cssflag, runId, randomService,
    zenithmin=0.0,
    zenithmax=180.0,
    elogmax=9.0,
    elogmin=1.0,
    gamma=2.0,
    injectionradius=1200.0,
    distanceentrance=1000.0,
    distanceexit=1000.0,
    CrustModel="PREM_mmc",           # density profile 
    PropagationWeightMode="AutoDetect",     
    stats={},
    **kwargs):
    """
    Tray Segment for basic NuGen simulation (generation)
    """

    # Configure IceTray modules 
    tray.AddModule("I3InfiniteSource", "source",
        prefix = gcdfile, 
        stream = icetray.I3Frame.DAQ ) 

    time = dataclasses.I3Time()
    time.set_mod_julian_time(int(mjd), 0, 0)

    tray.AddModule("I3MCEventHeaderGenerator","time-gen",
        Year    = time.utc_year,
        DAQTime = time.utc_daq_time,
        RunNumber = runId,
        IncrementEventID = True
    )
     
    tray.AddModule(DAQCounter,"counter3")(
        ("nevents",int(nevents)),
    )
   
    # Use NoWeight if Tau Regeneration.
    if PropagationWeightMode=="AutoDetect":
       if nuflavor == "NuTau":
          PropagationWeightMode="NoWeight"
       else: 
          PropagationWeightMode="NCGRWeighted"

    earthModelService = name+"_earthmodel"
    steering = name+"_steering"
    tray.AddService("I3EarthModelServiceFactory", earthModelService,
        ServiceName=earthModelService,
        DetectorDepth=1950*I3Units.m,       # old setting is 1945m
        EarthModel=CrustModel,              # density profile 
        )
    tray.AddService("I3NuGSteeringFactory", steering, ServiceName=steering,
        EarthModelServiceName=earthModelService,
        MuonRangeOpt=0.0*I3Units.m,         # use old Muon Range calculation
        )
          
    tray.AddModule("I3NeutrinoGenerator","neutrino",
        nevents=int(nevents),
        NeutrinoFlavor=nuflavor,
        zenithMin=zenithmin*I3Units.deg,
        zenithMax=zenithmax*I3Units.deg,
        azimuthMin=0.*I3Units.deg,
        azimuthMax=360.*I3Units.deg,
        EnergyMinLog=float(elogmin)*I3Units.GeV,
        EnergyMaxLog=float(elogmax)*I3Units.GeV,
        GammaIndex=float(gamma),
        injectionRadius=injectionradius*I3Units.m,
        distanceExit=distanceexit*I3Units.m,
        distanceEntrance=distanceentrance*I3Units.m,
        NoInEarthPropagation=False,
        PropagationWeightMode=PropagationWeightMode,
        AutomaticMuonFinalVolume=True,
        #TablesDir=expandvars("$PWD/anis.v1.8.2/data"),
        EarthModelServiceName = earthModelService,
    )

    if cssflag:
       tray.SetParameter("neutrino", "SigmaFileCCNu","/css/css_cc_nu.data")
       tray.SetParameter("neutrino", "SigmaFileCCNuBar","/css/css_cc_nubar.data")
       tray.SetParameter("neutrino", "SigmaFileNCNu","/css/css_nc_nu.data")
       tray.SetParameter("neutrino", "SigmaFileNCNuBar","/css/css_nc_nubar.data")

    def renameMCTree(frame):
        mctree = frame["I3MCTree"]
        del frame["I3MCTree"]
        frame["I3MCTree_preMuonProp"] = mctree
    tray.AddModule(renameMCTree, name+"_renameMCTree", Streams=[icetray.I3Frame.DAQ])
    tray.AddSegment(PropagateMuons,'propagator', RandomService=randomService,**kwargs) 

    tray.AddModule(BasicCounter,"count_g", 
                    Streams = [icetray.I3Frame.DAQ] ,
                    name="Generated Events",
                    Stats=stats)



def SetAtmosphereWeight(frame, atmosphere = 13):
      if frame.Has("CorsikaWeightMap"):
         weightmap = frame["CorsikaWeightMap"]
         weightmap["Atmosphere"] = atmosphere
         frame.Delete("CorsikaWeightMap")
         frame.Put("CorsikaWeightMap",weightmap)
      return True

def weighted_choice(weights,rng):
    total     = sum(weights)
    intervals = [ sum(weights[0:i]) for i in range(len(weights)+1) ]
    r = rng.uniform(0,total)
    for i in range(len(intervals)):
        if r >= intervals[i] and r < intervals[i+1]:
           return i
    return 0


def CorsikaReaderTraySegment(tray,name,
      gcdfile, mjd, runId, fluxsum, randomService,cors, oversampling, stats={},**kwargs):

      """
      Tray Segment for basic CORSIKA simulation (generation)
      """

      # Load libraries 
      from icecube import icetray, dataclasses, interfaces, phys_services
      from icecube import simclasses, sim_services, corsika_reader, PROPOSAL,ucr_icetray
      from icecube.simprod.generators.weights import Corsika5CompWeightModule, CorsikaWeightModule, PolygonatoWeightModule

      print "Readig: ",cors.outfile
      tray.AddModule('I3CORSIKAReader','corsika', 
                     FilenameList=[cors.outfile], Prefix=gcdfile, NEvents=cors.nevents)

      if cors.ranpri == 3: # 5-component case

        tray.AddModule(Corsika5CompWeightModule,"5compCorsikaweight",
          name                       = "CorsikaWeightMap",
          nevents                    = cors.nevents,
          spric                      = cors.spric == 'T',
          ThetaMin                   = cors.cthmin*I3Units.degree,
          ThetaMax                   = cors.cthmax*I3Units.degree,
          cylinderLength             = cors.length*I3Units.meter,
          cylinderRadius             = cors.radius*I3Units.meter,
          energyprimarymin           = cors.emin*I3Units.GeV,
          energyprimarymax           = cors.emax*I3Units.GeV,

          PrimaryNormalizationH      =  cors.pnormH,
          PrimaryNormalizationHe     =  cors.pnormHe,
          PrimaryNormalizationCNO    =  cors.pnormN,
          PrimaryNormalizationMgAlSi =  cors.pnormAl,
          PrimaryNormalizationFe     =  cors.pnormFe,

          PrimarySpectralIndexH      = -cors.pgamH,
          PrimarySpectralIndexHe     = -cors.pgamHe,
          PrimarySpectralIndexCNO    = -cors.pgamN,
          PrimarySpectralIndexMgAlSi = -cors.pgamAl,
          PrimarySpectralIndexFe     = -cors.pgamFe,

        )
        tray.AddModule(PolygonatoWeightModule,"polygonato")

      else:  # Standard CORSIKA
        tray.AddModule(CorsikaWeightModule,"corsikaweight",
        spric            = cors.spric == 'T',
        nevents          = cors.nevents,
        ThetaMin         = cors.cthmin*I3Units.degree,
        ThetaMax         = cors.cthmax*I3Units.degree,
#        OverSampling     = oversampling,  #need to add this to  weighting
        )
# Polygonato needs to be fixed to include other particle types besided
# the 5-comp
#      tray.AddModule(PolygonatoWeightModule,"polygonato")

      def renameMCTree(frame):
        mctree = frame["I3MCTree"]
        del frame["I3MCTree"]
        frame["I3MCTree_preMuonProp"] = mctree
      tray.AddModule(renameMCTree, name+"_renameMCTree", Streams=[icetray.I3Frame.DAQ])
      tray.AddSegment(PropagateMuons,'propagator', RandomService=randomService,**kwargs) 

      tray.AddModule(SetAtmosphereWeight,"atmo", atmosphere = cors.atmod,
            Streams = [icetray.I3Frame.DAQ] ) 



def configure_corsika(params):
    """
    Configure corsika from paramters
    """
    from icecube.simprod.generators import dcorsika

    atmod = params.atmospheres[params.procnum % len(params.atmospheres)]

    cors = dcorsika.dCorsika()

    cutoff_type  = params.cutofftype
    if cutoff_type == "EnergyPerNucleon" : spric = True
    elif cutoff_type == "EnergyPerParticle" : spric = False
    else: raise Exception, "Undefined CutoffType %s" % cutoff_type

    cors.version  = params.dcorsikaversion
    cors.platform = ''
    if cors.version.endswith('5comp'):
       cors.platform = os.uname()[4]

    if params.locut:
       cors.locut = "T %f" % params.locut
    else:
       cors.locut = "F 0"
    cors.ectus1  = params.ecuts1
    cors.ectus2  = params.ecuts2
    cors.ectus3  = params.ecuts3
    cors.ectus4  = params.ecuts4
    cors.emin    = params.eprimarymin
    cors.emax    = params.eprimarymax

    cors.cthmin  = params.cthmin
    cors.cthmax  = params.cthmax
    cors.seed    = params.seed
    cors.radius  = params.radius
    cors.length  = params.length
    cors.atmod   = atmod
    cors.runnum  = params.runid % 1000000
    cors.cache   = 1
    cors.nevents = params.nshowers
    cors.dslope  = params.dslope
    cors.eslope  = params.eslope
    cors.ranpri  = params.ranpri
    cors.donkg   = 0
    cors.doegs   = 0
    cors.outfile = "CORSIKA%06u.gz" % params.runid
    cors.outdir  = expandvars("$PWD/")
    cors.topdir  = expandvars("$PWD/")
    cors.logfile = cors.outfile+".log"
    cors.url     = params.repourl
    cors.depth   = 1950.
    cors.model   = params.model
    cors.spric   = spric
    cors.f2k     = 'T'
    if params.ranpri == 3:
       if len(params.pnorm) != 5:
          raise Exception, "corsika: list of pnorm indices is the wrong size: %s!!!!" % len(params.pnorm)
       if len(params.pgam) != 5:
          raise Exception, "corsika: list of pgam indices is the wrong size: %s!!!!" % len(params.pgam)

       cors.pnormH  = params.pnorm[0]
       cors.pnormHe = params.pnorm[1]
       cors.pnormN  = params.pnorm[2]
       cors.pnormAl = params.pnorm[3]
       cors.pnormFe = params.pnorm[4]

       cors.pgamH   = params.pgam[0]
       cors.pgamHe  = params.pgam[1]
       cors.pgamN   = params.pgam[2]
       cors.pgamAl  = params.pgam[3]
       cors.pgamFe  = params.pgam[4]

    return cors
    
                    
    
if __name__ == '__main__':
    n = NuGen()
    n.ExecuteOpts({})
