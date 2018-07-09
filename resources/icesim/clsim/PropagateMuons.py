#!/usr/bin/env python

from os.path import expandvars
from icecube import icetray, dataclasses
from I3Tray import I3Units

import math

def MakePropagator(
    radius=800*I3Units.m,
    length=1600*I3Units.m,
    particleType=dataclasses.I3Particle.ParticleType.MuMinus,
    impl='proposal',
    mediadef=None):
        """
        Create a muon propagator service.

        :param radius: radius of the target cylinder
        :param length: full height of the target cylinder
        :param impl: if "mmc", use MMC, otherwise use PROPOSAL
        :param mediadef: path to MMC media definition file
        """

        seed = 12345 # fixed seed, RNG will be re-set by the propagator module

        if impl.lower() == 'mmc':
            from icecube import sim_services, c2j_icetray, mmc_icetray, icetray
            if mediadef is None:
                mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/mediadef')
            jvmOpts = icetray.vector_string()    # fill this with parameters passed directly to the JavaVM
            jvmOpts.append(expandvars("-Djava.class.path=$I3_BUILD/lib/mmc.jar"))
            jvmOpts.append("-Xms256m")
            jvmOpts.append("-Xmx512m")
            jvmOpts.append("-XX:-HeapDumpOnOutOfMemoryError")

            # Now create the MMC propagators, but first *all* of the options must be set here.
            # There's no special options added behind the scenes.  This is much more flexible.
            #  Below are the standard options.  To interpret them see the MMC docs.
            mmcOpts = "-romb=5 -raw -user -sdec -time -lpm -bs=1 -ph=3 -bb=2 -sh=2 -frho -cont "
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
                mediadef=expandvars('$I3_BUILD/PROPOSAL/resources/config.json')
            return PROPOSAL.I3PropagatorServicePROPOSAL(
                config_file=mediadef)
        else:
            raise RuntimeError("unknown propagator: %s" % impl)

@icetray.traysegment
def PropagateMuons(tray, name,
    RandomService = None,
    CylinderRadius = 1200.*I3Units.m,
    CylinderLength = 1700.*I3Units.m,
    ):

    from I3Tray import I3Units

    from icecube import icetray, dataclasses, phys_services, sim_services, simclasses
    from icecube import cmc

    muPropagator      = MakePropagator(radius=CylinderRadius, length=CylinderLength, particleType=dataclasses.I3Particle.ParticleType.MuMinus)
    tauPropagator     = MakePropagator(radius=CylinderRadius, length=CylinderLength, particleType=dataclasses.I3Particle.ParticleType.TauMinus)
    cascadePropagator = cmc.I3CascadeMCService(phys_services.I3GSLRandomService(1)) # dummy RNG

    # set up propagators
    propagators = sim_services.I3ParticleTypePropagatorServiceMap()
    for pt in 'MuMinus', 'MuPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = muPropagator
    for pt in 'TauMinus', 'TauPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = tauPropagator
    for pt in 'DeltaE', 'Brems', 'PairProd', 'NuclInt', 'Hadrons', 'EMinus', 'EPlus':
        propagators[getattr(dataclasses.I3Particle.ParticleType, pt)] = cascadePropagator

    tray.AddModule('I3PropagatorModule', name+'_propagator',
        PropagatorServices=propagators,
        RandomService=RandomService,
        RNGStateName="RNGState",
        InputMCTreeName="I3MCTree_preMuonProp",
        OutputMCTreeName="I3MCTree")

    # add empty MMCTrackList objects for events that have none
    def addEmptyMMCTrackList(frame):
        if "MMCTrackList" not in frame:
            frame["MMCTrackList"] = simclasses.I3MMCTrackList()
    tray.AddModule(addEmptyMMCTrackList, name+'_addEmptyMMCTrackList',
        Streams=[icetray.I3Frame.DAQ])


if __name__=="__main__":
    from optparse import OptionParser
    from os.path import expandvars

    usage = "usage: %prog [options] inputfile"
    parser = OptionParser(usage)
    parser.add_option("-o", "--outfile",default="test_muons_propagated.i3",
                      dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
    parser.add_option("-s", "--seed",type="int",default=12345,
                      dest="SEED", help="Initial seed for the random number generator")
    parser.add_option("-r", "--runnumber", type="int", default=1,
                      dest="RUNNUMBER", help="The run number for this simulation")

    # parse cmd line args, bail out if anything is not understood
    (options,args) = parser.parse_args()
    if len(args) != 0:
            crap = "Got undefined options:"
            for a in args:
                    crap += a
                    crap += " "
            parser.error(crap)

    from I3Tray import *
    import os
    import sys

    from icecube import icetray, dataclasses, dataio, phys_services

    tray = I3Tray()

    # set up a random number generator
    randomService = phys_services.I3SPRNGRandomService(
        seed = options.SEED*2,
        nstreams = 10000,
        streamnum = options.RUNNUMBER)

    # re-use the same RNG for modules that need it on the context
    tray.context['I3RandomService'] = randomService

    tray.AddSegment(PropagateMuons, "PropagateMuons",
        RandomService = randomService)

    tray.AddModule("I3Writer","writer",
        Filename = options.OUTFILE)


    tray.Execute()



