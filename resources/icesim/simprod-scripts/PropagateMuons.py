"""Tray segments for muon propagation

"""
import os

import icecube
import icecube.icetray
import icecube.dataclasses
import icecube.phys_services
import icecube.sim_services
import icecube.simclasses
import icecube.cmc
import icecube.PROPOSAL

default_media_definition = os.path.expandvars(
    "$I3_BUILD/PROPOSAL/resources/config_icesim.json")


@icecube.icetray.traysegment
def PropagateMuons(tray, name,
                   RandomService=None,
                   CylinderRadius=800.,
                   CylinderLength=1600.,
                   SliceMuons=False,
                   SaveState=True,
                   InputMCTreeName="I3MCTree_preMuonProp",
                   OutputMCTreeName="I3MCTree",
                   **kwargs):
    r"""Propagate muons.

    This segment propagates muons through ice with ``PROPOSAL``; it
    simulates lepton decays and energy losses due to ionization,
    bremsstrahlung, photonuclear interactions, and pair production.

    :param I3RandomService RandomService:
        Random number generator service
    :param float CylinderRadius:
        Radius of the target volume in m
    :param float CylinderLength:
        Full height of the target volume in m
    :param bool SliceMuons:
        Emit constant-energy track slices in addition to stochastic losses
        (similar to the output of I3MuonSlicer)
    :param bool SaveState:
        If set to `True`, store the state of the supplied RNG.
    :param str InputMCTree:
        Name of input :cpp:class:`I3MCTree` frame object
    :param str OutputMCTree:
        Name of output :cpp:class:`I3MCTree` frame object
    :param \**kwargs:
        Additional keyword arguments are passed to
        :func:`icecube.simprod.segments.make_propagator`.

    """
    
    kwargs.update(dict(CylinderRadius=CylinderRadius, CylinderLength=CylinderLength, EmitTrackSegments=SliceMuons))
    propagators = make_standard_propagators(**kwargs)

    # Set up propagators.
    if "I3ParticleTypePropagatorServiceMap" in tray.context:
        propagator_map = tray.context["I3ParticleTypePropagatorServiceMap"]
        for k, v in propagators.items():
            propagator_map[k] = v
    else:
        propagator_map = propagators

    if SaveState:
        rng_state = "RNGState"
    else:
        rng_state = ""

    tray.AddModule("I3PropagatorModule", name+"_propagator",
                   PropagatorServices=propagator_map,
                   RandomService=RandomService,
                   InputMCTreeName=InputMCTreeName,
                   OutputMCTreeName=OutputMCTreeName,
                   RNGStateName=rng_state)

    # Add empty MMCTrackList objects for events that have none.
    def add_empty_tracklist(frame):
        if "MMCTrackList" not in frame:
            frame["MMCTrackList"] = icecube.simclasses.I3MMCTrackList()
        return True

    tray.AddModule(add_empty_tracklist, name+"_add_empty_tracklist",
                   Streams=[icecube.icetray.I3Frame.DAQ])

    return

def make_standard_propagators(CylinderRadius=800.,
                   CylinderLength=1600.,
                   EmitTrackSegments=False,
                   SplitSubPeVCascades=True,
                   MaxMuons=10,
                   **kwargs):
    """
    Set up standard propagators (PROPOSAL for muons and taus, CMC for cascades)
    Keyword arguments will be passed to I3PropagatorServicePROPOSAL
    """
    from icecube.icetray import I3Units
    
    muon_propagator = make_propagator(
        particle_type=icecube.dataclasses.I3Particle.ParticleType.MuMinus,
        cylinder_radius=CylinderRadius*icecube.icetray.I3Units.m,
        cylinder_length=CylinderLength*icecube.icetray.I3Units.m,
        emitTrackSegments=EmitTrackSegments,
        **kwargs)
    tau_propagator = make_propagator(
        particle_type=icecube.dataclasses.I3Particle.ParticleType.TauMinus,
        cylinder_radius=CylinderRadius*icecube.icetray.I3Units.m,
        cylinder_length=CylinderLength*icecube.icetray.I3Units.m,
        emitTrackSegments=EmitTrackSegments,
        **kwargs)

    cascade_propagator = icecube.cmc.I3CascadeMCService(
        icecube.phys_services.I3GSLRandomService(1))  # Dummy RNG
    cascade_propagator.SetEnergyThresholdSimulation(1*I3Units.PeV)
    if SplitSubPeVCascades:
        cascade_propagator.SetThresholdSplit(1*I3Units.TeV)
    else:
        cascade_propagator.SetThresholdSplit(1*I3Units.PeV)
    cascade_propagator.SetMaxMuons(MaxMuons)
    propagator_map =\
        icecube.sim_services.I3ParticleTypePropagatorServiceMap()

    for pt in "MuMinus", "MuPlus":
        key = getattr(icecube.dataclasses.I3Particle.ParticleType, pt)
        propagator_map[key] = muon_propagator

    for pt in "TauMinus", "TauPlus":
        key = getattr(icecube.dataclasses.I3Particle.ParticleType, pt)
        propagator_map[key] = tau_propagator

    for pt in "DeltaE", "Brems", "PairProd", "NuclInt", "Hadrons",\
              "EMinus", "EPlus":
        key = getattr(icecube.dataclasses.I3Particle.ParticleType, pt)
        propagator_map[key] = cascade_propagator
    
    return propagator_map

def make_propagator(
        media_definition=default_media_definition,
        particle_type=icecube.dataclasses.I3Particle.ParticleType.MuMinus,
        emitTrackSegments=False,
        cylinder_radius=800.,
        cylinder_length=1600.,
        bremsstrahlung="KelnerKokoulinPetrukhin",
        photonuclear_family="AbramowiczLevinLevyMaor",
        photonuclear="AbramowiczLevinLevyMaor97",
        nuclear_shadowing="Butkevich"):
    """Create a muon propagator service.

    :param str media_definition:
        Path the the media definition file
    :param ParticleType:
        Type of particle to propagate
    :param float cylinder_radius:
        Radius of the target volume in m
    :param float cylinder_length:
        Full height of the target volume in m
    :param str bremsstrahlung:
        Parametrization of the bremsstrahlung cross-section to use
    :param str photonuclear_family:
        Family of photonuclear cross-section parametrization to use
    :param str photonuclear:
        Specific edition of the photonuclear cross-section
        parametrization to use
    :param str nuclear_shadowing:
        Nuclear shadowing parametrization to use
    
    Note: the new PROPOSAL version only needs a configuration file
    everything (cross section, continuous loss output, etc.)
    can/should be defined there.

    """

    propagator_service = icecube.PROPOSAL.I3PropagatorServicePROPOSAL(
        config_file=media_definition)

    return propagator_service
