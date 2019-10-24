/**
 * class: I3PropagatorServicePROPOSAL
 *
 * Mainly a copy of:
 * Version $Id: I3PropagatorService.cxx 69179 2010-12-04 01:33:36Z dima $
 *
 *
 * Date 07 Jul 2012
 *
 * (c) 2003 IceCube Collaboration
 */


#include <sstream>
#include <unistd.h> // check for write permissions

#include <dataclasses/physics/I3Particle.h>

#include "PROPOSAL-icetray/Converter.h"
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"

using namespace std;
using namespace PROPOSAL;


// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL(std::string configfile,
    bool slice_tracks,
    I3Particle::ParticleType final_loss,
    double distance)
    : I3PropagatorService()
    , config_file_(PROPOSAL::Helper::ResolvePath(configfile))
    , proposal_service_()
    , final_stochastic_loss_(final_loss)
    , slice_tracks_(slice_tracks)
    , distance_to_propagate_(distance)
{
    if (config_file_.empty())
    {
        config_file_ = GetDefaultConfigFile();
    }

    log_info("Using configuration file: \"%s\"", config_file_.c_str());
    for (auto ptype : {I3Particle::MuMinus, I3Particle::MuPlus, I3Particle::TauMinus, I3Particle::TauPlus}) {
        RegisterParticleType(ptype);
    }
}

// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::~I3PropagatorServicePROPOSAL()
{
    if (rng_)
    {
        // rng_ = NULL;

        // Reset the random nuber generator to default at this place.
        // Reason: Clean up issues with boost::python objects and the python runtime
        // itself.
        // PROPOSAL::RandomGenerator::Get().SetDefaultRandomNumberGenerator();
    }
}

// ------------------------------------------------------------------------- //
std::string I3PropagatorServicePROPOSAL::GetDefaultConfigFile()
{
    const char* I3_BUILD = getenv("I3_BUILD");
    if (!I3_BUILD)
        log_fatal("$I3_BUILD is not set!");
    std::string s(I3_BUILD);

    return s + "/PROPOSAL/resources/config_icesim.json";
}

// ------------------------------------------------------------------------- //
void I3PropagatorServicePROPOSAL::SetRandomNumberGenerator(I3RandomServicePtr random)
{
    rng_                        = random;
    PROPOSAL::RandomGenerator::Get().SetI3RandomNumberGenerator(random);
}

// ------------------------------------------------------------------------- //
void I3PropagatorServicePROPOSAL::RegisterParticleType(I3Particle::ParticleType ptype)
{
    ParticleDef particle_def = particle_converter_.GeneratePROPOSALType(ptype);
    log_debug_stream(particle_def);

    // If ptype is not known an empty particle definition is returned (e.g. empty name)
    // Do not register unkown particles
    if (!particle_def.name.empty())
    {
        auto propagator = Propagator(particle_def, config_file_);
        // ensure that secondary output modes match settings
        for(auto sector : propagator.GetSectors()) {
            i3_assert(sector);
            sector->GetSectorDef().do_continuous_energy_loss_output = slice_tracks_;
            sector->GetSectorDef().only_loss_inside_detector = true;
        }
        proposal_service_.RegisterPropagator(propagator);
        log_debug_stream("particle registered");
    }
}

// ------------------------------------------------------------------------- //
std::vector<I3Particle> I3PropagatorServicePROPOSAL::Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr)
{
    // saying where we are
    log_debug("Entering I3PropagatorServicePROPOSAL::Propagate()");

    std::vector<I3Particle> daughters;

    log_trace("location type = %d", p.GetLocationType());
    if (p.GetLocationType() != I3Particle::InIce || 
        !proposal_service_.IsRegistered(particle_converter_.GeneratePROPOSALType(p.GetType())))
        return daughters;

    log_trace("particle to propagate:\n"
              "type/energy[GeV]/posx[m]/posy[m]/posz[m]/theta[deg]/phi[deg]/length[m]\n"
              "%d/%.2e/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f",
              p.GetType(),
              p.GetEnergy() / I3Units::GeV,
              p.GetPos().GetX() / I3Units::m,
              p.GetPos().GetY() / I3Units::m,
              p.GetPos().GetZ() / I3Units::m,
              p.GetZenith() / I3Units::deg,
              p.GetAzimuth() / I3Units::deg,
              p.GetLength() / I3Units::m);

    I3MMCTrackPtr mmcTrack = propagate(p, daughters);

    if (mmcTrack && frame)
    {
        if (!frame->Has("MMCTrackList"))
        {
            frame->Put("MMCTrackList", I3MMCTrackListPtr(new I3MMCTrackList));
        }

        I3MMCTrackListPtr trackList = frame->Get<I3MMCTrackListPtr>("MMCTrackList");
        i3_assert(trackList != NULL);

        trackList->push_back(*mmcTrack);
    }

    return daughters;
}

// ------------------------------------------------------------------------- //
I3MMCTrackPtr I3PropagatorServicePROPOSAL::propagate(I3Particle& p, std::vector<I3Particle>& daughters)
{
    Particle particle = particle_converter_.GeneratePROPOSALParticle(p);

    std::vector<DynamicData*> secondaries = proposal_service_.Propagate(particle, distance_to_propagate_);

    // get the propagated length of the particle
    double length = particle.GetPropagatedDistance();

    p.SetLength(length * I3Units::cm);
    log_trace(" length = %f cm ", length);

    int nParticles = secondaries.size();
    log_trace("nParticles = %d", nParticles);

    I3MMCTrackPtr mmcTrack;

    mmcTrack = GenerateMMCTrack(&particle);

    if (mmcTrack)
        mmcTrack->SetParticle(p);

    double currentEnergy = (mmcTrack && mmcTrack->GetEi() > 0) ?
        mmcTrack->GetEi() : p.GetEnergy();

    for (int i = 0; i < nParticles; i++)
    {
        // this should be a stochastic

        if (p.GetType() == I3Particle::EMinus || p.GetType() == I3Particle::EPlus)
        // || particle_def == ParticleType::Hadrons) //TODO(mario):  Wed 2017/10/25
        {
            if (p.GetShape() != I3Particle::TopShower)
            {
                log_fatal("The particle '%s' has no TopShower shape, but 'e-', 'e+' and 'Hadrons' need that. I don't "
                          "know why?",
                          p.GetTypeString().c_str());
            }
        }

        // Tomasz
        // in mmc the particle relationships are stored

        // this is not the particle you're looking for
        // move along...and add it to the daughter list
        daughters.push_back(particle_converter_.GenerateI3Particle(*secondaries.at(i)));
        I3Particle &segment = daughters.back();

        // I3PropagatorModule and related clients will not return secondaries
        // to the originating propagator by default. To force muon secondaries
        // (e.g. from pair production) to themselves be propagated, we unset
        // their length.
        if (secondaries.at(i)->GetTypeId() == PROPOSAL::DynamicData::Particle
            && proposal_service_.IsRegistered(static_cast<const PROPOSAL::Particle*>(secondaries.at(i))->GetParticleDef())
            && segment.GetLength() == 0)
        {
            segment.SetLength(NAN);
        }

        // Emit track segments with the current energy for use with
        // parameterizations such as https://arxiv.org/abs/1206.5530
        if (slice_tracks_) {
            double energyLoss = segment.GetEnergy();
            if (segment.GetType() == I3Particle::ContinuousEnergyLoss) {
                segment.SetType(p.GetType());
                segment.SetEnergy(currentEnergy);
            }
            currentEnergy -= std::min(energyLoss, currentEnergy);
        }
    }
    // Mark the parent track as Dark; its light yield is already accounted for
    // by the track segments.
    if (slice_tracks_) {
        p.SetShape(I3Particle::Dark);
    }

    if (final_stochastic_loss_ != I3Particle::unknown)
    {
        I3Particle i3_particle = particle_converter_.GenerateI3Particle(particle);
        i3_particle.SetType(final_stochastic_loss_);
        i3_particle.SetEnergy(particle.GetEnergy() - particle.GetParticleDef().mass);

        daughters.push_back(i3_particle);
    }

    Output::getInstance().ClearSecondaryVector(); // Tomasz

    return mmcTrack;
}

// ------------------------------------------------------------------------- //
I3MMCTrackPtr I3PropagatorServicePROPOSAL::GenerateMMCTrack(PROPOSAL::Particle* particle)
{

    // explicitly specifying the units from MMC
    double xi = particle->GetEntryPoint().GetX() * I3Units::cm / I3Units::m;
    double yi = particle->GetEntryPoint().GetY() * I3Units::cm / I3Units::m;
    double zi = particle->GetEntryPoint().GetZ() * I3Units::cm / I3Units::m;
    double ti = particle->GetEntryTime() * I3Units::second;
    double Ei = particle->GetEntryEnergy() * I3Units::MeV / I3Units::GeV;

    double xf = particle->GetExitPoint().GetX() * I3Units::cm / I3Units::m;
    double yf = particle->GetExitPoint().GetY() * I3Units::cm / I3Units::m;
    double zf = particle->GetExitPoint().GetZ() * I3Units::cm / I3Units::m;
    double tf = particle->GetExitTime() * I3Units::second;
    double Ef = particle->GetExitEnergy() * I3Units::MeV / I3Units::GeV;

    if (Ef == particle->GetMass() * I3Units::MeV / I3Units::GeV)
    {
        // If a particle decays, it gets its mass as the final energy
        // (in the new version of PROPOSAL).
        // However, in the old version (or in MMC),
        // it gets the negative length to the border of the detector in cm.
        // As this is not intuitive, this was changed in the core library.
        // But to make make the MMCTrack consistent with the old one,
        // this is changed here in the interface.
        PROPOSAL::Geometry& tmp_detector = proposal_service_.GetPropagatorToParticleDef(particle->GetParticleDef())->GetDetector();
        Ef = -tmp_detector.DistanceToBorder(particle->GetExitPoint(), particle->GetDirection()).first;
    }

    double xc = particle->GetClosestApproachPoint().GetX() * I3Units::cm / I3Units::m;
    double yc = particle->GetClosestApproachPoint().GetY() * I3Units::cm / I3Units::m;
    double zc = particle->GetClosestApproachPoint().GetZ() * I3Units::cm / I3Units::m;
    double tc = particle->GetClosestApproachTime() * I3Units::second;
    double Ec = particle->GetClosestApproachEnergy() * I3Units::MeV / I3Units::GeV;

    double Elost = particle->GetElost() * I3Units::MeV / I3Units::GeV;

    I3MMCTrackPtr mmcTrack(new I3MMCTrack);
    mmcTrack->SetEnter(xi, yi, zi, ti, Ei);
    mmcTrack->SetCenter(xc, yc, zc, tc, Ec);
    mmcTrack->SetExit(xf, yf, zf, tf, Ef);

    mmcTrack->SetDepositedEnergy(Elost);
    log_debug("Elost = %f", Elost);

    if (Elost > 0)
        return mmcTrack;
    // return null pointer if Elost <= 0
    return I3MMCTrackPtr();
}
