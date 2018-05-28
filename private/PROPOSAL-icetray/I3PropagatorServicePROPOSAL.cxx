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

// class header
#include <boost/assign.hpp>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>

#include <sstream>
#include <unistd.h> // check for write permissions

#include <dataclasses/physics/I3Particle.h>

#include "PROPOSAL-icetray/Converter.h"
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"

using namespace std;
using namespace PROPOSAL;


// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL(std::string configfile, I3Particle::ParticleType final_loss, double distance)
    : I3PropagatorService()
    , config_file_(PROPOSAL::Helper::ResolvePath(configfile))
    , proposal_service_()
    , final_stochastic_loss_(final_loss)
    , distance_to_propagate_(distance)
{
    if (config_file_.empty())
    {
        config_file_ = GetDefaultConfigFile();
    }

    log_info("Using configuration file: \"%s\"", config_file_.c_str());

    proposal_service_.RegisterPropagator(Propagator(MuMinusDef::Get(), config_file_));
    proposal_service_.RegisterPropagator(Propagator(MuPlusDef::Get(), config_file_));
    proposal_service_.RegisterPropagator(Propagator(TauMinusDef::Get(), config_file_));
    proposal_service_.RegisterPropagator(Propagator(TauPlusDef::Get(), config_file_));

    log_info("Propagator created for muons and taus");
}

// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::~I3PropagatorServicePROPOSAL()
{
    if (rng_)
    {
        rng_ = NULL;

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
    boost::function<double()> f = boost::bind(&I3RandomService::Uniform, random, 0, 1);
    PROPOSAL::RandomGenerator::Get().SetRandomNumberGenerator(f);
}

// ------------------------------------------------------------------------- //
void I3PropagatorServicePROPOSAL::RegisterParticleType(I3Particle::ParticleType ptype)
{
    ParticleDef particle_def = I3PROPOSALParticleConverter::GeneratePROPOSALType(ptype);
    proposal_service_.RegisterPropagator(Propagator(particle_def, config_file_));
}

// ------------------------------------------------------------------------- //
std::vector<I3Particle> I3PropagatorServicePROPOSAL::Propagate(I3Particle& p, DiagnosticMapPtr frame, I3FramePtr)
{
    // saying where we are
    log_debug("Entering I3PropagatorServicePROPOSAL::Propagate()");

    vector<I3Particle> daughters;

    log_trace("location type = %d", p.GetLocationType());
    if (p.GetLocationType() != I3Particle::InIce)
        return std::vector<I3Particle>();

    if ((p.GetType() == I3Particle::NuE) || (p.GetType() == I3Particle::NuEBar) || (p.GetType() == I3Particle::NuMu) ||
        (p.GetType() == I3Particle::NuMuBar) || (p.GetType() == I3Particle::NuTau) ||
        (p.GetType() == I3Particle::NuTauBar))
        return std::vector<I3Particle>();

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
I3MMCTrackPtr I3PropagatorServicePROPOSAL::propagate(I3Particle& p, vector<I3Particle>& daughters)
{
    Particle particle = I3PROPOSALParticleConverter::GeneratePROPOSALParticle(p);


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

    for (int i = 0; i < nParticles; i++)
    {
        // this should be a stochastic

        ParticleDef particle_def = I3PROPOSALParticleConverter::GeneratePROPOSALType(p.GetType());
        if (particle_def == EMinusDef::Get() || particle_def == EPlusDef::Get())
        // || particle_def == ParticleType::Hadrons) //TODO(mario):  Wed 2017/10/25
        {
            if (p.GetShape() != I3Particle::TopShower)
            {
                log_fatal("The particle '%s' has no TopShower shape, but 'e-', 'e+' and 'Hadrons' need that. I don't "
                          "know why?",
                          particle_def.name.c_str());
            }
        }

        // Tomasz
        // in mmc the particle relationships are stored

        // this is not the particle you're looking for
        // move along...and add it to the daughter list
        daughters.push_back(I3PROPOSALParticleConverter::GenerateI3Particle(*secondaries.at(i)));
    }

    if (final_stochastic_loss_ != I3Particle::unknown)
    {
        I3Particle i3_particle = I3PROPOSALParticleConverter::GenerateI3Particle(particle);
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
    double xi = particle->GetEntryPoint().GetX() * I3Units::m;
    double yi = particle->GetEntryPoint().GetY() * I3Units::m;
    double zi = particle->GetEntryPoint().GetZ() * I3Units::m;
    double ti = particle->GetEntryTime() * I3Units::second;
    double Ei = particle->GetEntryEnergy() * I3Units::GeV;

    double xf = particle->GetExitPoint().GetX() * I3Units::m;
    double yf = particle->GetExitPoint().GetY() * I3Units::m;
    double zf = particle->GetExitPoint().GetZ() * I3Units::m;
    double tf = particle->GetExitTime() * I3Units::second;
    double Ef = particle->GetExitEnergy() * I3Units::GeV;

    double xc = particle->GetClosestApproachPoint().GetX() * I3Units::m;
    double yc = particle->GetClosestApproachPoint().GetY() * I3Units::m;
    double zc = particle->GetClosestApproachPoint().GetZ() * I3Units::m;
    double tc = particle->GetClosestApproachTime() * I3Units::second;
    double Ec = particle->GetClosestApproachEnergy() * I3Units::GeV;

    I3MMCTrackPtr mmcTrack(new I3MMCTrack);
    mmcTrack->SetEnter(xi, yi, zi, ti, Ei);
    mmcTrack->SetCenter(xc, yc, zc, tc, Ec);
    mmcTrack->SetExit(xf, yf, zf, tf, Ef);

    double Elost = particle->GetElost() * I3Units::GeV;
    mmcTrack->SetDepositedEnergy(Elost);
    log_debug("Elost = %f", Elost);

    if (Elost > 0)
        return mmcTrack;
    // return null pointer if Elost <= 0
    return I3MMCTrackPtr();
}
