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
#include <boost/foreach.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/bimap.hpp>

#include <sstream>
#include <unistd.h> // check for write permissions

#include <dataclasses/physics/I3Particle.h>

// #include "PROPOSAL/src/PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"

using namespace std;
using namespace PROPOSAL;

class Output;

// ------------------------------------------------------------------------- //
bool IsWritable(std::string table_dir)
{
    bool writeable = false;

    if (access(table_dir.c_str(), F_OK) == 0)
    {
        if ((access(table_dir.c_str(), R_OK) == 0) && (access(table_dir.c_str(), W_OK) == 0))
        {
            writeable = true;
            log_info("Table directory does exist and has read and write permissions: %s", table_dir.c_str());
        }
        else
        {
            if (access(table_dir.c_str(), R_OK) != 0)
                log_info("Table directory is not readable: %s", table_dir.c_str());
            else
                log_info("Table directory is not writable: %s", table_dir.c_str());
        }
    }
    else
        log_info("Table directory does not exist: %s", table_dir.c_str());

    return writeable;
}

// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL( I3Particle& p, std::string configfile)
    : tearDownPerCall_(false)
    , particle_def_(GeneratePROPOSALParticleDef(p))
    , config_file_(configfile.empty() || !boost::filesystem::exists(configfile) ? GetDefaultConfigFile() : configfile)
    , proposal_(new Propagator(particle_def_, config_file_))
{
}

// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::~I3PropagatorServicePROPOSAL()
{
    delete proposal_;
}

// ------------------------------------------------------------------------- //
std::string I3PropagatorServicePROPOSAL::GetDefaultConfigFile()
{
	const char *I3_BUILD = getenv("I3_BUILD");
	if (!I3_BUILD)
		log_fatal("$I3_BUILD is not set!");
	std::string s(I3_BUILD);

    return s + "/PROPOSAL/resources/config.json";
}

// ------------------------------------------------------------------------- //
void I3PropagatorServicePROPOSAL::SetRandomNumberGenerator(I3RandomServicePtr random)
{
    rng_ = random;
    boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, random, 0, 1);

    PROPOSAL::RandomGenerator::Get().SetRandomNumberGenerator(f);
}

// ------------------------------------------------------------------------- //
std::vector<I3Particle> I3PropagatorServicePROPOSAL::Propagate(I3Particle& p, DiagnosticMapPtr frame){
    // saying where we are
    log_debug("Entering I3PropagatorServicePROPOSAL::Propagate()");

    vector<I3Particle> daughters;

    log_trace("location type = %d",p.GetLocationType());
    if(p.GetLocationType()!=I3Particle::InIce) return std::vector<I3Particle>();

    if((p.GetType() == I3Particle::NuE)   || (p.GetType() == I3Particle::NuEBar)   ||
       (p.GetType() == I3Particle::NuMu)  || (p.GetType() == I3Particle::NuMuBar)  ||
       (p.GetType() == I3Particle::NuTau) || (p.GetType() == I3Particle::NuTauBar))
        return std::vector<I3Particle>();

    log_trace(
        "particle to propagate:\n""type/energy[GeV]/posx[m]/posy[m]/posz[m]/theta[deg]/phi[deg]/length[m]\n"
        "%d/%.2e/%.2f/%.2f/%.2f/%.2f/%.2f/%.2f",
        p.GetType(),
        p.GetEnergy()/I3Units::GeV,
        p.GetPos().GetX()/I3Units::m,
        p.GetPos().GetY()/I3Units::m,
        p.GetPos().GetZ()/I3Units::m,
        p.GetZenith ()/I3Units::deg,
        p.GetAzimuth ()/I3Units::deg,
        p.GetLength()/I3Units::m
    );

    //TODO(mario):  Thu 2017/09/28
    // if (tearDownPerCall_)
    // {
    //     delete proposal;
    //     proposal = new Propagator(mediadef_,false);
    //     // Apply Settings
    //     Geometry* geo = new Geometry();
    //     geo->InitCylinder(0,0,0,cylinderRadius_,0,cylinderHeight_);
    //     proposal->SetDetector(geo);
    //     proposal->SetBrems(brems_param_);
    //     proposal->SetPhoto(photo_param_);
    //     proposal->SetPath_to_tables(tabledir_);
    //     // proposal->SetParticle(particle_type)
    //     proposal->ApplyOptions();
    //
    //     boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, rng_, 0, 1);
    //     proposal->SetRandomNumberGenerator(f);
    // }

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

typedef boost::bimap<I3Particle::ParticleType, ParticleType::Enum> bimap_ParticleType;
static const bimap_ParticleType I3_PROPOSAL_ParticleType_bimap = boost::assign::list_of<bimap_ParticleType::relation>
    (I3Particle::MuMinus,   ParticleType::MuMinus)
    (I3Particle::MuPlus,    ParticleType::MuPlus)
    (I3Particle::TauMinus,  ParticleType::TauMinus)
    (I3Particle::TauPlus,   ParticleType::TauPlus)
    (I3Particle::EMinus,    ParticleType::EMinus)
    (I3Particle::EPlus,     ParticleType::EPlus)
    (I3Particle::NuMu,      ParticleType::NuMu)
    (I3Particle::NuMuBar,   ParticleType::NuMuBar)
    (I3Particle::NuE,       ParticleType::NuE)
    (I3Particle::NuEBar,    ParticleType::NuEBar)
    (I3Particle::NuTau,     ParticleType::NuTau)
    (I3Particle::NuTauBar,  ParticleType::NuTauBar)
    (I3Particle::Brems,     ParticleType::Brems)
    (I3Particle::DeltaE,    ParticleType::DeltaE)
    (I3Particle::PairProd,  ParticleType::EPair)
    (I3Particle::NuclInt,   ParticleType::NuclInt)
    (I3Particle::MuPair,    ParticleType::MuPair)
    (I3Particle::Hadrons,   ParticleType::Hadrons)
    (I3Particle::Monopole,  ParticleType::Monopole)
    (I3Particle::STauMinus, ParticleType::STauMinus)
    (I3Particle::STauPlus,  ParticleType::STauPlus)
    // (I3Particle::SMP,       ParticleType::StableMassiveParticle)
    (I3Particle::Gamma,     ParticleType::Gamma)
    (I3Particle::Pi0,       ParticleType::Pi0)
    (I3Particle::PiPlus,    ParticleType::PiPlus)
    (I3Particle::PiMinus,   ParticleType::PiMinus)
    (I3Particle::KPlus,     ParticleType::KPlus)
    (I3Particle::KMinus,    ParticleType::KMinus)
    (I3Particle::PPlus,     ParticleType::PPlus)
    (I3Particle::PMinus,    ParticleType::PMinus);


// ------------------------------------------------------------------------- //
PROPOSAL::ParticleDef I3PropagatorServicePROPOSAL::GeneratePROPOSALParticleDef(const I3Particle& p)
{
    I3Particle::ParticleType ptype_I3 = p.GetType();

    if (ptype_I3 == I3Particle::MuMinus) return MuMinusDef::Get();
    else if (ptype_I3 == I3Particle::MuPlus) return MuPlusDef::Get();
    else if (ptype_I3 == I3Particle::EMinus) return EMinusDef::Get();
    else if (ptype_I3 == I3Particle::EPlus) return EPlusDef::Get();
    else if (ptype_I3 == I3Particle::TauMinus) return TauMinusDef::Get();
    else if (ptype_I3 == I3Particle::TauPlus) return TauPlusDef::Get();
    else
    {
        log_fatal("The I3Particle '%s' with type '%i' can not be converted to a PROPOSALParticle"
            , p.GetTypeString().c_str(), ptype_I3);
    }
}

I3Particle::ParticleType I3PropagatorServicePROPOSAL::GenerateI3Type(const PROPOSAL::DynamicData& secondary)
{
    DynamicData::Type type = secondary.GetTypeId();

    if (type == DynamicData::Particle)
    {
        const PROPOSAL::Particle& particle = static_cast<const PROPOSAL::Particle&>(secondary);
        ParticleDef particle_def = particle.GetParticleDef();
        if (particle_def == MuMinusDef::Get()) return I3Particle::MuMinus;
        else if (particle_def == MuPlusDef::Get()) return I3Particle::MuPlus;
        else if (particle_def == EMinusDef::Get()) return I3Particle::EMinus;
        else if (particle_def == EPlusDef::Get()) return I3Particle::EPlus;
        else if (particle_def == TauMinusDef::Get()) return I3Particle::TauMinus;
        else if (particle_def == TauPlusDef::Get()) return I3Particle::TauPlus;
        else
        {
            log_fatal("PROPOSALParticle '%s' can not be converted to a I3Particle" , particle_def.name.c_str());
        }
    }
    else if(type == DynamicData::Brems) return I3Particle::Brems;
    else if(type == DynamicData::Epair) return I3Particle::PairProd;
    else if(type == DynamicData::DeltaE) return I3Particle::DeltaE;
    else if(type == DynamicData::NuclInt) return I3Particle::NuclInt;
    else
    {
        log_fatal("PROPOSALParticle can not be converted to a I3Particle");
    }


    // bimap_ParticleType::right_const_iterator proposal_iterator = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL);
    // if (proposal_iterator == I3_PROPOSAL_ParticleType_bimap.right.end())
    // {
    //     log_fatal("The PROPOSALParticle '%s' with type '%i' can not be converted to a I3Particle"
    //         , PROPOSALParticle::GetName(ptype_PROPOSAL).c_str(), ptype_PROPOSAL);
    // }

    // ptype_I3 = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL) -> second;
}

// ------------------------------------------------------------------------- //
I3MMCTrackPtr I3PropagatorServicePROPOSAL::GenerateMMCTrack(PROPOSAL::Particle* particle){

    //explicitly specifying the units from MMC
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

    I3MMCTrackPtr mmcTrack( new I3MMCTrack);
    mmcTrack->SetEnter(xi,yi,zi,ti,Ei);
    mmcTrack->SetCenter(xc,yc,zc,tc,Ec);
    mmcTrack->SetExit(xf,yf,zf,tf,Ef);

    double Elost = particle->GetElost() * I3Units::GeV;
    mmcTrack->SetDepositedEnergy(Elost);
    log_debug("Elost = %f", Elost);

    if(Elost>0)
        return mmcTrack;
    //return null pointer if Elost <= 0
    return I3MMCTrackPtr();
}


// ------------------------------------------------------------------------- //
I3MMCTrackPtr I3PropagatorServicePROPOSAL::propagate( I3Particle& p, vector<I3Particle>& daughters){
    /**
     * Natural units of MMC is cm, deg, MeV, and s.
     * Therefore we need to convert explicitly to
     * MMC units before passing the propagate method
     */
    double x_0 = p.GetPos().GetX()/I3Units::cm;     // [cm]
    double y_0 = p.GetPos().GetY()/I3Units::cm;    // [cm]
    double z_0 = p.GetPos().GetZ()/I3Units::cm;     // [cm]
    double theta_0 = p.GetDir().CalcTheta()/I3Units::deg; // [deg]
    double phi_0 = p.GetDir().CalcPhi()/I3Units::deg;   // [deg]
    double e_0 = p.GetEnergy()/I3Units::MeV;  // [MeV]
    double t_0 = p.GetTime()/I3Units::s;     // [s]

    // log_debug("Name of particle to propagate: %s", PROPOSALParticle::GetName(GeneratePROPOSALType(p)).c_str());


    PROPOSAL::Particle& particle = proposal_->GetParticle();

    particle.SetPosition(PROPOSAL::Vector3D(x_0, y_0, z_0));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi_0, theta_0);
    direction.CalculateCartesianFromSpherical();
    particle.SetDirection(direction);
    particle.SetEnergy(e_0);
    particle.SetTime(t_0);

    proposal_->Propagate();

    std::vector<DynamicData*> secondaries = Output::getInstance().GetSecondarys();
    //get the propagated length of the particle
    double length = particle.GetPropagatedDistance();

    p.SetLength( length * I3Units::cm );
    log_trace(" length = %f cm ", length );

    int nParticles =  secondaries.size();
    log_trace("nParticles = %d", nParticles);

    I3MMCTrackPtr mmcTrack;

    mmcTrack = GenerateMMCTrack(&particle);

    if(mmcTrack)
        mmcTrack->SetParticle( p );

    for(int i=0; i < nParticles; i++)
    {
        //Tomasz
        //in mmc the particle relationships are stored
        double x = secondaries.at(i)->GetPosition().GetX() * I3Units::cm;
        double y = secondaries.at(i)->GetPosition().GetY() * I3Units::cm;
        double z = secondaries.at(i)->GetPosition().GetZ() * I3Units::cm;
        double theta = secondaries.at(i)->GetDirection().GetTheta() * I3Units::deg;
        double phi = secondaries.at(i)->GetDirection().GetPhi() * I3Units::deg;
        double t = secondaries.at(i)->GetTime() * I3Units::s;
        double e = secondaries.at(i)->GetEnergy() * I3Units::MeV;
        double l = secondaries.at(i)->GetPropagatedDistance() * I3Units::cm;

        log_trace("MMC DEBUG SEC  \n    pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g",
                  x, y, z, theta, phi, e, t, l);

        //this should be a stochastic
        I3Particle new_particle;

        ParticleDef particle_def = GeneratePROPOSALParticleDef(p);
        if (particle_def == EMinusDef::Get()
            || particle_def == EPlusDef::Get())
            // || particle_def == ParticleType::Hadrons) //TODO(mario):  Wed 2017/10/25
        {
            if (p.GetShape() != I3Particle::TopShower)
            {
                log_fatal("The particle '%s' has no TopShower shape, but 'e-', 'e+' and 'Hadrons' need that. I don't know why?", particle_def.name.c_str());
            }
        }

        new_particle.SetType(GenerateI3Type(*secondaries.at(i)));

        new_particle.SetLocationType(I3Particle::InIce);
        new_particle.SetPos(x, y, z);
        new_particle.SetTime(t);
        new_particle.SetLength(l);
        new_particle.SetThetaPhi(theta,phi);
        new_particle.SetEnergy(e);

        // this is not the particle you're looking for
        // move along...and add it to the daughter list
        daughters.push_back(new_particle);
    }

    Output::getInstance().ClearSecondaryVector(); //Tomasz

    return mmcTrack;
}

