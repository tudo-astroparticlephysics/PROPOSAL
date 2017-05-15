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

#include "dataclasses/physics/I3Particle.h"
#include "simclasses/I3MMCTrack.h"

#include "PROPOSAL-icetray/I3PropagatorServicePROPOSAL.h"
#include "PROPOSAL/Propagator.h"
#include "PROPOSAL/PROPOSALParticle.h"
#include "PROPOSAL/Output.h"

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
I3PropagatorServicePROPOSAL::I3PropagatorServicePROPOSAL(
                                                         std::string mediadef
                                                         , std::string tabledir
                                                         , double cylinderRadius
                                                         , double cylinderHeight
                                                         , I3Particle::ParticleType type
                                                         , double particleMass
                                                         , ParametrizationType::Enum brems_param
                                                         , ParametrizationType::Enum photo_param)
    : particleMass_(particleMass)
    , mediadef_(mediadef)
    , tabledir_(tabledir)
    , cylinderRadius_(cylinderRadius)
    , cylinderHeight_(cylinderHeight)
    , brems_param_(brems_param)
    , photo_param_(photo_param)
{

    I3Particle i3particle;
    i3particle.SetType(type);
    PROPOSALParticle* particle = new PROPOSALParticle(GeneratePROPOSALType(i3particle));
    if (isnan(particleMass_) == false)
    {
        particle->SetMass(particleMass_);
    }
    if (std::isinf(particleMass_))
    {
        log_fatal("You asked for Particle of non-finite mass %f", particleMass_);
    }

    // ----[ Check table dir and mediadef ]------------------ //

    if (mediadef_.empty())
        mediadef_ = GetDefaultMediaDef();
    if (tabledir_.empty())
        tabledir_ = GetDefaultTableDir();

    namespace fs = boost::filesystem;

    if (!fs::exists(mediadef_))
        log_fatal("The mediadef file '%s' can't be read!", mediadef_.c_str());
    if (!fs::is_directory(tabledir_))
        log_fatal("The table directory '%s' doesn't exist!", tabledir_.c_str());


    // ----[ Check, apply options ]-------------------------- //


    // Define propagator but do not apply option yet
    // proposal = new Propagator(mediadef_, false);
    proposal = new Propagator(mediadef_,particle, false);

    stringstream options;

    options << "You choose the following parameter by passing arguments:" << std::endl;
    options << "\tcylinderRadius = " << cylinderRadius_ << std::endl;
    options << "\tcylinderHeight = " << cylinderHeight_ << std::endl;
    options << "\tBremsstrahlungParametrization = " << brems_param_ << std::endl;
    options << "\tPhotonuclearParametrization = " << photo_param_ << std::endl;

    if (brems_param_ != proposal->GetBrems())
    {
        options << "\tChosen BremsstrahlungParametrization differs from parametrization in config file!" << std::endl;
        options << "Passed parametrization will be used: " << brems_param_ << std::endl;
    }

    if (photo_param_ != proposal->GetPhoto())
    {
        options << "\tChosen PhotonuclearParametrization differs from parametrization in config file!" << std::endl;
        options << "Passed parametrization will be used: " << photo_param_ << std::endl;
    }

    log_info("%s", options.str().c_str());


    Geometry* geo = new Geometry();
    geo->InitCylinder(0,0,0,cylinderRadius_,0,cylinderHeight_);
    proposal->SetDetector(geo);
    proposal->SetBrems(brems_param_);
    proposal->SetPhoto(photo_param_);
    proposal->SetPath_to_tables(tabledir_);

    // proposal->SetParticle(particle_type)

    proposal->ApplyOptions();

	tearDownPerCall_ = false;
}

// ------------------------------------------------------------------------- //
I3PropagatorServicePROPOSAL::~I3PropagatorServicePROPOSAL()
{
    delete proposal;
}

// ------------------------------------------------------------------------- //
std::string I3PropagatorServicePROPOSAL::GetDefaultMediaDef()
{
	const char *I3_BUILD = getenv("I3_BUILD");
	if (!I3_BUILD)
		log_fatal("$I3_BUILD is not set!");
	std::string s(I3_BUILD);

    return s + "/PROPOSAL/resources/configuration";
}

// ------------------------------------------------------------------------- //
std::string I3PropagatorServicePROPOSAL::GetDefaultTableDir()
{
    std::string append_string  = "/PROPOSAL/resources/tables";
    std::string append_string2 = "/PROPOSAL/tables";

    // --------------------------------------------------------------------- //
    // Environment variable set?
    // --------------------------------------------------------------------- //

    // Initializing a std::string with a NULL ptr is undefined behavior.
    // Why it doens't just return an empty string, I have no idea.
    std::string table_dir(getenv("PROPOSALTABLEDIR") ? getenv("PROPOSALTABLEDIR") : "");

    if (table_dir.empty())
    {
        log_info("$PROPOSALTABLEDIR is not set in env variables. Falling back to defaults.");
    }
    else
    {
        if(IsWritable(table_dir)) return table_dir;
    }

    // --------------------------------------------------------------------- //
    // Writable in I3_TESTDATA?
    // --------------------------------------------------------------------- //

    table_dir = std::string(getenv("I3_TESTDATA") ? getenv("I3_TESTDATA") : "");

    if (table_dir.empty())
    {
        log_warn("$I3_TESTDATA is not set, falling back to build folder!");
    }
    else
    {
        if (IsWritable(table_dir + append_string))
        {
            return table_dir + append_string;
        }
        else
        {
            log_warn("Falling back to build folder!");
        }
    }

    // --------------------------------------------------------------------- //
    // Fall back to I3_BUID
    // --------------------------------------------------------------------- //

    table_dir = std::string(getenv("I3_BUILD") ? getenv("I3_BUILD") : "");

    if (table_dir.empty())
    {
        log_fatal("$I3_BUILD is not set");
    }
    else
    {
        if (IsWritable(table_dir + append_string))
        {
            return table_dir + append_string;
        }
        else if (IsWritable(table_dir + append_string2))
        {
            return table_dir + append_string2;
        }
        else
        {
            log_fatal("No folder availble to fall back! Abort search for table directory.");
        }
    }
}

// ------------------------------------------------------------------------- //
void I3PropagatorServicePROPOSAL::SetRandomNumberGenerator(I3RandomServicePtr random)
{
    rng_ = random;
    boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, random, 0, 1);

    proposal->SetRandomNumberGenerator(f);
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

    if (tearDownPerCall_)
    {
        delete proposal;
        proposal = new Propagator(mediadef_,false);
        // Apply Settings
        Geometry* geo = new Geometry();
        geo->InitCylinder(0,0,0,cylinderRadius_,0,cylinderHeight_);
        proposal->SetDetector(geo);
        proposal->SetBrems(brems_param_);
        proposal->SetPhoto(photo_param_);
        proposal->SetPath_to_tables(tabledir_);
        // proposal->SetParticle(particle_type)
        proposal->ApplyOptions();

        boost::function<double ()> f = boost::bind(&I3RandomService::Uniform, rng_, 0, 1);
        proposal->SetRandomNumberGenerator(f);
    }

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
    (I3Particle::Gamma,     ParticleType::Gamma)
    (I3Particle::Pi0,       ParticleType::Pi0)
    (I3Particle::PiPlus,    ParticleType::PiPlus)
    (I3Particle::PiMinus,   ParticleType::PiMinus)
    (I3Particle::KPlus,     ParticleType::KPlus)
    (I3Particle::KMinus,    ParticleType::KMinus)
    (I3Particle::PPlus,     ParticleType::PPlus)
    (I3Particle::PMinus,    ParticleType::PMinus);


// ------------------------------------------------------------------------- //
ParticleType::Enum I3PropagatorServicePROPOSAL::GeneratePROPOSALType(const I3Particle& p)
{
    I3Particle::ParticleType ptype_I3 = p.GetType();
    ParticleType::Enum ptype_PROPOSAL;

    bimap_ParticleType::left_const_iterator i3_iterator = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3);
    if (i3_iterator == I3_PROPOSAL_ParticleType_bimap.left.end())
    {
        log_fatal("The I3Particle '%s' with type '%i' can not be converted to a PROPOSALParticle"
            , p.GetTypeString().c_str(), ptype_I3);
    }

    ptype_PROPOSAL = I3_PROPOSAL_ParticleType_bimap.left.find(ptype_I3) -> second;

    return ptype_PROPOSAL;
}

I3Particle::ParticleType I3PropagatorServicePROPOSAL::GenerateI3Type(ParticleType::Enum ptype_PROPOSAL)
{
    I3Particle::ParticleType ptype_I3;

    bimap_ParticleType::right_const_iterator proposal_iterator = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL);
    if (proposal_iterator == I3_PROPOSAL_ParticleType_bimap.right.end())
    {
        log_fatal("The PROPOSALParticle '%s' with type '%i' can not be converted to a I3Particle"
            , PROPOSALParticle::GetName(ptype_PROPOSAL).c_str(), ptype_PROPOSAL);
    }

    ptype_I3 = I3_PROPOSAL_ParticleType_bimap.right.find(ptype_PROPOSAL) -> second;

    return ptype_I3;
}

// ------------------------------------------------------------------------- //
I3MMCTrackPtr I3PropagatorServicePROPOSAL::GenerateMMCTrack(PROPOSALParticle* particle){

    //explicitly specifying the units from MMC
    double xi = particle->GetXi() * I3Units::m;
    double yi = particle->GetYi() * I3Units::m;
    double zi = particle->GetZi() * I3Units::m;
    double ti = particle->GetTi() * I3Units::second;
    double Ei = particle->GetEi() * I3Units::GeV;

    double xf = particle->GetXf() * I3Units::m;
    double yf = particle->GetYf() * I3Units::m;
    double zf = particle->GetZf() * I3Units::m;
    double tf = particle->GetTf() * I3Units::second;
    double Ef = particle->GetEf() * I3Units::GeV;

    double xc = particle->GetXc() * I3Units::m;
    double yc = particle->GetYc() * I3Units::m;
    double zc = particle->GetZc() * I3Units::m;
    double tc = particle->GetTc() * I3Units::second;
    double Ec = particle->GetEc() * I3Units::GeV;

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

    log_debug("Name of particle to propagate: %s", PROPOSALParticle::GetName(GeneratePROPOSALType(p)).c_str());


    // PROPOSALParticle* particle = new PROPOSALParticle(GeneratePROPOSALType(p), x_0, y_0, z_0, theta_0, phi_0, e_0, t_0);
    proposal->ResetParticle();
    PROPOSALParticle* particle = proposal->GetParticle();

    if (particle == NULL)
        log_fatal("Error calling the Particle constructor");

    particle->SetX(x_0);
    particle->SetY(y_0);
    particle->SetZ(z_0);
    particle->SetTheta(theta_0);
    particle->SetPhi(phi_0);
    particle->SetEnergy(e_0);
    particle->SetT(t_0);

    proposal->propagate();
    // proposal->Propagate(particle);

    vector<PROPOSALParticle*> aobj_l = Output::getInstance().GetSecondarys();
    //get the propagated length of the particle
    double length = particle->GetPropagatedDistance();

    p.SetLength( length * I3Units::cm );
    log_trace(" length = %f cm ", length );

    int nParticles =  aobj_l.size();
    log_trace("nParticles = %d", nParticles);

    I3MMCTrackPtr mmcTrack;

    mmcTrack = GenerateMMCTrack(particle);

    if(mmcTrack)
        mmcTrack->SetParticle( p );

    for(int i=0; i < nParticles; i++)
    {
        //Tomasz
        //in mmc the particle relationships are stored
        ParticleType::Enum type = aobj_l.at(i)->GetType();
        double x = aobj_l.at(i)->GetX() * I3Units::cm;
        double y = aobj_l.at(i)->GetY() * I3Units::cm;
        double z = aobj_l.at(i)->GetZ() * I3Units::cm;
        double theta = aobj_l.at(i)->GetTheta() * I3Units::deg;
        double phi = aobj_l.at(i)->GetPhi() * I3Units::deg;
        double t = aobj_l.at(i)->GetT() * I3Units::s;
        double e = aobj_l.at(i)->GetEnergy() * I3Units::MeV;
        double l = aobj_l.at(i)->GetPropagatedDistance() * I3Units::cm;

        log_trace("MMC DEBUG SEC  \n    type=%d pos=(%g,%g,%g) ang=(%g,%g)  e=%g t=%g  l=%g",
                  type, x, y, z, theta, phi, e, t, l);

        //this should be a stochastic
        I3Particle new_particle;


        if (GeneratePROPOSALType(p) == ParticleType::EMinus
            || GeneratePROPOSALType(p) == ParticleType::EPlus
            || GeneratePROPOSALType(p) == ParticleType::Hadrons)
        {
            if (p.GetShape() != I3Particle::TopShower)
            {
                log_fatal("The particle '%s' has no TopShower shape, but 'e-', 'e+' and 'Hadrons' need that. I don't know why?",
                    PROPOSALParticle::GetName(GeneratePROPOSALType(p)).c_str());
            }
        }


        new_particle.SetType(GenerateI3Type(type));

        new_particle.SetLocationType(I3Particle::InIce);
        new_particle.SetPos(x, y, z);
        new_particle.SetTime(t);
        new_particle.SetLength(l);
        new_particle.SetThetaPhi(theta,phi);
        new_particle.SetEnergy(e);

        // this is not the particle you're looking for
        // move along...and add it to the daughter list
        daughters.push_back(new_particle);

        //we're done with eobj
        //delete aobj_l.at(i); //Tomasz
    }

    Output::getInstance().ClearSecondaryVector(); //Tomasz
    // delete particle; //Mario

    return mmcTrack;
}

