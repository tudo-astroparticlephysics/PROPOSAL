/** $Id: SimplePropagator.cxx 150208 2016-09-21 11:32:47Z tfuchs $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 150208 $
 * $Date: 2016-09-21 13:32:47 +0200 (Mi, 21. Sep 2016) $
 */

#include <boost/assign.hpp>
#include <boost/bimap.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include "PROPOSAL-icetray/Converter.h"
#include "PROPOSAL-icetray/SimplePropagator.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;


SimplePropagator::SimplePropagator(I3Particle::ParticleType pt,
                                   const std::string& medium,
                                   double ecut,
                                   double vcut,
                                   double rho)
{
    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    // Sector definition

    PROPOSAL::Sector::Definition sec_def;

    sec_def.stopping_decay              = true;
    sec_def.scattering_model            = PROPOSAL::ScatteringFactory::Moliere;
    sec_def.do_exact_time_calculation   = true;
    sec_def.do_continuous_randomization = ecut < 0;

    // Medium

    sec_def.SetMedium(*PROPOSAL::MediumFactory::Get().CreateMedium(medium, 1.0));

    // Geometry

    PROPOSAL::Vector3D position(0.0, 0.0, 0.0);
    PROPOSAL::Sphere detector(position, 0.0, 1e18);

    sec_def.SetGeometry(detector);

    // Cuts

    sec_def.cut_settings.SetEcut(ecut);
    sec_def.cut_settings.SetVcut(vcut);

    // Parametrizations

    // Using defaults:
    // lpm effect = true
    // Bremsstrahlung = KelnerKokoulinPetrukhin
    // Photonuclear = AbramowiczLevinLevyMaor97

    // Interpolation

    PROPOSAL::InterpolationDef interpolation_def;
    interpolation_def.path_to_tables = prefix.str();

    // Init new propagator

    // TODO(mario): Check for muon, tau Mon 2017/11/06
    propagator_ = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(),
                                           boost::assign::list_of<PROPOSAL::Sector::Definition>(sec_def),
                                           detector,
                                           interpolation_def);
}

SimplePropagator::~SimplePropagator()
{
    delete propagator_;
}

void SimplePropagator::SetSeed(int seed)
{
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
}

void SimplePropagator::SetRandomNumberGenerator(I3RandomServicePtr rng)
{
    boost::function<double()> f = boost::bind(&I3RandomService::Uniform, rng, 0, 1);
    PROPOSAL::RandomGenerator::Get().SetRandomNumberGenerator(f);
}

I3Particle SimplePropagator::to_I3Particle(const PROPOSAL::DynamicData& pp)
{
    I3Particle p;

    double x = pp.GetPosition().GetX() * I3Units::cm;
    double y = pp.GetPosition().GetY() * I3Units::cm;
    double z = pp.GetPosition().GetZ() * I3Units::cm;

    double theta = pp.GetDirection().GetTheta() * I3Units::degree;
    double phi   = pp.GetDirection().GetPhi() * I3Units::degree;

    p.SetPos(x, y, z);
    p.SetThetaPhi(theta, phi);
    p.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    p.SetTime(pp.GetTime() * I3Units::second);

    p.SetType(I3PROPOSALParticleConverter::GenerateI3Type(pp));
    p.SetLocationType(I3Particle::InIce);
    p.SetEnergy(pp.GetEnergy() * I3Units::MeV);

    return p;
}

I3Particle SimplePropagator::propagate(const I3Particle& p,
                                       double distance,
                                       boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x, y, z, theta, phi = 0.0;

    PROPOSAL::Particle pp = propagator_->GetParticle();
    pp.SetParentParticleId(0);
    pp.SetParticleId(0);
    pp.SetTime(p.GetTime() / I3Units::second);

    x = p.GetPos().GetX() / I3Units::cm;
    y = p.GetPos().GetY() / I3Units::cm;
    z = p.GetPos().GetZ() / I3Units::cm;

    pp.SetPosition(PROPOSAL::Vector3D(x, y, z));
    pp.SetEnergy(p.GetEnergy() / I3Units::MeV);

    propagator_->Propagate(distance / I3Units::cm);

    endpoint.SetEnergy(pp.GetEnergy() * I3Units::MeV);

    x = pp.GetPosition().GetX() * I3Units::cm;
    y = pp.GetPosition().GetY() * I3Units::cm;
    z = pp.GetPosition().GetZ() * I3Units::cm;

    theta = pp.GetDirection().GetTheta() * I3Units::degree;
    phi   = pp.GetDirection().GetPhi() * I3Units::degree;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);
    endpoint.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    endpoint.SetTime(pp.GetTime() * I3Units::second);

    // Tomasz
    if (losses)
    {
        std::vector<PROPOSAL::DynamicData*> history = Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(to_I3Particle(*history[i]));
        }
    }

    Output::getInstance().ClearSecondaryVector();

    return endpoint;
}
