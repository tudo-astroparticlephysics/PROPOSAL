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
                                   double rho,
                                   I3Particle::ParticleType final_loss)
    : propagator_(NULL)
    , final_stochastic_loss_(final_loss)
{
    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    // Sector definition

    PROPOSAL::Sector::Definition sec_def;

    sec_def.stopping_decay              = final_stochastic_loss_ == I3Particle::unknown;
    sec_def.scattering_model            = PROPOSAL::ScatteringFactory::Moliere;
    sec_def.do_exact_time_calculation   = true;
    sec_def.do_continuous_randomization = ecut < 0;

    // Medium

    sec_def.SetMedium(*PROPOSAL::MediumFactory::Get().CreateMedium(medium, rho));

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

    propagator_ = new PROPOSAL::Propagator(I3PROPOSALParticleConverter::GeneratePROPOSALType(pt),
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

I3Particle SimplePropagator::propagate(const I3Particle& p,
                                       double distance,
                                       boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x, y, z, theta, phi;

    PROPOSAL::Particle& pp = propagator_->GetParticle();
    pp.SetParentParticleId(0);
    pp.SetParticleId(0);

    x = p.GetPos().GetX() / I3Units::cm;
    y = p.GetPos().GetY() / I3Units::cm;
    z = p.GetPos().GetZ() / I3Units::cm;

    theta = p.GetDir().CalcTheta() / I3Units::radian;
    phi   = p.GetDir().CalcPhi() / I3Units::radian;

    pp.SetPosition(PROPOSAL::Vector3D(x, y, z));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi, theta);
    direction.CalculateCartesianFromSpherical();
    pp.SetDirection(direction);

    pp.SetEnergy(p.GetEnergy() / I3Units::MeV);
    pp.SetTime(p.GetTime() / I3Units::second);
    pp.SetPropagatedDistance(p.GetLength() / I3Units::cm);

    std::vector<PROPOSAL::DynamicData*> history = propagator_->Propagate(distance / I3Units::cm);

    x = pp.GetPosition().GetX() * I3Units::cm;
    y = pp.GetPosition().GetY() * I3Units::cm;
    z = pp.GetPosition().GetZ() * I3Units::cm;

    theta = pp.GetDirection().GetTheta() * I3Units::radian;
    phi   = pp.GetDirection().GetPhi() * I3Units::radian;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);

    endpoint.SetEnergy(pp.GetEnergy() * I3Units::MeV);
    endpoint.SetLength(pp.GetPropagatedDistance() * I3Units::cm);
    endpoint.SetTime(pp.GetTime() * I3Units::second);

    // Tomasz
    if (losses)
    {
        // std::vector<PROPOSAL::DynamicData*> history = Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(I3PROPOSALParticleConverter::GenerateI3Particle(*history[i]));
        }

        if (final_stochastic_loss_ != I3Particle::unknown)
        {
            I3Particle i3_particle = I3PROPOSALParticleConverter::GenerateI3Particle(pp);
            i3_particle.SetType(final_stochastic_loss_);
            i3_particle.SetEnergy(pp.GetEnergy() - pp.GetParticleDef().mass);

            losses->push_back(i3_particle);
        }
    }

    Output::getInstance().ClearSecondaryVector();

    return endpoint;
}
