/** $Id: SimplePropagator.cxx 150208 2016-09-21 11:32:47Z tfuchs $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 150208 $
 * $Date: 2016-09-21 13:32:47 +0200 (Mi, 21. Sep 2016) $
 */

#include <boost/bind.hpp>
#include <vector>

#include "PROPOSAL-icetray/Converter.h"
#include "PROPOSAL-icetray/SimplePropagator.h"

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
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables";

    // Sector definition

    PROPOSAL::Sector::Definition sec_def;

    sec_def.stopping_decay              = final_stochastic_loss_ == I3Particle::unknown;
    sec_def.scattering_model            = PROPOSAL::ScatteringFactory::Highland;
    sec_def.do_exact_time_calculation   = true;
    sec_def.do_continuous_randomization = ecut < 0;

    // Medium

    sec_def.SetMedium(*PROPOSAL::MediumFactory::Get().CreateMedium(medium, rho));

    // Geometry

    PROPOSAL::Vector3D position(0.0, 0.0, 0.0);
    PROPOSAL::Sphere detector(position, 1e18, 0.0);

    sec_def.SetGeometry(detector);

    // Cuts

    sec_def.cut_settings.SetEcut(ecut);
    sec_def.cut_settings.SetVcut(vcut);

    std::vector<Sector::Definition> sec_def_vec;
    sec_def_vec.push_back(sec_def);

    // Parametrizations

    // Using defaults:
    // lpm effect = true
    // Bremsstrahlung = KelnerKokoulinPetrukhin
    // Photonuclear = AbramowiczLevinLevyMaor97

    // Interpolation

    PROPOSAL::InterpolationDef interpolation_def;
    interpolation_def.path_to_tables = prefix.str();
    interpolation_def.path_to_tables_readonly = prefix.str();

    // Init new propagator

    propagator_ = new PROPOSAL::Propagator(I3PROPOSALParticleConverter::GeneratePROPOSALType(pt),
                                           sec_def_vec,
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
    std::function<double()> f = boost::bind(&I3RandomService::Uniform, rng, 0, 1);
    PROPOSAL::RandomGenerator::Get().SetRandomNumberGenerator(f);
}

I3Particle SimplePropagator::propagate(const I3Particle& p,
                                       double distance,
                                       boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x = p.GetPos().GetX() / I3Units::cm;
    double y = p.GetPos().GetY() / I3Units::cm;
    double z = p.GetPos().GetZ() / I3Units::cm;

    double theta = p.GetDir().CalcTheta() / I3Units::radian;
    double phi   = p.GetDir().CalcPhi() / I3Units::radian;

    double energy = p.GetEnergy() / I3Units::MeV;
    double time   = p.GetTime() / I3Units::second;
    double length = p.GetLength() / I3Units::cm;

    // The Muons from NuGen have NaN as default propagated length.
    // So this has to be corrected.
    if (std::isnan(length))
    {
        length = 0.0;
    }
    else if (std::isinf(length))
    {
        log_fatal("the propagated length is Inf, should be finite or NaN.");
    }

    PROPOSAL::Particle& particle = propagator_->GetParticle();
    particle.SetParentParticleId(0);
    particle.SetParticleId(0);

    particle.SetPosition(PROPOSAL::Vector3D(x, y, z));

    PROPOSAL::Vector3D direction;
    direction.SetSphericalCoordinates(1.0, phi, theta);
    direction.CalculateCartesianFromSpherical();
    particle.SetDirection(direction);

    particle.SetEnergy(energy);
    particle.SetTime(time);
    particle.SetPropagatedDistance(length);

    std::vector<PROPOSAL::DynamicData*> history = propagator_->Propagate(distance / I3Units::cm);

    x = particle.GetPosition().GetX() * I3Units::cm;
    y = particle.GetPosition().GetY() * I3Units::cm;
    z = particle.GetPosition().GetZ() * I3Units::cm;

    theta = particle.GetDirection().GetTheta() * I3Units::radian;
    phi   = particle.GetDirection().GetPhi() * I3Units::radian;

    energy = particle.GetEnergy() * I3Units::MeV;
    time   = particle.GetTime() * I3Units::second;
    length = particle.GetPropagatedDistance() * I3Units::cm;

    endpoint.SetPos(x, y, z);
    endpoint.SetThetaPhi(theta, phi);

    endpoint.SetEnergy(energy);
    endpoint.SetLength(length);
    endpoint.SetTime(time);

    // Tomasz
    if (losses)
    {
        // std::vector<PROPOSAL::DynamicData*> history = Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(particle_converter_.GenerateI3Particle(*history[i]));
        }

        if (final_stochastic_loss_ != I3Particle::unknown)
        {
            I3Particle i3_particle = particle_converter_.GenerateI3Particle(particle);
            i3_particle.SetType(final_stochastic_loss_);
            i3_particle.SetEnergy(particle.GetEnergy() - particle.GetParticleDef().mass);

            losses->push_back(i3_particle);
        }
    }

    Output::getInstance().ClearSecondaryVector();

    return endpoint;
}
