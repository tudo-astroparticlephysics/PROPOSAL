/** $Id: MuonPropagator.cxx 137064 2015-08-31 18:24:47Z jvansanten $
 * @file
 * @author Jakob van Santen <vansanten@wisc.edu>
 *
 * $Revision: 137064 $
 * $Date: 2015-08-31 20:24:47 +0200 (Mo, 31. Aug 2015) $
 */

#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include "MuonGun/MuonPropagator.h"
#include "PROPOSAL-icetray/Converter.h"

namespace I3MuonGun {

MuonPropagator::MuonPropagator(const std::string& medium, double ecut, double vcut, double rho)
    : propagator_(NULL)
{
    PROPOSAL::Sector::Definition sec_def;

    sec_def.cut_settings.SetEcut(ecut);
    sec_def.cut_settings.SetVcut(vcut);

    sec_def.SetMedium(*PROPOSAL::MediumFactory::Get().CreateMedium(medium, 1.0));

    PROPOSAL::Sphere detector(PROPOSAL::Vector3D(0.0, 0.0, 0.0), 0.0, 1e18);
    sec_def.SetGeometry(detector);

    std::vector<PROPOSAL::Sector::Definition> sector_definitions;
    sector_definitions.push_back(sec_def);

    PROPOSAL::InterpolationDef interpol_def;

    std::ostringstream prefix;
    prefix << getenv("I3_BUILD") << "/MuonGun/resources/tables/icecube";

    interpol_def.path_to_tables = prefix.str();

    propagator_ = new PROPOSAL::Propagator(PROPOSAL::MuMinusDef::Get(),
                                           sector_definitions,
                                           detector,
                                           interpol_def);
}

MuonPropagator::~MuonPropagator()
{
    delete propagator_;
}

void MuonPropagator::SetSeed(int seed)
{
    PROPOSAL::RandomGenerator::Get().SetSeed(seed);
}

inline std::string GetMMCName(I3Particle::ParticleType pt)
{
    std::string name;

    if (pt == I3Particle::MuMinus)
        name = "mu-";
    else if (pt == I3Particle::MuPlus)
        name = "mu+";

    return name;
}

std::string MuonPropagator::GetName(const I3Particle& p)
{
    return GetMMCName(p.GetType());
}


/** Differential stochastic rate: d^2N/dv/dx [1/m] */
double MuonPropagator::GetStochasticRate(double energy, double fraction, I3Particle::ParticleType type) const
{
    // propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
    //
    // // Check kinematics
    // if (fraction <= 0 || energy*(1-fraction) <= propagator_->get_particle()->m*I3Units::MeV)
    // 	return 0.;
    // propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
    // propagator_->get_cros()->get_ionization()->setEnergy();
    //
    // double contrib;
    // double rate = 0.;
    // // Separate contributions from each element for brems/epair/photonuclear interactions
    // for (int i=0; i < propagator_->get_cros()->get_medium()->get_numCompontents(); i++) {
    // 	propagator_->get_cros()->set_component(i);
    // 	if (std::isfinite(contrib = propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->function(fraction))
    // && contrib > 0) 		rate += contrib; 	if (std::isfinite(contrib =
    // propagator_->get_cros()->get_epairproduction()->get_Stochastic()->function(fraction)) && contrib > 0)
    // rate
    // += contrib; 	if (std::isfinite(contrib =
    // propagator_->get_cros()->get_photonuclear()->get_Stochastic()->function(fraction)) && contrib > 0)
    // rate
    // += contrib;
    // }
    // // Only one bulk ionization contribution
    // if (std::isfinite(contrib = propagator_->get_cros()->get_ionization()->get_Stochastic()->function(fraction)) &&
    // contrib > 0) 	rate += contrib;
    // // printf("brems dN/dx: %e\n", propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx());
    // // printf("epair dN/dx: %e\n", propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx());
    // // printf("photo dN/dx: %e\n", propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx());
    // // printf("ioniz dN/dx: %e\n", propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx());
    //
    // return rate*(I3Units::m/I3Units::cm);
}

/** total stochastic rate: dN/dx [1/m] */
double MuonPropagator::GetTotalStochasticRate(double energy, I3Particle::ParticleType type) const
{
    // propagator_->get_output()->initDefault(0, 0, GetMMCName(type), 0, 0, 0, 0, 0, 0);
    // propagator_->get_particle()->setEnergy(energy/I3Units::MeV);
    // propagator_->get_cros()->get_ionization()->setEnergy();
    //
    // double rate = 0;
    // rate += propagator_->get_cros()->get_bremsstrahlung()->get_Stochastic()->dNdx();
    // rate += propagator_->get_cros()->get_epairproduction()->get_Stochastic()->dNdx();
    // rate += propagator_->get_cros()->get_photonuclear()->get_Stochastic()->dNdx();
    // rate += propagator_->get_cros()->get_ionization()->get_Stochastic()->dNdx();
    //
    // return rate*(I3Units::m/I3Units::cm);
}

I3Particle MuonPropagator::propagate(const I3Particle& p,
                                     double distance,
                                     boost::shared_ptr<std::vector<I3Particle> > losses)
{
    I3Particle endpoint(p);

    double x, y, z, theta, phi;

    PROPOSAL::Particle pp = propagator_->GetParticle();
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

    propagator_->Propagate(distance / I3Units::cm);

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

    if (losses)
    {
        std::vector<PROPOSAL::DynamicData*> history = PROPOSAL::Output::getInstance().GetSecondarys();
        for (unsigned int i = 0; i < history.size(); i++)
        {
            losses->push_back(I3PROPOSALParticleConverter::GenerateI3Particle(*history[i]));
        }
    }

    PROPOSAL::Output::getInstance().ClearSecondaryVector();

    return endpoint;
}

void Crust::AddLayer(I3Surfaces::SurfacePtr s, boost::shared_ptr<MuonPropagator> p)
{
    boundaries_.push_back(s);
    propagators_.push_back(p);
}

I3Particle Crust::Ingest(const I3Particle& p)
{
    I3Particle propped(p);
    double l = 0;
    for (unsigned i = 0; (propped.GetEnergy() > 0) && (i < boundaries_.size()); i++)
    {
        double dx = boundaries_[i]->GetIntersection(propped.GetPos(), propped.GetDir()).first;
        if (dx > 0)
            propped = (i > 0 ? propagators_[i - 1] : defaultPropagator_)->propagate(propped, dx);
        // Force lengths to measure the distance back to the outermost surface
        if (i > 0)
            l += std::min(dx, propped.GetLength());
    }
    propped.SetLength(l);

    return propped;
}

} // namespace I3MuonGun
