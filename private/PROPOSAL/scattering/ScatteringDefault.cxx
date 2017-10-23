
#include <boost/bind.hpp>

#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/scattering/ScatteringDefault.h"
// #include "PROPOSAL/Output.h"
#include "PROPOSAL/methods.h"
#include "PROPOSAL/Output.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/particle/Particle.h"
#include "PROPOSAL/medium/Medium.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"


using namespace PROPOSAL;
using namespace std;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ScatteringDefault::ScatteringDefault(Particle& particle, Utility& utility)
    : Scattering(particle)
    , scatter_(new UtilityIntegralScattering(utility))
{
    if (particle.GetParticleDef() != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility paricle definition!");
    }
}

ScatteringDefault::ScatteringDefault(Particle& particle, Utility& utility, InterpolationDef interpolation_def)
    : Scattering(particle)
    , scatter_(new UtilityInterpolantScattering(utility, interpolation_def))
{
    if (particle.GetParticleDef() != utility.GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the utility paricle definition!");
    }
}

ScatteringDefault::ScatteringDefault(const ScatteringDefault& scattering)
    : Scattering(scattering)
    , scatter_(scattering.scatter_->clone(scattering.scatter_->GetUtility()))
{
}

ScatteringDefault::ScatteringDefault(Particle& particle, const Utility& utility, const ScatteringDefault& scattering)
    : Scattering(particle)
    , scatter_(scattering.scatter_->clone(utility))
{
    if (particle.GetParticleDef() != scattering.GetParticle().GetParticleDef())
    {
        log_fatal("Particle definition should be equal to the scattering paricle definition!");
    }
    if (utility != scattering.scatter_->GetUtility())
    {
        log_fatal("Utilities of the ScatteringDefault should have same values!");
    }
}

ScatteringDefault::~ScatteringDefault()
{
    delete scatter_;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//-------------------------operators and swap function------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// ScatteringDefault& ScatteringDefault::operator=(const ScatteringDefault &scattering){
//     if (this != &scattering)
//     {
//       ScatteringDefault tmp(scattering);
//       swap(tmp);
//     }
//     return *this;
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

bool ScatteringDefault::compare(const Scattering& scattering) const
{
    const ScatteringDefault* scatteringDefault = dynamic_cast<const ScatteringDefault*>(&scattering);

    if (!scatteringDefault)
        return false;
    else if (*scatter_ != *scatteringDefault->scatter_)
        return false;
    else
        return true;
}

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


// void ScatteringDefault::swap(ScatteringDefault &scattering)
// {
//     using std::swap;
//
//     swap(do_interpolation_,scattering.do_interpolation_);
//
//     if(scattering.interpolant_ != NULL)
//     {
//         interpolant_->swap(*scattering.interpolant_);
//     }
//     else
//     {
//         interpolant_ = NULL;
//     }
//
//     if(scattering.interpolant_diff_ != NULL)
//     {
//         interpolant_diff_->swap(*scattering.interpolant_diff_) ;
//     }
//     else
//     {
//         interpolant_diff_ = NULL;
//     }
//
//     integral_.swap(scattering.integral_);
// }


//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//------------------------private member functions----------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//


//----------------------------------------------------------------------------//
long double ScatteringDefault::CalculateTheta0(double dr, double ei, double ef)
{
    double aux              = scatter_->Calculate(ei, ef, 0.0);
    double cutoff           = 1;
    double radiation_lenght = scatter_->GetUtility().GetMedium().GetRadiationLength();

    // TODO: check if one has to take the absolute value of the particle charge
    aux = sqrt(max(aux, 0.0) / radiation_lenght) * particle_.GetCharge();
    aux *= max(1 + 0.038 * log(dr / radiation_lenght), 0.0);

    return min(aux, cutoff);
}

//----------------------------------------------------------------------------//

Scattering::RandomAngles ScatteringDefault::CalculateRandomAngle(double dr, double ei, double ef)
{
    double Theta0, rnd1, rnd2;
    Scattering::RandomAngles random_angles;

    Theta0 = CalculateTheta0(dr, ei, ef);

    rnd1 = SQRT2 * Theta0 * erfInv(2. * (RandomGenerator::Get().RandomDouble() - 0.5));
    rnd2 = SQRT2 * Theta0 * erfInv(2. * (RandomGenerator::Get().RandomDouble() - 0.5));

    random_angles.sx = (rnd1 / SQRT3 + rnd2) / 2;
    random_angles.tx = rnd2;

    rnd1 = SQRT2 * Theta0 * erfInv(2 * (RandomGenerator::Get().RandomDouble() - 0.5));
    rnd2 = SQRT2 * Theta0 * erfInv(2 * (RandomGenerator::Get().RandomDouble() - 0.5));

    random_angles.sy = (rnd1 / SQRT3 + rnd2) / 2;
    random_angles.ty = rnd2;

    return random_angles;
}
