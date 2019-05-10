
#include "PROPOSAL/propagation_utility/ContinuousRandomizer.h"

#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityInterpolant.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//
//--------------------------------constructors--------------------------------//
//----------------------------------------------------------------------------//
//----------------------------------------------------------------------------//

ContinuousRandomizer::ContinuousRandomizer(const Utility& utility)
{
    DE2de = new UtilityIntegralContRand(utility);
}

ContinuousRandomizer::ContinuousRandomizer(const Utility& utility, const InterpolationDef interpolation_def)
{
    DE2de = new UtilityInterpolantContRand(utility, interpolation_def);
}

ContinuousRandomizer::ContinuousRandomizer(const Utility& utility, const ContinuousRandomizer& randomizer)
    : DE2de(randomizer.DE2de->clone(utility))
{
    if (utility != randomizer.DE2de->GetUtility())
    {
        log_fatal("Utilities of the ContinuousRandomizer should have same values!");
    }
}

ContinuousRandomizer::ContinuousRandomizer(const ContinuousRandomizer& randomizer)
    : DE2de(randomizer.DE2de->clone(randomizer.DE2de->GetUtility()))
{
}

ContinuousRandomizer::~ContinuousRandomizer()
{
    delete DE2de;
}

bool ContinuousRandomizer::operator==(const ContinuousRandomizer& randomizer) const
{
    return *DE2de == *randomizer.DE2de;
}

bool ContinuousRandomizer::operator!=(const ContinuousRandomizer& randomizer) const
{
    return !(*this == randomizer);
}
// ------------------------------------------------------------------------- //
// Public member function
// ------------------------------------------------------------------------- //

double ContinuousRandomizer::Randomize(double initial_energy, double final_energy, double rnd)
{
    double sigma, xhi, xlo, rndtmp;

    // this happens if small distances are propagated and the
    // energy loss is so small that it is smaller than the precision
    // which is checked for in during the calculation.
    if (initial_energy == final_energy)
    {
        return final_energy;
    }

    sigma = std::sqrt(DE2de->Calculate(initial_energy, final_energy, 0.0));

    // It is not drawn from the real gaus distribution but rather from the
    // area which is possible due to the limits of the initial energy and the
    // particle mass. Another possibility would be to draw again but that would be
    // more expensive.
    //
    // calculate the allowed region
    xhi = 0.5 + std::erf((initial_energy - final_energy) / (SQRT2 * sigma)) / 2;
    xlo = 0.5 + std::erf((DE2de->GetUtility().GetParticleDef().low - final_energy) / (SQRT2 * sigma)) / 2;

    // draw random number from the allowed region.
    rndtmp = xlo + (xhi - xlo) * rnd;

    // Calculate and return the needed value.
    return sigma * inverseErrorFunction(rndtmp) + final_energy;
}
