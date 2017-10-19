
#include <boost/algorithm/string.hpp> // case insensitive string compare for configuration file

#include "PROPOSAL/crossection/EpairIntegral.h"
#include "PROPOSAL/crossection/EpairInterpolant.h"
#include "PROPOSAL/crossection/parametrization/EpairProduction.h"
#include "PROPOSAL/crossection/factories/EpairProductionFactory.h"

#include "PROPOSAL/Output.h"

using namespace PROPOSAL;

EpairProductionFactory::EpairProductionFactory()
{
}

EpairProductionFactory::~EpairProductionFactory()
{
}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def) const
{
    return new EpairIntegral(EpairProductionRhoIntegral(particle_def, medium, cuts, def.lpm_effect, def.multiplier));
}

// ------------------------------------------------------------------------- //
CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    return new EpairInterpolant(
        EpairProductionRhoInterpolant(particle_def, medium, cuts, def.lpm_effect, def.multiplier, interpolation_def),
        interpolation_def);
}
