
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
// PhotoQ2
// ------------------------------------------------------------------------- //

CrossSection* EpairProductionFactory::CreateEpairIntegral(const ParticleDef& particle_def,
                                                          const Medium& medium,
                                                          const EnergyCutSettings& cuts,
                                                          bool lpm,
                                                          double multiplier) const
{
    return new EpairIntegral(EpairProductionRhoIntegral(particle_def, medium, cuts, lpm, multiplier));
}

CrossSection* EpairProductionFactory::CreateEpairInterpolant(const ParticleDef& particle_def,
                                                             const Medium& medium,
                                                             const EnergyCutSettings& cuts,
                                                             bool lpm,
                                                             double multiplier,
                                                             InterpolationDef def) const
{
    return new EpairInterpolant(EpairProductionRhoInterpolant(particle_def, medium, cuts, lpm, multiplier, def), def);
}


// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def) const
{
    return CreateEpairIntegral(particle_def, medium, cuts, def.lpm_effect, def.multiplier);
}

CrossSection* EpairProductionFactory::CreateEpairProduction(const ParticleDef& particle_def,
                                                            const Medium& medium,
                                                            const EnergyCutSettings& cuts,
                                                            const Definition& def,
                                                            InterpolationDef interpolation_def) const
{
    return CreateEpairInterpolant(particle_def, medium, cuts, def.lpm_effect, def.multiplier, interpolation_def);
}
