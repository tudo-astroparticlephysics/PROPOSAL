
#include "PROPOSAL/crossection/IonizIntegral.h"
#include "PROPOSAL/crossection/IonizInterpolant.h"
#include "PROPOSAL/crossection/factories/IonizationFactory.h"
#include "PROPOSAL/crossection/parametrization/Ionization.h"

#include "PROPOSAL/Logging.h"

using namespace PROPOSAL;

IonizationFactory::IonizationFactory() {}

IonizationFactory::~IonizationFactory() {}

// ------------------------------------------------------------------------- //
// Most general creator
// ------------------------------------------------------------------------- //

// ------------------------------------------------------------------------- //
CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  const Medium& medium,
                                                  const EnergyCutSettings& cuts,
                                                  const Definition& def) const
{
    return new IonizIntegral(Ionization(particle_def, medium, cuts, def.multiplier));
}

// ------------------------------------------------------------------------- //
CrossSection* IonizationFactory::CreateIonization(const ParticleDef& particle_def,
                                                  const Medium& medium,
                                                  const EnergyCutSettings& cuts,
                                                  const Definition& def,
                                                  InterpolationDef interpolation_def) const
{
    return new IonizInterpolant(Ionization(particle_def, medium, cuts, def.multiplier), interpolation_def);
}
